/*!
 * \file output_structure_legacy.cpp
 * \brief Main subroutines for output solver information
 * \author F. Palacios, T. Economon
 * \version 7.3.1 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2022, SU2 Contributors (cf. AUTHORS.md)
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


#include "../../include/output/COutputLegacy.hpp"

#include "../../../Common/include/geometry/CGeometry.hpp"
#include "../../include/solvers/CBaselineSolver.hpp"

COutputLegacy::COutputLegacy(CConfig *config) {

  rank = SU2_MPI::GetRank();
  size = SU2_MPI::GetSize();

  unsigned short iDim, iSpan, iMarker;

  /*--- Initialize point and connectivity counters to zero. ---*/

  /*--- Initialize turbo flag ---*/
  turbo = config->GetBoolTurbomachinery();

  if(turbo){
    /*--- Initializate quantities for turboperformace ---*/
    nSpanWiseSections = config->GetnSpanMaxAllZones();
    nMarkerTurboPerf  = config->GetnMarker_TurboPerformance();


    TotalStaticEfficiency         = new su2double*[nMarkerTurboPerf];
    TotalTotalEfficiency          = new su2double*[nMarkerTurboPerf];
    KineticEnergyLoss             = new su2double*[nMarkerTurboPerf];
    TRadius                       = new su2double*[nMarkerTurboPerf];
    TotalPressureLoss             = new su2double*[nMarkerTurboPerf];
    MassFlowIn                    = new su2double*[nMarkerTurboPerf];
    MassFlowOut                   = new su2double*[nMarkerTurboPerf];
    FlowAngleIn                   = new su2double*[nMarkerTurboPerf];
    FlowAngleIn_BC                = new su2double*[nMarkerTurboPerf];
    FlowAngleOut                  = new su2double*[nMarkerTurboPerf];
    EulerianWork                  = new su2double*[nMarkerTurboPerf];
    TotalEnthalpyIn               = new su2double*[nMarkerTurboPerf];
    TotalEnthalpyIn_BC            = new su2double*[nMarkerTurboPerf];
    EntropyIn                     = new su2double*[nMarkerTurboPerf];
    EntropyOut                    = new su2double*[nMarkerTurboPerf];
    EntropyIn_BC                  = new su2double*[nMarkerTurboPerf];
    PressureRatio                 = new su2double*[nMarkerTurboPerf];
    TotalTemperatureIn            = new su2double*[nMarkerTurboPerf];
    EnthalpyOut                   = new su2double*[nMarkerTurboPerf];
    MachIn                        = new su2double**[nMarkerTurboPerf];
    MachOut                       = new su2double**[nMarkerTurboPerf];
    VelocityOutIs                 = new su2double*[nMarkerTurboPerf];
    DensityIn                     = new su2double*[nMarkerTurboPerf];
    PressureIn                    = new su2double*[nMarkerTurboPerf];
    TurboVelocityIn               = new su2double**[nMarkerTurboPerf];
    DensityOut                    = new su2double*[nMarkerTurboPerf];
    PressureOut                   = new su2double*[nMarkerTurboPerf];
    TurboVelocityOut              = new su2double**[nMarkerTurboPerf];
    EnthalpyOutIs                 = new su2double*[nMarkerTurboPerf];
    EntropyGen                    = new su2double*[nMarkerTurboPerf];
    AbsFlowAngleIn                = new su2double*[nMarkerTurboPerf];
    TotalEnthalpyOut              = new su2double*[nMarkerTurboPerf];
    TotalEnthalpyOutIs            = new su2double*[nMarkerTurboPerf];
    RothalpyIn                    = new su2double*[nMarkerTurboPerf];
    RothalpyOut                   = new su2double*[nMarkerTurboPerf];
    AbsFlowAngleOut               = new su2double*[nMarkerTurboPerf];
    PressureOut_BC                = new su2double*[nMarkerTurboPerf];
    TemperatureIn                 = new su2double*[nMarkerTurboPerf];
    TemperatureOut                = new su2double*[nMarkerTurboPerf];
    TotalPressureIn               = new su2double*[nMarkerTurboPerf];
    TotalPressureOut              = new su2double*[nMarkerTurboPerf];
    TotalTemperatureOut           = new su2double*[nMarkerTurboPerf];
    EnthalpyIn                    = new su2double*[nMarkerTurboPerf];
    TurbIntensityIn               = new su2double*[nMarkerTurboPerf];
    Turb2LamViscRatioIn           = new su2double*[nMarkerTurboPerf];
    TurbIntensityOut              = new su2double*[nMarkerTurboPerf];
    Turb2LamViscRatioOut          = new su2double*[nMarkerTurboPerf];
    NuFactorIn                    = new su2double*[nMarkerTurboPerf];
    NuFactorOut                   = new su2double*[nMarkerTurboPerf];

    for (iMarker = 0; iMarker < nMarkerTurboPerf; iMarker++){
      TotalStaticEfficiency   [iMarker] = new su2double [nSpanWiseSections + 1];
      TotalTotalEfficiency    [iMarker] = new su2double [nSpanWiseSections + 1];
      KineticEnergyLoss       [iMarker] = new su2double [nSpanWiseSections + 1];
      TRadius                 [iMarker] = new su2double [nSpanWiseSections + 1];
      TotalPressureLoss       [iMarker] = new su2double [nSpanWiseSections + 1];
      MassFlowIn              [iMarker] = new su2double [nSpanWiseSections + 1];
      MassFlowOut             [iMarker] = new su2double [nSpanWiseSections + 1];
      FlowAngleIn             [iMarker] = new su2double [nSpanWiseSections + 1];
      FlowAngleIn_BC          [iMarker] = new su2double [nSpanWiseSections + 1];
      FlowAngleOut            [iMarker] = new su2double [nSpanWiseSections + 1];
      EulerianWork            [iMarker] = new su2double [nSpanWiseSections + 1];
      TotalEnthalpyIn         [iMarker] = new su2double [nSpanWiseSections + 1];
      TotalEnthalpyIn_BC      [iMarker] = new su2double [nSpanWiseSections + 1];
      EntropyIn               [iMarker] = new su2double [nSpanWiseSections + 1];
      EntropyOut              [iMarker] = new su2double [nSpanWiseSections + 1];
      EntropyIn_BC            [iMarker] = new su2double [nSpanWiseSections + 1];
      PressureRatio           [iMarker] = new su2double [nSpanWiseSections + 1];
      TotalTemperatureIn      [iMarker] = new su2double [nSpanWiseSections + 1];
      EnthalpyOut             [iMarker] = new su2double [nSpanWiseSections + 1];
      MachIn                  [iMarker] = new su2double*[nSpanWiseSections + 1];
      MachOut                 [iMarker] = new su2double*[nSpanWiseSections + 1];
      VelocityOutIs           [iMarker] = new su2double [nSpanWiseSections + 1];
      DensityIn               [iMarker] = new su2double [nSpanWiseSections + 1];
      PressureIn              [iMarker] = new su2double [nSpanWiseSections + 1];
      TurboVelocityIn         [iMarker] = new su2double*[nSpanWiseSections + 1];
      DensityOut              [iMarker] = new su2double [nSpanWiseSections + 1];
      PressureOut             [iMarker] = new su2double [nSpanWiseSections + 1];
      TurboVelocityOut        [iMarker] = new su2double*[nSpanWiseSections + 1];
      EnthalpyOutIs           [iMarker] = new su2double [nSpanWiseSections + 1];
      EntropyGen              [iMarker] = new su2double [nSpanWiseSections + 1];
      AbsFlowAngleIn          [iMarker] = new su2double [nSpanWiseSections + 1];
      TotalEnthalpyOut        [iMarker] = new su2double [nSpanWiseSections + 1];
      TotalEnthalpyOutIs      [iMarker] = new su2double [nSpanWiseSections + 1];
      RothalpyIn              [iMarker] = new su2double [nSpanWiseSections + 1];
      RothalpyOut             [iMarker] = new su2double [nSpanWiseSections + 1];
      AbsFlowAngleOut         [iMarker] = new su2double [nSpanWiseSections + 1];
      PressureOut_BC          [iMarker] = new su2double [nSpanWiseSections + 1];
      TemperatureIn           [iMarker] = new su2double [nSpanWiseSections + 1];
      TemperatureOut          [iMarker] = new su2double [nSpanWiseSections + 1];
      TotalPressureIn         [iMarker] = new su2double [nSpanWiseSections + 1];
      TotalPressureOut        [iMarker] = new su2double [nSpanWiseSections + 1];
      TotalTemperatureOut     [iMarker] = new su2double [nSpanWiseSections + 1];
      EnthalpyIn              [iMarker] = new su2double [nSpanWiseSections + 1];
      TurbIntensityIn         [iMarker] = new su2double [nSpanWiseSections + 1];
      Turb2LamViscRatioIn     [iMarker] = new su2double [nSpanWiseSections + 1];
      TurbIntensityOut        [iMarker] = new su2double [nSpanWiseSections + 1];
      Turb2LamViscRatioOut    [iMarker] = new su2double [nSpanWiseSections + 1];
      NuFactorIn              [iMarker] = new su2double [nSpanWiseSections + 1];
      NuFactorOut             [iMarker] = new su2double [nSpanWiseSections + 1];


      for (iSpan = 0; iSpan < nSpanWiseSections + 1; iSpan++){
        TotalStaticEfficiency   [iMarker][iSpan] = 0.0;
        TotalTotalEfficiency    [iMarker][iSpan] = 0.0;
        KineticEnergyLoss       [iMarker][iSpan] = 0.0;
        TRadius                 [iMarker][iSpan] = 0.0;
        TotalPressureLoss       [iMarker][iSpan] = 0.0;
        MassFlowIn              [iMarker][iSpan] = 0.0;
        MassFlowOut             [iMarker][iSpan] = 0.0;
        FlowAngleIn             [iMarker][iSpan] = 0.0;
        FlowAngleIn_BC          [iMarker][iSpan] = config->GetFlowAngleIn_BC();
        FlowAngleOut            [iMarker][iSpan] = 0.0;
        EulerianWork            [iMarker][iSpan] = 0.0;
        TotalEnthalpyIn         [iMarker][iSpan] = 0.0;
        TotalEnthalpyIn_BC      [iMarker][iSpan] = 0.0;
        EntropyIn               [iMarker][iSpan] = 0.0;
        EntropyOut              [iMarker][iSpan] = 0.0;
        EntropyIn_BC            [iMarker][iSpan] = 0.0;
        PressureRatio           [iMarker][iSpan] = 0.0;
        TotalTemperatureIn      [iMarker][iSpan] = 0.0;
        EnthalpyOut             [iMarker][iSpan] = 0.0;


        VelocityOutIs           [iMarker][iSpan] = 0.0;
        DensityIn               [iMarker][iSpan] = 0.0;
        PressureIn              [iMarker][iSpan] = 0.0;

        DensityOut              [iMarker][iSpan] = 0.0;
        PressureOut             [iMarker][iSpan] = 0.0;

        EnthalpyOutIs           [iMarker][iSpan] = 0.0;
        EntropyGen              [iMarker][iSpan] = 0.0;
        AbsFlowAngleIn          [iMarker][iSpan] = 0.0;
        TotalEnthalpyOut        [iMarker][iSpan] = 0.0;
        TotalEnthalpyOutIs      [iMarker][iSpan] = 0.0;
        RothalpyIn              [iMarker][iSpan] = 0.0;
        RothalpyOut             [iMarker][iSpan] = 0.0;
        AbsFlowAngleOut         [iMarker][iSpan] = 0.0;
        PressureOut_BC          [iMarker][iSpan] = config->GetPressureOut_BC();

        TemperatureIn           [iMarker][iSpan] = 0.0;
        TemperatureOut          [iMarker][iSpan] = 0.0;
        TotalPressureIn         [iMarker][iSpan] = 0.0;
        TotalPressureOut        [iMarker][iSpan] = 0.0;
        TotalTemperatureOut     [iMarker][iSpan] = 0.0;
        EnthalpyIn              [iMarker][iSpan] = 0.0;
        TurbIntensityIn         [iMarker][iSpan] = 0.0;
        Turb2LamViscRatioIn     [iMarker][iSpan] = 0.0;
        TurbIntensityOut        [iMarker][iSpan] = 0.0;
        Turb2LamViscRatioOut    [iMarker][iSpan] = 0.0;
        NuFactorIn              [iMarker][iSpan] = 0.0;
        NuFactorOut             [iMarker][iSpan] = 0.0;
        MachIn                  [iMarker][iSpan] = new su2double[4];
        MachOut                 [iMarker][iSpan] = new su2double[4];
        TurboVelocityIn         [iMarker][iSpan] = new su2double[4];
        TurboVelocityOut        [iMarker][iSpan] = new su2double[4];

        for (iDim = 0; iDim < 4; iDim++){
          MachIn           [iMarker][iSpan][iDim]   = 0.0;
          MachOut          [iMarker][iSpan][iDim]   = 0.0;
          TurboVelocityIn  [iMarker][iSpan][iDim]   = 0.0;
          TurboVelocityOut [iMarker][iSpan][iDim]   = 0.0;
        }
      }
    }
  }
}

COutputLegacy::~COutputLegacy(void) {
  /* delete pointers initialized at construction*/
  /* Coords and Conn_*(Connectivity) have their own dealloc functions */
  /* Data is taken care of in DeallocateSolution function */

  /*--- Delete turboperformance pointers initiliazed at constrction  ---*/
  unsigned short iMarker, iSpan;
  if(turbo){
    for(iMarker = 0; iMarker< nMarkerTurboPerf; iMarker++){
      for(iSpan=0; iSpan<nSpanWiseSections+1; iSpan++){
        delete [] MachIn          [iMarker][iSpan];
        delete [] MachOut         [iMarker][iSpan];
        delete [] TurboVelocityIn [iMarker][iSpan];
        delete [] TurboVelocityOut[iMarker][iSpan];
      }
    }
    for(iMarker = 0; iMarker< nMarkerTurboPerf; iMarker++){
      delete [] TotalStaticEfficiency[iMarker];
      delete [] TotalTotalEfficiency [iMarker];
      delete [] KineticEnergyLoss    [iMarker];
      delete [] TRadius              [iMarker];
      delete [] TotalPressureLoss    [iMarker];
      delete [] MassFlowIn           [iMarker];
      delete [] MassFlowOut          [iMarker];
      delete [] FlowAngleIn          [iMarker];
      delete [] FlowAngleOut         [iMarker];
      delete [] EulerianWork         [iMarker];
      delete [] TotalEnthalpyIn      [iMarker];
      delete [] TotalEnthalpyOut     [iMarker];
      delete [] TotalEnthalpyOutIs   [iMarker];
      delete [] PressureRatio        [iMarker];
      delete [] EnthalpyOut          [iMarker];
      delete [] VelocityOutIs        [iMarker];
      delete [] TotalTemperatureIn   [iMarker];
      delete [] FlowAngleIn_BC       [iMarker];
      delete [] EntropyIn            [iMarker];
      delete [] EntropyIn_BC         [iMarker];
      delete [] EntropyOut           [iMarker];
      delete [] TotalEnthalpyIn_BC   [iMarker];
      delete [] DensityIn            [iMarker];
      delete [] PressureIn           [iMarker];
      delete [] DensityOut           [iMarker];
      delete [] PressureOut          [iMarker];
      delete [] EnthalpyOutIs        [iMarker];
      delete [] EntropyGen           [iMarker];
      delete [] AbsFlowAngleIn       [iMarker];
      delete [] RothalpyIn           [iMarker];
      delete [] RothalpyOut          [iMarker];
      delete [] AbsFlowAngleOut      [iMarker];
      delete [] PressureOut_BC       [iMarker];
      delete [] MachIn               [iMarker];
      delete [] MachOut              [iMarker];
      delete [] TurboVelocityIn      [iMarker];
      delete [] TurboVelocityOut     [iMarker];
      delete [] TemperatureIn        [iMarker];
      delete [] TemperatureOut       [iMarker];
      delete [] TotalPressureIn      [iMarker];
      delete [] TotalPressureOut     [iMarker];
      delete [] TotalTemperatureOut  [iMarker];
      delete [] EnthalpyIn           [iMarker];
      delete [] TurbIntensityIn      [iMarker];
      delete [] Turb2LamViscRatioIn  [iMarker];
      delete [] TurbIntensityOut     [iMarker];
      delete [] Turb2LamViscRatioOut [iMarker];
      delete [] NuFactorIn           [iMarker];
      delete [] NuFactorOut          [iMarker];


    }
    delete [] TotalStaticEfficiency;
    delete [] TotalTotalEfficiency;
    delete [] KineticEnergyLoss;
    delete [] TRadius;
    delete [] TotalPressureLoss;
    delete [] MassFlowIn;
    delete [] MassFlowOut;
    delete [] FlowAngleIn;
    delete [] FlowAngleOut;
    delete [] EulerianWork;
    delete [] TotalEnthalpyIn;
    delete [] TotalEnthalpyOut;
    delete [] TotalEnthalpyOutIs;
    delete [] PressureRatio;
    delete [] EnthalpyOut;
    delete [] VelocityOutIs;
    delete [] TotalTemperatureIn;
    delete [] FlowAngleIn_BC;
    delete [] EntropyIn;
    delete [] EntropyIn_BC;
    delete [] EntropyOut;
    delete [] TotalEnthalpyIn_BC;
    delete [] DensityIn;
    delete [] PressureIn;
    delete [] DensityOut;
    delete [] PressureOut;
    delete [] EnthalpyOutIs;
    delete [] EntropyGen;
    delete [] AbsFlowAngleIn;
    delete [] RothalpyIn;
    delete [] RothalpyOut;
    delete [] AbsFlowAngleOut;
    delete [] PressureOut_BC;
    delete [] MachIn;
    delete [] MachOut;
    delete [] TurboVelocityIn;
    delete [] TurboVelocityOut;
    delete [] TemperatureIn;
    delete [] TemperatureOut;
    delete [] TotalPressureIn;
    delete [] TotalPressureOut;
    delete [] TotalTemperatureOut;
    delete [] EnthalpyIn;
    delete [] TurbIntensityIn;
    delete [] Turb2LamViscRatioIn;
    delete [] TurbIntensityOut;
    delete [] Turb2LamViscRatioOut;
    delete [] NuFactorIn;
    delete [] NuFactorOut;
  }
}

void COutputLegacy::SetConvHistory_Header(ofstream *ConvHist_file, CConfig *config, unsigned short val_iZone, unsigned short val_iInst) {
  char cstr[200], turb_resid[1000], adj_turb_resid[1000];
  unsigned short iMarker_Monitoring;
  string Monitoring_Tag, monitoring_coeff, aeroelastic_coeff, turbo_coeff;

  bool rotating_frame = config->GetRotating_Frame();
  bool aeroelastic = config->GetAeroelastic_Simulation();
  bool equiv_area = config->GetEquivArea();
  bool buffet = (config->GetViscous() || config->GetKind_Regime() == ENUM_REGIME::COMPRESSIBLE);
  bool engine        = ((config->GetnMarker_EngineInflow() != 0) || (config->GetnMarker_EngineExhaust() != 0));
  bool actuator_disk = ((config->GetnMarker_ActDiskInlet() != 0) || (config->GetnMarker_ActDiskOutlet() != 0));
  bool turbulent = ((config->GetKind_Solver() == MAIN_SOLVER::RANS) || (config->GetKind_Solver() == MAIN_SOLVER::ADJ_RANS) ||
                    (config->GetKind_Solver() == MAIN_SOLVER::DISC_ADJ_RANS) || (config->GetKind_Solver() == MAIN_SOLVER::DISC_ADJ_INC_RANS) ||
                    (config->GetKind_Solver() == MAIN_SOLVER::INC_RANS));
  bool cont_adj = config->GetContinuous_Adjoint();
  bool disc_adj = config->GetDiscrete_Adjoint();
  bool frozen_visc = (cont_adj && config->GetFrozen_Visc_Cont()) ||( disc_adj && config->GetFrozen_Visc_Disc());
  bool inv_design = (config->GetInvDesign_Cp() || config->GetInvDesign_HeatFlux());

  bool output_surface = (config->GetnMarker_Analyze() != 0);
  bool output_comboObj = (config->GetnObj() > 1);
  bool output_per_surface = true;
  bool turbo = config->GetBoolTurbomachinery();
  unsigned short direct_diff = config->GetDirectDiff();

  bool compressible = (config->GetKind_Regime() == ENUM_REGIME::COMPRESSIBLE);
  bool incompressible = (config->GetKind_Regime() == ENUM_REGIME::INCOMPRESSIBLE);
  bool incload = config->GetIncrementalLoad();

  bool thermal = false; /* Flag for whether to print heat flux values */
  bool weakly_coupled_heat = config->GetWeakly_Coupled_Heat();
  bool radiation            = config->AddRadiation();

  if (config->GetKind_Solver() == MAIN_SOLVER::RANS || config->GetKind_Solver()  == MAIN_SOLVER::NAVIER_STOKES ||
      config->GetKind_Solver() == MAIN_SOLVER::INC_RANS || config->GetKind_Solver()  == MAIN_SOLVER::INC_NAVIER_STOKES) {
    thermal = true;
  }

  /*--- Write file name with extension ---*/
  string filename = config->GetConv_FileName();
  string hist_ext = ".csv";
  if (config->GetTabular_FileFormat() == TAB_OUTPUT::TAB_TECPLOT) hist_ext = ".dat";

  if(config->GetnZone() > 1){
    filename = config->GetMultizone_HistoryFileName(filename, val_iZone, hist_ext);
  }
  if(config->GetnTimeInstances() > 1){
    filename = config->GetMultiInstance_HistoryFileName(filename, val_iInst);
  }

  if (config->GetTime_Domain() && config->GetRestart()) {
    filename = config->GetUnsteady_FileName(filename, config->GetRestart_Iter(), hist_ext);
  }

  strcpy (cstr, filename.data());

  ConvHist_file->open(cstr, ios::out);
  ConvHist_file->precision(15);

  /*--- Begin of the header ---*/

  char begin[]= "\"Iteration\"";

  /*--- Header for the coefficients ---*/

  char flow_coeff[]= ",\"CL\",\"CD\",\"CSF\",\"CMx\",\"CMy\",\"CMz\",\"CFx\",\"CFy\",\"CFz\",\"CL/CD\",\"AoA\",\"Custom_ObjFunc\"";
  char heat_coeff[]= ",\"HeatFlux_Total\",\"HeatFlux_Maximum\",\"Temperature_Total\"";
  char equivalent_area_coeff[]= ",\"CEquivArea\",\"CNearFieldOF\"";
  char buffet_coeff[]= ",\"Buffet_Metric\"";
  char engine_coeff[]= ",\"NetThrust\",\"Power\",\"AeroCDrag\",\"SolidCDrag\",\"Radial_Distortion\",\"Circumferential_Distortion\"";
  char rotating_frame_coeff[]= ",\"CMerit\",\"CT\",\"CQ\"";
  char fem_coeff[]= ",\"VM_Stress\",\"Force_Coeff\"";
  char fem_incload[]= ",\"IncLoad\"";
  char adj_coeff[]= ",\"Sens_Geo\",\"Sens_Mach\",\"Sens_AoA\",\"Sens_Press\",\"Sens_Temp\",\"Sens_AoS\"";
  char adj_inc_coeff[]=",\"Sens_Geo\",\"Sens_Vin\",\"Sens_Pout\",\"Sens_Temp\"";
  char adj_turbo_coeff[]=",\"Sens_Geo\",\"Sens_PressOut\",\"Sens_TotTempIn\"";
  char surface_outputs[]= ",\"Avg_MassFlow\",\"Avg_Mach\",\"Avg_Temp\",\"Avg_Press\",\"Avg_Density\",\"Avg_Enthalpy\",\"Avg_NormalVel\",\"Uniformity\",\"Secondary_Strength\",\"Momentum_Distortion\",\"Secondary_Over_Uniformity\",\"Avg_TotalTemp\",\"Avg_TotalPress\",\"Pressure_Drop\"";
  char Cp_inverse_design[]= ",\"Cp_Diff\"";
  char Heat_inverse_design[]= ",\"HeatFlux_Diff\"";
  char d_flow_coeff[] = ",\"D(CL)\",\"D(CD)\",\"D(CSF)\",\"D(CMx)\",\"D(CMy)\",\"D(CMz)\",\"D(CFx)\",\"D(CFy)\",\"D(CFz)\",\"D(CL/CD)\",\"D(Custom_ObjFunc)\"";
  char d_thermal_coeff[] = ",\"D(HeatFlux_Total)\",\"D(HeatFlux_Maximum)\"";
  char d_engine[] = ",\"D(NetThrust)\",\"D(Power)\",\"D(AeroCDrag)\",\"D(SolidCDrag)\",\"D(Radial_Distortion)\",\"D(Circumferential_Distortion)\"";
  char d_turbo_coeff[] = ",\"D(TotalPressureLoss_0)\",\"D(FlowAngleOut_0)\",\"D(TotalEfficency)\",\"D(TotalStaticEfficiency)\", \"D(EntropyGen)\"";
  char d_surface_outputs[]= ",\"D(Uniformity)\",\"D(Secondary_Strength)\",\"D(Momentum_Distortion)\",\"D(Secondary_Over_Uniformity)\",\"D(Pressure_Drop)\"";

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
    if(buffet) monitoring_coeff += ",\"Buffet_Metric_"    + Monitoring_Tag + "\"";
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
    case TURB_MODEL::SA:case TURB_MODEL::SA_NEG:case TURB_MODEL::SA_E: case TURB_MODEL::SA_COMP: case TURB_MODEL::SA_E_COMP:
      SPRINTF (turb_resid, ",\"Res_Turb[0]\"");
      break;
    case TURB_MODEL::SST:case TURB_MODEL::SST_SUST:
      SPRINTF (turb_resid, ",\"Res_Turb[0]\",\"Res_Turb[1]\"");
      break;
    default: break;
  }
  switch (config->GetKind_Turb_Model()) {
    case TURB_MODEL::SA:case TURB_MODEL::SA_NEG:case TURB_MODEL::SA_E: case TURB_MODEL::SA_COMP: case TURB_MODEL::SA_E_COMP:
      SPRINTF (adj_turb_resid, ",\"Res_AdjTurb[0]\"");
      break;
    case TURB_MODEL::SST:case TURB_MODEL::SST_SUST:
      SPRINTF (adj_turb_resid, ",\"Res_AdjTurb[0]\",\"Res_AdjTurb[1]\"");
      break;
    default: break;
  }
  char fem_resid[]= ",\"Res_FEM[0]\",\"Res_FEM[1]\",\"Res_FEM[2]\"";
  char heat_resid[]= ",\"Res_Heat\"";
  char rad_resid[]= ",\"Res_P1-rad\"";

  /*--- End of the header ---*/

  char end[]= ",\"Linear_Solver_Iterations\",\"CFL_Number\",\"Time(min)\"\n";
  char endfea[]= ",\"Linear_Solver_Iterations\",\"Time(min)\"\n";

  if ((config->GetTabular_FileFormat() == TAB_OUTPUT::TAB_TECPLOT)) {
    ConvHist_file[0] << "TITLE = \"SU2 Simulation\"" << endl;
    ConvHist_file[0] << "VARIABLES = ";
  }

  /*--- Write the header, case depending ---*/

  switch (config->GetKind_Solver()) {

    case MAIN_SOLVER::EULER : case MAIN_SOLVER::NAVIER_STOKES: case MAIN_SOLVER::RANS :
    case MAIN_SOLVER::INC_EULER : case MAIN_SOLVER::INC_NAVIER_STOKES: case MAIN_SOLVER::INC_RANS :
    case MAIN_SOLVER::FEM_EULER : case MAIN_SOLVER::FEM_NAVIER_STOKES: case MAIN_SOLVER::FEM_RANS : case MAIN_SOLVER::FEM_LES:
      ConvHist_file[0] << begin;
      if (!turbo) ConvHist_file[0] << flow_coeff;
      if (buffet) ConvHist_file[0] << buffet_coeff;
      if (turbo) ConvHist_file[0] << turbo_coeff;
      if (thermal && !turbo) ConvHist_file[0] << heat_coeff;
      if (equiv_area) ConvHist_file[0] << equivalent_area_coeff;
      if (engine || actuator_disk) ConvHist_file[0] << engine_coeff;
      if (inv_design) {
        ConvHist_file[0] << Cp_inverse_design;
        if (thermal && !turbo) ConvHist_file[0] << Heat_inverse_design;
      }
      if (rotating_frame && !turbo) ConvHist_file[0] << rotating_frame_coeff;

      ConvHist_file[0] << flow_resid;
      if (turbulent) ConvHist_file[0] << turb_resid;
      if (weakly_coupled_heat) ConvHist_file[0] << heat_resid;
      if (radiation) ConvHist_file[0] << rad_resid;
      if (aeroelastic) ConvHist_file[0] << aeroelastic_coeff;
      if (output_per_surface) ConvHist_file[0] << monitoring_coeff;
      if (output_surface) ConvHist_file[0] << surface_outputs;
      if (direct_diff != NO_DERIVATIVE) {
        if (!turbo) ConvHist_file[0] << d_flow_coeff;
        else        ConvHist_file[0] << d_turbo_coeff;
        if (engine || actuator_disk) ConvHist_file[0] << d_engine;
        if (thermal) ConvHist_file[0] << d_thermal_coeff;
        if (output_surface) ConvHist_file[0] << d_surface_outputs;
      }
      if (output_comboObj) ConvHist_file[0] << combo_obj;
      ConvHist_file[0] << end;

      break;

    case MAIN_SOLVER::ADJ_EULER      : case MAIN_SOLVER::ADJ_NAVIER_STOKES      : case MAIN_SOLVER::ADJ_RANS:
    case MAIN_SOLVER::DISC_ADJ_EULER: case MAIN_SOLVER::DISC_ADJ_NAVIER_STOKES: case MAIN_SOLVER::DISC_ADJ_RANS:
    case MAIN_SOLVER::DISC_ADJ_INC_EULER: case MAIN_SOLVER::DISC_ADJ_INC_NAVIER_STOKES: case MAIN_SOLVER::DISC_ADJ_INC_RANS:
      if (!turbo) {
        if (compressible) {
          ConvHist_file[0] << begin << adj_coeff << adj_flow_resid;
        }
        if (incompressible) {
          ConvHist_file[0] << begin << adj_inc_coeff << adj_flow_resid;
        }
      }
      else ConvHist_file[0] << begin << adj_turbo_coeff << adj_flow_resid;
      if ((turbulent) && (!frozen_visc)) ConvHist_file[0] << adj_turb_resid;
      ConvHist_file[0] << end;
      break;

    case MAIN_SOLVER::HEAT_EQUATION:
      ConvHist_file[0] << begin << heat_coeff;
      ConvHist_file[0] << heat_resid << end;
      break;

    case MAIN_SOLVER::FEM_ELASTICITY:
      ConvHist_file[0] << begin << fem_coeff;
      if (incload) ConvHist_file[0] << fem_incload;
      ConvHist_file[0] << fem_resid << endfea;
      break;

    case MAIN_SOLVER::DISC_ADJ_FEM:
      ConvHist_file[0] << begin << fem_coeff;
      ConvHist_file[0] << fem_resid << endfea;
      break;

    default:
      break;
  }

  if (config->GetTabular_FileFormat() == TAB_OUTPUT::TAB_TECPLOT) {
    ConvHist_file[0] << "ZONE T= \"Convergence history\"" << endl;
  }

}


void COutputLegacy::SetConvHistory_Body(ofstream *ConvHist_file,
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
  bool fea                  = ((config[val_iZone]->GetKind_Solver()== MAIN_SOLVER::FEM_ELASTICITY)||(config[val_iZone]->GetKind_Solver()== MAIN_SOLVER::DISC_ADJ_FEM));
  unsigned long iIntIter    = config[val_iZone]->GetInnerIter();
  unsigned long iExtIter    = config[val_iZone]->GetInnerIter();
  unsigned short FinestMesh = config[val_iZone]->GetFinestMesh();
  unsigned short nZone      = config[val_iZone]->GetnZone();
  unsigned short nInst      = config[val_iZone]->GetnTimeInstances();
  bool cont_adj             = config[val_iZone]->GetContinuous_Adjoint();
  bool disc_adj             = config[val_iZone]->GetDiscrete_Adjoint();
  bool energy               = config[val_iZone]->GetEnergy_Equation();
  bool incload              = config[val_iZone]->GetIncrementalLoad();
  bool fixed_cl             = config[val_iZone]->GetFixed_CL_Mode();
  bool output_files         = true;

  bool radiation            = config[val_iZone]->AddRadiation();

  bool compressible = (config[val_iZone]->GetKind_Regime() == ENUM_REGIME::COMPRESSIBLE);
  bool incompressible = (config[val_iZone]->GetKind_Regime() == ENUM_REGIME::INCOMPRESSIBLE);

  if (!disc_adj && !cont_adj && !DualTime_Iteration) {

    if (fixed_cl &&
        (solver_container[val_iZone][val_iInst][MESH_0][FLOW_SOL]->GetStart_AoA_FD()) &&
        (iExtIter != solver_container[val_iZone][val_iInst][MESH_0][FLOW_SOL]->GetIter_Update_AoA())) {
      output_files = false;
    }

    if (fea || fluid_structure) output_files = false;

    /*--- We need to evaluate some of the objective functions to write the value on the history file ---*/

    if (((iExtIter % (config[val_iZone]->GetVolumeOutputFrequency(0) )) == 0) ||
        (!fixed_cl && (iExtIter == (config[val_iZone]->GetnInner_Iter()-1))) ||
        /*--- If CL mode we need to compute the complete solution at two very particular iterations ---*/
        (fixed_cl && (iExtIter == (config[val_iZone]->GetnInner_Iter()-2) ||
          (solver_container[val_iZone][val_iInst][MESH_0][FLOW_SOL]->GetStart_AoA_FD() &&
          iExtIter == solver_container[val_iZone][val_iInst][MESH_0][FLOW_SOL]->GetIter_Update_AoA())))) {


      if ((rank == MASTER_NODE) && output_files) cout << endl << "------------------------ Evaluate Special Output ------------------------";

      switch (config[val_iZone]->GetKind_Solver()) {

        case MAIN_SOLVER::EULER: case MAIN_SOLVER::NAVIER_STOKES: case MAIN_SOLVER::RANS:
        case MAIN_SOLVER::INC_EULER: case MAIN_SOLVER::INC_NAVIER_STOKES: case MAIN_SOLVER::INC_RANS:
          /*--- For specific applications, evaluate and plot the surface. ---*/

          if (config[val_iZone]->GetnMarker_Analyze() != 0) {
            SpecialOutput_AnalyzeSurface(solver_container[val_iZone][val_iInst][MESH_0][FLOW_SOL],
                                         geometry[val_iZone][val_iInst][MESH_0], config[val_iZone], output_files);
          }

          /*--- For specific applications, evaluate and plot the surface. ---*/

          if ((config[val_iZone]->GetnMarker_Analyze() != 0) && compressible) {
            SpecialOutput_Distortion(solver_container[val_iZone][val_iInst][MESH_0][FLOW_SOL],
                                     geometry[val_iZone][val_iInst][MESH_0], config[val_iZone], output_files);
          }

          /*--- For specific applications, evaluate and plot the cp coefficent at different stations. ---*/

          if (config[val_iZone]->GetPlot_Section_Forces()) {
            SpecialOutput_SpanLoad(solver_container[val_iZone][val_iInst][MESH_0][FLOW_SOL],
                                   geometry[val_iZone][val_iInst][MESH_0], config[val_iZone], output_files);
          }

          break;

        default:
          break;
      }

      /*--- Output a file with the forces breakdown. ---*/

      if (config[val_iZone]->GetTime_Marching() == TIME_MARCHING::HARMONIC_BALANCE) {
        SpecialOutput_HarmonicBalance(solver_container, geometry, config, val_iInst, nInst, output_files);
      }

      /*--- Compute span-wise values file for turbomachinery. ---*/

      if (config[val_iZone]->GetBoolTurbomachinery()) {
        SpecialOutput_Turbo(solver_container, geometry, config, val_iZone, output_files);
      }

      /*--- Output a file with the forces breakdown. ---*/

      SpecialOutput_ForcesBreakdown(solver_container, geometry, config, val_iZone, output_files);

      if ((rank == MASTER_NODE) && output_files) cout << "-------------------------------------------------------------------------" << endl << endl;

    }

  }

  /*--- Output using only the master node ---*/

  if (rank == MASTER_NODE) {

    /*-- Compute the total objective if a "combo" objective is used ---*/

    if (output_comboObj) {
      switch (config[val_iZone]->GetKind_Solver()) {
      case MAIN_SOLVER::EULER:                   case MAIN_SOLVER::NAVIER_STOKES:                   case MAIN_SOLVER::RANS:
      case MAIN_SOLVER::INC_EULER:               case MAIN_SOLVER::INC_NAVIER_STOKES:               case MAIN_SOLVER::INC_RANS:
        solver_container[val_iZone][val_iInst][FinestMesh][FLOW_SOL]->Evaluate_ObjFunc(config[val_iZone],
                                                                                       solver_container[val_iZone][val_iInst][FinestMesh]);
        break;
      default:
        break;
      }
    }

    unsigned long ExtIter_OffSet = config[val_iZone]->GetExtIter_OffSet();
    if (config[val_iZone]->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_1ST ||
        config[val_iZone]->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_2ND)
      ExtIter_OffSet = 0;

    /*--- WARNING: These buffers have hard-coded lengths. Note that you
     may have to adjust them to be larger if adding more entries. ---*/

    char begin[1000], direct_coeff[1000], heat_coeff[1000], equivalent_area_coeff[1000], engine_coeff[1000], rotating_frame_coeff[1000], Cp_inverse_design[1000], Heat_inverse_design[1000], surface_coeff[1000], aeroelastic_coeff[1000], monitoring_coeff[10000], buffet_coeff[1000],
    adjoint_coeff[1000], flow_resid[1000], adj_flow_resid[1000], turb_resid[1000], trans_resid[1000],
    adj_turb_resid[1000], rad_resid[1000],
    begin_fem[1000], fem_coeff[1000], heat_resid[1000], combo_obj[1000],
    fem_resid[1000], end[1000], end_fem[1000], surface_outputs[1000], d_surface_outputs[1000], d_direct_coeff[1000], turbo_coeff[10000];


    su2double dummy = 0.0;
    const su2double *Coord = nullptr;
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
    bool transition = (config[val_iZone]->GetKind_Trans_Model() == TURB_TRANS_MODEL::LM);
    bool thermal = (config[val_iZone]->GetKind_Solver() == MAIN_SOLVER::RANS || config[val_iZone]->GetKind_Solver()  == MAIN_SOLVER::NAVIER_STOKES ||
                    config[val_iZone]->GetKind_Solver() == MAIN_SOLVER::INC_NAVIER_STOKES || config[val_iZone]->GetKind_Solver() == MAIN_SOLVER::INC_RANS );
    bool turbulent = ((config[val_iZone]->GetKind_Solver() == MAIN_SOLVER::RANS) || (config[val_iZone]->GetKind_Solver() == MAIN_SOLVER::ADJ_RANS) ||
                      (config[val_iZone]->GetKind_Solver() == MAIN_SOLVER::DISC_ADJ_RANS) || config[val_iZone]->GetKind_Solver()  == MAIN_SOLVER::INC_RANS ||
                      config[val_iZone]->GetKind_Solver()  == MAIN_SOLVER::DISC_ADJ_INC_RANS);
    bool adjoint =  cont_adj || disc_adj;
    bool frozen_visc = (cont_adj && config[val_iZone]->GetFrozen_Visc_Cont()) ||( disc_adj && config[val_iZone]->GetFrozen_Visc_Disc());
    bool heat =  ((config[val_iZone]->GetKind_Solver() == MAIN_SOLVER::HEAT_EQUATION) || (config[val_iZone]->GetWeakly_Coupled_Heat()));
    bool weakly_coupled_heat = config[val_iZone]->GetWeakly_Coupled_Heat();
    bool flow = (config[val_iZone]->GetKind_Solver() == MAIN_SOLVER::EULER) || (config[val_iZone]->GetKind_Solver() == MAIN_SOLVER::NAVIER_STOKES) ||
                (config[val_iZone]->GetKind_Solver() == MAIN_SOLVER::RANS) || (config[val_iZone]->GetKind_Solver() == MAIN_SOLVER::FEM_EULER) ||
                (config[val_iZone]->GetKind_Solver() == MAIN_SOLVER::FEM_NAVIER_STOKES) || (config[val_iZone]->GetKind_Solver() == MAIN_SOLVER::FEM_RANS) ||
                (config[val_iZone]->GetKind_Solver() == MAIN_SOLVER::FEM_LES) || (config[val_iZone]->GetKind_Solver() == MAIN_SOLVER::ADJ_EULER) ||
                (config[val_iZone]->GetKind_Solver() == MAIN_SOLVER::ADJ_NAVIER_STOKES) || (config[val_iZone]->GetKind_Solver() == MAIN_SOLVER::ADJ_RANS) ||
                (config[val_iZone]->GetKind_Solver() == MAIN_SOLVER::INC_EULER) || (config[val_iZone]->GetKind_Solver() == MAIN_SOLVER::INC_NAVIER_STOKES) ||
                (config[val_iZone]->GetKind_Solver() == MAIN_SOLVER::INC_RANS);
    bool buffet = (config[val_iZone]->GetViscous() || config[val_iZone]->GetKind_Regime() == ENUM_REGIME::COMPRESSIBLE);

    bool fem = ((config[val_iZone]->GetKind_Solver() == MAIN_SOLVER::FEM_ELASTICITY) ||          // FEM structural solver.
                (config[val_iZone]->GetKind_Solver() == MAIN_SOLVER::DISC_ADJ_FEM));
    bool linear_analysis = (config[val_iZone]->GetGeometricConditions() == STRUCT_DEFORMATION::SMALL);
    bool nonlinear_analysis = (config[val_iZone]->GetGeometricConditions() == STRUCT_DEFORMATION::LARGE);
    bool fsi = (config[val_iZone]->GetFSI_Simulation());
    bool discadj_fem = (config[val_iZone]->GetKind_Solver() == MAIN_SOLVER::DISC_ADJ_FEM);

    bool turbo = config[val_iZone]->GetBoolTurbomachinery();

    unsigned short nTurboPerf  = config[val_iZone]->GetnMarker_TurboPerformance();

    bool output_per_surface = true;

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
      else if ((config[ExtraHeatOutputZone]->GetKind_Solver() != MAIN_SOLVER::HEAT_EQUATION)) {
        SU2_MPI::Error("Error in output routine: No heat solver in extra output zone.", CURRENT_FUNCTION);
      }
      else {
        extra_heat_output = true;
      }
    }

    /*--- Initialize variables to store information from all domains (direct solution) ---*/

    su2double Total_CL = 0.0, Total_CD = 0.0, Total_CSF = 0.0, Total_CMx = 0.0, Total_CMy = 0.0, Total_CMz = 0.0, Total_CEff = 0.0,
    Total_CEquivArea = 0.0, Total_CNearFieldOF = 0.0, Total_CFx = 0.0, Total_CFy = 0.0, Total_CFz = 0.0, Total_CMerit = 0.0,
    Total_CT = 0.0, Total_CQ = 0.0,
    Total_Heat = 0.0, Total_MaxHeat = 0.0, Total_Temperature = 0.0, Total_Custom_ObjFunc = 0.0,
    Total_ComboObj = 0.0, Total_NetThrust = 0.0, Total_Power = 0.0, Total_AeroCD = 0.0, Total_SolidCD = 0.0, Total_IDR = 0.0, Total_IDC = 0.0,
    Total_AoA = 0.0, Total_Buffet_Metric = 0.0;
    su2double Surface_MassFlow = 0.0, Surface_Mach = 0.0, Surface_Temperature = 0.0, Surface_Pressure = 0.0, Surface_Density = 0.0, Surface_Enthalpy = 0.0, Surface_NormalVelocity = 0.0, Surface_TotalTemperature = 0.0, Surface_TotalPressure = 0.0, Surface_Uniformity = 0.0, Surface_SecondaryStrength = 0.0,Surface_MomentumDistortion = 0.0, Surface_SecondOverUniform = 0.0, Surface_PressureDrop = 0.0;

    su2double Total_ForceCoeff = 0.0, Total_VMStress = 0.0, Total_IncLoad = 0.0;
    su2double Total_SensE = 0.0, Total_SensNu = 0.0;

    unsigned short iSpan;

    /*--- Initialize variables to store information from all domains (adjoint solution) ---*/
    su2double Total_Sens_Geo = 0.0, Total_Sens_Mach = 0.0, Total_Sens_AoA = 0.0;
    su2double Total_Sens_Press = 0.0, Total_Sens_Temp = 0.0;

    su2double Total_Sens_BPressure = 0.0;
    su2double Total_Sens_ModVel = 0.0;

    /*--- Initialize variables to store information from all domains (direct differentiation) ---*/
    su2double D_Total_CL = 0.0, D_Total_CD = 0.0, D_Total_CSF = 0.0, D_Total_CMx = 0.0, D_Total_CMy = 0.0, D_Total_CMz = 0.0, D_Total_CEff = 0.0, D_Total_CFx = 0.0,
        D_Total_CFy = 0.0, D_Total_CFz = 0.0, D_Total_NetThrust = 0.0, D_Total_Power = 0.0, D_Total_AeroCD = 0.0, D_Total_SolidCD = 0.0, D_Total_IDR = 0.0, D_Total_IDC = 0.0, D_Total_Custom_ObjFunc = 0.0, D_Total_Heat = 0.0, D_Total_MaxHeat = 0.0,
        D_TotalPressure_Loss = 0.0, D_FlowAngle_Out = 0.0, D_TotalStaticEfficiency = 0.0,
        D_TotalTotalEfficiency = 0.0, D_EntropyGen = 0.0,
        D_Surface_Uniformity = 0.0, D_Surface_SecondaryStrength = 0.0, D_Surface_MomentumDistortion = 0.0, D_Surface_SecondOverUniform = 0.0, D_Surface_PressureDrop = 0.0;

    /*--- Residual arrays ---*/
    su2double *residual_flow         = nullptr,
    *residual_turbulent    = nullptr,
    *residual_transition   = nullptr;
    su2double *residual_adjflow      = nullptr,
    *residual_adjturbulent = nullptr,
    *residual_adjheat = nullptr;
    su2double *residual_fea          = nullptr;
    su2double *residual_fem          = nullptr;
    su2double *residual_heat         = nullptr;
    su2double *residual_rad          = nullptr;

    /*--- Coefficients Monitored arrays ---*/
    su2double *aeroelastic_plunge = nullptr,
    *aeroelastic_pitch     = nullptr,
    *Surface_CL            = nullptr,
    *Surface_CD            = nullptr,
    *Surface_CSF           = nullptr,
    *Surface_CEff          = nullptr,
    *Surface_CFx           = nullptr,
    *Surface_CFy           = nullptr,
    *Surface_CFz           = nullptr,
    *Surface_CMx           = nullptr,
    *Surface_CMy           = nullptr,
    *Surface_CMz           = nullptr,
    *Surface_Buffet_Metric = nullptr;

    /*--- Initialize number of variables ---*/
    unsigned short nVar_Flow = 0, nVar_Turb = 0,
    nVar_Trans = 0, nVar_Heat = 0,
    nVar_AdjFlow = 0, nVar_AdjTurb = 0, nVar_AdjHeat = 0,
    nVar_FEM = 0, nVar_Rad = 0;

    /*--- Direct problem variables ---*/
    if (compressible) nVar_Flow = nDim+2; else nVar_Flow = nDim+2;
    if (turbulent) {
      switch (config[val_iZone]->GetKind_Turb_Model()) {
        case TURB_MODEL::SA: case TURB_MODEL::SA_NEG: case TURB_MODEL::SA_E: case TURB_MODEL::SA_E_COMP: case TURB_MODEL::SA_COMP: nVar_Turb = 1; break;
        case TURB_MODEL::SST: case TURB_MODEL::SST_SUST: nVar_Turb = 2; break;
        default: break;
      }
    }
    if (transition) nVar_Trans = 2;
    if (heat) nVar_Heat = 1;

    if (fem) {
      if (linear_analysis) nVar_FEM = nDim;
      if (nonlinear_analysis) nVar_FEM = 3;

      if (config[val_iZone]->GetKind_Solver() == MAIN_SOLVER::DISC_ADJ_FEM) nVar_FEM = nDim;

    }

    if (radiation) nVar_Rad = 1;

    /*--- Adjoint problem variables ---*/
    if (compressible) nVar_AdjFlow = nDim+2; else nVar_AdjFlow = nDim+2;
    if (turbulent) {
      switch (config[val_iZone]->GetKind_Turb_Model()) {
        case TURB_MODEL::SA: case TURB_MODEL::SA_NEG: case TURB_MODEL::SA_E: case TURB_MODEL::SA_E_COMP: case TURB_MODEL::SA_COMP: nVar_AdjTurb = 1; break;
        case TURB_MODEL::SST: case TURB_MODEL::SST_SUST: nVar_AdjTurb = 2; break;
        default: break;
      }
    }
    if (weakly_coupled_heat) nVar_AdjHeat = 1;

    /*--- Allocate memory for the residual ---*/
    residual_flow       = new su2double[nVar_Flow];
    residual_turbulent  = new su2double[nVar_Turb];
    residual_transition = new su2double[nVar_Trans];
    residual_heat       = new su2double[nVar_Heat];
    residual_fem        = new su2double[nVar_FEM];
    residual_rad        = new su2double[nVar_Rad];
    residual_adjflow      = new su2double[nVar_AdjFlow];
    residual_adjturbulent = new su2double[nVar_AdjTurb];
    residual_adjheat      = new su2double[nVar_AdjHeat];

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
    if(buffet) Surface_Buffet_Metric = new su2double[config[ZONE_0]->GetnMarker_Monitoring()];

    /*--- Write information from nodes ---*/

    switch (config[val_iZone]->GetKind_Solver()) {

      case MAIN_SOLVER::EULER:                   case MAIN_SOLVER::NAVIER_STOKES:                   case MAIN_SOLVER::RANS:
      case MAIN_SOLVER::INC_EULER:               case MAIN_SOLVER::INC_NAVIER_STOKES:               case MAIN_SOLVER::INC_RANS:
      case MAIN_SOLVER::FEM_EULER:               case MAIN_SOLVER::FEM_NAVIER_STOKES:               case MAIN_SOLVER::FEM_RANS:      case MAIN_SOLVER::FEM_LES:
      case MAIN_SOLVER::ADJ_EULER:               case MAIN_SOLVER::ADJ_NAVIER_STOKES:               case MAIN_SOLVER::ADJ_RANS:
      case MAIN_SOLVER::DISC_ADJ_EULER:          case MAIN_SOLVER::DISC_ADJ_NAVIER_STOKES:          case MAIN_SOLVER::DISC_ADJ_RANS:
      case MAIN_SOLVER::DISC_ADJ_INC_EULER:      case MAIN_SOLVER::DISC_ADJ_INC_NAVIER_STOKES:      case MAIN_SOLVER::DISC_ADJ_INC_RANS:

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

        if (thermal) {
          Total_Heat         = solver_container[val_iZone][val_iInst][FinestMesh][FLOW_SOL]->GetTotal_HeatFlux();
          Total_MaxHeat      = solver_container[val_iZone][val_iInst][FinestMesh][FLOW_SOL]->GetTotal_MaxHeatFlux();
          Total_Temperature  = solver_container[val_iZone][val_iInst][FinestMesh][FLOW_SOL]->GetTotal_AvgTemperature();

          if(weakly_coupled_heat) {
            Total_Heat         = solver_container[val_iZone][val_iInst][FinestMesh][HEAT_SOL]->GetTotal_HeatFlux();
            Total_MaxHeat      = solver_container[val_iZone][val_iInst][FinestMesh][HEAT_SOL]->GetTotal_MaxHeatFlux();
            Total_Temperature  = solver_container[val_iZone][val_iInst][FinestMesh][HEAT_SOL]->GetTotal_AvgTemperature();
          }
        }

        if(buffet){
          Total_Buffet_Metric = solver_container[val_iZone][val_iInst][FinestMesh][FLOW_SOL]->GetTotal_Buffet_Metric();
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
            D_Total_NetThrust = SU2_TYPE::GetDerivative(Total_NetThrust);
            D_Total_Power     = SU2_TYPE::GetDerivative(Total_Power);
            D_Total_AeroCD    = SU2_TYPE::GetDerivative(Total_AeroCD);
            D_Total_SolidCD   = SU2_TYPE::GetDerivative(Total_SolidCD);
            D_Total_IDR       = SU2_TYPE::GetDerivative(Total_IDR);
            D_Total_IDC       = SU2_TYPE::GetDerivative(Total_IDC);
          }

        }

        if (equiv_area) {
          Total_CEquivArea    = solver_container[val_iZone][val_iInst][FinestMesh][FLOW_SOL]->GetTotal_CEquivArea();
          Total_CNearFieldOF  = solver_container[val_iZone][val_iInst][FinestMesh][FLOW_SOL]->GetTotal_CNearFieldOF();

          Total_CEquivArea    = config[val_iZone]->GetWeightCd()*Total_CD + (1.0-config[val_iZone]->GetWeightCd())*Total_CEquivArea;
          Total_CNearFieldOF  = config[val_iZone]->GetWeightCd()*Total_CD + (1.0-config[val_iZone]->GetWeightCd())*Total_CNearFieldOF;
        }

        if (engine || actuator_disk) {
          Total_NetThrust = solver_container[val_iZone][val_iInst][FinestMesh][FLOW_SOL]->GetTotal_NetThrust();
          Total_Power     = solver_container[val_iZone][val_iInst][FinestMesh][FLOW_SOL]->GetTotal_Power();
          Total_AeroCD    = solver_container[val_iZone][val_iInst][FinestMesh][FLOW_SOL]->GetTotal_AeroCD();
          Total_SolidCD   = solver_container[val_iZone][val_iInst][FinestMesh][FLOW_SOL]->GetTotal_SolidCD();
          Total_IDR       = solver_container[val_iZone][val_iInst][FinestMesh][FLOW_SOL]->GetTotal_IDR();
          Total_IDC       = solver_container[val_iZone][val_iInst][FinestMesh][FLOW_SOL]->GetTotal_IDC();
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

            if(buffet) Surface_Buffet_Metric[iMarker_Monitoring] = solver_container[val_iZone][val_iInst][FinestMesh][FLOW_SOL]->GetSurface_Buffet_Metric(iMarker_Monitoring);
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

          if (direct_diff != NO_DERIVATIVE){
            D_Surface_Uniformity = SU2_TYPE::GetDerivative(Surface_Uniformity);
            D_Surface_SecondaryStrength = SU2_TYPE::GetDerivative(Surface_SecondaryStrength);
            D_Surface_MomentumDistortion = SU2_TYPE::GetDerivative(Surface_MomentumDistortion);
            D_Surface_SecondOverUniform = SU2_TYPE::GetDerivative(Surface_SecondOverUniform);
            D_Surface_PressureDrop = SU2_TYPE::GetDerivative(Surface_PressureDrop);
          }
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


        /*--- FEA residual ---*/
        //        if (fluid_structure) {
        //          for (iVar = 0; iVar < nVar_FEA; iVar++)
        //            residual_fea[iVar] = solver_container[ZONE_0][FinestMesh][FEA_SOL]->GetRes_RMS(iVar);
        //        }

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
          Total_Sens_ModVel    = solver_container[val_iZone][val_iInst][FinestMesh][ADJFLOW_SOL]->GetTotal_Sens_ModVel();

          if (radiation){
            Total_Sens_Temp      += solver_container[val_iZone][val_iInst][FinestMesh][ADJRAD_SOL]->GetTotal_Sens_Temp();
          }

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

          if (weakly_coupled_heat) {
            for (iVar = 0; iVar < nVar_Heat; iVar++) {
              residual_adjheat[iVar] = solver_container[val_iZone][val_iInst][FinestMesh][ADJHEAT_SOL]->GetRes_RMS(iVar);
            }
          }

        }

        break;


      case MAIN_SOLVER::HEAT_EQUATION:

        /*--- Heat coefficients  ---*/

        Total_Heat         = solver_container[val_iZone][val_iInst][FinestMesh][HEAT_SOL]->GetTotal_HeatFlux();
        Total_MaxHeat      = solver_container[val_iZone][val_iInst][FinestMesh][HEAT_SOL]->GetTotal_MaxHeatFlux();
        Total_Temperature  = solver_container[val_iZone][val_iInst][FinestMesh][HEAT_SOL]->GetTotal_AvgTemperature();

        /*--- Heat Residuals ---*/

        for (iVar = 0; iVar < nVar_Heat; iVar++) {
          residual_heat[iVar] = solver_container[val_iZone][val_iInst][FinestMesh][HEAT_SOL]->GetRes_RMS(iVar);
        }

        break;

      case MAIN_SOLVER::FEM_ELASTICITY:

        /*--- FEM coefficients -- As of now, this is the Von Mises Stress ---*/

        Total_VMStress = solver_container[val_iZone][val_iInst][FinestMesh][FEA_SOL]->GetTotal_CFEA();

        Total_ForceCoeff = solver_container[val_iZone][val_iInst][FinestMesh][FEA_SOL]->GetForceCoeff();

        Total_IncLoad = solver_container[val_iZone][val_iInst][FinestMesh][FEA_SOL]->GetLoad_Increment();

        LinSolvIter = (unsigned long) solver_container[val_iZone][val_iInst][FinestMesh][FEA_SOL]->GetIterLinSolver();

        /*--- Residuals: ---*/
        /*--- Linear analysis: RMS of the displacements in the nDim coordinates ---*/
        /*--- Nonlinear analysis: UTOL, RTOL and DTOL (defined in the Postprocessing function) ---*/

        if (linear_analysis) {
          for (iVar = 0; iVar < nVar_FEM; iVar++) {
            residual_fem[iVar] = solver_container[val_iZone][val_iInst][FinestMesh][FEA_SOL]->GetRes_RMS(iVar);
          }
        }
        else if (nonlinear_analysis) {
          for (iVar = 0; iVar < nVar_FEM; iVar++) {
            residual_fem[iVar] = solver_container[val_iZone][val_iInst][FinestMesh][FEA_SOL]->GetRes_FEM(iVar);
          }
        }

        break;

      case MAIN_SOLVER::DISC_ADJ_FEM:

        /*--- FEM coefficients -- As of now, this is the Von Mises Stress ---*/

        Total_VMStress = solver_container[val_iZone][val_iInst][FinestMesh][FEA_SOL]->GetTotal_CFEA();

        /*--- Residuals: ---*/
        /*--- Linear analysis: RMS of the displacements in the nDim coordinates ---*/
        /*--- Nonlinear analysis: UTOL, RTOL and DTOL (defined in the Postprocessing function) ---*/
         for (iVar = 0; iVar < nVar_FEM; iVar++) {
           residual_fem[iVar] = solver_container[val_iZone][val_iInst][FinestMesh][ADJFEA_SOL]->GetRes_RMS(iVar);
         }

        break;

      default:
        break;

    }

    if (extra_heat_output) {
      Extra_Total_Heat      = solver_container[ExtraHeatOutputZone][val_iInst][FinestMesh][HEAT_SOL]->GetTotal_HeatFlux();
      //Extra_Total_Temperature   = solver_container[ExtraHeatOutputZone][val_iInst][FinestMesh][HEAT_SOL]->GetTotal_Temperature();
      Extra_Heat_Residual   = log10(solver_container[ExtraHeatOutputZone][val_iInst][FinestMesh][HEAT_SOL]->GetRes_RMS(0));
    }

    if (radiation){
      if (disc_adj){
        for (iVar = 0; iVar < nVar_Rad; iVar++) {
          residual_rad[iVar] = solver_container[val_iZone][val_iInst][FinestMesh][ADJRAD_SOL]->GetRes_RMS(iVar);
        }
      }
      else{
        for (iVar = 0; iVar < nVar_Rad; iVar++) {
          residual_rad[iVar] = solver_container[val_iZone][val_iInst][FinestMesh][RAD_SOL]->GetRes_RMS(iVar);
        }
      }
    }

    /*--- Header frequency ---*/

    bool Unsteady = ((config[val_iZone]->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_1ST) ||
                     (config[val_iZone]->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_2ND));
    bool In_NoDualTime = (!DualTime_Iteration && (iExtIter % config[val_iZone]->GetScreen_Wrt_Freq(2) == 0));
    bool In_DualTime_0 = (DualTime_Iteration && (iIntIter % config[val_iZone]->GetScreen_Wrt_Freq(0) == 0));
    bool In_DualTime_1 = (!DualTime_Iteration && Unsteady);
    bool In_DualTime_2 = (Unsteady && DualTime_Iteration && (iExtIter % config[val_iZone]->GetScreen_Wrt_Freq(2) == 0));
    bool In_DualTime_3 = (Unsteady && !DualTime_Iteration && (iExtIter % config[val_iZone]->GetScreen_Wrt_Freq(2) == 0));

    /*--- Header frequency: analogy for dynamic structural analysis ---*/
    /*--- DualTime_Iteration is a bool we receive, which is true if it comes from FEM_StructuralIteration and false from SU2_CFD ---*/
    /*--- We maintain the name, as it is an input of the function ---*/
    /*--- dynamic determines if the problem is, or not, time dependent ---*/
    bool dynamic = (config[val_iZone]->GetTime_Domain());              // Dynamic simulations.
    bool In_NoDynamic = (!DualTime_Iteration && (iExtIter % config[val_iZone]->GetScreen_Wrt_Freq(2) == 0));
    bool In_Dynamic_0 = (DualTime_Iteration && (iIntIter % config[val_iZone]->GetScreen_Wrt_Freq(0) == 0));
    bool In_Dynamic_1 = (!DualTime_Iteration && nonlinear_analysis);
    bool In_Dynamic_2 = (nonlinear_analysis && DualTime_Iteration && (iExtIter % config[val_iZone]->GetScreen_Wrt_Freq(2) == 0));
    bool In_Dynamic_3 = (nonlinear_analysis && !DualTime_Iteration && (iExtIter % config[val_iZone]->GetScreen_Wrt_Freq(2) == 0));

    bool write_heads;
    if (Unsteady) write_heads = (iIntIter == 0);
    else write_heads = (((iExtIter % (config[val_iZone]->GetScreen_Wrt_Freq(2)*40)) == 0));

    bool write_turbo = (((iExtIter % (config[val_iZone]->GetScreen_Wrt_Freq(2)*40)) == 0) || (iExtIter == (config[val_iZone]->GetnInner_Iter() -1)));

    /*--- Analogous for dynamic problems (as of now I separate the problems, it may be worthy to do all together later on ---*/
    bool write_heads_FEM;
    if (nonlinear_analysis) write_heads_FEM = (iIntIter == 0);
    else write_heads_FEM = (((iExtIter % (config[val_iZone]->GetScreen_Wrt_Freq(2)*40)) == 0));

    if (  (!fem && ((In_NoDualTime || In_DualTime_0 || In_DualTime_1) && (In_NoDualTime || In_DualTime_2 || In_DualTime_3))) ||
        (fem  && ( (In_NoDynamic || In_Dynamic_0 || In_Dynamic_1) && (In_NoDynamic || In_Dynamic_2 || In_Dynamic_3)))
        ) {


      /*--- Prepare the history file output, note that the dual
       time output don't write to the history file ---*/
      if (!DualTime_Iteration) {

        /*--- Write the begining of the history file ---*/
        SPRINTF(begin, "%12d", SU2_TYPE::Int(iExtIter+ExtIter_OffSet));

        /*--- Write the end of the history file ---*/
        SPRINTF (end, ", %12.10f, %12.10f, %12.10f\n", su2double(LinSolvIter), config[val_iZone]->GetCFL(MESH_0), timeused/60.0);

        /*--- Write the solution and residual of the history file ---*/
        switch (config[val_iZone]->GetKind_Solver()) {

          case MAIN_SOLVER::EULER : case MAIN_SOLVER::NAVIER_STOKES: case MAIN_SOLVER::RANS:
          case MAIN_SOLVER::INC_EULER : case MAIN_SOLVER::INC_NAVIER_STOKES: case MAIN_SOLVER::INC_RANS:
          case MAIN_SOLVER::FEM_EULER : case MAIN_SOLVER::FEM_NAVIER_STOKES: case MAIN_SOLVER::FEM_RANS: case MAIN_SOLVER::FEM_LES:
          case MAIN_SOLVER::ADJ_EULER: case MAIN_SOLVER::ADJ_NAVIER_STOKES: case MAIN_SOLVER::ADJ_RANS: case MAIN_SOLVER::DISC_ADJ_EULER:
          case MAIN_SOLVER::DISC_ADJ_NAVIER_STOKES: case MAIN_SOLVER::DISC_ADJ_RANS:
          case MAIN_SOLVER::DISC_ADJ_INC_EULER: case MAIN_SOLVER::DISC_ADJ_INC_NAVIER_STOKES: case MAIN_SOLVER::DISC_ADJ_INC_RANS:

            /*--- Direct coefficients ---*/
            SPRINTF (direct_coeff, ", %14.8e, %14.8e, %14.8e, %14.8e, %14.8e, %14.8e, %14.8e, %14.8e, %14.8e, %14.8e, %14.8e, %14.8e",
                     Total_CL, Total_CD, Total_CSF, Total_CMx, Total_CMy, Total_CMz, Total_CFx, Total_CFy,
                     Total_CFz, Total_CEff, Total_AoA, Total_Custom_ObjFunc);
            if (buffet) SPRINTF (buffet_coeff, ", %14.8e",  Total_Buffet_Metric);
            if (thermal || heat) SPRINTF (heat_coeff, ", %14.8e, %14.8e, %14.8e",  Total_Heat, Total_MaxHeat, Total_Temperature);
            if (equiv_area) SPRINTF (equivalent_area_coeff, ", %14.8e, %14.8e", Total_CEquivArea, Total_CNearFieldOF);
            if (engine || actuator_disk) SPRINTF (engine_coeff, ", %14.8e, %14.8e, %14.8e, %14.8e, %14.8e, %14.8e", Total_NetThrust, Total_Power, Total_AeroCD, Total_SolidCD, Total_IDR, Total_IDC);
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
              SPRINTF (d_direct_coeff, ", %14.8e, %14.8e, %14.8e, %14.8e, %14.8e, %14.8e, %14.8e, %14.8e, %14.8e, %14.8e, %14.8e, %14.8e, %14.8e, %14.8e, %14.8e, %14.8e, %14.8e",
                       D_Total_CL, D_Total_CD, D_Total_CSF, D_Total_CMx, D_Total_CMy, D_Total_CMz, D_Total_CFx, D_Total_CFy,
                       D_Total_CFz, D_Total_CEff, D_Total_Custom_ObjFunc, D_Total_NetThrust, D_Total_Power, D_Total_AeroCD, D_Total_SolidCD, D_Total_IDR, D_Total_IDC);
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

                if(buffet){
                  SPRINTF(surface_coeff, ", %12.10f", Surface_Buffet_Metric[iMarker_Monitoring]);
                  strcat(monitoring_coeff, surface_coeff);
                }
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

              if (direct_diff != NO_DERIVATIVE) {
                SPRINTF( d_surface_outputs, ", %14.8e, %14.8e, %14.8e, %14.8e, %14.8e",
                        D_Surface_Uniformity, D_Surface_SecondaryStrength, D_Surface_MomentumDistortion, D_Surface_SecondOverUniform, D_Surface_PressureDrop);
              }
            }

            /*--- Transition residual ---*/
            if (transition) {
              SPRINTF (trans_resid, ", %12.10f, %12.10f", log10(residual_transition[0]), log10(residual_transition[1]));
            }

            /*--- Combo objective ---*/
            if (output_comboObj) {
              SPRINTF(combo_obj,", %12.10f", Total_ComboObj);
            }

            /*--- Fluid structure residual ---*/
            //            if (fluid_structure) {
            //              if (nDim == 2) SPRINTF (levelset_resid, ", %12.10f, %12.10f, 0.0", log10 (residual_fea[0]), log10 (residual_fea[1]));
            //              else SPRINTF (levelset_resid, ", %12.10f, %12.10f, %12.10f", log10 (residual_fea[0]), log10 (residual_fea[1]), log10 (residual_fea[2]));
            //            }

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

            if (radiation){
              SPRINTF (rad_resid, ", %14.8e", log10 (residual_rad[0]));
            }

            break;

          case MAIN_SOLVER::HEAT_EQUATION:

            SPRINTF (direct_coeff, ", %14.8e, %14.8e, %14.8e", Total_Heat, Total_MaxHeat, Total_Temperature);
            SPRINTF (heat_resid, ", %14.8e", log10 (residual_heat[0]));

            break;

          case MAIN_SOLVER::FEM_ELASTICITY:

            SPRINTF (begin_fem, ", %14.8e", 0.0);

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

            break;

          case MAIN_SOLVER::DISC_ADJ_FEM:

            SPRINTF (direct_coeff, ", %12.10f", Total_VMStress);
            if (nDim == 2) {
              SPRINTF (fem_resid, ", %14.8e, %14.8e, %14.8e, %14.8e, %14.8e", log10 (residual_fem[0]), log10 (residual_fem[1]), dummy, dummy, dummy);
            }
            else {
              SPRINTF (fem_resid, ", %14.8e, %14.8e, %14.8e, %14.8e, %14.8e", log10 (residual_fem[0]), log10 (residual_fem[1]), log10 (residual_fem[2]), dummy, dummy);
            }

            break;
          default:
            break;
        }
      }
      if ((val_iZone == 0 && val_iInst == 0)|| fluid_structure){
        /*--- Write the screen header---*/
        if (  (!fem && ((write_heads) && !(!DualTime_Iteration && Unsteady))) ||
            (fem && ((write_heads_FEM) && !(!DualTime_Iteration && nonlinear_analysis)))
        ) {

          if (!fem) {
            if (!Unsteady && (config[val_iZone]->GetTime_Marching() != TIME_MARCHING::TIME_STEPPING)) {
              switch (config[val_iZone]->GetKind_Solver()) {
              case MAIN_SOLVER::EULER : case MAIN_SOLVER::NAVIER_STOKES: case MAIN_SOLVER::RANS:
              case MAIN_SOLVER::INC_EULER : case MAIN_SOLVER::INC_NAVIER_STOKES: case MAIN_SOLVER::INC_RANS:
              case MAIN_SOLVER::FEM_EULER : case MAIN_SOLVER::FEM_NAVIER_STOKES: case MAIN_SOLVER::FEM_RANS: case MAIN_SOLVER::FEM_LES:
              case MAIN_SOLVER::ADJ_EULER : case MAIN_SOLVER::ADJ_NAVIER_STOKES: case MAIN_SOLVER::ADJ_RANS:

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
                break;

              case MAIN_SOLVER::DISC_ADJ_EULER: case MAIN_SOLVER::DISC_ADJ_NAVIER_STOKES: case MAIN_SOLVER::DISC_ADJ_RANS:
              case MAIN_SOLVER::DISC_ADJ_INC_EULER: case MAIN_SOLVER::DISC_ADJ_INC_NAVIER_STOKES: case MAIN_SOLVER::DISC_ADJ_INC_RANS:

                cout << endl;
                cout << "------------------------ Discrete Adjoint Summary -----------------------" << endl;
                cout << "Total Geometry Sensitivity (updated every "  << config[val_iZone]->GetVolumeOutputFrequency(0) << " iterations): ";
                cout.precision(4);
                cout.setf(ios::scientific, ios::floatfield);
                cout << Total_Sens_Geo;
                cout << endl << "-------------------------------------------------------------------------" << endl;
                break;
              default:
                break;
              }
            }
            else {
              if (flow) {
                if ((config[val_iZone]->GetTime_Marching() == TIME_MARCHING::TIME_STEPPING) && (config[val_iZone]->GetUnst_CFL()== 0.0))
                {
                  cout << endl << "Min DT: " << solver_container[val_iZone][val_iInst][FinestMesh][FLOW_SOL]->GetMin_Delta_Time()<< ".Max DT: " << solver_container[val_iZone][val_iInst][FinestMesh][FLOW_SOL]->GetMax_Delta_Time() << ".Time step: " << config[val_iZone]->GetDelta_UnstTimeND() << ".";
                } else if ((config[val_iZone]->GetTime_Marching() == TIME_MARCHING::TIME_STEPPING) && (config[val_iZone]->GetUnst_CFL()!= 0.0)) {
                  cout << endl << "Min DT: " << solver_container[val_iZone][val_iInst][FinestMesh][FLOW_SOL]->GetMin_Delta_Time()<< ".Max DT: " << solver_container[val_iZone][val_iInst][FinestMesh][FLOW_SOL]->GetMax_Delta_Time() << ". Time step: " << solver_container[val_iZone][val_iInst][config[val_iZone]->GetFinestMesh()][FLOW_SOL]->GetMin_Delta_Time() << ". CFL: " << config[val_iZone]->GetUnst_CFL()<<".";
                } else {
                  cout << endl << "Min DT: " << solver_container[val_iZone][val_iInst][FinestMesh][FLOW_SOL]->GetMin_Delta_Time()<< ".Max DT: " << solver_container[val_iZone][val_iInst][FinestMesh][FLOW_SOL]->GetMax_Delta_Time() << ".Dual Time step: " << config[val_iZone]->GetDelta_UnstTimeND() << ".";
                }
              } else {
                cout << endl << "Dual Time step: " << config[val_iZone]->GetDelta_UnstTimeND() << ".";
              }
            }
          }
          else if (fem && !fsi) {
            if (dynamic) {
              cout << endl << "Simulation time: " << config[val_iZone]->GetCurrent_DynTime() << ". Time step: " << config[val_iZone]->GetDelta_DynTime() << ".";
            }
          }

          switch (config[val_iZone]->GetKind_Solver()) {
          case MAIN_SOLVER::EULER :                  case MAIN_SOLVER::NAVIER_STOKES:
          case MAIN_SOLVER::INC_EULER :              case MAIN_SOLVER::INC_NAVIER_STOKES:
          case MAIN_SOLVER::FEM_EULER : case MAIN_SOLVER::FEM_NAVIER_STOKES: case MAIN_SOLVER::FEM_LES:

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

            if (radiation){
              cout << "     Res[P1-RAD]";
            }

            cout << endl;

            break;

          case MAIN_SOLVER::RANS : case MAIN_SOLVER::INC_RANS:

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
              case TURB_MODEL::SA: case TURB_MODEL::SA_NEG: case TURB_MODEL::SA_E: case TURB_MODEL::SA_E_COMP: case TURB_MODEL::SA_COMP:        cout << "       Res[nu]"; break;
              case TURB_MODEL::SST: case TURB_MODEL::SST_SUST: cout << "     Res[kine]" << "    Res[omega]"; break;
              default: break;
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

            if (radiation){
              cout << "     Res[P1-RAD]";
            }

            cout << endl;

            break;

            case MAIN_SOLVER::HEAT_EQUATION :
              if (!Unsteady) cout << endl << " Iter" << "    Time(s)";
              else cout << endl << " IntIter" << "  ExtIter";

              cout <<  "      Res[Heat]" << "   HFlux(Total)";
              break;

            case MAIN_SOLVER::FEM_ELASTICITY :
              if (!nonlinear_analysis) cout << endl << " Iter" << "    Time(s)";
              else cout << endl << " IntIter" << " ExtIter";

              if (linear_analysis) {
                if (nDim == 2) cout << "    Res[Displx]" << "    Res[Disply]" << "      VMS(Max)"<<  endl;
                if (nDim == 3) cout << "    Res[Displx]" << "    Res[Disply]" << "    Res[Displz]" << "      VMS(Max)"<<  endl;
              }
              else if (nonlinear_analysis) {
                cout << "     Res[UTOL]" << "     Res[RTOL]" << "     Res[ETOL]"  << "      VMS(Max)"<<  endl;
              }
            break;

            case MAIN_SOLVER::ADJ_EULER :              case MAIN_SOLVER::ADJ_NAVIER_STOKES :
            case MAIN_SOLVER::DISC_ADJ_EULER:          case MAIN_SOLVER::DISC_ADJ_NAVIER_STOKES:
            case MAIN_SOLVER::DISC_ADJ_INC_EULER:      case MAIN_SOLVER::DISC_ADJ_INC_NAVIER_STOKES:

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

              if (radiation){
                cout << "     Res[Psi_P1]";
              }

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

            case MAIN_SOLVER::ADJ_RANS : case MAIN_SOLVER::DISC_ADJ_RANS: case MAIN_SOLVER::DISC_ADJ_INC_RANS:

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
                cout << " Res[Psi_Turb[0]]";
              }
              else {
                if (incompressible) {if (energy) {cout << "   Res[Psi_Temp]";}
                  else {cout << "   Res[Psi_Velx]";}}
                else cout << "     Res[Psi_E]";
              }
              if (weakly_coupled_heat) {
                cout << "     Res[Psi_E]";
              }
              if (radiation){
                cout << "    Res[Psi_P1]";
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

            case MAIN_SOLVER::DISC_ADJ_FEM :
              cout << endl << " IntIter" << " ExtIter";

              if (nDim == 2) cout << "    Res[Ux_adj]" << "    Res[Uy_adj]" << "       Sens[E]" << "      Sens[Nu]"<<  endl;
              if (nDim == 3) cout << "    Res[Ux_adj]" << "    Res[Uy_adj]" << "    Res[Uz_adj]" << "       Sens[E]" << "      Sens[Nu]"<<  endl;

           break;

          default:
            break;
          }
        }
      }

      /*--- Write the solution on the screen ---*/

      if ((val_iZone == 0 && val_iInst == 0)|| fluid_structure){
        cout.precision(6);
        cout.setf(ios::fixed, ios::floatfield);
        if (!fem) {
          if (!Unsteady) {
            cout.width(5); cout << iExtIter + ExtIter_OffSet;
            cout.width(11); cout << timeiter;

          } else if (Unsteady && DualTime_Iteration) {
            cout.width(8); cout << iIntIter;
            cout.width(8); cout << iExtIter;
          }
        }
        else if (fem) {
          if (!DualTime_Iteration) {
            if (!nonlinear_analysis) {
              cout.width(5); cout << iExtIter;
              cout.width(11); cout << timeiter;

            } else {
              cout.width(8); cout << iIntIter;
              cout.width(8); cout << iExtIter;
            }
          }
          else if (discadj_fem){
              cout.width(8); cout << iIntIter;
              cout.width(8); cout << iExtIter;
          }
        }
      }

      switch (config[val_iZone]->GetKind_Solver()) {
        case MAIN_SOLVER::EULER : case MAIN_SOLVER::NAVIER_STOKES:
        case MAIN_SOLVER::INC_EULER : case MAIN_SOLVER::INC_NAVIER_STOKES:
        case MAIN_SOLVER::FEM_EULER : case MAIN_SOLVER::FEM_NAVIER_STOKES: case MAIN_SOLVER::FEM_LES:

          /*--- Write history file ---*/

          if ((!DualTime_Iteration) && (output_files)) {
            if (!turbo) {
              ConvHist_file[0] << begin << direct_coeff;
              if (buffet) ConvHist_file[0] << buffet_coeff;
              if (thermal) ConvHist_file[0] << heat_coeff;
              if (equiv_area) ConvHist_file[0] << equivalent_area_coeff;
              if (engine || actuator_disk) ConvHist_file[0] << engine_coeff;
              if (inv_design) {
                ConvHist_file[0] << Cp_inverse_design;
                if (thermal) ConvHist_file[0] << Heat_inverse_design;
              }
              if (rotating_frame && !turbo) ConvHist_file[0] << rotating_frame_coeff;
              ConvHist_file[0] << flow_resid;
              if (weakly_coupled_heat) ConvHist_file[0] << heat_resid;
              if (radiation) ConvHist_file[0] << rad_resid;
            }
            else {
              ConvHist_file[0] << begin << turbo_coeff << flow_resid;
            }

            if (aeroelastic) ConvHist_file[0] << aeroelastic_coeff;
            if (output_per_surface) ConvHist_file[0] << monitoring_coeff;
            if (output_surface) ConvHist_file[0] << surface_outputs;
            if (direct_diff != NO_DERIVATIVE) {
              ConvHist_file[0] << d_direct_coeff;
              if (output_surface) ConvHist_file[0] << d_surface_outputs;
            }
            if (output_comboObj) ConvHist_file[0] << combo_obj;
            ConvHist_file[0] << end;
            ConvHist_file[0].flush();
          }

          /*--- Write screen output ---*/
          if ((val_iZone == 0 && val_iInst == 0)|| fluid_structure) {
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
              if (radiation) { cout.width(15); cout << log10(residual_rad[0]); }
            }
            cout << endl;
          }
          break;

        case MAIN_SOLVER::RANS : case MAIN_SOLVER::INC_RANS:

          /*--- Write history file ---*/

          if ((!DualTime_Iteration) && (output_files)) {

            if (!turbo) {
              ConvHist_file[0] << begin << direct_coeff;
              if (buffet) ConvHist_file[0] << buffet_coeff;
              if (thermal) ConvHist_file[0] << heat_coeff;
              if (equiv_area) ConvHist_file[0] << equivalent_area_coeff;
              if (engine || actuator_disk) ConvHist_file[0] << engine_coeff;
              if (inv_design) {
                ConvHist_file[0] << Cp_inverse_design;
                if (thermal) ConvHist_file[0] << Heat_inverse_design;
              }
              if (rotating_frame && !turbo) ConvHist_file[0] << rotating_frame_coeff;
              ConvHist_file[0] << flow_resid << turb_resid;
              if (weakly_coupled_heat) ConvHist_file[0] << heat_resid;
            }
            else {
              ConvHist_file[0] << begin << turbo_coeff << flow_resid << turb_resid;
            }

            if (aeroelastic) ConvHist_file[0] << aeroelastic_coeff;
            if (output_per_surface) ConvHist_file[0] << monitoring_coeff;
            if (output_surface) ConvHist_file[0] << surface_outputs;
            if (direct_diff != NO_DERIVATIVE) {
              ConvHist_file[0] << d_direct_coeff;
              if (output_surface) ConvHist_file[0] << d_surface_outputs;
            }
            if (output_comboObj) ConvHist_file[0] << combo_obj;
            ConvHist_file[0] << end;
            ConvHist_file[0].flush();
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


        case MAIN_SOLVER::HEAT_EQUATION:

          if (!DualTime_Iteration) {
            ConvHist_file[0] << begin << direct_coeff << heat_resid << end;
            ConvHist_file[0].flush();
          }
          break;

        case MAIN_SOLVER::FEM_ELASTICITY:

          if (!DualTime_Iteration) {
            ConvHist_file[0] << begin << fem_coeff << fem_resid << end_fem;
            ConvHist_file[0].flush();

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
          break;

        case MAIN_SOLVER::DISC_ADJ_FEM:

          cout.precision(6);
          cout.setf(ios::fixed, ios::floatfield);

          cout.width(15); cout << log10(residual_fem[0]);
          cout.width(15); cout << log10(residual_fem[1]);
          if (nDim == 3) { cout.width(15); cout << log10(residual_fem[2]); }

          cout.precision(4);
          cout.setf(ios::scientific, ios::floatfield);


          if (config[val_iZone]->GetnElasticityMod() == 1){
            cout.width(14); cout << solver_container[val_iZone][val_iInst][FinestMesh][ADJFEA_SOL]->GetGlobal_Sens_E(0);
            cout.width(14); cout << solver_container[val_iZone][val_iInst][FinestMesh][ADJFEA_SOL]->GetGlobal_Sens_Nu(0);
          }
          else{
            Total_SensE = 0.0; Total_SensNu = 0.0;
            for (unsigned short iVar = 0; iVar < config[val_iZone]->GetnElasticityMod(); iVar++){
                Total_SensE += solver_container[val_iZone][val_iInst][FinestMesh][ADJFEA_SOL]->GetGlobal_Sens_E(0)
                    *solver_container[val_iZone][val_iInst][FinestMesh][ADJFEA_SOL]->GetGlobal_Sens_E(0);
                Total_SensNu += solver_container[val_iZone][val_iInst][FinestMesh][ADJFEA_SOL]->GetGlobal_Sens_Nu(0)
                    *solver_container[val_iZone][val_iInst][FinestMesh][ADJFEA_SOL]->GetGlobal_Sens_Nu(0);
            }
            Total_SensE = sqrt(Total_SensE);
            Total_SensNu = sqrt(Total_SensNu);
            cout.width(14); cout << Total_SensE;
            cout.width(14); cout << Total_SensNu;
          }

          cout << endl;
          break;

        case MAIN_SOLVER::ADJ_EULER :              case MAIN_SOLVER::ADJ_NAVIER_STOKES :
        case MAIN_SOLVER::DISC_ADJ_EULER:          case MAIN_SOLVER::DISC_ADJ_NAVIER_STOKES:
        case MAIN_SOLVER::DISC_ADJ_INC_EULER:      case MAIN_SOLVER::DISC_ADJ_INC_NAVIER_STOKES:

          if (!DualTime_Iteration) {
            ConvHist_file[0] << begin << adjoint_coeff << adj_flow_resid << end;
            ConvHist_file[0].flush();
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

              if (radiation) { cout.width(15); cout << log10(residual_rad[0]); }

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

        case MAIN_SOLVER::ADJ_RANS : case MAIN_SOLVER::DISC_ADJ_RANS: case MAIN_SOLVER::DISC_ADJ_INC_RANS:

          if (!DualTime_Iteration) {
            ConvHist_file[0] << begin << adjoint_coeff << adj_flow_resid;
            if (!frozen_visc)
              ConvHist_file[0] << adj_turb_resid;
            ConvHist_file[0] << end;
            ConvHist_file[0].flush();
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
            if (weakly_coupled_heat) {
              cout.width(17); cout << log10(residual_adjheat[0]);
            }

            if (radiation) { cout.width(15); cout << log10(residual_rad[0]); }

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

        default:
          break;
      }
      cout.unsetf(ios::fixed);

    }


    delete [] residual_flow;
    delete [] residual_turbulent;
    delete [] residual_transition;
    delete [] residual_fea;
    delete [] residual_fem;
    delete [] residual_heat;
    delete [] residual_rad;

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

void COutputLegacy::SpecialOutput_ForcesBreakdown(CSolver *****solver, CGeometry ****geometry, CConfig **config, unsigned short val_iZone, bool output) const {

  char cstr[200];
  unsigned short iDim, iMarker_Monitoring;
  ofstream Breakdown_file;

  bool compressible       = (config[val_iZone]->GetKind_Regime() == ENUM_REGIME::COMPRESSIBLE);
  bool incompressible     = (config[val_iZone]->GetKind_Regime() == ENUM_REGIME::INCOMPRESSIBLE);
  bool unsteady           = (config[val_iZone]->GetTime_Marching() != TIME_MARCHING::STEADY);
  bool viscous            = config[val_iZone]->GetViscous();
  bool dynamic_grid       = config[val_iZone]->GetDynamic_Grid();
  bool gravity            = config[val_iZone]->GetGravityForce();
  bool turbulent          = config[val_iZone]->GetKind_Solver() == MAIN_SOLVER::RANS;
  bool fixed_cl           = config[val_iZone]->GetFixed_CL_Mode();
  MAIN_SOLVER Kind_Solver = config[val_iZone]->GetKind_Solver();
  TURB_MODEL Kind_Turb_Model = config[val_iZone]->GetKind_Turb_Model();
  unsigned short Ref_NonDim = config[val_iZone]->GetRef_NonDim();

  unsigned short FinestMesh = config[val_iZone]->GetFinestMesh();
  unsigned short nDim = geometry[val_iZone][INST_0][FinestMesh]->GetnDim();
  bool flow = ((config[val_iZone]->GetKind_Solver() == MAIN_SOLVER::EULER) ||
               (config[val_iZone]->GetKind_Solver() == MAIN_SOLVER::NAVIER_STOKES) ||
               (config[val_iZone]->GetKind_Solver() == MAIN_SOLVER::RANS) ||
               (config[val_iZone]->GetKind_Solver() == MAIN_SOLVER::INC_EULER) ||
               (config[val_iZone]->GetKind_Solver() == MAIN_SOLVER::INC_NAVIER_STOKES) ||
               (config[val_iZone]->GetKind_Solver() == MAIN_SOLVER::INC_RANS));

  /*--- Output the mean flow solution using only the master node ---*/

  if ((rank == MASTER_NODE) && (flow) && (output)) {

    cout << endl << "Writing the forces breakdown file ("<< config[val_iZone]->GetBreakdown_FileName() << ")." << endl;

    /*--- Initialize variables to store information from all domains (direct solution) ---*/

    su2double Total_CL = 0.0, Total_CD = 0.0, Total_CSF = 0.0,
    Total_CMx = 0.0, Total_CMy = 0.0, Total_CMz = 0.0, Total_CEff = 0.0,
    Total_CoPx = 0.0, Total_CoPy = 0.0, Total_CoPz = 0.0,
    Total_CFx = 0.0, Total_CFy = 0.0, Total_CFz = 0.0, Inv_CL = 0.0,
    Inv_CD = 0.0, Inv_CSF = 0.0, Inv_CMx = 0.0, Inv_CMy = 0.0,
    Inv_CMz = 0.0, Inv_CEff = 0.0, Inv_CFx = 0.0, Inv_CFy = 0.0, Inv_CFz =
    0.0,      Mnt_CL = 0.0,
    Mnt_CD = 0.0, Mnt_CSF = 0.0, Mnt_CMx = 0.0, Mnt_CMy = 0.0,
    Mnt_CMz = 0.0, Mnt_CEff = 0.0, Mnt_CFx = 0.0, Mnt_CFy = 0.0, Mnt_CFz =
    0.0, Visc_CL = 0.0,
    Visc_CD = 0.0, Visc_CSF = 0.0, Visc_CMx = 0.0, Visc_CMy = 0.0,
    Visc_CMz = 0.0, Visc_CEff = 0.0, Visc_CFx = 0.0, Visc_CFy = 0.0, Visc_CFz =
    0.0, *Surface_CL = nullptr, *Surface_CD = nullptr,
    *Surface_CSF = nullptr, *Surface_CEff = nullptr, *Surface_CFx = nullptr,
    *Surface_CFy = nullptr, *Surface_CFz = nullptr,
    *Surface_CMx = nullptr, *Surface_CMy = nullptr, *Surface_CMz = nullptr,
    *Surface_CL_Inv = nullptr,
    *Surface_CD_Inv = nullptr, *Surface_CSF_Inv = nullptr,
    *Surface_CEff_Inv = nullptr, *Surface_CFx_Inv = nullptr, *Surface_CFy_Inv =
    nullptr, *Surface_CFz_Inv = nullptr, *Surface_CMx_Inv = nullptr,
    *Surface_CMy_Inv = nullptr, *Surface_CMz_Inv = nullptr,
    *Surface_CL_Visc = nullptr,
    *Surface_CD_Visc = nullptr, *Surface_CSF_Visc = nullptr,
    *Surface_CEff_Visc = nullptr, *Surface_CFx_Visc = nullptr, *Surface_CFy_Visc =
    nullptr, *Surface_CFz_Visc = nullptr, *Surface_CMx_Visc = nullptr,
    *Surface_CMy_Visc = nullptr, *Surface_CMz_Visc = nullptr,
    *Surface_CL_Mnt = nullptr,
    *Surface_CD_Mnt = nullptr, *Surface_CSF_Mnt = nullptr,
    *Surface_CEff_Mnt = nullptr, *Surface_CFx_Mnt = nullptr, *Surface_CFy_Mnt =
    nullptr, *Surface_CFz_Mnt = nullptr, *Surface_CMx_Mnt = nullptr,
    *Surface_CMy_Mnt = nullptr, *Surface_CMz_Mnt = nullptr;

    /*--- Allocate memory for the coefficients being monitored ---*/

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

    Surface_CL_Inv      = new su2double[config[ZONE_0]->GetnMarker_Monitoring()];
    Surface_CD_Inv      = new su2double[config[ZONE_0]->GetnMarker_Monitoring()];
    Surface_CSF_Inv = new su2double[config[ZONE_0]->GetnMarker_Monitoring()];
    Surface_CEff_Inv       = new su2double[config[ZONE_0]->GetnMarker_Monitoring()];
    Surface_CFx_Inv        = new su2double[config[ZONE_0]->GetnMarker_Monitoring()];
    Surface_CFy_Inv        = new su2double[config[ZONE_0]->GetnMarker_Monitoring()];
    Surface_CFz_Inv        = new su2double[config[ZONE_0]->GetnMarker_Monitoring()];
    Surface_CMx_Inv        = new su2double[config[ZONE_0]->GetnMarker_Monitoring()];
    Surface_CMy_Inv        = new su2double[config[ZONE_0]->GetnMarker_Monitoring()];
    Surface_CMz_Inv        = new su2double[config[ZONE_0]->GetnMarker_Monitoring()];

    Surface_CL_Visc = new su2double[config[ZONE_0]->GetnMarker_Monitoring()];
    Surface_CD_Visc = new su2double[config[ZONE_0]->GetnMarker_Monitoring()];
    Surface_CSF_Visc =
    new su2double[config[ZONE_0]->GetnMarker_Monitoring()];
    Surface_CEff_Visc = new su2double[config[ZONE_0]->GetnMarker_Monitoring()];
    Surface_CFx_Visc = new su2double[config[ZONE_0]->GetnMarker_Monitoring()];
    Surface_CFy_Visc = new su2double[config[ZONE_0]->GetnMarker_Monitoring()];
    Surface_CFz_Visc = new su2double[config[ZONE_0]->GetnMarker_Monitoring()];
    Surface_CMx_Visc = new su2double[config[ZONE_0]->GetnMarker_Monitoring()];
    Surface_CMy_Visc = new su2double[config[ZONE_0]->GetnMarker_Monitoring()];
    Surface_CMz_Visc = new su2double[config[ZONE_0]->GetnMarker_Monitoring()];


    Surface_CL_Mnt = new su2double[config[ZONE_0]->GetnMarker_Monitoring()];
    Surface_CD_Mnt = new su2double[config[ZONE_0]->GetnMarker_Monitoring()];
    Surface_CSF_Mnt =
    new su2double[config[ZONE_0]->GetnMarker_Monitoring()];
    Surface_CEff_Mnt = new su2double[config[ZONE_0]->GetnMarker_Monitoring()];
    Surface_CFx_Mnt = new su2double[config[ZONE_0]->GetnMarker_Monitoring()];
    Surface_CFy_Mnt = new su2double[config[ZONE_0]->GetnMarker_Monitoring()];
    Surface_CFz_Mnt = new su2double[config[ZONE_0]->GetnMarker_Monitoring()];
    Surface_CMx_Mnt = new su2double[config[ZONE_0]->GetnMarker_Monitoring()];
    Surface_CMy_Mnt = new su2double[config[ZONE_0]->GetnMarker_Monitoring()];
    Surface_CMz_Mnt = new su2double[config[ZONE_0]->GetnMarker_Monitoring()];

    /*--- Flow solution coefficients ---*/

    Total_CL       = solver[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetTotal_CL();
    Total_CD       = solver[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetTotal_CD();
    Total_CSF      = solver[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetTotal_CSF();
    Total_CEff        = solver[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetTotal_CEff();
    Total_CMx         = solver[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetTotal_CMx();
    Total_CMy         = solver[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetTotal_CMy();
    Total_CMz         = solver[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetTotal_CMz();
    Total_CFx         = solver[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetTotal_CFx();
    Total_CFy         = solver[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetTotal_CFy();
    Total_CFz         = solver[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetTotal_CFz();

    if (nDim == 2) {
      Total_CoPx = solver[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetTotal_CoPx() / solver[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetTotal_CFy();
      Total_CoPy = solver[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetTotal_CoPy() / solver[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetTotal_CFx();
      Total_CoPz = 0.0;
    }
    if (nDim == 3) {
      Total_CoPx = solver[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetTotal_CoPx() / solver[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetTotal_CFz();
      Total_CoPy = 0.0;
      Total_CoPz = solver[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetTotal_CoPz() / solver[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetTotal_CFx();
    }

    if (config[ZONE_0]->GetSystemMeasurements() == US) { Total_CoPx *= 12.0; Total_CoPy *= 12.0; Total_CoPz *= 12.0; }

    /*--- Flow inviscid solution coefficients ---*/

    Inv_CL =
    solver[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetAllBound_CL_Inv();
    Inv_CD =
    solver[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetAllBound_CD_Inv();
    Inv_CSF =
    solver[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetAllBound_CSF_Inv();
    Inv_CEff =
    solver[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetAllBound_CEff_Inv();
    Inv_CMx =
    solver[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetAllBound_CMx_Inv();
    Inv_CMy =
    solver[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetAllBound_CMy_Inv();
    Inv_CMz =
    solver[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetAllBound_CMz_Inv();
    Inv_CFx =
    solver[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetAllBound_CFx_Inv();
    Inv_CFy =
    solver[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetAllBound_CFy_Inv();
    Inv_CFz =
    solver[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetAllBound_CFz_Inv();

    /*--- Flow viscous solution coefficients ---*/

    Visc_CL =
    solver[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetAllBound_CL_Visc();
    Visc_CD =
    solver[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetAllBound_CD_Visc();
    Visc_CSF =
    solver[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetAllBound_CSF_Visc();
    Visc_CEff =
    solver[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetAllBound_CEff_Visc();
    Visc_CMx =
    solver[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetAllBound_CMx_Visc();
    Visc_CMy =
    solver[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetAllBound_CMy_Visc();
    Visc_CMz =
    solver[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetAllBound_CMz_Visc();
    Visc_CFx =
    solver[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetAllBound_CFx_Visc();
    Visc_CFy =
    solver[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetAllBound_CFy_Visc();
    Visc_CFz =
    solver[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetAllBound_CFz_Visc();

    /*--- Flow momentum solution coefficients ---*/

    Mnt_CL =
    solver[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetAllBound_CL_Mnt();
    Mnt_CD =
    solver[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetAllBound_CD_Mnt();
    Mnt_CSF =
    solver[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetAllBound_CSF_Mnt();
    Mnt_CEff =
    solver[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetAllBound_CEff_Mnt();
    Mnt_CMx =
    solver[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetAllBound_CMx_Mnt();
    Mnt_CMy =
    solver[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetAllBound_CMy_Mnt();
    Mnt_CMz =
    solver[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetAllBound_CMz_Mnt();
    Mnt_CFx =
    solver[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetAllBound_CFx_Mnt();
    Mnt_CFy =
    solver[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetAllBound_CFy_Mnt();
    Mnt_CFz =
    solver[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetAllBound_CFz_Mnt();


    /*--- Look over the markers being monitored and get the desired values ---*/

    for (iMarker_Monitoring = 0;
         iMarker_Monitoring < config[ZONE_0]->GetnMarker_Monitoring();
         iMarker_Monitoring++) {
      Surface_CL[iMarker_Monitoring] =
      solver[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetSurface_CL(
                                                             iMarker_Monitoring);
      Surface_CD[iMarker_Monitoring] =
      solver[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetSurface_CD(
                                                             iMarker_Monitoring);
      Surface_CSF[iMarker_Monitoring] =
      solver[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetSurface_CSF(
                                                              iMarker_Monitoring);
      Surface_CEff[iMarker_Monitoring] =
      solver[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetSurface_CEff(
                                                               iMarker_Monitoring);
      Surface_CMx[iMarker_Monitoring] =
      solver[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetSurface_CMx(
                                                              iMarker_Monitoring);
      Surface_CMy[iMarker_Monitoring] =
      solver[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetSurface_CMy(
                                                              iMarker_Monitoring);
      Surface_CMz[iMarker_Monitoring] =
      solver[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetSurface_CMz(
                                                              iMarker_Monitoring);
      Surface_CFx[iMarker_Monitoring] =
      solver[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetSurface_CFx(
                                                              iMarker_Monitoring);
      Surface_CFy[iMarker_Monitoring] =
      solver[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetSurface_CFy(
                                                              iMarker_Monitoring);
      Surface_CFz[iMarker_Monitoring] =
      solver[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetSurface_CFz(
                                                              iMarker_Monitoring);

      Surface_CL_Inv[iMarker_Monitoring] =
      solver[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetSurface_CL_Inv(
                                                                 iMarker_Monitoring);
      Surface_CD_Inv[iMarker_Monitoring] =
      solver[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetSurface_CD_Inv(
                                                                 iMarker_Monitoring);
      Surface_CSF_Inv[iMarker_Monitoring] =
      solver[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetSurface_CSF_Inv(
                                                                  iMarker_Monitoring);
      Surface_CEff_Inv[iMarker_Monitoring] =
      solver[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetSurface_CEff_Inv(
                                                                   iMarker_Monitoring);
      Surface_CMx_Inv[iMarker_Monitoring] =
      solver[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetSurface_CMx_Inv(
                                                                  iMarker_Monitoring);
      Surface_CMy_Inv[iMarker_Monitoring] =
      solver[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetSurface_CMy_Inv(
                                                                  iMarker_Monitoring);
      Surface_CMz_Inv[iMarker_Monitoring] =
      solver[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetSurface_CMz_Inv(
                                                                  iMarker_Monitoring);
      Surface_CFx_Inv[iMarker_Monitoring] =
      solver[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetSurface_CFx_Inv(
                                                                  iMarker_Monitoring);
      Surface_CFy_Inv[iMarker_Monitoring] =
      solver[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetSurface_CFy_Inv(
                                                                  iMarker_Monitoring);
      Surface_CFz_Inv[iMarker_Monitoring] =
      solver[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetSurface_CFz_Inv(
                                                                  iMarker_Monitoring);
      Surface_CL_Visc[iMarker_Monitoring] =
      solver[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetSurface_CL_Visc(
                                                                  iMarker_Monitoring);
      Surface_CD_Visc[iMarker_Monitoring] =
      solver[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetSurface_CD_Visc(
                                                                  iMarker_Monitoring);
      Surface_CSF_Visc[iMarker_Monitoring] =
      solver[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetSurface_CSF_Visc(
                                                                   iMarker_Monitoring);
      Surface_CEff_Visc[iMarker_Monitoring] =
      solver[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetSurface_CEff_Visc(
                                                                    iMarker_Monitoring);
      Surface_CMx_Visc[iMarker_Monitoring] =
      solver[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetSurface_CMx_Visc(
                                                                   iMarker_Monitoring);
      Surface_CMy_Visc[iMarker_Monitoring] =
      solver[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetSurface_CMy_Visc(
                                                                   iMarker_Monitoring);
      Surface_CMz_Visc[iMarker_Monitoring] =
      solver[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetSurface_CMz_Visc(
                                                                   iMarker_Monitoring);
      Surface_CFx_Visc[iMarker_Monitoring] =
      solver[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetSurface_CFx_Visc(
                                                                   iMarker_Monitoring);
      Surface_CFy_Visc[iMarker_Monitoring] =
      solver[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetSurface_CFy_Visc(
                                                                   iMarker_Monitoring);
      Surface_CFz_Visc[iMarker_Monitoring] =
      solver[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetSurface_CFz_Visc(
                                                                   iMarker_Monitoring);

      Surface_CL_Mnt[iMarker_Monitoring] =
      solver[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetSurface_CL_Mnt(
                                                                 iMarker_Monitoring);
      Surface_CD_Mnt[iMarker_Monitoring] =
      solver[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetSurface_CD_Mnt(
                                                                 iMarker_Monitoring);
      Surface_CSF_Mnt[iMarker_Monitoring] =
      solver[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetSurface_CSF_Mnt(
                                                                  iMarker_Monitoring);
      Surface_CEff_Mnt[iMarker_Monitoring] =
      solver[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetSurface_CEff_Mnt(
                                                                   iMarker_Monitoring);
      Surface_CMx_Mnt[iMarker_Monitoring] =
      solver[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetSurface_CMx_Mnt(
                                                                  iMarker_Monitoring);
      Surface_CMy_Mnt[iMarker_Monitoring] =
      solver[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetSurface_CMy_Mnt(
                                                                  iMarker_Monitoring);
      Surface_CMz_Mnt[iMarker_Monitoring] =
      solver[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetSurface_CMz_Mnt(
                                                                  iMarker_Monitoring);
      Surface_CFx_Mnt[iMarker_Monitoring] =
      solver[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetSurface_CFx_Mnt(
                                                                  iMarker_Monitoring);
      Surface_CFy_Mnt[iMarker_Monitoring] =
      solver[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetSurface_CFy_Mnt(
                                                                  iMarker_Monitoring);
      Surface_CFz_Mnt[iMarker_Monitoring] =
      solver[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetSurface_CFz_Mnt(
                                                                  iMarker_Monitoring);

    }


    /*--- Write file name with extension ---*/

    string filename = config[val_iZone]->GetBreakdown_FileName();
    strcpy (cstr, filename.data());

    Breakdown_file.open(cstr, ios::out);

    Breakdown_file << "\n" <<"-------------------------------------------------------------------------" << "\n";
    Breakdown_file <<"|    ___ _   _ ___                                                      |" << "\n";
    Breakdown_file <<"|   / __| | | |_  )   Release 7.3.0  \"Blackbird\"                        |" << "\n";
    Breakdown_file <<"|   \\__ \\ |_| |/ /                                                      |" << "\n";
    Breakdown_file <<"|   |___/\\___//___|   Suite (Computational Fluid Dynamics Code)         |" << "\n";
    Breakdown_file << "|                                                                       |" << "\n";
    Breakdown_file <<"-------------------------------------------------------------------------" << "\n";
    Breakdown_file << "| SU2 Project Website: https://su2code.github.io                        |" << "\n";
    Breakdown_file << "|                                                                       |" << "\n";
    Breakdown_file << "| The SU2 Project is maintained by the SU2 Foundation                   |" << "\n";
    Breakdown_file << "| (http://su2foundation.org)                                            |" << "\n";
    Breakdown_file <<"-------------------------------------------------------------------------" << "\n";
    Breakdown_file << "| Copyright 2012-2022, SU2 Contributors (cf. AUTHORS.md)                |" << "\n";
    Breakdown_file << "|                                                                       |" << "\n";
    Breakdown_file << "| SU2 is free software; you can redistribute it and/or                  |" << "\n";
    Breakdown_file << "| modify it under the terms of the GNU Lesser General Public            |" << "\n";
    Breakdown_file << "| License as published by the Free Software Foundation; either          |" << "\n";
    Breakdown_file << "| version 2.1 of the License, or (at your option) any later version.    |" << "\n";
    Breakdown_file << "|                                                                       |" << "\n";
    Breakdown_file << "| SU2 is distributed in the hope that it will be useful,                |" << "\n";
    Breakdown_file << "| but WITHOUT ANY WARRANTY; without even the implied warranty of        |" << "\n";
    Breakdown_file << "| MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU      |" << "\n";
    Breakdown_file << "| Lesser General Public License for more details.                       |" << "\n";
    Breakdown_file << "|                                                                       |" << "\n";
    Breakdown_file << "| You should have received a copy of the GNU Lesser General Public      |" << "\n";
    Breakdown_file << "| License along with SU2. If not, see <http://www.gnu.org/licenses/>.   |" << "\n";
    Breakdown_file <<"-------------------------------------------------------------------------" << "\n";

    Breakdown_file.precision(6);

    Breakdown_file << "\n" << "\n" << "Problem definition:" << "\n" << "\n";

    switch (Kind_Solver) {
      case MAIN_SOLVER::EULER: case MAIN_SOLVER::INC_EULER:
        if (compressible) Breakdown_file << "Compressible Euler equations." << "\n";
        if (incompressible) Breakdown_file << "Incompressible Euler equations." << "\n";
        break;
      case MAIN_SOLVER::NAVIER_STOKES: case MAIN_SOLVER::INC_NAVIER_STOKES:
        if (compressible) Breakdown_file << "Compressible Laminar Navier-Stokes' equations." << "\n";
        if (incompressible) Breakdown_file << "Incompressible Laminar Navier-Stokes' equations." << "\n";
        break;
      case MAIN_SOLVER::RANS: case MAIN_SOLVER::INC_RANS:
        if (compressible) Breakdown_file << "Compressible RANS equations." << "\n";
        if (incompressible) Breakdown_file << "Incompressible RANS equations." << "\n";
        Breakdown_file << "Turbulence model: ";
        switch (Kind_Turb_Model) {
          case TURB_MODEL::SA:        Breakdown_file << "Spalart Allmaras" << "\n"; break;
          case TURB_MODEL::SA_NEG:    Breakdown_file << "Negative Spalart Allmaras" << "\n"; break;
          case TURB_MODEL::SA_E:      Breakdown_file << "Edwards Spalart Allmaras" << "\n"; break;
          case TURB_MODEL::SA_COMP:   Breakdown_file << "Compressibility Correction Spalart Allmaras" << "\n"; break;
          case TURB_MODEL::SA_E_COMP: Breakdown_file << "Compressibility Correction Edwards Spalart Allmaras" << "\n"; break;
          case TURB_MODEL::SST:       Breakdown_file << "Menter's TURB_MODEL::SST"     << "\n"; break;
          case TURB_MODEL::SST_SUST:  Breakdown_file << "Menter's TURB_MODEL::SST with sustaining terms" << "\n"; break;
          default: break;
        }
        break;
      default:
        break;
    }


    /*--- Compressible version of console output ---*/

    if (compressible) {


    if ((compressible) && (Kind_Solver != MAIN_SOLVER::FEM_ELASTICITY)) {
      Breakdown_file << "Mach number: " << config[val_iZone]->GetMach() <<"."<< "\n";
      Breakdown_file << "Angle of attack (AoA): " << config[val_iZone]->GetAoA() <<" deg, and angle of sideslip (AoS): " << config[val_iZone]->GetAoS() <<" deg."<< "\n";
      if ((Kind_Solver == MAIN_SOLVER::NAVIER_STOKES) || (Kind_Solver == MAIN_SOLVER::ADJ_NAVIER_STOKES) ||
          (Kind_Solver == MAIN_SOLVER::RANS) || (Kind_Solver == MAIN_SOLVER::ADJ_RANS))
        Breakdown_file << "Reynolds number: " << config[val_iZone]->GetReynolds() <<"."<< "\n";
    }

    if (fixed_cl) {
      Breakdown_file << "Simulation at a cte. CL: " << config[val_iZone]->GetTarget_CL() << ".\n";
      Breakdown_file << "Approx. Delta CL / Delta AoA: " << config[val_iZone]->GetdCL_dAlpha() << " (1/deg).\n";
      Breakdown_file << "Approx. Delta CD / Delta CL: " << config[val_iZone]->GetdCD_dCL() << ".\n";
      if (nDim == 3 ) {
        Breakdown_file << "Approx. Delta CMx / Delta CL: " << config[val_iZone]->GetdCMx_dCL() << ".\n";
        Breakdown_file << "Approx. Delta CMy / Delta CL: " << config[val_iZone]->GetdCMy_dCL() << ".\n";
      }
      Breakdown_file << "Approx. Delta CMz / Delta CL: " << config[val_iZone]->GetdCMz_dCL() << ".\n";
    }

    if (Ref_NonDim == DIMENSIONAL) { Breakdown_file << "Dimensional simulation." << "\n"; }
    else if (Ref_NonDim == FREESTREAM_PRESS_EQ_ONE) { Breakdown_file << "Non-Dimensional simulation (P=1.0, Rho=1.0, T=1.0 at the farfield)." << "\n"; }
    else if (Ref_NonDim == FREESTREAM_VEL_EQ_MACH) { Breakdown_file << "Non-Dimensional simulation (V=Mach, Rho=1.0, T=1.0 at the farfield)." << "\n"; }
    else if (Ref_NonDim == FREESTREAM_VEL_EQ_ONE) { Breakdown_file << "Non-Dimensional simulation (V=1.0, Rho=1.0, T=1.0 at the farfield)." << "\n"; }

    if (config[val_iZone]->GetSystemMeasurements() == SI) {
      Breakdown_file << "The reference area is " << config[val_iZone]->GetRefArea() << " m^2." << "\n";
      Breakdown_file << "The reference length is " << config[val_iZone]->GetRefLength() << " m." << "\n";
    }

    if (config[val_iZone]->GetSystemMeasurements() == US) {
      Breakdown_file << "The reference area is " << config[val_iZone]->GetRefArea()*12.0*12.0 << " in^2." << "\n";
      Breakdown_file << "The reference length is " << config[val_iZone]->GetRefLength()*12.0 << " in." << "\n";
    }
    Breakdown_file << "\n" << "\n" <<"Problem definition:" << "\n" << "\n";
    if (compressible) {
      if (viscous) {
        Breakdown_file << "Viscous flow: Computing pressure using the ideal gas law" << "\n";
        Breakdown_file << "based on the free-stream temperature and a density computed" << "\n";
        Breakdown_file << "from the Reynolds number." << "\n";
      } else {
        Breakdown_file << "Inviscid flow: Computing density based on free-stream" << "\n";
        Breakdown_file << "temperature and pressure using the ideal gas law." << "\n";
      }
    }

    if (dynamic_grid) Breakdown_file << "Force coefficients computed using MACH_MOTION." << "\n";
    else Breakdown_file << "Force coefficients computed using free-stream values." << "\n";

    if (incompressible) {
      Breakdown_file << "Viscous and Inviscid flow: rho_ref, and vel_ref" << "\n";
      Breakdown_file << "are based on the free-stream values, p_ref = rho_ref*vel_ref^2." << "\n";
      Breakdown_file << "The free-stream value of the pressure is 0." << "\n";
      Breakdown_file << "Mach number: "<< config[val_iZone]->GetMach() << ", computed using the Bulk modulus." << "\n";
      Breakdown_file << "Angle of attack (deg): "<< config[val_iZone]->GetAoA() << ", computed using the the free-stream velocity." << "\n";
      Breakdown_file << "Side slip angle (deg): "<< config[val_iZone]->GetAoS() << ", computed using the the free-stream velocity." << "\n";
      if (viscous) Breakdown_file << "Reynolds number: " << config[val_iZone]->GetReynolds() << ", computed using free-stream values."<< "\n";
      Breakdown_file << "Only dimensional computation, the grid should be dimensional." << "\n";
    }

    Breakdown_file <<"-- Input conditions:"<< "\n";

    if (compressible) {
      switch (config[val_iZone]->GetKind_FluidModel()) {

        case STANDARD_AIR:
          Breakdown_file << "Fluid Model: STANDARD_AIR "<< "\n";
          Breakdown_file << "Specific gas constant: " << config[val_iZone]->GetGas_Constant();
          if (config[val_iZone]->GetSystemMeasurements() == SI) Breakdown_file << " N.m/kg.K." << "\n";
          else if (config[val_iZone]->GetSystemMeasurements() == US) Breakdown_file << " lbf.ft/slug.R." << "\n";
          Breakdown_file << "Specific gas constant (non-dim): " << config[val_iZone]->GetGas_ConstantND()<< "\n";
          Breakdown_file << "Specific Heat Ratio: 1.4000 "<< "\n";
          break;

        case IDEAL_GAS:
          Breakdown_file << "Fluid Model: IDEAL_GAS "<< "\n";
          Breakdown_file << "Specific gas constant: " << config[val_iZone]->GetGas_Constant() << " N.m/kg.K." << "\n";
          Breakdown_file << "Specific gas constant (non-dim): " << config[val_iZone]->GetGas_ConstantND()<< "\n";
          Breakdown_file << "Specific Heat Ratio: "<< config[val_iZone]->GetGamma() << "\n";
          break;

        case VW_GAS:
          Breakdown_file << "Fluid Model: Van der Waals "<< "\n";
          Breakdown_file << "Specific gas constant: " << config[val_iZone]->GetGas_Constant() << " N.m/kg.K." << "\n";
          Breakdown_file << "Specific gas constant (non-dim): " << config[val_iZone]->GetGas_ConstantND()<< "\n";
          Breakdown_file << "Specific Heat Ratio: "<< config[val_iZone]->GetGamma() << "\n";
          Breakdown_file << "Critical Pressure:   " << config[val_iZone]->GetPressure_Critical()  << " Pa." << "\n";
          Breakdown_file << "Critical Temperature:  " << config[val_iZone]->GetTemperature_Critical() << " K." << "\n";
          Breakdown_file << "Critical Pressure (non-dim):   " << config[val_iZone]->GetPressure_Critical() /config[val_iZone]->GetPressure_Ref() << "\n";
          Breakdown_file << "Critical Temperature (non-dim) :  " << config[val_iZone]->GetTemperature_Critical() /config[val_iZone]->GetTemperature_Ref() << "\n";
          break;

        case PR_GAS:
          Breakdown_file << "Fluid Model: Peng-Robinson "<< "\n";
          Breakdown_file << "Specific gas constant: " << config[val_iZone]->GetGas_Constant() << " N.m/kg.K." << "\n";
          Breakdown_file << "Specific gas constant(non-dim): " << config[val_iZone]->GetGas_ConstantND()<< "\n";
          Breakdown_file << "Specific Heat Ratio: "<< config[val_iZone]->GetGamma() << "\n";
          Breakdown_file << "Critical Pressure:   " << config[val_iZone]->GetPressure_Critical()  << " Pa." << "\n";
          Breakdown_file << "Critical Temperature:  " << config[val_iZone]->GetTemperature_Critical() << " K." << "\n";
          Breakdown_file << "Critical Pressure (non-dim):   " << config[val_iZone]->GetPressure_Critical() /config[val_iZone]->GetPressure_Ref() << "\n";
          Breakdown_file << "Critical Temperature (non-dim) :  " << config[val_iZone]->GetTemperature_Critical() /config[val_iZone]->GetTemperature_Ref() << "\n";
          break;
      }

      if (viscous) {

        switch (config[val_iZone]->GetKind_ViscosityModel()) {

          case VISCOSITYMODEL::CONSTANT:
            Breakdown_file << "Viscosity Model: CONSTANT_VISCOSITY  "<< "\n";
            Breakdown_file << "Laminar Viscosity: " << config[val_iZone]->GetMu_Constant();
            if (config[val_iZone]->GetSystemMeasurements() == SI) Breakdown_file << " N.s/m^2." << "\n";
            else if (config[val_iZone]->GetSystemMeasurements() == US) Breakdown_file << " lbf.s/ft^2." << "\n";
            Breakdown_file << "Laminar Viscosity (non-dim): " << config[val_iZone]->GetMu_ConstantND()<< "\n";
            break;

          case VISCOSITYMODEL::SUTHERLAND:
            Breakdown_file << "Viscosity Model: SUTHERLAND "<< "\n";
            Breakdown_file << "Ref. Laminar Viscosity: " << config[val_iZone]->GetMu_Ref();
            if (config[val_iZone]->GetSystemMeasurements() == SI) Breakdown_file << " N.s/m^2." << "\n";
            else if (config[val_iZone]->GetSystemMeasurements() == US) Breakdown_file << " lbf.s/ft^2." << "\n";
            Breakdown_file << "Ref. Temperature: " << config[val_iZone]->GetMu_Temperature_Ref();
            if (config[val_iZone]->GetSystemMeasurements() == SI) Breakdown_file << " K." << "\n";
            else if (config[val_iZone]->GetSystemMeasurements() == US) Breakdown_file << " R." << "\n";
            Breakdown_file << "Sutherland Constant: "<< config[val_iZone]->GetMu_S();
            if (config[val_iZone]->GetSystemMeasurements() == SI) Breakdown_file << " K." << "\n";
            else if (config[val_iZone]->GetSystemMeasurements() == US) Breakdown_file << " R." << "\n";
            Breakdown_file << "Laminar Viscosity (non-dim): " << config[val_iZone]->GetMu_ConstantND()<< "\n";
            Breakdown_file << "Ref. Temperature (non-dim): " << config[val_iZone]->GetMu_Temperature_RefND()<< "\n";
            Breakdown_file << "Sutherland constant (non-dim): "<< config[val_iZone]->GetMu_SND()<< "\n";
            break;

          default:
            break;

        }
        switch (config[val_iZone]->GetKind_ConductivityModel()) {

          case CONDUCTIVITYMODEL::CONSTANT_PRANDTL:
            Breakdown_file << "Conductivity Model: CONSTANT_PRANDTL  "<< "\n";
            Breakdown_file << "Prandtl: " << config[val_iZone]->GetPrandtl_Lam()<< "\n";
            break;

          case CONDUCTIVITYMODEL::CONSTANT:
            Breakdown_file << "Conductivity Model: CONSTANT "<< "\n";
            Breakdown_file << "Molecular Conductivity: " << config[val_iZone]->GetThermal_Conductivity_Constant()<< " W/m^2.K." << "\n";
            Breakdown_file << "Molecular Conductivity (non-dim): " << config[val_iZone]->GetThermal_Conductivity_ConstantND()<< "\n";
            break;

          default:
            break;

        }

        if ((Kind_Solver == MAIN_SOLVER::RANS) || (Kind_Solver == MAIN_SOLVER::ADJ_RANS) || (Kind_Solver == MAIN_SOLVER::DISC_ADJ_RANS)) {
          switch (config[val_iZone]->GetKind_ConductivityModel_Turb()) {
            case CONDUCTIVITYMODEL_TURB::CONSTANT_PRANDTL:
              Breakdown_file << "Turbulent Conductivity Model: CONSTANT_PRANDTL  "<< "\n";
              Breakdown_file << "Turbulent Prandtl: " << config[val_iZone]->GetPrandtl_Turb()<< "\n";
              break;
            case CONDUCTIVITYMODEL_TURB::NONE:
              Breakdown_file << "Turbulent Conductivity Model: NONE "<< "\n";
              Breakdown_file << "No turbulent component in effective thermal conductivity." << "\n";
              break;
          }
        }

      }
    }

    if (incompressible) {
      Breakdown_file << "Bulk modulus: " << config[val_iZone]->GetBulk_Modulus();
      if (config[val_iZone]->GetSystemMeasurements() == SI) Breakdown_file << " Pa." << "\n";
      else if (config[val_iZone]->GetSystemMeasurements() == US) Breakdown_file << " psf." << "\n";
      Breakdown_file << "Epsilon^2 multiplier of Beta for incompressible preconditioner: " << config[val_iZone]->GetBeta_Factor();
      if (config[val_iZone]->GetSystemMeasurements() == SI) Breakdown_file << " Pa." << "\n";
      else if (config[val_iZone]->GetSystemMeasurements() == US) Breakdown_file << " psf." << "\n";
    }

    Breakdown_file << "Free-stream static pressure: " << config[val_iZone]->GetPressure_FreeStream();
    if (config[val_iZone]->GetSystemMeasurements() == SI) Breakdown_file << " Pa." << "\n";
    else if (config[val_iZone]->GetSystemMeasurements() == US) Breakdown_file << " psf." << "\n";

    Breakdown_file << "Free-stream total pressure: " << config[val_iZone]->GetPressure_FreeStream() * pow( 1.0+config[val_iZone]->GetMach()*config[val_iZone]->GetMach()*0.5*(config[val_iZone]->GetGamma()-1.0), config[val_iZone]->GetGamma()/(config[val_iZone]->GetGamma()-1.0) );
    if (config[val_iZone]->GetSystemMeasurements() == SI) Breakdown_file << " Pa." << "\n";
    else if (config[val_iZone]->GetSystemMeasurements() == US) Breakdown_file << " psf." << "\n";

    if (compressible) {
      Breakdown_file << "Free-stream temperature: " << config[val_iZone]->GetTemperature_FreeStream();
      if (config[val_iZone]->GetSystemMeasurements() == SI) Breakdown_file << " K." << "\n";
      else if (config[val_iZone]->GetSystemMeasurements() == US) Breakdown_file << " R." << "\n";

      Breakdown_file << "Free-stream total temperature: " << config[val_iZone]->GetTemperature_FreeStream() * (1.0 + config[val_iZone]->GetMach() * config[val_iZone]->GetMach() * 0.5 * (config[val_iZone]->GetGamma() - 1.0));
      if (config[val_iZone]->GetSystemMeasurements() == SI) Breakdown_file << " K." << "\n";
      else if (config[val_iZone]->GetSystemMeasurements() == US) Breakdown_file << " R." << "\n";
    }

    Breakdown_file << "Free-stream density: " << config[val_iZone]->GetDensity_FreeStream();
    if (config[val_iZone]->GetSystemMeasurements() == SI) Breakdown_file << " kg/m^3." << "\n";
    else if (config[val_iZone]->GetSystemMeasurements() == US) Breakdown_file << " slug/ft^3." << "\n";

    if (nDim == 2) {
      Breakdown_file << "Free-stream velocity: (" << config[val_iZone]->GetVelocity_FreeStream()[0] << ", ";
      Breakdown_file << config[val_iZone]->GetVelocity_FreeStream()[1] << ")";
    }
    if (nDim == 3) {
      Breakdown_file << "Free-stream velocity: (" << config[val_iZone]->GetVelocity_FreeStream()[0] << ", ";
      Breakdown_file << config[val_iZone]->GetVelocity_FreeStream()[1] << ", " << config[val_iZone]->GetVelocity_FreeStream()[2] << ")";
    }
    if (config[val_iZone]->GetSystemMeasurements() == SI) Breakdown_file << " m/s. ";
    else if (config[val_iZone]->GetSystemMeasurements() == US) Breakdown_file << " ft/s. ";

    Breakdown_file << "Magnitude: "  << config[val_iZone]->GetModVel_FreeStream();
    if (config[val_iZone]->GetSystemMeasurements() == SI) Breakdown_file << " m/s." << "\n";
    else if (config[val_iZone]->GetSystemMeasurements() == US) Breakdown_file << " ft/s." << "\n";

    if (compressible) {
      Breakdown_file << "Free-stream total energy per unit mass: " << config[val_iZone]->GetEnergy_FreeStream();
      if (config[val_iZone]->GetSystemMeasurements() == SI) Breakdown_file << " m^2/s^2." << "\n";
      else if (config[val_iZone]->GetSystemMeasurements() == US) Breakdown_file << " ft^2/s^2." << "\n";
    }

    if (viscous) {
      Breakdown_file << "Free-stream viscosity: " << config[val_iZone]->GetViscosity_FreeStream();
      if (config[val_iZone]->GetSystemMeasurements() == SI) Breakdown_file << " N.s/m^2." << "\n";
      else if (config[val_iZone]->GetSystemMeasurements() == US) Breakdown_file << " lbf.s/ft^2." << "\n";
      if (turbulent) {
        Breakdown_file << "Free-stream turb. kinetic energy per unit mass: " << config[val_iZone]->GetTke_FreeStream();
        if (config[val_iZone]->GetSystemMeasurements() == SI) Breakdown_file << " m^2/s^2." << "\n";
        else if (config[val_iZone]->GetSystemMeasurements() == US) Breakdown_file << " ft^2/s^2." << "\n";
        Breakdown_file << "Free-stream specific dissipation: " << config[val_iZone]->GetOmega_FreeStream();
        if (config[val_iZone]->GetSystemMeasurements() == SI) Breakdown_file << " 1/s." << "\n";
        else if (config[val_iZone]->GetSystemMeasurements() == US) Breakdown_file << " 1/s." << "\n";
      }
    }

    if (unsteady) { Breakdown_file << "Total time: " << config[val_iZone]->GetTotal_UnstTime() << " s. Time step: " << config[val_iZone]->GetDelta_UnstTime() << " s." << "\n"; }

    /*--- Print out reference values. ---*/

    Breakdown_file <<"-- Reference values:"<< "\n";

    if (compressible) {
      Breakdown_file << "Reference specific gas constant: " << config[val_iZone]->GetGas_Constant_Ref();
      if (config[val_iZone]->GetSystemMeasurements() == SI) Breakdown_file << " N.m/kg.K." << "\n";
      else if (config[val_iZone]->GetSystemMeasurements() == US) Breakdown_file << " lbf.ft/slug.R." << "\n";
    }

    Breakdown_file << "Reference pressure: " << config[val_iZone]->GetPressure_Ref();
    if (config[val_iZone]->GetSystemMeasurements() == SI) Breakdown_file << " Pa." << "\n";
    else if (config[val_iZone]->GetSystemMeasurements() == US) Breakdown_file << " psf." << "\n";

    if (compressible) {
      Breakdown_file << "Reference temperature: " << config[val_iZone]->GetTemperature_Ref();
      if (config[val_iZone]->GetSystemMeasurements() == SI) Breakdown_file << " K." << "\n";
      else if (config[val_iZone]->GetSystemMeasurements() == US) Breakdown_file << " R." << "\n";
    }

    Breakdown_file << "Reference density: " << config[val_iZone]->GetDensity_Ref();
    if (config[val_iZone]->GetSystemMeasurements() == SI) Breakdown_file << " kg/m^3." << "\n";
    else if (config[val_iZone]->GetSystemMeasurements() == US) Breakdown_file << " slug/ft^3." << "\n";

    Breakdown_file << "Reference velocity: " << config[val_iZone]->GetVelocity_Ref();
    if (config[val_iZone]->GetSystemMeasurements() == SI) Breakdown_file << " m/s." << "\n";
    else if (config[val_iZone]->GetSystemMeasurements() == US) Breakdown_file << " ft/s." << "\n";

    if (compressible) {
      Breakdown_file << "Reference energy per unit mass: " << config[val_iZone]->GetEnergy_Ref();
      if (config[val_iZone]->GetSystemMeasurements() == SI) Breakdown_file << " m^2/s^2." << "\n";
      else if (config[val_iZone]->GetSystemMeasurements() == US) Breakdown_file << " ft^2/s^2." << "\n";
    }

    if (incompressible) {
      Breakdown_file << "Reference length: " << config[val_iZone]->GetLength_Ref();
      if (config[val_iZone]->GetSystemMeasurements() == SI) Breakdown_file << " m." << "\n";
      else if (config[val_iZone]->GetSystemMeasurements() == US) Breakdown_file << " in." << "\n";
    }

    if (viscous) {
      Breakdown_file << "Reference viscosity: " << config[val_iZone]->GetViscosity_Ref();
      if (config[val_iZone]->GetSystemMeasurements() == SI) Breakdown_file << " N.s/m^2." << "\n";
      else if (config[val_iZone]->GetSystemMeasurements() == US) Breakdown_file << " lbf.s/ft^2." << "\n";
      if (compressible){
        Breakdown_file << "Reference conductivity: " << config[val_iZone]->GetThermal_Conductivity_Ref();
        if (config[val_iZone]->GetSystemMeasurements() == SI) Breakdown_file << " W/m^2.K." << "\n";
        else if (config[val_iZone]->GetSystemMeasurements() == US) Breakdown_file << " lbf/ft.s.R." << "\n";
      }
    }


    if (unsteady) Breakdown_file << "Reference time: " << config[val_iZone]->GetTime_Ref() <<" s." << "\n";

    /*--- Print out resulting non-dim values here. ---*/

    Breakdown_file << "-- Resulting non-dimensional state:" << "\n";
    Breakdown_file << "Mach number (non-dim): " << config[val_iZone]->GetMach() << "\n";
    if (viscous) {
      Breakdown_file << "Reynolds number (non-dim): " << config[val_iZone]->GetReynolds() <<". Re length: " << config[val_iZone]->GetLength_Reynolds();
      if (config[val_iZone]->GetSystemMeasurements() == SI) Breakdown_file << " m." << "\n";
      else if (config[val_iZone]->GetSystemMeasurements() == US) Breakdown_file << " ft." << "\n";
    }
    if (gravity) {
      Breakdown_file << "Froude number (non-dim): " << config[val_iZone]->GetFroude() << "\n";
      Breakdown_file << "Lenght of the baseline wave (non-dim): " << 2.0*PI_NUMBER*config[val_iZone]->GetFroude()*config[val_iZone]->GetFroude() << "\n";
    }

    if (compressible) {
      Breakdown_file << "Specific gas constant (non-dim): " << config[val_iZone]->GetGas_ConstantND() << "\n";
      Breakdown_file << "Free-stream temperature (non-dim): " << config[val_iZone]->GetTemperature_FreeStreamND() << "\n";
    }

    Breakdown_file << "Free-stream pressure (non-dim): " << config[val_iZone]->GetPressure_FreeStreamND() << "\n";

    Breakdown_file << "Free-stream density (non-dim): " << config[val_iZone]->GetDensity_FreeStreamND() << "\n";

    if (nDim == 2) {
      Breakdown_file << "Free-stream velocity (non-dim): (" << config[val_iZone]->GetVelocity_FreeStreamND()[0] << ", ";
      Breakdown_file << config[val_iZone]->GetVelocity_FreeStreamND()[1] << "). ";
    } else {
      Breakdown_file << "Free-stream velocity (non-dim): (" << config[val_iZone]->GetVelocity_FreeStreamND()[0] << ", ";
      Breakdown_file << config[val_iZone]->GetVelocity_FreeStreamND()[1] << ", " << config[val_iZone]->GetVelocity_FreeStreamND()[2] << "). ";
    }
    Breakdown_file << "Magnitude: "   << config[val_iZone]->GetModVel_FreeStreamND() << "\n";

    if (compressible)
      Breakdown_file << "Free-stream total energy per unit mass (non-dim): " << config[val_iZone]->GetEnergy_FreeStreamND() << "\n";

    if (viscous) {
      Breakdown_file << "Free-stream viscosity (non-dim): " << config[val_iZone]->GetViscosity_FreeStreamND() << "\n";
      if (turbulent) {
        Breakdown_file << "Free-stream turb. kinetic energy (non-dim): " << config[val_iZone]->GetTke_FreeStreamND() << "\n";
        Breakdown_file << "Free-stream specific dissipation (non-dim): " << config[val_iZone]->GetOmega_FreeStreamND() << "\n";
      }
    }

    if (unsteady) {
      Breakdown_file << "Total time (non-dim): " << config[val_iZone]->GetTotal_UnstTimeND() << "\n";
      Breakdown_file << "Time step (non-dim): " << config[val_iZone]->GetDelta_UnstTimeND() << "\n";
    }

    } else {

    /*--- Incompressible version of the console output ---*/

      bool energy     = config[val_iZone]->GetEnergy_Equation();
      bool boussinesq = (config[val_iZone]->GetKind_DensityModel() == INC_DENSITYMODEL::BOUSSINESQ);

      if (config[val_iZone]->GetRef_Inc_NonDim() == DIMENSIONAL) {
        Breakdown_file << "Viscous and Inviscid flow: rho_ref, vel_ref, temp_ref, p_ref" << "\n";
        Breakdown_file << "are set to 1.0 in order to perform a dimensional calculation." << "\n";
        if (dynamic_grid) Breakdown_file << "Force coefficients computed using MACH_MOTION." << "\n";
        else Breakdown_file << "Force coefficients computed using initial values." << "\n";
      }
      else if (config[val_iZone]->GetRef_Inc_NonDim() == INITIAL_VALUES) {
        Breakdown_file << "Viscous and Inviscid flow: rho_ref, vel_ref, and temp_ref" << "\n";
        Breakdown_file << "are based on the initial values, p_ref = rho_ref*vel_ref^2." << "\n";
        if (dynamic_grid) Breakdown_file << "Force coefficients computed using MACH_MOTION." << "\n";
        else Breakdown_file << "Force coefficients computed using initial values." << "\n";
      }
      else if (config[val_iZone]->GetRef_Inc_NonDim() == REFERENCE_VALUES) {
        Breakdown_file << "Viscous and Inviscid flow: rho_ref, vel_ref, and temp_ref" << "\n";
        Breakdown_file << "are user-provided reference values, p_ref = rho_ref*vel_ref^2." << "\n";
        if (dynamic_grid) Breakdown_file << "Force coefficients computed using MACH_MOTION." << "\n";
        else Breakdown_file << "Force coefficients computed using reference values." << "\n";
      }
      Breakdown_file << "The reference area for force coeffs. is " << config[val_iZone]->GetRefArea() << " m^2." << "\n";
      Breakdown_file << "The reference length for force coeffs. is " << config[val_iZone]->GetRefLength() << " m." << "\n";

      Breakdown_file << "The pressure is decomposed into thermodynamic and dynamic components." << "\n";
      Breakdown_file << "The initial value of the dynamic pressure is 0." << "\n";

      Breakdown_file << "Mach number: "<< config[val_iZone]->GetMach();
      if (config[val_iZone]->GetKind_FluidModel() == CONSTANT_DENSITY) {
        Breakdown_file << ", computed using the Bulk modulus." << "\n";
      } else {
        Breakdown_file << ", computed using fluid speed of sound." << "\n";
      }

      Breakdown_file << "For external flows, the initial state is imposed at the far-field." << "\n";
      Breakdown_file << "Angle of attack (deg): "<< config[val_iZone]->GetAoA() << ", computed using the initial velocity." << "\n";
      Breakdown_file << "Side slip angle (deg): "<< config[val_iZone]->GetAoS() << ", computed using the initial velocity." << "\n";

      if (viscous) {
        Breakdown_file << "Reynolds number per meter: " << config[val_iZone]->GetReynolds() << ", computed using initial values."<< "\n";
        Breakdown_file << "Reynolds number is a byproduct of inputs only (not used internally)." << "\n";
      }
      Breakdown_file << "SI units only. The grid should be dimensional (meters)." << "\n";

      switch (config[val_iZone]->GetKind_DensityModel()) {

        case INC_DENSITYMODEL::CONSTANT:
          if (energy) Breakdown_file << "Energy equation is active and decoupled." << "\n";
          else Breakdown_file << "No energy equation." << "\n";
          break;

        case INC_DENSITYMODEL::BOUSSINESQ:
          if (energy) Breakdown_file << "Energy equation is active and coupled through Boussinesq approx." << "\n";
          break;

        case INC_DENSITYMODEL::VARIABLE:
          if (energy) Breakdown_file << "Energy equation is active and coupled for variable density." << "\n";
          break;

      }

      Breakdown_file <<"-- Input conditions:"<< "\n";

      switch (config[val_iZone]->GetKind_FluidModel()) {

        case CONSTANT_DENSITY:
          Breakdown_file << "Fluid Model: CONSTANT_DENSITY "<< "\n";
          if (energy) {
            Breakdown_file << "Specific heat at constant pressure (Cp): " << config[val_iZone]->GetSpecific_Heat_Cp() << " N.m/kg.K." << "\n";
          }
          if (boussinesq) Breakdown_file << "Thermal expansion coefficient: " << config[val_iZone]->GetThermal_Expansion_Coeff() << " K^-1." << "\n";
          Breakdown_file << "Thermodynamic pressure not required." << "\n";
          break;

        case INC_IDEAL_GAS:
          Breakdown_file << "Fluid Model: INC_IDEAL_GAS "<< endl;
          Breakdown_file << "Variable density incompressible flow using ideal gas law." << endl;
          Breakdown_file << "Density is a function of temperature (constant thermodynamic pressure)." << endl;
          Breakdown_file << "Specific heat at constant pressure (Cp): " << config[val_iZone]->GetSpecific_Heat_Cp() << " N.m/kg.K." << endl;
          Breakdown_file << "Molecular weight : "<< config[val_iZone]->GetMolecular_Weight() << " g/mol" << endl;
          Breakdown_file << "Specific gas constant: " << config[val_iZone]->GetGas_Constant() << " N.m/kg.K." << endl;
          Breakdown_file << "Thermodynamic pressure: " << config[val_iZone]->GetPressure_Thermodynamic();
          if (config[val_iZone]->GetSystemMeasurements() == SI) Breakdown_file << " Pa." << endl;
          else if (config[val_iZone]->GetSystemMeasurements() == US) Breakdown_file << " psf." << endl;
          break;

        case INC_IDEAL_GAS_POLY:
          Breakdown_file << "Fluid Model: INC_IDEAL_GAS_POLY "<< endl;
          Breakdown_file << "Variable density incompressible flow using ideal gas law." << endl;
          Breakdown_file << "Density is a function of temperature (constant thermodynamic pressure)." << endl;
          Breakdown_file << "Molecular weight: " << config[val_iZone]->GetMolecular_Weight() << " g/mol." << endl;
          Breakdown_file << "Specific gas constant: " << config[val_iZone]->GetGas_Constant() << " N.m/kg.K." << endl;
          Breakdown_file << "Specific gas constant (non-dim): " << config[val_iZone]->GetGas_ConstantND() << endl;
          Breakdown_file << "Thermodynamic pressure: " << config[val_iZone]->GetPressure_Thermodynamic();
          if (config[val_iZone]->GetSystemMeasurements() == SI) Breakdown_file << " Pa." << endl;
          else if (config[val_iZone]->GetSystemMeasurements() == US) Breakdown_file << " psf." << endl;
          Breakdown_file << "Cp(T) polynomial coefficients: \n  (";
          for (unsigned short iVar = 0; iVar < config[val_iZone]->GetnPolyCoeffs(); iVar++) {
            Breakdown_file << config[val_iZone]->GetCp_PolyCoeff(iVar);
            if (iVar < config[val_iZone]->GetnPolyCoeffs()-1) Breakdown_file << ", ";
          }
          Breakdown_file << ")." << endl;
          Breakdown_file << "Cp(T) polynomial coefficients (non-dim.): \n  (";
          for (unsigned short iVar = 0; iVar < config[val_iZone]->GetnPolyCoeffs(); iVar++) {
            Breakdown_file << config[val_iZone]->GetCp_PolyCoeffND(iVar);
            if (iVar < config[val_iZone]->GetnPolyCoeffs()-1) Breakdown_file << ", ";
          }
          Breakdown_file << ")." << endl;
          break;

      }
      if (viscous) {
        switch (config[val_iZone]->GetKind_ViscosityModel()) {

          case VISCOSITYMODEL::CONSTANT:
            Breakdown_file << "Viscosity Model: CONSTANT_VISCOSITY  "<< "\n";
            Breakdown_file << "Constant Laminar Viscosity: " << config[val_iZone]->GetMu_Constant();
            if (config[val_iZone]->GetSystemMeasurements() == SI) Breakdown_file << " N.s/m^2." << "\n";
            else if (config[val_iZone]->GetSystemMeasurements() == US) Breakdown_file << " lbf.s/ft^2." << "\n";
            Breakdown_file << "Laminar Viscosity (non-dim): " << config[val_iZone]->GetMu_ConstantND()<< "\n";
            break;

          case VISCOSITYMODEL::SUTHERLAND:
            Breakdown_file << "Viscosity Model: SUTHERLAND "<< "\n";
            Breakdown_file << "Ref. Laminar Viscosity: " << config[val_iZone]->GetMu_Ref();
            if (config[val_iZone]->GetSystemMeasurements() == SI) Breakdown_file << " N.s/m^2." << "\n";
            else if (config[val_iZone]->GetSystemMeasurements() == US) Breakdown_file << " lbf.s/ft^2." << "\n";
            Breakdown_file << "Ref. Temperature: " << config[val_iZone]->GetMu_Temperature_Ref();
            if (config[val_iZone]->GetSystemMeasurements() == SI) Breakdown_file << " K." << "\n";
            else if (config[val_iZone]->GetSystemMeasurements() == US) Breakdown_file << " R." << "\n";
            Breakdown_file << "Sutherland Constant: "<< config[val_iZone]->GetMu_S();
            if (config[val_iZone]->GetSystemMeasurements() == SI) Breakdown_file << " K." << "\n";
            else if (config[val_iZone]->GetSystemMeasurements() == US) Breakdown_file << " R." << "\n";
            Breakdown_file << "Laminar Viscosity (non-dim): " << config[val_iZone]->GetMu_ConstantND()<< "\n";
            Breakdown_file << "Ref. Temperature (non-dim): " << config[val_iZone]->GetMu_Temperature_RefND()<< "\n";
            Breakdown_file << "Sutherland constant (non-dim): "<< config[val_iZone]->GetMu_SND()<< "\n";
            break;

          case VISCOSITYMODEL::POLYNOMIAL:
            Breakdown_file << "Viscosity Model: POLYNOMIAL_VISCOSITY  "<< endl;
            Breakdown_file << "Mu(T) polynomial coefficients: \n  (";
            for (unsigned short iVar = 0; iVar < config[val_iZone]->GetnPolyCoeffs(); iVar++) {
              Breakdown_file << config[val_iZone]->GetMu_PolyCoeff(iVar);
              if (iVar < config[val_iZone]->GetnPolyCoeffs()-1) Breakdown_file << ", ";
            }
            Breakdown_file << ")." << endl;
            Breakdown_file << "Mu(T) polynomial coefficients (non-dim.): \n  (";
            for (unsigned short iVar = 0; iVar < config[val_iZone]->GetnPolyCoeffs(); iVar++) {
              Breakdown_file << config[val_iZone]->GetMu_PolyCoeffND(iVar);
              if (iVar < config[val_iZone]->GetnPolyCoeffs()-1) Breakdown_file << ", ";
            }
            Breakdown_file << ")." << endl;
            break;

        }

        if (energy) {
          switch (config[val_iZone]->GetKind_ConductivityModel()) {

            case CONDUCTIVITYMODEL::CONSTANT_PRANDTL:
              Breakdown_file << "Conductivity Model: CONSTANT_PRANDTL  "<< "\n";
              Breakdown_file << "Prandtl (Laminar): " << config[val_iZone]->GetPrandtl_Lam()<< "\n";
              break;

            case CONDUCTIVITYMODEL::CONSTANT:
              Breakdown_file << "Conductivity Model: CONSTANT "<< "\n";
              Breakdown_file << "Molecular Conductivity: " << config[val_iZone]->GetThermal_Conductivity_Constant()<< " W/m^2.K." << "\n";
              Breakdown_file << "Molecular Conductivity (non-dim): " << config[val_iZone]->GetThermal_Conductivity_ConstantND()<< "\n";
              break;

            case CONDUCTIVITYMODEL::POLYNOMIAL:
              Breakdown_file << "Viscosity Model: POLYNOMIAL "<< endl;
              Breakdown_file << "Kt(T) polynomial coefficients: \n  (";
              for (unsigned short iVar = 0; iVar < config[val_iZone]->GetnPolyCoeffs(); iVar++) {
                Breakdown_file << config[val_iZone]->GetKt_PolyCoeff(iVar);
                if (iVar < config[val_iZone]->GetnPolyCoeffs()-1) Breakdown_file << ", ";
              }
              Breakdown_file << ")." << endl;
              Breakdown_file << "Kt(T) polynomial coefficients (non-dim.): \n  (";
              for (unsigned short iVar = 0; iVar < config[val_iZone]->GetnPolyCoeffs(); iVar++) {
                Breakdown_file << config[val_iZone]->GetKt_PolyCoeffND(iVar);
                if (iVar < config[val_iZone]->GetnPolyCoeffs()-1) Breakdown_file << ", ";
              }
              Breakdown_file << ")." << endl;
              break;

          }

          if ((Kind_Solver == MAIN_SOLVER::RANS) || (Kind_Solver == MAIN_SOLVER::ADJ_RANS) || (Kind_Solver == MAIN_SOLVER::DISC_ADJ_RANS)) {
            switch (config[val_iZone]->GetKind_ConductivityModel_Turb()) {
              case CONDUCTIVITYMODEL_TURB::CONSTANT_PRANDTL:
                Breakdown_file << "Turbulent Conductivity Model: CONSTANT_PRANDTL  "<< "\n";
                Breakdown_file << "Turbulent Prandtl: " << config[val_iZone]->GetPrandtl_Turb()<< "\n";
                break;
              case CONDUCTIVITYMODEL_TURB::NONE:
                Breakdown_file << "Turbulent Conductivity Model: NONE "<< "\n";
                Breakdown_file << "No turbulent component in effective thermal conductivity." << "\n";
                break;
            }
          }

        }

      }

      if (config[val_iZone]->GetKind_FluidModel() == CONSTANT_DENSITY) {
        Breakdown_file << "Bulk modulus: " << config[val_iZone]->GetBulk_Modulus();
        if (config[val_iZone]->GetSystemMeasurements() == SI) Breakdown_file << " Pa." << "\n";
        else if (config[val_iZone]->GetSystemMeasurements() == US) Breakdown_file << " psf." << "\n";
      }

      Breakdown_file << "Initial dynamic pressure: " << config[val_iZone]->GetPressure_FreeStream();
      if (config[val_iZone]->GetSystemMeasurements() == SI) Breakdown_file << " Pa." << "\n";
      else if (config[val_iZone]->GetSystemMeasurements() == US) Breakdown_file << " psf." << "\n";

      Breakdown_file << "Initial total pressure: " << config[val_iZone]->GetPressure_FreeStream() + 0.5*config[val_iZone]->GetDensity_FreeStream()*config[val_iZone]->GetModVel_FreeStream()*config[val_iZone]->GetModVel_FreeStream();
      if (config[val_iZone]->GetSystemMeasurements() == SI) Breakdown_file << " Pa." << "\n";
      else if (config[val_iZone]->GetSystemMeasurements() == US) Breakdown_file << " psf." << "\n";

      if (energy) {
        Breakdown_file << "Initial temperature: " << config[val_iZone]->GetTemperature_FreeStream();
        if (config[val_iZone]->GetSystemMeasurements() == SI) Breakdown_file << " K." << "\n";
        else if (config[val_iZone]->GetSystemMeasurements() == US) Breakdown_file << " R." << "\n";
      }

      Breakdown_file << "Initial density: " << config[val_iZone]->GetDensity_FreeStream();
      if (config[val_iZone]->GetSystemMeasurements() == SI) Breakdown_file << " kg/m^3." << "\n";
      else if (config[val_iZone]->GetSystemMeasurements() == US) Breakdown_file << " slug/ft^3." << "\n";

      if (nDim == 2) {
        Breakdown_file << "Initial velocity: (" << config[val_iZone]->GetVelocity_FreeStream()[0] << ", ";
        Breakdown_file << config[val_iZone]->GetVelocity_FreeStream()[1] << ")";
      }
      if (nDim == 3) {
        Breakdown_file << "Initial velocity: (" << config[val_iZone]->GetVelocity_FreeStream()[0] << ", ";
        Breakdown_file << config[val_iZone]->GetVelocity_FreeStream()[1] << ", " << config[val_iZone]->GetVelocity_FreeStream()[2] << ")";
      }
      if (config[val_iZone]->GetSystemMeasurements() == SI) Breakdown_file << " m/s. ";
      else if (config[val_iZone]->GetSystemMeasurements() == US) Breakdown_file << " ft/s. ";

      Breakdown_file << "Magnitude: "  << config[val_iZone]->GetModVel_FreeStream();
      if (config[val_iZone]->GetSystemMeasurements() == SI) Breakdown_file << " m/s." << "\n";
      else if (config[val_iZone]->GetSystemMeasurements() == US) Breakdown_file << " ft/s." << "\n";

      if (viscous) {
        Breakdown_file << "Initial laminar viscosity: " << config[val_iZone]->GetViscosity_FreeStream();
        if (config[val_iZone]->GetSystemMeasurements() == SI) Breakdown_file << " N.s/m^2." << "\n";
        else if (config[val_iZone]->GetSystemMeasurements() == US) Breakdown_file << " lbf.s/ft^2." << "\n";
        if (turbulent) {
          Breakdown_file << "Initial turb. kinetic energy per unit mass: " << config[val_iZone]->GetTke_FreeStream();
          if (config[val_iZone]->GetSystemMeasurements() == SI) Breakdown_file << " m^2/s^2." << "\n";
          else if (config[val_iZone]->GetSystemMeasurements() == US) Breakdown_file << " ft^2/s^2." << "\n";
          Breakdown_file << "Initial specific dissipation: " << config[val_iZone]->GetOmega_FreeStream();
          if (config[val_iZone]->GetSystemMeasurements() == SI) Breakdown_file << " 1/s." << "\n";
          else if (config[val_iZone]->GetSystemMeasurements() == US) Breakdown_file << " 1/s." << "\n";
        }
      }

      if (unsteady) { Breakdown_file << "Total time: " << config[val_iZone]->GetTotal_UnstTime() << " s. Time step: " << config[val_iZone]->GetDelta_UnstTime() << " s." << "\n"; }

      /*--- Print out reference values. ---*/

      Breakdown_file <<"-- Reference values:"<< "\n";

      if (config[val_iZone]->GetKind_FluidModel() != CONSTANT_DENSITY) {
        Breakdown_file << "Reference specific gas constant: " << config[val_iZone]->GetGas_Constant_Ref();
        if (config[val_iZone]->GetSystemMeasurements() == SI) Breakdown_file << " N.m/kg.K." << "\n";
        else if (config[val_iZone]->GetSystemMeasurements() == US) Breakdown_file << " lbf.ft/slug.R." << "\n";
      } else {
        if (energy) {
          Breakdown_file << "Reference specific heat: " << config[val_iZone]->GetGas_Constant_Ref();
          if (config[val_iZone]->GetSystemMeasurements() == SI) Breakdown_file << " N.m/kg.K." << "\n";
          else if (config[val_iZone]->GetSystemMeasurements() == US) Breakdown_file << " lbf.ft/slug.R." << "\n";
        }
      }

      Breakdown_file << "Reference pressure: " << config[val_iZone]->GetPressure_Ref();
      if (config[val_iZone]->GetSystemMeasurements() == SI) Breakdown_file << " Pa." << "\n";
      else if (config[val_iZone]->GetSystemMeasurements() == US) Breakdown_file << " psf." << "\n";

      if (energy) {
        Breakdown_file << "Reference temperature: " << config[val_iZone]->GetTemperature_Ref();
        if (config[val_iZone]->GetSystemMeasurements() == SI) Breakdown_file << " K." << "\n";
        else if (config[val_iZone]->GetSystemMeasurements() == US) Breakdown_file << " R." << "\n";
      }

      Breakdown_file << "Reference density: " << config[val_iZone]->GetDensity_Ref();
      if (config[val_iZone]->GetSystemMeasurements() == SI) Breakdown_file << " kg/m^3." << "\n";
      else if (config[val_iZone]->GetSystemMeasurements() == US) Breakdown_file << " slug/ft^3." << "\n";

      Breakdown_file << "Reference velocity: " << config[val_iZone]->GetVelocity_Ref();
      if (config[val_iZone]->GetSystemMeasurements() == SI) Breakdown_file << " m/s." << "\n";
      else if (config[val_iZone]->GetSystemMeasurements() == US) Breakdown_file << " ft/s." << "\n";

      Breakdown_file << "Reference length: " << config[val_iZone]->GetLength_Ref();
      if (config[val_iZone]->GetSystemMeasurements() == SI) Breakdown_file << " m." << "\n";
      else if (config[val_iZone]->GetSystemMeasurements() == US) Breakdown_file << " in." << "\n";

      if (viscous) {
        Breakdown_file << "Reference viscosity: " << config[val_iZone]->GetViscosity_Ref();
        if (config[val_iZone]->GetSystemMeasurements() == SI) Breakdown_file << " N.s/m^2." << "\n";
        else if (config[val_iZone]->GetSystemMeasurements() == US) Breakdown_file << " lbf.s/ft^2." << "\n";
      }

      if (unsteady) Breakdown_file << "Reference time: " << config[val_iZone]->GetTime_Ref() <<" s." << "\n";

      /*--- Print out resulting non-dim values here. ---*/

      Breakdown_file << "-- Resulting non-dimensional state:" << "\n";
      Breakdown_file << "Mach number (non-dim): " << config[val_iZone]->GetMach() << "\n";
      if (viscous) {
        Breakdown_file << "Reynolds number (per m): " << config[val_iZone]->GetReynolds() << "\n";
      }

      if (config[val_iZone]->GetKind_FluidModel() != CONSTANT_DENSITY) {
        Breakdown_file << "Specific gas constant (non-dim): " << config[val_iZone]->GetGas_ConstantND() << "\n";
        Breakdown_file << "Initial thermodynamic pressure (non-dim): " << config[val_iZone]->GetPressure_ThermodynamicND() << "\n";
      } else {
        if (energy) {
          Breakdown_file << "Specific heat at constant pressure (non-dim): " << config[val_iZone]->GetSpecific_Heat_CpND() << "\n";
          if (boussinesq) Breakdown_file << "Thermal expansion coefficient (non-dim.): " << config[val_iZone]->GetThermal_Expansion_CoeffND() << " K^-1." << "\n";
        }
      }

      if (energy) Breakdown_file << "Initial temperature (non-dim): " << config[val_iZone]->GetTemperature_FreeStreamND() << "\n";
      Breakdown_file << "Initial pressure (non-dim): " << config[val_iZone]->GetPressure_FreeStreamND() << "\n";
      Breakdown_file << "Initial density (non-dim): " << config[val_iZone]->GetDensity_FreeStreamND() << "\n";

      if (nDim == 2) {
        Breakdown_file << "Initial velocity (non-dim): (" << config[val_iZone]->GetVelocity_FreeStreamND()[0] << ", ";
        Breakdown_file << config[val_iZone]->GetVelocity_FreeStreamND()[1] << "). ";
      } else {
        Breakdown_file << "Initial velocity (non-dim): (" << config[val_iZone]->GetVelocity_FreeStreamND()[0] << ", ";
        Breakdown_file << config[val_iZone]->GetVelocity_FreeStreamND()[1] << ", " << config[val_iZone]->GetVelocity_FreeStreamND()[2] << "). ";
      }
      Breakdown_file << "Magnitude: "   << config[val_iZone]->GetModVel_FreeStreamND() << "\n";

      if (viscous) {
        Breakdown_file << "Initial viscosity (non-dim): " << config[val_iZone]->GetViscosity_FreeStreamND() << "\n";
        if (turbulent) {
          Breakdown_file << "Initial turb. kinetic energy (non-dim): " << config[val_iZone]->GetTke_FreeStreamND() << "\n";
          Breakdown_file << "Initial specific dissipation (non-dim): " << config[val_iZone]->GetOmega_FreeStreamND() << "\n";
        }
      }

      if (unsteady) {
        Breakdown_file << "Total time (non-dim): " << config[val_iZone]->GetTotal_UnstTimeND() << "\n";
        Breakdown_file << "Time step (non-dim): " << config[val_iZone]->GetDelta_UnstTimeND() << "\n";
      }

    }

    /*--- Begin forces breakdown info. ---*/

    Breakdown_file << fixed;
    Breakdown_file << "\n" << "\n" <<"Forces breakdown:" << "\n" << "\n";

    if (nDim == 3) {
      su2double m = solver[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetTotal_CFz()/solver[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetTotal_CFx();
      su2double term = (Total_CoPz/m)-Total_CoPx;

      if (term > 0) Breakdown_file << "Center of Pressure: X="  << 1/m <<"Z-"<< term << "." << "\n\n";
      else Breakdown_file << "Center of Pressure: X="  << 1/m <<"Z+"<< fabs(term);
      if (config[val_iZone]->GetSystemMeasurements() == SI) Breakdown_file << " m." << "\n\n";
      else Breakdown_file << " in." << "\n\n";
    }
    else {
      su2double m = solver[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetTotal_CFy()/solver[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetTotal_CFx();
      su2double term = (Total_CoPy/m)-Total_CoPx;
      if (term > 0) Breakdown_file << "Center of Pressure: X="  << 1/m <<"Y-"<< term << "." << "\n\n";
      else Breakdown_file << "Center of Pressure: X="  << 1/m <<"Y+"<< fabs(term);
      if (config[val_iZone]->GetSystemMeasurements() == SI) Breakdown_file << " m." << "\n\n";
      else Breakdown_file << " in." << "\n\n";
    }

    /*--- Reference area and force factors. ---*/

    su2double RefDensity, RefArea, RefVel, Factor, Ref;
    RefArea     = config[val_iZone]->GetRefArea();
    if (compressible) {
      RefDensity  = solver[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetDensity_Inf();
      RefVel = solver[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetModVelocity_Inf();
    } else {
      if ((config[val_iZone]->GetRef_Inc_NonDim() == DIMENSIONAL) ||
          (config[val_iZone]->GetRef_Inc_NonDim() == INITIAL_VALUES)) {
        RefDensity  = solver[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetDensity_Inf();
        RefVel = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          RefVel  += solver[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetVelocity_Inf(iDim)*solver[val_iZone][INST_0][FinestMesh][FLOW_SOL]->GetVelocity_Inf(iDim);
        RefVel = sqrt(RefVel);
      } else {
        RefDensity = config[val_iZone]->GetInc_Density_Ref();
        RefVel    = config[val_iZone]->GetInc_Velocity_Ref();
      }
    }
    Factor = (0.5*RefDensity*RefArea*RefVel*RefVel);
    Ref = config[val_iZone]->GetDensity_Ref() * config[val_iZone]->GetVelocity_Ref() * config[val_iZone]->GetVelocity_Ref() * 1.0 * 1.0;

    Breakdown_file << "NOTE: Multiply forces by the non-dimensional factor: " << Factor << ", and the reference factor: " << Ref  << "\n";
    Breakdown_file << "to obtain the dimensional force."  << "\n" << "\n";

    Breakdown_file << "Total CL:    ";
    Breakdown_file.width(11);
    Breakdown_file << Total_CL;
    Breakdown_file << " | Pressure (";
    Breakdown_file.width(5);
    Breakdown_file << SU2_TYPE::Int((Inv_CL * 100.0) / (Total_CL + EPS));
    Breakdown_file << "%): ";
    Breakdown_file.width(11);
    Breakdown_file << Inv_CL;
    Breakdown_file << " | Friction (";
    Breakdown_file.width(5);
    Breakdown_file << SU2_TYPE::Int((Visc_CL * 100.0) / (Total_CL + EPS));
    Breakdown_file << "%): ";
    Breakdown_file.width(11);
    Breakdown_file << Visc_CL;
    Breakdown_file << " | Momentum (";
    Breakdown_file.width(5);
    Breakdown_file << SU2_TYPE::Int((Mnt_CL * 100.0) / (Total_CL + EPS));
    Breakdown_file << "%): ";
    Breakdown_file.width(11);
    Breakdown_file << Mnt_CL << "\n";

    Breakdown_file << "Total CD:    ";
    Breakdown_file.width(11);
    Breakdown_file << Total_CD;
    Breakdown_file << " | Pressure (";
    Breakdown_file.width(5);
    Breakdown_file << SU2_TYPE::Int((Inv_CD * 100.0) / (Total_CD + EPS)) << "%): ";
    Breakdown_file.width(11);
    Breakdown_file << Inv_CD;
    Breakdown_file << " | Friction (";
    Breakdown_file.width(5);
    Breakdown_file << SU2_TYPE::Int((Visc_CD * 100.0) / (Total_CD + EPS)) << "%): ";
    Breakdown_file.width(11);
    Breakdown_file << Visc_CD;
    Breakdown_file << " | Momentum (";
    Breakdown_file.width(5);
    Breakdown_file << SU2_TYPE::Int((Mnt_CD * 100.0) / (Total_CD + EPS)) << "%): ";
    Breakdown_file.width(11);
    Breakdown_file << Mnt_CD << "\n";

    if (nDim == 3) {
      Breakdown_file << "Total CSF:   ";
      Breakdown_file.width(11);
      Breakdown_file << Total_CSF;
      Breakdown_file << " | Pressure (";
      Breakdown_file.width(5);
      Breakdown_file << SU2_TYPE::Int((Inv_CSF * 100.0) / (Total_CSF + EPS));
      Breakdown_file << "%): ";
      Breakdown_file.width(11);
      Breakdown_file << Inv_CSF;
      Breakdown_file << " | Friction (";
      Breakdown_file.width(5);
      Breakdown_file <<  SU2_TYPE::Int((Visc_CSF * 100.0) / (Total_CSF + EPS));
      Breakdown_file << "%): ";
      Breakdown_file.width(11);
      Breakdown_file << Visc_CSF;
      Breakdown_file << " | Momentum (";
      Breakdown_file.width(5);
      Breakdown_file << SU2_TYPE::Int((Mnt_CSF * 100.0) / (Total_CSF + EPS));
      Breakdown_file << "%): ";
      Breakdown_file.width(11);
      Breakdown_file << Mnt_CSF << "\n";
    }

    Breakdown_file << "Total CL/CD: ";
    Breakdown_file.width(11);
    Breakdown_file << Total_CEff;
    Breakdown_file << " | Pressure (";
    Breakdown_file.width(5);
    Breakdown_file << SU2_TYPE::Int((Inv_CEff * 100.0) / (Total_CEff + EPS));
    Breakdown_file << "%): ";
    Breakdown_file.width(11);
    Breakdown_file << Inv_CEff;
    Breakdown_file << " | Friction (";
    Breakdown_file.width(5);
    Breakdown_file <<  SU2_TYPE::Int((Visc_CEff * 100.0) / (Total_CEff + EPS));
    Breakdown_file << "%): ";
    Breakdown_file.width(11);
    Breakdown_file << Visc_CEff;
    Breakdown_file << " | Momentum (";
    Breakdown_file.width(5);
    Breakdown_file << SU2_TYPE::Int((Mnt_CEff * 100.0) / (Total_CEff + EPS));
    Breakdown_file << "%): ";
    Breakdown_file.width(11);
    Breakdown_file << Mnt_CEff << "\n";

    if (nDim == 3) {
      Breakdown_file << "Total CMx:   ";
      Breakdown_file.width(11);
      Breakdown_file << Total_CMx;
      Breakdown_file << " | Pressure (";
      Breakdown_file.width(5);
      Breakdown_file << SU2_TYPE::Int((Inv_CMx * 100.0) / (Total_CMx + EPS));
      Breakdown_file << "%): ";
      Breakdown_file.width(11);
      Breakdown_file << Inv_CMx;
      Breakdown_file << " | Friction (";
      Breakdown_file.width(5);
      Breakdown_file << SU2_TYPE::Int((Visc_CMx * 100.0) / (Total_CMx + EPS));
      Breakdown_file << "%): ";
      Breakdown_file.width(11);
      Breakdown_file << Visc_CMx;
      Breakdown_file << " | Momentum (";
      Breakdown_file.width(5);
      Breakdown_file << SU2_TYPE::Int((Mnt_CMx * 100.0) / (Total_CMx + EPS));
      Breakdown_file << "%): ";
      Breakdown_file.width(11);
      Breakdown_file << Mnt_CMx << "\n";

      Breakdown_file << "Total CMy:   ";
      Breakdown_file.width(11);
      Breakdown_file << Total_CMy;
      Breakdown_file << " | Pressure (";
      Breakdown_file.width(5);
      Breakdown_file << SU2_TYPE::Int((Inv_CMy * 100.0) / (Total_CMy + EPS));
      Breakdown_file << "%): ";
      Breakdown_file.width(11);
      Breakdown_file << Inv_CMy;
      Breakdown_file << " | Friction (";
      Breakdown_file.width(5);
      Breakdown_file << SU2_TYPE::Int((Visc_CMy * 100.0) / (Total_CMy + EPS));
      Breakdown_file << "%): ";
      Breakdown_file.width(11);
      Breakdown_file << Visc_CMy;
      Breakdown_file << " | Momentum (";
      Breakdown_file.width(5);
      Breakdown_file << SU2_TYPE::Int((Mnt_CMz * 100.0) / (Total_CMz + EPS));
      Breakdown_file << "%): ";
      Breakdown_file.width(11);
      Breakdown_file << Mnt_CMy << "\n";
    }

    Breakdown_file << "Total CMz:   ";
    Breakdown_file.width(11);
    Breakdown_file << Total_CMz;
    Breakdown_file << " | Pressure (";
    Breakdown_file.width(5);
    Breakdown_file << SU2_TYPE::Int((Inv_CMz * 100.0) / (Total_CMz + EPS));
    Breakdown_file << "%): ";
    Breakdown_file.width(11);
    Breakdown_file << Inv_CMz;
    Breakdown_file << " | Friction (";
    Breakdown_file.width(5);
    Breakdown_file << SU2_TYPE::Int((Visc_CMz * 100.0) / (Total_CMz + EPS));
    Breakdown_file << "%): ";
    Breakdown_file.width(11);
    Breakdown_file << Visc_CMz;
    Breakdown_file << " | Momentum (";
    Breakdown_file.width(5);
    Breakdown_file << SU2_TYPE::Int((Mnt_CMz * 100.0) / (Total_CMz + EPS));
    Breakdown_file << "%): ";
    Breakdown_file.width(11);
    Breakdown_file << Mnt_CMz << "\n";

    Breakdown_file << "Total CFx:   ";
    Breakdown_file.width(11);
    Breakdown_file << Total_CFx;
    Breakdown_file << " | Pressure (";
    Breakdown_file.width(5);
    Breakdown_file << SU2_TYPE::Int((Inv_CFx * 100.0) / (Total_CFx + EPS));
    Breakdown_file << "%): ";
    Breakdown_file.width(11);
    Breakdown_file << Inv_CFx;
    Breakdown_file << " | Friction (";
    Breakdown_file.width(5);
    Breakdown_file << SU2_TYPE::Int((Visc_CFx * 100.0) / (Total_CFx + EPS));
    Breakdown_file << "%): ";
    Breakdown_file.width(11);
    Breakdown_file << Visc_CFx;
    Breakdown_file << " | Momentum (";
    Breakdown_file.width(5);
    Breakdown_file << SU2_TYPE::Int((Mnt_CFx * 100.0) / (Total_CFx + EPS));
    Breakdown_file << "%): ";
    Breakdown_file.width(11);
    Breakdown_file << Mnt_CFx << "\n";

    Breakdown_file << "Total CFy:   ";
    Breakdown_file.width(11);
    Breakdown_file << Total_CFy;
    Breakdown_file << " | Pressure (";
    Breakdown_file.width(5);
    Breakdown_file << SU2_TYPE::Int((Inv_CFy * 100.0) / (Total_CFy + EPS));
    Breakdown_file << "%): ";
    Breakdown_file.width(11);
    Breakdown_file << Inv_CFy;
    Breakdown_file << " | Friction (";
    Breakdown_file.width(5);
    Breakdown_file << SU2_TYPE::Int((Visc_CFy * 100.0) / (Total_CFy + EPS));
    Breakdown_file << "%): ";
    Breakdown_file.width(11);
    Breakdown_file << Visc_CFy;
    Breakdown_file << " | Momentum (";
    Breakdown_file.width(5);
    Breakdown_file << SU2_TYPE::Int((Mnt_CFy * 100.0) / (Total_CFy + EPS));
    Breakdown_file << "%): ";
    Breakdown_file.width(11);
    Breakdown_file << Mnt_CFy << "\n";

    if (nDim == 3) {
      Breakdown_file << "Total CFz:   ";
      Breakdown_file.width(11);
      Breakdown_file << Total_CFz;
      Breakdown_file << " | Pressure (";
      Breakdown_file.width(5);
      Breakdown_file << SU2_TYPE::Int((Inv_CFz * 100.0) / (Total_CFz + EPS));
      Breakdown_file << "%): ";
      Breakdown_file.width(11);
      Breakdown_file << Inv_CFz;
      Breakdown_file << " | Friction (";
      Breakdown_file.width(5);
      Breakdown_file << SU2_TYPE::Int((Visc_CFz * 100.0) / (Total_CFz + EPS));
      Breakdown_file << "%): ";
      Breakdown_file.width(11);
      Breakdown_file << Visc_CFz;
      Breakdown_file << " | Momentum (";
      Breakdown_file.width(5);
      Breakdown_file << SU2_TYPE::Int((Mnt_CFz * 100.0) / (Total_CFz + EPS));
      Breakdown_file << "%): ";
      Breakdown_file.width(11);
      Breakdown_file << Mnt_CFz << "\n";
    }

    Breakdown_file << "\n" << "\n";

    for (iMarker_Monitoring = 0;
         iMarker_Monitoring < config[val_iZone]->GetnMarker_Monitoring();
         iMarker_Monitoring++) {

      Breakdown_file << "Surface name: "
      << config[val_iZone]->GetMarker_Monitoring_TagBound(
                                                          iMarker_Monitoring) << "\n" << "\n";

      Breakdown_file << "Total CL    (";
      Breakdown_file.width(5);
      Breakdown_file
      << SU2_TYPE::Int(
                       (Surface_CL[iMarker_Monitoring] * 100.0)
                       / (Total_CL + EPS));
      Breakdown_file << "%): ";
      Breakdown_file.width(11);
      Breakdown_file << Surface_CL[iMarker_Monitoring];
      Breakdown_file << " | Pressure (";
      Breakdown_file.width(5);
      Breakdown_file
      << SU2_TYPE::Int(
                       (Surface_CL_Inv[iMarker_Monitoring] * 100.0)
                       / (Surface_CL[iMarker_Monitoring] + EPS));
      Breakdown_file << "%): ";
      Breakdown_file.width(11);
      Breakdown_file << Surface_CL_Inv[iMarker_Monitoring];
      Breakdown_file << " | Friction (";
      Breakdown_file.width(5);
      Breakdown_file
      << SU2_TYPE::Int(
                       (Surface_CL_Visc[iMarker_Monitoring] * 100.0)
                       / (Surface_CL[iMarker_Monitoring] + EPS));
      Breakdown_file << "%): ";
      Breakdown_file.width(11);
      Breakdown_file << Surface_CL_Visc[iMarker_Monitoring];
      Breakdown_file << " | Momentum (";
      Breakdown_file.width(5);
      Breakdown_file
      << SU2_TYPE::Int(
                       (Surface_CL_Mnt[iMarker_Monitoring] * 100.0)
                       / (Surface_CL[iMarker_Monitoring] + EPS));
      Breakdown_file << "%): ";
      Breakdown_file.width(11);
      Breakdown_file << Surface_CL_Mnt[iMarker_Monitoring] << "\n";

      Breakdown_file << "Total CD    (";
      Breakdown_file.width(5);
      Breakdown_file
      << SU2_TYPE::Int(
                       (Surface_CD[iMarker_Monitoring] * 100.0)
                       / (Total_CD + EPS));
      Breakdown_file << "%): ";
      Breakdown_file.width(11);
      Breakdown_file << Surface_CD[iMarker_Monitoring];
      Breakdown_file << " | Pressure (";
      Breakdown_file.width(5);
      Breakdown_file
      << SU2_TYPE::Int(
                       (Surface_CD_Inv[iMarker_Monitoring] * 100.0)
                       / (Surface_CD[iMarker_Monitoring] + EPS));
      Breakdown_file << "%): ";
      Breakdown_file.width(11);
      Breakdown_file << Surface_CD_Inv[iMarker_Monitoring];
      Breakdown_file << " | Friction (";
      Breakdown_file.width(5);
      Breakdown_file
      << SU2_TYPE::Int(
                       (Surface_CD_Visc[iMarker_Monitoring] * 100.0)
                       / (Surface_CD[iMarker_Monitoring] + EPS));
      Breakdown_file << "%): ";
      Breakdown_file.width(11);
      Breakdown_file << Surface_CD_Visc[iMarker_Monitoring];
      Breakdown_file << " | Momentum (";
      Breakdown_file.width(5);
      Breakdown_file
      << SU2_TYPE::Int(
                       (Surface_CD_Mnt[iMarker_Monitoring] * 100.0)
                       / (Surface_CD[iMarker_Monitoring] + EPS));
      Breakdown_file << "%): ";
      Breakdown_file.width(11);
      Breakdown_file << Surface_CD_Mnt[iMarker_Monitoring] << "\n";

      if (nDim == 3) {
        Breakdown_file << "Total CSF   (";
        Breakdown_file.width(5);
        Breakdown_file
        << SU2_TYPE::Int(
                         (Surface_CSF[iMarker_Monitoring] * 100.0)
                         / (Total_CSF + EPS));
        Breakdown_file << "%): ";
        Breakdown_file.width(11);
        Breakdown_file << Surface_CSF[iMarker_Monitoring];
        Breakdown_file << " | Pressure (";
        Breakdown_file.width(5);
        Breakdown_file
        << SU2_TYPE::Int(
                         (Surface_CSF_Inv[iMarker_Monitoring] * 100.0)
                         / (Surface_CSF[iMarker_Monitoring] + EPS));
        Breakdown_file << "%): ";
        Breakdown_file.width(11);
        Breakdown_file << Surface_CSF_Inv[iMarker_Monitoring];
        Breakdown_file << " | Friction (";
        Breakdown_file.width(5);
        Breakdown_file
        << SU2_TYPE::Int(
                         (Surface_CSF_Visc[iMarker_Monitoring] * 100.0)
                         / (Surface_CSF[iMarker_Monitoring] + EPS));
        Breakdown_file << "%): ";
        Breakdown_file.width(11);
        Breakdown_file
        << Surface_CSF_Visc[iMarker_Monitoring];
        Breakdown_file << " | Momentum (";
        Breakdown_file.width(5);
        Breakdown_file
        << SU2_TYPE::Int(
                         (Surface_CSF_Mnt[iMarker_Monitoring] * 100.0)
                         / (Surface_CSF[iMarker_Monitoring] + EPS));
        Breakdown_file << "%): ";
        Breakdown_file.width(11);
        Breakdown_file
        << Surface_CSF_Mnt[iMarker_Monitoring] << "\n";
      }

      Breakdown_file << "Total CL/CD (";
      Breakdown_file.width(5);
      Breakdown_file
      << SU2_TYPE::Int(
                       (Surface_CEff[iMarker_Monitoring] * 100.0) / (Total_CEff + EPS));
      Breakdown_file << "%): ";
      Breakdown_file.width(11);
      Breakdown_file << Surface_CEff[iMarker_Monitoring];
      Breakdown_file << " | Pressure (";
      Breakdown_file.width(5);
      Breakdown_file
      << SU2_TYPE::Int(
                       (Surface_CEff_Inv[iMarker_Monitoring] * 100.0)
                       / (Surface_CEff[iMarker_Monitoring] + EPS));
      Breakdown_file << "%): ";
      Breakdown_file.width(11);
      Breakdown_file << Surface_CEff_Inv[iMarker_Monitoring];
      Breakdown_file << " | Friction (";
      Breakdown_file.width(5);
      Breakdown_file
      << SU2_TYPE::Int(
                       (Surface_CEff_Visc[iMarker_Monitoring] * 100.0)
                       / (Surface_CEff[iMarker_Monitoring] + EPS));
      Breakdown_file << "%): ";
      Breakdown_file.width(11);
      Breakdown_file
      << Surface_CEff_Visc[iMarker_Monitoring];
      Breakdown_file << " | Momentum (";
      Breakdown_file.width(5);
      Breakdown_file
      << SU2_TYPE::Int(
                       (Surface_CEff_Mnt[iMarker_Monitoring] * 100.0)
                       / (Surface_CEff[iMarker_Monitoring] + EPS));
      Breakdown_file << "%): ";
      Breakdown_file.width(11);
      Breakdown_file
      << Surface_CEff_Mnt[iMarker_Monitoring] << "\n";

      if (nDim == 3) {

        Breakdown_file << "Total CMx   (";
        Breakdown_file.width(5);
        Breakdown_file
        << SU2_TYPE::Int(
                         (Surface_CMx[iMarker_Monitoring] * 100.0) / (Total_CMx + EPS));
        Breakdown_file << "%): ";
        Breakdown_file.width(11);
        Breakdown_file << Surface_CMx[iMarker_Monitoring];
        Breakdown_file << " | Pressure (";
        Breakdown_file.width(5);
        Breakdown_file
        << SU2_TYPE::Int(
                         (Surface_CMx_Inv[iMarker_Monitoring] * 100.0)
                         / (Surface_CMx[iMarker_Monitoring] + EPS));
        Breakdown_file << "%): ";
        Breakdown_file.width(11);
        Breakdown_file << Surface_CMx_Inv[iMarker_Monitoring];
        Breakdown_file << " | Friction (";
        Breakdown_file.width(5);
        Breakdown_file
        << SU2_TYPE::Int(
                         (Surface_CMx_Visc[iMarker_Monitoring] * 100.0)
                         / (Surface_CMx[iMarker_Monitoring] + EPS));
        Breakdown_file << "%): ";
        Breakdown_file.width(11);
        Breakdown_file
        << Surface_CMx_Visc[iMarker_Monitoring];
        Breakdown_file << " | Momentum (";
        Breakdown_file.width(5);
        Breakdown_file
        << SU2_TYPE::Int(
                         (Surface_CMx_Mnt[iMarker_Monitoring] * 100.0)
                         / (Surface_CMx[iMarker_Monitoring] + EPS));
        Breakdown_file << "%): ";
        Breakdown_file.width(11);
        Breakdown_file
        << Surface_CMx_Mnt[iMarker_Monitoring] << "\n";

        Breakdown_file << "Total CMy   (";
        Breakdown_file.width(5);
        Breakdown_file
        << SU2_TYPE::Int(
                         (Surface_CMy[iMarker_Monitoring] * 100.0) / (Total_CMy + EPS));
        Breakdown_file << "%): ";
        Breakdown_file.width(11);
        Breakdown_file << Surface_CMy[iMarker_Monitoring];
        Breakdown_file << " | Pressure (";
        Breakdown_file.width(5);
        Breakdown_file
        << SU2_TYPE::Int(
                         (Surface_CMy_Inv[iMarker_Monitoring] * 100.0)
                         / (Surface_CMy[iMarker_Monitoring] + EPS));
        Breakdown_file << "%): ";
        Breakdown_file.width(11);
        Breakdown_file << Surface_CMy_Inv[iMarker_Monitoring];
        Breakdown_file << " | Friction (";
        Breakdown_file.width(5);
        Breakdown_file
        << SU2_TYPE::Int(
                         (Surface_CMy_Visc[iMarker_Monitoring] * 100.0)
                         / (Surface_CMy[iMarker_Monitoring] + EPS));
        Breakdown_file << "%): ";
        Breakdown_file.width(11);
        Breakdown_file
        << Surface_CMy_Visc[iMarker_Monitoring];
        Breakdown_file << " | Momentum (";
        Breakdown_file.width(5);
        Breakdown_file
        << SU2_TYPE::Int(
                         (Surface_CMy_Mnt[iMarker_Monitoring] * 100.0)
                         / (Surface_CMy[iMarker_Monitoring] + EPS));
        Breakdown_file << "%): ";
        Breakdown_file.width(11);
        Breakdown_file
        << Surface_CMy_Mnt[iMarker_Monitoring] << "\n";
      }

      Breakdown_file << "Total CMz   (";
      Breakdown_file.width(5);
      Breakdown_file
      << SU2_TYPE::Int((Surface_CMz[iMarker_Monitoring] * 100.0) / (Total_CMz + EPS));
      Breakdown_file << "%): ";
      Breakdown_file.width(11);
      Breakdown_file << Surface_CMz[iMarker_Monitoring];
      Breakdown_file << " | Pressure (";
      Breakdown_file.width(5);
      Breakdown_file
      << SU2_TYPE::Int(
                       (Surface_CMz_Inv[iMarker_Monitoring] * 100.0)
                       / (Surface_CMz[iMarker_Monitoring] + EPS));
      Breakdown_file << "%): ";
      Breakdown_file.width(11);
      Breakdown_file << Surface_CMz_Inv[iMarker_Monitoring];
      Breakdown_file << " | Friction (";
      Breakdown_file.width(5);
      Breakdown_file
      << SU2_TYPE::Int(
                       (Surface_CMz_Visc[iMarker_Monitoring] * 100.0)
                       / (Surface_CMz[iMarker_Monitoring] + EPS));
      Breakdown_file << "%): ";
      Breakdown_file.width(11);
      Breakdown_file
      << Surface_CMz_Visc[iMarker_Monitoring];
      Breakdown_file << " | Momentum (";
      Breakdown_file.width(5);
      Breakdown_file
      << SU2_TYPE::Int(
                       (Surface_CMz_Mnt[iMarker_Monitoring] * 100.0)
                       / (Surface_CMz[iMarker_Monitoring] + EPS));
      Breakdown_file << "%): ";
      Breakdown_file.width(11);
      Breakdown_file
      << Surface_CMz_Mnt[iMarker_Monitoring] << "\n";

      Breakdown_file << "Total CFx   (";
      Breakdown_file.width(5);
      Breakdown_file
      << SU2_TYPE::Int((Surface_CFx[iMarker_Monitoring] * 100.0) / (Total_CFx + EPS));
      Breakdown_file << "%): ";
      Breakdown_file.width(11);
      Breakdown_file << Surface_CFx[iMarker_Monitoring];
      Breakdown_file << " | Pressure (";
      Breakdown_file.width(5);
      Breakdown_file
      << SU2_TYPE::Int(
                       (Surface_CFx_Inv[iMarker_Monitoring] * 100.0)
                       / (Surface_CFx[iMarker_Monitoring] + EPS));
      Breakdown_file << "%): ";
      Breakdown_file.width(11);
      Breakdown_file << Surface_CFx_Inv[iMarker_Monitoring];
      Breakdown_file << " | Friction (";
      Breakdown_file.width(5);
      Breakdown_file
      << SU2_TYPE::Int(
                       (Surface_CFx_Visc[iMarker_Monitoring] * 100.0)
                       / (Surface_CFx[iMarker_Monitoring] + EPS));
      Breakdown_file << "%): ";
      Breakdown_file.width(11);
      Breakdown_file
      << Surface_CFx_Visc[iMarker_Monitoring];
      Breakdown_file << " | Momentum (";
      Breakdown_file.width(5);
      Breakdown_file
      << SU2_TYPE::Int(
                       (Surface_CFx_Mnt[iMarker_Monitoring] * 100.0)
                       / (Surface_CFx[iMarker_Monitoring] + EPS));
      Breakdown_file << "%): ";
      Breakdown_file.width(11);
      Breakdown_file
      << Surface_CFx_Mnt[iMarker_Monitoring] << "\n";

      Breakdown_file << "Total CFy   (";
      Breakdown_file.width(5);
      Breakdown_file
      << SU2_TYPE::Int((Surface_CFy[iMarker_Monitoring] * 100.0) / (Total_CFy + EPS));
      Breakdown_file << "%): ";
      Breakdown_file.width(11);
      Breakdown_file << Surface_CFy[iMarker_Monitoring];
      Breakdown_file << " | Pressure (";
      Breakdown_file.width(5);
      Breakdown_file
      << SU2_TYPE::Int(
                       (Surface_CFy_Inv[iMarker_Monitoring] * 100.0)
                       / (Surface_CFy[iMarker_Monitoring] + EPS));
      Breakdown_file << "%): ";
      Breakdown_file.width(11);
      Breakdown_file << Surface_CFy_Inv[iMarker_Monitoring];
      Breakdown_file << " | Friction (";
      Breakdown_file.width(5);
      Breakdown_file
      << SU2_TYPE::Int(
                       (Surface_CFy_Visc[iMarker_Monitoring] * 100.0)
                       / (Surface_CFy[iMarker_Monitoring] + EPS));
      Breakdown_file << "%): ";
      Breakdown_file.width(11);
      Breakdown_file
      << Surface_CFy_Visc[iMarker_Monitoring];
      Breakdown_file << " | Momentum (";
      Breakdown_file.width(5);
      Breakdown_file
      << SU2_TYPE::Int(
                       (Surface_CFy_Mnt[iMarker_Monitoring] * 100.0)
                       / (Surface_CFy[iMarker_Monitoring] + EPS));
      Breakdown_file << "%): ";
      Breakdown_file.width(11);
      Breakdown_file
      << Surface_CFy_Mnt[iMarker_Monitoring] << "\n";

      if (nDim == 3) {
        Breakdown_file << "Total CFz   (";
        Breakdown_file.width(5);
        Breakdown_file
        << SU2_TYPE::Int(
                         (Surface_CFz[iMarker_Monitoring] * 100.0) / (Total_CFz + EPS));
        Breakdown_file << "%): ";
        Breakdown_file.width(11);
        Breakdown_file << Surface_CFz[iMarker_Monitoring];
        Breakdown_file << " | Pressure (";
        Breakdown_file.width(5);
        Breakdown_file
        << SU2_TYPE::Int(
                         (Surface_CFz_Inv[iMarker_Monitoring] * 100.0)
                         / (Surface_CFz[iMarker_Monitoring] + EPS));
        Breakdown_file << "%): ";
        Breakdown_file.width(11);
        Breakdown_file << Surface_CFz_Inv[iMarker_Monitoring];
        Breakdown_file << " | Friction (";
        Breakdown_file.width(5);
        Breakdown_file
        << SU2_TYPE::Int(
                         (Surface_CFz_Visc[iMarker_Monitoring] * 100.0)
                         / (Surface_CFz[iMarker_Monitoring] + EPS));
        Breakdown_file << "%): ";
        Breakdown_file.width(11);
        Breakdown_file
        << Surface_CFz_Visc[iMarker_Monitoring];
        Breakdown_file << " | Momentum (";
        Breakdown_file.width(5);
        Breakdown_file
        << SU2_TYPE::Int(
                         (Surface_CFz_Mnt[iMarker_Monitoring] * 100.0)
                         / (Surface_CFz[iMarker_Monitoring] + EPS));
        Breakdown_file << "%): ";
        Breakdown_file.width(11);
        Breakdown_file
        << Surface_CFz_Mnt[iMarker_Monitoring] << "\n";

      }

      Breakdown_file << "\n";


    }

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

    delete [] Surface_CL_Inv;
    delete [] Surface_CD_Inv;
    delete [] Surface_CSF_Inv;
    delete [] Surface_CEff_Inv;
    delete [] Surface_CFx_Inv;
    delete [] Surface_CFy_Inv;
    delete [] Surface_CFz_Inv;
    delete [] Surface_CMx_Inv;
    delete [] Surface_CMy_Inv;
    delete [] Surface_CMz_Inv;

    delete [] Surface_CL_Visc;
    delete [] Surface_CD_Visc;
    delete [] Surface_CSF_Visc;
    delete [] Surface_CEff_Visc;
    delete [] Surface_CFx_Visc;
    delete [] Surface_CFy_Visc;
    delete [] Surface_CFz_Visc;
    delete [] Surface_CMx_Visc;
    delete [] Surface_CMy_Visc;
    delete [] Surface_CMz_Visc;

    delete [] Surface_CL_Mnt;
    delete [] Surface_CD_Mnt;
    delete [] Surface_CSF_Mnt;
    delete [] Surface_CEff_Mnt;
    delete [] Surface_CFx_Mnt;
    delete [] Surface_CFy_Mnt;
    delete [] Surface_CFz_Mnt;
    delete [] Surface_CMx_Mnt;
    delete [] Surface_CMy_Mnt;
    delete [] Surface_CMz_Mnt;

    Breakdown_file.close();

  }

}

void COutputLegacy::SpecialOutput_SpanLoad(CSolver *solver, CGeometry *geometry, CConfig *config, bool output) const {

  short iSection, nSection;
  unsigned long iVertex, iPoint, Trailing_Point;
  su2double *Plane_P0, *Plane_P0_, *Plane_Normal, *Plane_Normal_, *CPressure,
  Force[3], ForceInviscid[3], MomentInviscid[3] =
  { 0.0, 0.0, 0.0 }, MomentDist[3] = { 0.0, 0.0, 0.0 }, RefDensity,
  RefPressure, RefArea, *Velocity_Inf, Gas_Constant, Mach2Vel,
  Mach_Motion, Gamma, RefVel2 = 0.0, factor, NDPressure,
  RefLength, Alpha, CL_Inv,
  Xcoord_LeadingEdge = 0.0, Ycoord_LeadingEdge = 0.0, Zcoord_LeadingEdge = 0.0,
  Xcoord_TrailingEdge = 0.0, Ycoord_TrailingEdge = 0.0, Zcoord_TrailingEdge = 0.0,
  Xcoord_LeadingEdge_ = 0.0,
  Xcoord_TrailingEdge_ = 0.0, Ycoord_TrailingEdge_ = 0.0, Zcoord_TrailingEdge_ = 0.0,
  MaxDistance, Distance, Chord, Aux, Dihedral_Trailing;

  su2double B, Y, C_L, C_L0, Elliptic_Spanload;

  vector<su2double> Xcoord_Airfoil, Ycoord_Airfoil, Zcoord_Airfoil,
  CPressure_Airfoil;
  vector<su2double> Xcoord_Airfoil_, Ycoord_Airfoil_, Zcoord_Airfoil_,
  CPressure_Airfoil_;
  string Marker_Tag, Slice_Filename, Slice_Ext;
  ofstream Cp_File;
  unsigned short iDim;

  bool dynamic_grid = config->GetDynamic_Grid();

  Plane_P0 = new su2double[3];
  Plane_P0_ = new su2double[3];
  Plane_Normal = new su2double[3];
  Plane_Normal_ = new su2double[3];
  CPressure = new su2double[geometry->GetnPoint()];

  if ((rank == MASTER_NODE) && (output)) {
    cout << endl << "Writing the spanload file (load_distribution.dat).";
  }

  /*--- Compute some reference quantities and necessary values ---*/

  RefDensity = solver->GetDensity_Inf();
  RefPressure = solver->GetPressure_Inf();
  RefArea = config->GetRefArea();
  Velocity_Inf = solver->GetVelocity_Inf();
  Gamma = config->GetGamma();
  const auto Origin = config->GetRefOriginMoment(0);
  RefLength = config->GetRefLength();
  Alpha = config->GetAoA() * PI_NUMBER / 180.0;

  if (dynamic_grid) {
    Gas_Constant = config->GetGas_ConstantND();
    Mach2Vel = sqrt(
                    Gamma * Gas_Constant * config->GetTemperature_FreeStreamND());
    Mach_Motion = config->GetMach_Motion();
    RefVel2 = (Mach_Motion * Mach2Vel) * (Mach_Motion * Mach2Vel);
  } else {
    RefVel2 = 0.0;
    for (iDim = 0; iDim < geometry->GetnDim(); iDim++)
      RefVel2 += Velocity_Inf[iDim] * Velocity_Inf[iDim];
  }
  factor = 1.0 / (0.5 * RefDensity * RefArea * RefVel2);

  if (geometry->GetnDim() == 3) {

    /*--- Copy the pressure to an auxiliar structure ---*/

    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
      CPressure[iPoint] = (solver->GetNodes()->GetPressure(iPoint)
                           - RefPressure) * factor * RefArea;
    }

    nSection = config->GetnLocationStations();

    for (iSection = 0; iSection < nSection; iSection++) {

      /*--- Read the values from the config file ---*/

      Plane_Normal[0] = 0.0; Plane_P0[0] = 0.0;
      Plane_Normal[1] = 0.0; Plane_P0[1] = 0.0;
      Plane_Normal[2] = 0.0; Plane_P0[2] = 0.0;

      if (config->GetGeo_Description() == FUSELAGE) {
        Plane_Normal[0] = 1.0;
        Plane_P0[0] = config->GetLocationStations(iSection);
      }
      else if (config->GetGeo_Description() == NACELLE) {
        Plane_Normal[0] = 0.0;
        Plane_Normal[1] = -sin(config->GetLocationStations(iSection)*PI_NUMBER/180.0);
        Plane_Normal[2] = cos(config->GetLocationStations(iSection)*PI_NUMBER/180.0);

        /*--- Apply tilt angle to the plane ---*/

        su2double Tilt_Angle = config->GetNacelleLocation(3)*PI_NUMBER/180;
        su2double Plane_NormalX_Tilt = Plane_Normal[0]*cos(Tilt_Angle) + Plane_Normal[2]*sin(Tilt_Angle);
        su2double Plane_NormalY_Tilt = Plane_Normal[1];
        su2double Plane_NormalZ_Tilt = Plane_Normal[2]*cos(Tilt_Angle) - Plane_Normal[0]*sin(Tilt_Angle);

        /*--- Apply toe angle to the plane ---*/

        su2double Toe_Angle = config->GetNacelleLocation(4)*PI_NUMBER/180;
        su2double Plane_NormalX_Tilt_Toe = Plane_NormalX_Tilt*cos(Toe_Angle) - Plane_NormalY_Tilt*sin(Toe_Angle);
        su2double Plane_NormalY_Tilt_Toe = Plane_NormalX_Tilt*sin(Toe_Angle) + Plane_NormalY_Tilt*cos(Toe_Angle);
        su2double Plane_NormalZ_Tilt_Toe = Plane_NormalZ_Tilt;

        /*--- Update normal vector ---*/

        Plane_Normal[0] = Plane_NormalX_Tilt_Toe;
        Plane_Normal[1] = Plane_NormalY_Tilt_Toe;
        Plane_Normal[2] = Plane_NormalZ_Tilt_Toe;

      }
      else {
        Plane_Normal[1] = 1.0;
        Plane_P0[1] = config->GetLocationStations(iSection);
      }

      /*--- Compute the airfoil sections (note that we feed in the Cp) ---*/

      geometry->ComputeAirfoil_Section(Plane_P0, Plane_Normal, -1E6, 1E6, -1E6, 1E6, -1E6, 1E6,
                                       CPressure, Xcoord_Airfoil, Ycoord_Airfoil, Zcoord_Airfoil,
                                       CPressure_Airfoil, true, config);

      if ((rank == MASTER_NODE) && (Xcoord_Airfoil.size() == 0)) {
        if ((config->GetGeo_Description() == FUSELAGE) || (config->GetGeo_Description() == WING))
          cout << endl << "Please check the config file, the section (" << Plane_P0[0] <<", " << Plane_P0[1] <<", " << Plane_P0[2] << ") has not been detected." << endl;
        if (config->GetGeo_Description() == NACELLE)
          cout << endl << "Please check the config file, the section (" << Plane_Normal[0] <<", " << Plane_Normal[1] <<", " << Plane_Normal[2] << ") has not been detected." << endl;
      }


      /*--- Compute dihedral using a step in the station value ---*/

      Plane_P0_[0] = 0.0; Plane_Normal_[0] = 0.0;
      Plane_P0_[1] = 0.0; Plane_Normal_[1] = 0.0;
      Plane_P0_[2] = 0.0; Plane_Normal_[2] = 0.0;

      if (config->GetGeo_Description() == FUSELAGE) {
        Plane_Normal_[0] = 1.0;
        if (iSection == 0) Plane_P0_[0] = config->GetLocationStations(iSection) + 0.01;
        else Plane_P0_[0] = config->GetLocationStations(iSection) - 0.01;
      }
      else if (config->GetGeo_Description() == NACELLE) {
        if (iSection == 0) {
          Plane_Normal_[0] = 0.0;
          Plane_Normal_[1] = -sin((config->GetLocationStations(iSection) + 0.01)*PI_NUMBER/180.0);
          Plane_Normal_[2] = cos((config->GetLocationStations(iSection) + 0.01)*PI_NUMBER/180.0);

          /*--- Apply tilt angle to the plane ---*/

          su2double Tilt_Angle = config->GetNacelleLocation(3)*PI_NUMBER/180;
          su2double Plane_NormalX_Tilt = Plane_Normal[0]*cos(Tilt_Angle) + Plane_Normal[2]*sin(Tilt_Angle);
          su2double Plane_NormalY_Tilt = Plane_Normal[1];
          su2double Plane_NormalZ_Tilt = Plane_Normal[2]*cos(Tilt_Angle) - Plane_Normal[0]*sin(Tilt_Angle);

          /*--- Apply toe angle to the plane ---*/

          su2double Toe_Angle = config->GetNacelleLocation(4)*PI_NUMBER/180;
          su2double Plane_NormalX_Tilt_Toe = Plane_NormalX_Tilt*cos(Toe_Angle) - Plane_NormalY_Tilt*sin(Toe_Angle);
          su2double Plane_NormalY_Tilt_Toe = Plane_NormalX_Tilt*sin(Toe_Angle) + Plane_NormalY_Tilt*cos(Toe_Angle);
          su2double Plane_NormalZ_Tilt_Toe = Plane_NormalZ_Tilt;

          /*--- Update normal vector ---*/

          Plane_Normal[0] = Plane_NormalX_Tilt_Toe;
          Plane_Normal[1] = Plane_NormalY_Tilt_Toe;
          Plane_Normal[2] = Plane_NormalZ_Tilt_Toe;

        }
        else {
          Plane_Normal_[0] = 0.0;
          Plane_Normal_[1] = -sin((config->GetLocationStations(iSection) - 0.01)*PI_NUMBER/180.0);
          Plane_Normal_[2] = cos((config->GetLocationStations(iSection) - 0.01)*PI_NUMBER/180.0);

          /*--- Apply tilt angle to the plane ---*/

          su2double Tilt_Angle = config->GetNacelleLocation(3)*PI_NUMBER/180;
          su2double Plane_NormalX_Tilt = Plane_Normal[0]*cos(Tilt_Angle) + Plane_Normal[2]*sin(Tilt_Angle);
          su2double Plane_NormalY_Tilt = Plane_Normal[1];
          su2double Plane_NormalZ_Tilt = Plane_Normal[2]*cos(Tilt_Angle) - Plane_Normal[0]*sin(Tilt_Angle);

          /*--- Apply toe angle to the plane ---*/

          su2double Toe_Angle = config->GetNacelleLocation(4)*PI_NUMBER/180;
          su2double Plane_NormalX_Tilt_Toe = Plane_NormalX_Tilt*cos(Toe_Angle) - Plane_NormalY_Tilt*sin(Toe_Angle);
          su2double Plane_NormalY_Tilt_Toe = Plane_NormalX_Tilt*sin(Toe_Angle) + Plane_NormalY_Tilt*cos(Toe_Angle);
          su2double Plane_NormalZ_Tilt_Toe = Plane_NormalZ_Tilt;

          /*--- Update normal vector ---*/

          Plane_Normal[0] = Plane_NormalX_Tilt_Toe;
          Plane_Normal[1] = Plane_NormalY_Tilt_Toe;
          Plane_Normal[2] = Plane_NormalZ_Tilt_Toe;

        }
      }
      else {
        Plane_Normal_[1] = 1.0;
        if (iSection == 0) Plane_P0_[1] = config->GetLocationStations(iSection) + 0.01;
        else Plane_P0_[1] = config->GetLocationStations(iSection) - 0.01;
      }

      geometry->ComputeAirfoil_Section(Plane_P0_, Plane_Normal_, -1E6, 1E6, -1E6, 1E6, -1E6, 1E6,
                                       CPressure, Xcoord_Airfoil_, Ycoord_Airfoil_, Zcoord_Airfoil_,
                                       CPressure_Airfoil_, true, config);

      /*--- Output the pressure on each section (tecplot format) ---*/

      if ((rank == MASTER_NODE) && (Xcoord_Airfoil.size() != 0)) {

        /*--- Find leading and trailing edge ---*/

        Xcoord_LeadingEdge = 1E6; Xcoord_TrailingEdge = -1E6;
        for (iVertex = 0; iVertex < Xcoord_Airfoil.size(); iVertex++) {
          if (Xcoord_Airfoil[iVertex] < Xcoord_LeadingEdge) {
            Xcoord_LeadingEdge = Xcoord_Airfoil[iVertex];
            Ycoord_LeadingEdge = Ycoord_Airfoil[iVertex];
            Zcoord_LeadingEdge = Zcoord_Airfoil[iVertex];
          }
          if (Xcoord_Airfoil[iVertex] > Xcoord_TrailingEdge) {
            Xcoord_TrailingEdge = Xcoord_Airfoil[iVertex];
            Ycoord_TrailingEdge = Ycoord_Airfoil[iVertex];
            Zcoord_TrailingEdge = Zcoord_Airfoil[iVertex];
          }
        }

        Chord = (Xcoord_TrailingEdge-Xcoord_LeadingEdge);

        /*--- Compute dihedral ---*/

        Xcoord_LeadingEdge_ = 1E6; Xcoord_TrailingEdge_ = -1E6;
        for (iVertex = 0; iVertex < Xcoord_Airfoil_.size(); iVertex++) {
          if (Xcoord_Airfoil_[iVertex] < Xcoord_LeadingEdge_) {
            Xcoord_LeadingEdge_ = Xcoord_Airfoil_[iVertex];
          }
          if (Xcoord_Airfoil_[iVertex] > Xcoord_TrailingEdge_) {
            Xcoord_TrailingEdge_ = Xcoord_Airfoil_[iVertex];
            Ycoord_TrailingEdge_ = Ycoord_Airfoil_[iVertex];
            Zcoord_TrailingEdge_ = Zcoord_Airfoil_[iVertex];
          }
        }

        if (iSection == 0) {
          Dihedral_Trailing = atan((Zcoord_TrailingEdge_ - Zcoord_TrailingEdge) / (Ycoord_TrailingEdge_ - Ycoord_TrailingEdge))*180/PI_NUMBER;
        }
        else {
          Dihedral_Trailing = atan((Zcoord_TrailingEdge - Zcoord_TrailingEdge_) / (Ycoord_TrailingEdge - Ycoord_TrailingEdge_))*180/PI_NUMBER;
        }

        /*--- Write Cp at each section (tecplot format) ---*/

        if (output) {

          ofstream Cp_File;

          if (iSection == 0) {
            Cp_File.open("cp_sections.dat", ios::out);
            Cp_File << "TITLE = \"Airfoil sections\"" << endl;
            Cp_File << "VARIABLES = \"x/c\",\"C<sub>p</sub>\",\"x\",\"y\",\"z\",\"y/c\",\"z/c\"" << endl;
          } else
            Cp_File.open("cp_sections.dat", ios::app);

          if (config->GetGeo_Description() == NACELLE) {
            su2double theta_deg = atan2(Plane_Normal[1], -Plane_Normal[2])/PI_NUMBER*180 + 180;
            Cp_File << "ZONE T=\"Theta = " << theta_deg << " deg\", I= " << Xcoord_Airfoil.size() << ", F=POINT" << "\n";
          }
          else {
            if (config->GetSystemMeasurements() == SI) Cp_File << "ZONE T=\"y = " << Plane_P0[1] << " m\", I= "
              << Xcoord_Airfoil.size() << ", F=POINT" << "\n";

            if (config->GetSystemMeasurements() == US) Cp_File << "ZONE T=\"y = " << Plane_P0[1]*12.0 << " in\", I= "
              << Xcoord_Airfoil.size() << ", F=POINT" << "\n";
          }



          /*--- Coordinates and pressure value ---*/

          for (iVertex = 0; iVertex < Xcoord_Airfoil.size(); iVertex++) {

            su2double XCoord = Xcoord_Airfoil[iVertex];
            su2double YCoord = Ycoord_Airfoil[iVertex];
            su2double ZCoord = Zcoord_Airfoil[iVertex];

            /*--- Undo the transformation based on the Theta angle ---*/

            if (config->GetGeo_Description() == NACELLE) {
              su2double theta_deg = atan2(Plane_Normal[1],-Plane_Normal[2])/PI_NUMBER*180 + 180;
              su2double Angle = theta_deg*PI_NUMBER/180 - 0.5*PI_NUMBER;

              XCoord = Xcoord_Airfoil[iVertex] + config->GetNacelleLocation(0);
              YCoord = (Ycoord_Airfoil[iVertex]*cos(Angle) - Zcoord_Airfoil[iVertex]*sin(Angle)) + config->GetNacelleLocation(1);
              ZCoord = (Zcoord_Airfoil[iVertex]*cos(Angle) + Ycoord_Airfoil[iVertex]*sin(Angle)) + config->GetNacelleLocation(2);

            }

            if (config->GetSystemMeasurements() == US) {
              Cp_File <<  (Xcoord_Airfoil[iVertex] - Xcoord_LeadingEdge) / Chord  << " " << CPressure_Airfoil[iVertex]
              << " " << XCoord * 12.0 << " " << YCoord * 12.0 << " " << ZCoord * 12.0
              << " " << (Ycoord_Airfoil[iVertex] - Ycoord_LeadingEdge) / Chord << " " << (Zcoord_Airfoil[iVertex] - Zcoord_LeadingEdge)  / Chord << "\n";
            }
            else {
              Cp_File << (Xcoord_Airfoil[iVertex] - Xcoord_LeadingEdge) / Chord << " " << CPressure_Airfoil[iVertex]
              << " " << XCoord << " " << YCoord << " " << ZCoord
              << " " << (Ycoord_Airfoil[iVertex] - Ycoord_LeadingEdge) / Chord << " " << (Zcoord_Airfoil[iVertex] - Zcoord_LeadingEdge)  / Chord << "\n";
            }

          }

          Cp_File.close();

        }

        /*--- Compute load distribution ---*/

        ForceInviscid[0] = 0.0; ForceInviscid[1] = 0.0; ForceInviscid[2] = 0.0; MomentInviscid[1] = 0.0;

        for (iVertex = 0; iVertex < Xcoord_Airfoil.size() - 1; iVertex++) {

          NDPressure = 0.5 * (CPressure_Airfoil[iVertex] + CPressure_Airfoil[iVertex + 1]);

          Force[0] = -(Zcoord_Airfoil[iVertex + 1] - Zcoord_Airfoil[iVertex]) * NDPressure;
          Force[1] = 0.0;
          Force[2] = (Xcoord_Airfoil[iVertex + 1] - Xcoord_Airfoil[iVertex]) * NDPressure;

          ForceInviscid[0] += Force[0];
          ForceInviscid[1] += Force[1];
          ForceInviscid[2] += Force[2];

          MomentDist[0] = 0.5 * (Xcoord_Airfoil[iVertex] + Xcoord_Airfoil[iVertex + 1]) - Origin[0];
          MomentDist[1] = 0.5 * (Ycoord_Airfoil[iVertex] + Ycoord_Airfoil[iVertex + 1]) - Origin[1];
          MomentDist[2] = 0.5 * (Zcoord_Airfoil[iVertex] + Zcoord_Airfoil[iVertex + 1]) - Origin[3];

          MomentInviscid[1] += (Force[0] * MomentDist[2] - Force[2] * MomentDist[0]) / RefLength;

        }

        /*--- Compute local chord, for the nondimensionalization  ---*/

        MaxDistance = 0.0; Trailing_Point = 0;

        for (iVertex = 1; iVertex < Xcoord_Airfoil.size(); iVertex++) {

          Distance = sqrt(pow(Xcoord_Airfoil[iVertex] - Xcoord_Airfoil[Trailing_Point], 2.0) +
                          pow(Ycoord_Airfoil[iVertex] - Ycoord_Airfoil[Trailing_Point], 2.0) +
                          pow(Zcoord_Airfoil[iVertex] - Zcoord_Airfoil[Trailing_Point], 2.0));

          if (MaxDistance < Distance) { MaxDistance = Distance; }

        }

        Chord = MaxDistance;

        CL_Inv = cos(Dihedral_Trailing * PI_NUMBER / 180.0) * fabs( -ForceInviscid[0] * sin(Alpha) + ForceInviscid[2] * cos(Alpha) )/ Chord;

        /*--- Compute sectional lift at the root ---*/

        B                  = 2.0*config->GetSemiSpan();
        RefArea            = config->GetRefArea();
        C_L                = solver->GetTotal_CL();
        C_L0               = 8.0*C_L*RefArea/(B*PI_NUMBER);
        Y                  = Ycoord_Airfoil[0];
        Aux                = Y/(0.5*B);
        Elliptic_Spanload  = (C_L0 / RefLength) * sqrt(fabs(1.0-Aux*Aux));


        /*--- Write load distribution ---*/

        if (output) {

          ofstream Load_File;
          if (iSection == 0) {

            if ((config->GetTabular_FileFormat() == TAB_OUTPUT::TAB_CSV)) {
              Load_File.open("load_distribution.csv", ios::out);
              Load_File << "\"Percent Semispan\",\"Sectional C_L\",\"Spanload (c C_L / c_ref) \",\"Elliptic Spanload\"" << endl;
            }
            else {
              Load_File.open("load_distribution.dat", ios::out);
              Load_File << "TITLE = \"Load distribution\"" << endl;
              Load_File << "VARIABLES = \"Percent Semispan\",\"Sectional C<sub>L</sub>\",\"Spanload (c C<sub>L</sub> / c<sub>ref</sub>) \",\"Elliptic Spanload\"" << endl;
              Load_File << "ZONE T=\"Wing load distribution\"" << endl;
            }
          } else {
            if ((config->GetTabular_FileFormat() == TAB_OUTPUT::TAB_CSV)) Load_File.open("load_distribution.csv", ios::app);
            else Load_File.open("load_distribution.dat", ios::app);
          }


          /*--- CL and spanload ---*/

          if ((config->GetTabular_FileFormat() == TAB_OUTPUT::TAB_CSV))
            Load_File << 100.0*Ycoord_Airfoil[0]/(0.5*B) << ", " << CL_Inv  << ", " << Chord*CL_Inv / RefLength <<", " << Elliptic_Spanload   << endl;
          else
            Load_File << 100.0*Ycoord_Airfoil[0]/(0.5*B) << " " << CL_Inv  << " " << Chord*CL_Inv / RefLength <<" " << Elliptic_Spanload   << endl;

          Load_File.close();

        }

      }

    }

  }

  /*--- Delete dynamically allocated memory ---*/

  delete[] Plane_P0;
  delete[] Plane_P0_;
  delete[] Plane_Normal;
  delete[] Plane_Normal_;
  delete[] CPressure;

}


void COutputLegacy::SetCp_InverseDesign(CSolver *solver_container, CGeometry *geometry, CConfig *config, unsigned long iExtIter) {

  unsigned short iMarker, icommas, Boundary, iDim;
  unsigned long iVertex, iPoint, (*Point2Vertex)[2], nPointLocal = 0, nPointGlobal = 0;
  su2double XCoord, YCoord, ZCoord, Pressure, PressureCoeff = 0, Cp, CpTarget, *Normal = nullptr, Area, PressDiff;
  bool *PointInDomain;
  string text_line, surfCp_filename;
  ifstream Surface_file;
  char buffer[50], cstr[200];


  nPointLocal = geometry->GetnPoint();
#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&nPointLocal, &nPointGlobal, 1, MPI_UNSIGNED_LONG, MPI_SUM, SU2_MPI::GetComm());
#else
  nPointGlobal = nPointLocal;
#endif

  Point2Vertex = new unsigned long[nPointGlobal][2];
  PointInDomain = new bool[nPointGlobal];

  for (iPoint = 0; iPoint < nPointGlobal; iPoint ++)
    PointInDomain[iPoint] = false;

  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    Boundary   = config->GetMarker_All_KindBC(iMarker);

    if ((Boundary == EULER_WALL             ) ||
        (Boundary == HEAT_FLUX              ) ||
        (Boundary == ISOTHERMAL             ) ||
        (Boundary == NEARFIELD_BOUNDARY)) {
      for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {

        /*--- The Pressure file uses the global numbering ---*/

#ifndef HAVE_MPI
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
#else
        iPoint = geometry->nodes->GetGlobalIndex(geometry->vertex[iMarker][iVertex]->GetNode());
#endif

        if (geometry->vertex[iMarker][iVertex]->GetNode() < geometry->GetnPointDomain()) {
          Point2Vertex[iPoint][0] = iMarker;
          Point2Vertex[iPoint][1] = iVertex;
          PointInDomain[iPoint] = true;
          solver_container->SetCPressureTarget(iMarker, iVertex, 0.0);
        }

      }
    }
  }

  /*--- Prepare to read the surface pressure files (CSV) ---*/

  surfCp_filename = "TargetCp";
  strcpy (cstr, surfCp_filename.c_str());

  /*--- Write file name with extension if unsteady or steady ---*/

  if (((config->GetTime_Marching() != TIME_MARCHING::STEADY) && config->GetTime_Domain()) ||
       (config->GetTime_Marching() == TIME_MARCHING::HARMONIC_BALANCE)) {
    if ((SU2_TYPE::Int(iExtIter) >= 0)    && (SU2_TYPE::Int(iExtIter) < 10))    SPRINTF (buffer, "_0000%d.dat", SU2_TYPE::Int(iExtIter));
    if ((SU2_TYPE::Int(iExtIter) >= 10)   && (SU2_TYPE::Int(iExtIter) < 100))   SPRINTF (buffer, "_000%d.dat",  SU2_TYPE::Int(iExtIter));
    if ((SU2_TYPE::Int(iExtIter) >= 100)  && (SU2_TYPE::Int(iExtIter) < 1000))  SPRINTF (buffer, "_00%d.dat",   SU2_TYPE::Int(iExtIter));
    if ((SU2_TYPE::Int(iExtIter) >= 1000) && (SU2_TYPE::Int(iExtIter) < 10000)) SPRINTF (buffer, "_0%d.dat",    SU2_TYPE::Int(iExtIter));
    if (SU2_TYPE::Int(iExtIter) >= 10000) SPRINTF (buffer, "_%d.dat", SU2_TYPE::Int(iExtIter));
  }
  else
    SPRINTF (buffer, ".dat");

  strcat (cstr, buffer);

  /*--- Read the surface pressure file ---*/

  string::size_type position;

  Surface_file.open(cstr, ios::in);

  if (!(Surface_file.fail())) {

    getline(Surface_file, text_line);

    while (getline(Surface_file, text_line)) {
      for (icommas = 0; icommas < 50; icommas++) {
        position = text_line.find( ",", 0 );
        if (position!=string::npos) text_line.erase (position,1);
      }
      stringstream  point_line(text_line);

      if (geometry->GetnDim() == 2) point_line >> iPoint >> XCoord >> YCoord >> Pressure >> PressureCoeff;
      if (geometry->GetnDim() == 3) point_line >> iPoint >> XCoord >> YCoord >> ZCoord >> Pressure >> PressureCoeff;

      if (PointInDomain[iPoint]) {

        /*--- Find the vertex for the Point and Marker ---*/

        iMarker = Point2Vertex[iPoint][0];
        iVertex = Point2Vertex[iPoint][1];

        solver_container->SetCPressureTarget(iMarker, iVertex, PressureCoeff);

      }

    }

    Surface_file.close();

  }

  /*--- Compute the pressure difference ---*/

  PressDiff = 0.0;
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    Boundary   = config->GetMarker_All_KindBC(iMarker);

    if ((Boundary == EULER_WALL             ) ||
        (Boundary == HEAT_FLUX              ) ||
        (Boundary == ISOTHERMAL             ) ||
        (Boundary == NEARFIELD_BOUNDARY)) {
      for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {

        Normal = geometry->vertex[iMarker][iVertex]->GetNormal();

        Cp = solver_container->GetCPressure(iMarker, iVertex);
        CpTarget = solver_container->GetCPressureTarget(iMarker, iVertex);

        Area = 0.0;
        for (iDim = 0; iDim < geometry->GetnDim(); iDim++)
          Area += Normal[iDim]*Normal[iDim];
        Area = sqrt(Area);

        PressDiff += Area * (CpTarget - Cp) * (CpTarget - Cp);
      }

    }
  }

#ifdef HAVE_MPI
  su2double MyPressDiff = PressDiff;   PressDiff = 0.0;
  SU2_MPI::Allreduce(&MyPressDiff, &PressDiff, 1, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());
#endif

  /*--- Update the total Cp difference coeffient ---*/

  solver_container->SetTotal_CpDiff(PressDiff);

  delete [] Point2Vertex;
  delete [] PointInDomain;

}

void COutputLegacy::SetHeatFlux_InverseDesign(CSolver *solver_container, CGeometry *geometry, CConfig *config, unsigned long iExtIter) {

  unsigned short iMarker, icommas, Boundary, iDim;
  unsigned long iVertex, iPoint, (*Point2Vertex)[2], nPointLocal = 0, nPointGlobal = 0;
  su2double XCoord, YCoord, ZCoord, PressureCoeff, HeatFlux = 0.0, HeatFluxDiff, HeatFluxTarget, *Normal = nullptr, Area,
  Pressure, Cf;
  bool *PointInDomain;
  string text_line, surfHeatFlux_filename;
  ifstream Surface_file;
  char buffer[50], cstr[200];


  nPointLocal = geometry->GetnPoint();
#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&nPointLocal, &nPointGlobal, 1, MPI_UNSIGNED_LONG, MPI_SUM, SU2_MPI::GetComm());
#else
  nPointGlobal = nPointLocal;
#endif

  Point2Vertex = new unsigned long[nPointGlobal][2];
  PointInDomain = new bool[nPointGlobal];

  for (iPoint = 0; iPoint < nPointGlobal; iPoint ++)
    PointInDomain[iPoint] = false;

  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    Boundary   = config->GetMarker_All_KindBC(iMarker);

    if ((Boundary == EULER_WALL             ) ||
        (Boundary == HEAT_FLUX              ) ||
        (Boundary == ISOTHERMAL             ) ||
        (Boundary == NEARFIELD_BOUNDARY)) {
      for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {

        /*--- The Pressure file uses the global numbering ---*/

#ifndef HAVE_MPI
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
#else
        iPoint = geometry->nodes->GetGlobalIndex(geometry->vertex[iMarker][iVertex]->GetNode());
#endif

        if (geometry->vertex[iMarker][iVertex]->GetNode() < geometry->GetnPointDomain()) {
          Point2Vertex[iPoint][0] = iMarker;
          Point2Vertex[iPoint][1] = iVertex;
          PointInDomain[iPoint] = true;
          solver_container->SetHeatFluxTarget(iMarker, iVertex, 0.0);
        }
      }
    }
  }

  /*--- Prepare to read the surface pressure files (CSV) ---*/

  surfHeatFlux_filename = "TargetHeatFlux";
  strcpy (cstr, surfHeatFlux_filename.c_str());

  /*--- Write file name with extension if unsteady or steady ---*/

  if (((config->GetTime_Marching() != TIME_MARCHING::STEADY) && config->GetTime_Domain()) ||
       (config->GetTime_Marching() == TIME_MARCHING::HARMONIC_BALANCE)) {
    if ((SU2_TYPE::Int(iExtIter) >= 0)    && (SU2_TYPE::Int(iExtIter) < 10))    SPRINTF (buffer, "_0000%d.dat", SU2_TYPE::Int(iExtIter));
    if ((SU2_TYPE::Int(iExtIter) >= 10)   && (SU2_TYPE::Int(iExtIter) < 100))   SPRINTF (buffer, "_000%d.dat",  SU2_TYPE::Int(iExtIter));
    if ((SU2_TYPE::Int(iExtIter) >= 100)  && (SU2_TYPE::Int(iExtIter) < 1000))  SPRINTF (buffer, "_00%d.dat",   SU2_TYPE::Int(iExtIter));
    if ((SU2_TYPE::Int(iExtIter) >= 1000) && (SU2_TYPE::Int(iExtIter) < 10000)) SPRINTF (buffer, "_0%d.dat",    SU2_TYPE::Int(iExtIter));
    if (SU2_TYPE::Int(iExtIter) >= 10000) SPRINTF (buffer, "_%d.dat", SU2_TYPE::Int(iExtIter));
  }
  else
    SPRINTF (buffer, ".dat");

  strcat (cstr, buffer);

  /*--- Read the surface pressure file ---*/

  string::size_type position;

  Surface_file.open(cstr, ios::in);

  if (!(Surface_file.fail())) {

    getline(Surface_file, text_line);

    while (getline(Surface_file, text_line)) {
      for (icommas = 0; icommas < 50; icommas++) {
        position = text_line.find( ",", 0 );
        if (position!=string::npos) text_line.erase (position,1);
      }
      stringstream  point_line(text_line);

      if (geometry->GetnDim() == 2) point_line >> iPoint >> XCoord >> YCoord >> Pressure >> PressureCoeff >> Cf >> HeatFlux;
      if (geometry->GetnDim() == 3) point_line >> iPoint >> XCoord >> YCoord >> ZCoord >> Pressure >> PressureCoeff >> Cf >> HeatFlux;

      if (PointInDomain[iPoint]) {

        /*--- Find the vertex for the Point and Marker ---*/

        iMarker = Point2Vertex[iPoint][0];
        iVertex = Point2Vertex[iPoint][1];

        solver_container->SetHeatFluxTarget(iMarker, iVertex, HeatFlux);

      }

    }

    Surface_file.close();
  }

  /*--- Compute the pressure difference ---*/

  HeatFluxDiff = 0.0;
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    Boundary   = config->GetMarker_All_KindBC(iMarker);

    if ((Boundary == EULER_WALL             ) ||
        (Boundary == HEAT_FLUX              ) ||
        (Boundary == ISOTHERMAL             ) ||
        (Boundary == NEARFIELD_BOUNDARY)) {
      for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {

        Normal = geometry->vertex[iMarker][iVertex]->GetNormal();

        HeatFlux = solver_container->GetHeatFlux(iMarker, iVertex);
        HeatFluxTarget = solver_container->GetHeatFluxTarget(iMarker, iVertex);

        Area = 0.0;
        for (iDim = 0; iDim < geometry->GetnDim(); iDim++)
          Area += Normal[iDim]*Normal[iDim];
        Area = sqrt(Area);

        HeatFluxDiff += Area * (HeatFluxTarget - HeatFlux) * (HeatFluxTarget - HeatFlux);

      }

    }
  }

#ifdef HAVE_MPI
  su2double MyHeatFluxDiff = HeatFluxDiff;   HeatFluxDiff = 0.0;
  SU2_MPI::Allreduce(&MyHeatFluxDiff, &HeatFluxDiff, 1, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());
#endif

  /*--- Update the total HeatFlux difference coeffient ---*/

  solver_container->SetTotal_HeatFluxDiff(HeatFluxDiff);

  delete [] Point2Vertex;
  delete [] PointInDomain;

}



void COutputLegacy::SpecialOutput_Distortion(CSolver *solver, CGeometry *geometry, CConfig *config, bool output) const {

  unsigned short iMarker, iDim, iMarker_Analyze;
  unsigned long iPoint, iVertex;
  su2double xCoord = 0.0, yCoord = 0.0, zCoord = 0.0, Area = 0.0, *Vector, TotalArea = 0.0;
  su2double xCoord_CG = 0.0, yCoord_CG = 0.0, zCoord_CG = 0.0, TipRadius, HubRadius, Distance = 0.0, Distance_Mirror = 0.0;
  su2double *r, MinDistance, xCoord_ = 0.0, yCoord_ = 0.0, zCoord_ = 0;
  unsigned short iStation, iAngle, nAngle;
  char cstr[200];
  su2double *** ProbeArray, dx = 0.0, dy = 0.0, dz = 0.0, dx_ = 0.0, dy_ = 0.0, dz_ = 0.0, UpVector[3], radians, RotatedVector[3];
  su2double Pressure, SoundSpeed, Velocity2, Mach,  Gamma, TotalPressure, Mach_Inf, TotalPressure_Inf,
  Temperature, TotalTemperature, Pressure_Inf, Temperature_Inf, TotalTemperature_Inf, Velocity_Inf, Density;
  unsigned short nDim = geometry->GetnDim();
  unsigned short Theta, nStation;
  unsigned long nVertex_Surface, nLocalVertex_Surface, MaxLocalVertex_Surface;
  unsigned long Buffer_Send_nVertex[1], *Buffer_Recv_nVertex = nullptr;
  unsigned long Total_Index;
  unsigned short Theta_DC60 = 60, nStation_DC60 = 5;
  su2double PT_Mean, Mach_Mean, q_Mean, PT, q, *PT_Sector, PT_Sector_Min, DC60, *PT_Station, *PT_Station_Min, *Mach_Station, *Mach_Station_Min, IDR, IDC, IDC_Mach;


  bool Engine_HalfModel = config->GetEngine_HalfModel();
  su2double SignFlip = 1.0;
  su2double Beta, Alpha;
  su2double Mach_ij, Mach_ip1j, Mach_im1j, Mach_ijp1, Mach_ijm1, Filtered_Mach;
  su2double Alpha_ij, Alpha_ip1j, Alpha_im1j, Alpha_ijp1, Alpha_ijm1, Filtered_Alpha;
  su2double Beta_ij, Beta_ip1j, Beta_im1j, Beta_ijp1, Beta_ijm1, Filtered_Beta;
  su2double a, b, c, d;

  int iProcessor, nProcessor;
  nProcessor = size;


  if (rank == MASTER_NODE && !config->GetDiscrete_Adjoint()) cout << endl << "Writing Surface Analysis file (surface_analysis.dat).";

  /*--- Open and rrite file name with extension if unsteady ---*/

  ofstream SurfFlow_file;

  if (output && (rank == MASTER_NODE)) {

    if ((config->GetTabular_FileFormat() == TAB_OUTPUT::TAB_CSV)) strcpy (cstr, "surface_analysis.vtk");
    else strcpy (cstr, "surface_analysis.dat");

    SurfFlow_file.precision(15);

    SurfFlow_file.open(cstr, ios::out);

    if ((config->GetTabular_FileFormat() == TAB_OUTPUT::TAB_CSV)) {
      SurfFlow_file << "# vtk DataFile Version 3.0" << endl;
      SurfFlow_file << "vtk output" << endl;
      SurfFlow_file << "ASCII" << endl;
    }
    else {
      SurfFlow_file <<"TITLE = \"Surface Analysis\"" <<endl;
      SurfFlow_file <<"VARIABLES = \"y(in)\", \"z(in)\", \"PT/PT<sub>inf</sub>\", \"TT/TT<sub>inf</sub>\", \"P/P<sub>inf</sub>\", \"T/T<sub>inf</sub>\", \"v<sub>x</sub>/v<sub>inf</sub>\", \"v<sub>y</sub>/v<sub>inf</sub>\", \"v<sub>z</sub>/v<sub>inf</sub>\", \"<greek>a</greek> (deg)\", \"<greek>b</greek> (deg)\", \"Mach\", \"Filtered <greek>a</greek> (deg)\", \"Filtered <greek>b</greek> (deg)\", \"Filtered Mach\"" << endl;
    }

  }

  /*--- Loop over all the markers to analyze ---*/

  for (iMarker_Analyze = 0; iMarker_Analyze < config->GetnMarker_Analyze(); iMarker_Analyze++) {

    string Analyze_TagBound = config->GetMarker_Analyze_TagBound(iMarker_Analyze);

    nVertex_Surface = 0; nLocalVertex_Surface = 0; MaxLocalVertex_Surface = 0;

    /*--- Find the max number of surface vertices among all
     partitions and set up buffers. The master node will handle the
     writing of the CSV file after gathering all of the data. ---*/

    nLocalVertex_Surface = 0;
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      string Marker_TagBound = config->GetMarker_All_TagBound(iMarker);
      if (Marker_TagBound == Analyze_TagBound) {
        for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          if (geometry->nodes->GetDomain(iPoint)) nLocalVertex_Surface++;
        }
      }
    }

    /*--- Communicate the number of local vertices on each partition
     to the master node ---*/

    Buffer_Send_nVertex[0] = nLocalVertex_Surface;
    if (rank == MASTER_NODE) Buffer_Recv_nVertex = new unsigned long [nProcessor];

#ifdef HAVE_MPI
    SU2_MPI::Allreduce(&nLocalVertex_Surface, &MaxLocalVertex_Surface, 1, MPI_UNSIGNED_LONG, MPI_MAX, SU2_MPI::GetComm());
    SU2_MPI::Gather(&Buffer_Send_nVertex, 1, MPI_UNSIGNED_LONG, Buffer_Recv_nVertex, 1, MPI_UNSIGNED_LONG, MASTER_NODE, SU2_MPI::GetComm());
#else
    MaxLocalVertex_Surface = nLocalVertex_Surface;
    Buffer_Recv_nVertex[MASTER_NODE] = Buffer_Send_nVertex[MASTER_NODE];
#endif

    /*--- Send and Recv buffers ---*/

    su2double *Buffer_Send_Coord_x = nullptr, *Buffer_Recv_Coord_x = nullptr;
    Buffer_Send_Coord_x = new su2double [MaxLocalVertex_Surface];

    su2double *Buffer_Send_Coord_y = nullptr, *Buffer_Recv_Coord_y = nullptr;
    Buffer_Send_Coord_y = new su2double [MaxLocalVertex_Surface];

    su2double *Buffer_Send_Coord_z = nullptr, *Buffer_Recv_Coord_z = nullptr;
    if (nDim == 3)  Buffer_Send_Coord_z = new su2double [MaxLocalVertex_Surface];

    su2double *Buffer_Send_PT = nullptr, *Buffer_Recv_PT = nullptr;
    Buffer_Send_PT = new su2double [MaxLocalVertex_Surface];

    su2double *Buffer_Send_TT = nullptr, *Buffer_Recv_TT = nullptr;
    Buffer_Send_TT = new su2double [MaxLocalVertex_Surface];

    su2double *Buffer_Send_P = nullptr, *Buffer_Recv_P = nullptr;
    Buffer_Send_P = new su2double [MaxLocalVertex_Surface];

    su2double *Buffer_Send_T = nullptr, *Buffer_Recv_T = nullptr;
    Buffer_Send_T = new su2double [MaxLocalVertex_Surface];

    su2double *Buffer_Send_Mach = nullptr, *Buffer_Recv_Mach = nullptr;
    Buffer_Send_Mach = new su2double [MaxLocalVertex_Surface];

    su2double *Buffer_Send_Vel_x = nullptr, *Buffer_Recv_Vel_x = nullptr;
    Buffer_Send_Vel_x = new su2double [MaxLocalVertex_Surface];

    su2double *Buffer_Send_Vel_y = nullptr, *Buffer_Recv_Vel_y = nullptr;
    Buffer_Send_Vel_y = new su2double [MaxLocalVertex_Surface];

    su2double *Buffer_Send_Vel_z = nullptr, *Buffer_Recv_Vel_z = nullptr;
    if (nDim == 3) Buffer_Send_Vel_z = new su2double [MaxLocalVertex_Surface];

    su2double *Buffer_Send_q = nullptr, *Buffer_Recv_q = nullptr;
    Buffer_Send_q = new su2double [MaxLocalVertex_Surface];

    su2double *Buffer_Send_Area = nullptr, *Buffer_Recv_Area = nullptr;
    Buffer_Send_Area = new su2double [MaxLocalVertex_Surface];

    /*--- Prepare the receive buffers on the master node only. ---*/

    if (rank == MASTER_NODE) {
      Buffer_Recv_Coord_x = new su2double [nProcessor*MaxLocalVertex_Surface];
      Buffer_Recv_Coord_y = new su2double [nProcessor*MaxLocalVertex_Surface];
      if (nDim == 3) Buffer_Recv_Coord_z = new su2double [nProcessor*MaxLocalVertex_Surface];
      Buffer_Recv_PT = new su2double [nProcessor*MaxLocalVertex_Surface];
      Buffer_Recv_TT = new su2double [nProcessor*MaxLocalVertex_Surface];
      Buffer_Recv_P = new su2double [nProcessor*MaxLocalVertex_Surface];
      Buffer_Recv_T = new su2double [nProcessor*MaxLocalVertex_Surface];
      Buffer_Recv_Mach = new su2double [nProcessor*MaxLocalVertex_Surface];
      Buffer_Recv_Vel_x = new su2double [nProcessor*MaxLocalVertex_Surface];
      Buffer_Recv_Vel_y = new su2double [nProcessor*MaxLocalVertex_Surface];
      if (nDim == 3) {
        Buffer_Recv_Vel_z = new su2double [nProcessor*MaxLocalVertex_Surface];
      }
      Buffer_Recv_q = new su2double [nProcessor*MaxLocalVertex_Surface];
      Buffer_Recv_Area = new su2double [nProcessor*MaxLocalVertex_Surface];
    }

    /*--- Loop over all vertices in this partition and load the
     data of the specified type into the buffer to be sent to
     the master node. ---*/

    nVertex_Surface = 0;
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      string Marker_TagBound = config->GetMarker_All_TagBound(iMarker);
      if (Marker_TagBound == Analyze_TagBound) {

        for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

          if (geometry->nodes->GetDomain(iPoint)) {

            Buffer_Send_Coord_x[nVertex_Surface] = geometry->nodes->GetCoord(iPoint, 0);
            Buffer_Send_Coord_y[nVertex_Surface] = geometry->nodes->GetCoord(iPoint, 1);
            if (nDim == 3) { Buffer_Send_Coord_z[nVertex_Surface] = geometry->nodes->GetCoord(iPoint, 2); }

            Pressure         = solver->GetNodes()->GetPressure(iPoint);
            Density          = solver->GetNodes()->GetDensity(iPoint);
            Temperature      = solver->GetNodes()->GetTemperature(iPoint);
            SoundSpeed       = solver->GetNodes()->GetSoundSpeed(iPoint);
            Velocity2        = solver->GetNodes()->GetVelocity2(iPoint);
            Mach             = sqrt(Velocity2)/SoundSpeed;
            Gamma            = config->GetGamma();

            Mach_Inf         = config->GetMach();
            Pressure_Inf     = config->GetPressure_FreeStreamND();
            Temperature_Inf  = config->GetTemperature_FreeStreamND();
            Velocity_Inf     = sqrt(config->GetVelocity_FreeStreamND()[0]*config->GetVelocity_FreeStreamND()[0]
                                    + config->GetVelocity_FreeStreamND()[1]*config->GetVelocity_FreeStreamND()[1]
                                    + config->GetVelocity_FreeStreamND()[2]*config->GetVelocity_FreeStreamND()[2]);

            Buffer_Send_P[nVertex_Surface] = Pressure / Pressure_Inf;
            Buffer_Send_T[nVertex_Surface]     = Temperature / Temperature_Inf;
            Buffer_Send_Mach[nVertex_Surface] = Mach;

            TotalPressure    = Pressure * pow( 1.0 + Mach * Mach * 0.5 * (Gamma - 1.0), Gamma / (Gamma - 1.0));
            TotalPressure_Inf  = Pressure_Inf * pow( 1.0 + Mach_Inf * Mach_Inf * 0.5 * (Gamma - 1.0), Gamma / (Gamma - 1.0));
            Buffer_Send_PT[nVertex_Surface] = TotalPressure / TotalPressure_Inf;

            TotalTemperature = Temperature * (1.0 + Mach * Mach  * 0.5 * (Gamma - 1.0));
            TotalTemperature_Inf  = Temperature_Inf * (1.0 + Mach * Mach  * 0.5 * (Gamma - 1.0));
            Buffer_Send_TT[nVertex_Surface] = TotalTemperature / TotalTemperature_Inf;

            Buffer_Send_Vel_x[nVertex_Surface] = solver->GetNodes()->GetVelocity(iPoint,0) / Velocity_Inf;
            Buffer_Send_Vel_y[nVertex_Surface] = solver->GetNodes()->GetVelocity(iPoint,1) / Velocity_Inf;
            if (nDim == 3) {
              Buffer_Send_Vel_z[nVertex_Surface] = solver->GetNodes()->GetVelocity(iPoint,2) / Velocity_Inf;
            }

            Buffer_Send_q[nVertex_Surface] = 0.5*Density*Velocity2;

            Vector = geometry->vertex[iMarker][iVertex]->GetNormal();
            Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) { Area += Vector[iDim]*Vector[iDim]; } Area = sqrt(Area);
            Buffer_Send_Area[nVertex_Surface] = Area;

            /*--- If US system, the output should be in inches ---*/

            if (config->GetSystemMeasurements() == US) {

              Buffer_Send_Coord_x[nVertex_Surface] *= 12.0;
              Buffer_Send_Coord_y[nVertex_Surface] *= 12.0;
              if (nDim == 3) Buffer_Send_Coord_z[nVertex_Surface] *= 12.0;
              Buffer_Send_Area[nVertex_Surface] *= 144.0;

            }

            nVertex_Surface++;

          }
        }
        break;
      }
    }

    /*--- Send the information to the master node ---*/

#ifdef HAVE_MPI

    SU2_MPI::Gather(Buffer_Send_Coord_x, MaxLocalVertex_Surface, MPI_DOUBLE, Buffer_Recv_Coord_x, MaxLocalVertex_Surface, MPI_DOUBLE, MASTER_NODE, SU2_MPI::GetComm());
    SU2_MPI::Gather(Buffer_Send_Coord_y, MaxLocalVertex_Surface, MPI_DOUBLE, Buffer_Recv_Coord_y, MaxLocalVertex_Surface, MPI_DOUBLE, MASTER_NODE, SU2_MPI::GetComm());
    if (nDim == 3) SU2_MPI::Gather(Buffer_Send_Coord_z, MaxLocalVertex_Surface, MPI_DOUBLE, Buffer_Recv_Coord_z, MaxLocalVertex_Surface, MPI_DOUBLE, MASTER_NODE, SU2_MPI::GetComm());
    SU2_MPI::Gather(Buffer_Send_PT, MaxLocalVertex_Surface, MPI_DOUBLE, Buffer_Recv_PT, MaxLocalVertex_Surface, MPI_DOUBLE, MASTER_NODE, SU2_MPI::GetComm());
    SU2_MPI::Gather(Buffer_Send_TT, MaxLocalVertex_Surface, MPI_DOUBLE, Buffer_Recv_TT, MaxLocalVertex_Surface, MPI_DOUBLE, MASTER_NODE, SU2_MPI::GetComm());
    SU2_MPI::Gather(Buffer_Send_P, MaxLocalVertex_Surface, MPI_DOUBLE, Buffer_Recv_P, MaxLocalVertex_Surface, MPI_DOUBLE, MASTER_NODE, SU2_MPI::GetComm());
    SU2_MPI::Gather(Buffer_Send_T, MaxLocalVertex_Surface, MPI_DOUBLE, Buffer_Recv_T, MaxLocalVertex_Surface, MPI_DOUBLE, MASTER_NODE, SU2_MPI::GetComm());
    SU2_MPI::Gather(Buffer_Send_Mach, MaxLocalVertex_Surface, MPI_DOUBLE, Buffer_Recv_Mach, MaxLocalVertex_Surface, MPI_DOUBLE, MASTER_NODE, SU2_MPI::GetComm());
    SU2_MPI::Gather(Buffer_Send_Vel_x, MaxLocalVertex_Surface, MPI_DOUBLE, Buffer_Recv_Vel_x, MaxLocalVertex_Surface, MPI_DOUBLE, MASTER_NODE, SU2_MPI::GetComm());
    SU2_MPI::Gather(Buffer_Send_Vel_y, MaxLocalVertex_Surface, MPI_DOUBLE, Buffer_Recv_Vel_y, MaxLocalVertex_Surface, MPI_DOUBLE, MASTER_NODE, SU2_MPI::GetComm());
    if (nDim == 3) SU2_MPI::Gather(Buffer_Send_Vel_z, MaxLocalVertex_Surface, MPI_DOUBLE, Buffer_Recv_Vel_z, MaxLocalVertex_Surface, MPI_DOUBLE, MASTER_NODE, SU2_MPI::GetComm());
    SU2_MPI::Gather(Buffer_Send_q, MaxLocalVertex_Surface, MPI_DOUBLE, Buffer_Recv_q, MaxLocalVertex_Surface, MPI_DOUBLE, MASTER_NODE, SU2_MPI::GetComm());
    SU2_MPI::Gather(Buffer_Send_Area, MaxLocalVertex_Surface, MPI_DOUBLE, Buffer_Recv_Area, MaxLocalVertex_Surface, MPI_DOUBLE, MASTER_NODE, SU2_MPI::GetComm());

#else

    for (iVertex = 0; iVertex < MaxLocalVertex_Surface; iVertex++) {
      Buffer_Recv_Coord_x[iVertex] = Buffer_Send_Coord_x[iVertex];
      Buffer_Recv_Coord_y[iVertex] = Buffer_Send_Coord_y[iVertex];
      if (nDim == 3) Buffer_Recv_Coord_z[iVertex] = Buffer_Send_Coord_z[iVertex];
      Buffer_Recv_PT[iVertex] = Buffer_Send_PT[iVertex];
      Buffer_Recv_TT[iVertex] = Buffer_Send_TT[iVertex];
      Buffer_Recv_P[iVertex] = Buffer_Send_P[iVertex];
      Buffer_Recv_T[iVertex] = Buffer_Send_T[iVertex];
      Buffer_Recv_Mach[iVertex] = Buffer_Send_Mach[iVertex];
      Buffer_Recv_Vel_x[iVertex] = Buffer_Send_Vel_x[iVertex];
      Buffer_Recv_Vel_y[iVertex] = Buffer_Send_Vel_y[iVertex];
      if (nDim == 3) Buffer_Recv_Vel_z[iVertex] = Buffer_Send_Vel_z[iVertex];
      Buffer_Recv_q[iVertex] = Buffer_Send_q[iVertex];
      Buffer_Recv_Area[iVertex] = Buffer_Send_Area[iVertex];
    }

#endif

    if (rank == MASTER_NODE) {

      /*--- Compute the location of the critical points of the distortion measure, and center of gravity ---*/

      TotalArea = 0.0; xCoord_CG = 0.0; yCoord_CG = 0.0; zCoord_CG = 0.0; PT_Mean = 0.0; Mach_Mean = 0.0;  q_Mean = 0.0;

      for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
        for (iVertex = 0; iVertex < Buffer_Recv_nVertex[iProcessor]; iVertex++) {

          /*--- Current index position and global index ---*/

          Total_Index = iProcessor*MaxLocalVertex_Surface+iVertex;

          /*--- Retrieve the merged data for this node ---*/

          xCoord = Buffer_Recv_Coord_x[Total_Index];
          yCoord = Buffer_Recv_Coord_y[Total_Index];
          if (nDim == 3) zCoord = Buffer_Recv_Coord_z[Total_Index];
          PT   = Buffer_Recv_PT[Total_Index];
          Mach = Buffer_Recv_Mach[Total_Index];
          q    = Buffer_Recv_q[Total_Index];

          Area       = Buffer_Recv_Area[Total_Index];
          TotalArea += Area;
          xCoord_CG += xCoord*Area;
          yCoord_CG += yCoord*Area;
          zCoord_CG += zCoord*Area;
          PT_Mean   += PT*Area;
          Mach_Mean += PT*Area;
          q_Mean    += q*Area;

        }
      }

      /*--- Evaluate the area averaged pressure and CG ---*/

      xCoord_CG = xCoord_CG / TotalArea;
      yCoord_CG = yCoord_CG / TotalArea;
      zCoord_CG = zCoord_CG / TotalArea;
      PT_Mean   /= TotalArea;
      Mach_Mean /= TotalArea;
      q_Mean    /=  TotalArea;

      /*--- If it is a half model, CGy = 0 ---*/

      if (Engine_HalfModel) { yCoord_CG = 0.0; }

      /*--- Compute hub and tip radius ---*/

      TipRadius = 1E-6; HubRadius = 1E6;
      for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
        for (iVertex = 0; iVertex < Buffer_Recv_nVertex[iProcessor]; iVertex++) {

          /*--- Current index position and global index ---*/

          Total_Index = iProcessor*MaxLocalVertex_Surface+iVertex;

          /*--- Retrieve the merged data for this node ---*/

          xCoord = Buffer_Recv_Coord_x[Total_Index];
          yCoord = Buffer_Recv_Coord_y[Total_Index];
          if (nDim == 3) zCoord = Buffer_Recv_Coord_z[Total_Index];

          if (nDim == 2)
            Distance = sqrt((xCoord_CG-xCoord)*(xCoord_CG-xCoord) +
                            (yCoord_CG-yCoord)*(yCoord_CG-yCoord));

          if (nDim == 3)
            Distance = sqrt((xCoord_CG-xCoord)*(xCoord_CG-xCoord) +
                            (yCoord_CG-yCoord)*(yCoord_CG-yCoord) +
                            (zCoord_CG-zCoord)*(zCoord_CG-zCoord));

          if (Distance > TipRadius) TipRadius = Distance;
          if (Distance < HubRadius) HubRadius = Distance;

        }
      }

      if (HubRadius/TipRadius < 0.05) HubRadius = 0.0;

      /*--- Evaluate the DC60 parameter ---*/

      Theta = Theta_DC60;
      nStation = nStation_DC60;

      nAngle = SU2_TYPE::Int(360/float(Theta));
      r = new su2double [nStation+1];

      /*--- Allocate memory ---*/

      PT_Sector = new su2double [nAngle];
      ProbeArray = new su2double ** [nAngle];
      for (iAngle = 0; iAngle < nAngle; iAngle++) {
        ProbeArray[iAngle] = new su2double * [nStation];
        for (iStation = 0; iStation < nStation; iStation++) {
          ProbeArray[iAngle][iStation] = new su2double [5];
        }
      }

      /*--- Define the radius for each probe ---*/

      r[0] = HubRadius; r[nStation] = TipRadius;
      for (iStation = 1; iStation < nStation; iStation++) {
        r[iStation] = sqrt(  r[iStation-1]*r[iStation-1] + (r[nStation]*r[nStation] - r[0]*r[0])/float(nStation) );
      }

      /*--- Define the probe rack ---*/

      UpVector[0] = 0.0; UpVector[1] = 0.0; UpVector[2] = 1.0;

      for (iAngle = 0; iAngle < nAngle; iAngle++) {

        radians = -iAngle*Theta*2.0*PI_NUMBER/360;
        RotatedVector[0] =  UpVector[0];
        RotatedVector[1] =  UpVector[1] * cos(radians) - UpVector[2] * sin(radians);
        RotatedVector[2] =  UpVector[1] * sin(radians) + UpVector[2] * cos(radians);

        for (iStation = 1; iStation <= nStation; iStation++) {
          ProbeArray[iAngle][iStation-1][0] = xCoord_CG+RotatedVector[0]*sqrt(0.5*(r[iStation]*r[iStation]+r[iStation-1]*r[iStation-1]));
          ProbeArray[iAngle][iStation-1][1] = yCoord_CG+RotatedVector[1]*sqrt(0.5*(r[iStation]*r[iStation]+r[iStation-1]*r[iStation-1]));
          ProbeArray[iAngle][iStation-1][2] = zCoord_CG+RotatedVector[2]*sqrt(0.5*(r[iStation]*r[iStation]+r[iStation-1]*r[iStation-1]));
        }

      }

      /*--- Compute the Total pressure at each probe, closes grid point to the location ---*/

      for (iAngle = 0; iAngle < nAngle; iAngle++) {

        for (iStation = 0; iStation < nStation; iStation++) {
          xCoord_ = ProbeArray[iAngle][iStation][0];
          yCoord_ = ProbeArray[iAngle][iStation][1];
          zCoord_ = ProbeArray[iAngle][iStation][2];

          MinDistance = 1E6;

          for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
            for (iVertex = 0; iVertex < Buffer_Recv_nVertex[iProcessor]; iVertex++) {

              Total_Index = iProcessor*MaxLocalVertex_Surface+iVertex;
              xCoord = Buffer_Recv_Coord_x[Total_Index];
              yCoord = Buffer_Recv_Coord_y[Total_Index];
              if (nDim == 3) zCoord = Buffer_Recv_Coord_z[Total_Index];

              dx = (xCoord_ - xCoord); dy = (yCoord_ - yCoord);
              if (nDim == 3) dz = (zCoord_ - zCoord);

              Distance = dx*dx + dy*dy; if (nDim == 3) Distance += dz*dz; Distance = sqrt(Distance);

              if (Engine_HalfModel) {

                yCoord = -yCoord;

                dx_ = (xCoord_ - xCoord); dy_ = (yCoord_ - yCoord);
                if (nDim == 3) dz_ = (zCoord_ - zCoord);

                Distance_Mirror = dx_*dx_ + dy_*dy_;
                if (nDim == 3) Distance_Mirror += dz_*dz_;
                Distance_Mirror = sqrt(Distance_Mirror);

                if (Distance_Mirror < Distance) {
                  Distance = Distance_Mirror;
                  dx = dx_; dy = dy_;
                  if (nDim == 3) dz = dz_;
                }

              }

              if (Distance <= MinDistance) {
                MinDistance = Distance;
                ProbeArray[iAngle][iStation][3] = Buffer_Recv_PT[Total_Index];
                ProbeArray[iAngle][iStation][4] = Buffer_Recv_q[Total_Index];
              }

            }
          }

        }

      }

      /*--- Evaluate the average pressure at each sector, fan face and dynamic pressure ---*/

      PT_Mean = 0.0; q_Mean = 0.0;
      for (iAngle = 0; iAngle < nAngle; iAngle++) {
        PT_Sector[iAngle] = 0.0;
        for (iStation = 0; iStation < nStation; iStation++) {
          PT_Sector[iAngle] += ProbeArray[iAngle][iStation][3]/float(nStation);
          PT_Mean           += ProbeArray[iAngle][iStation][3]/float(nStation*nAngle);
          q_Mean            += ProbeArray[iAngle][iStation][4]/float(nStation*nAngle);
        }
      }

      /*--- Compute the min value of the averaged pressure at each sector ---*/

      PT_Sector_Min = PT_Sector[0];
      for (iAngle = 1; iAngle < nAngle; iAngle++) {
        if (PT_Sector[iAngle] <= PT_Sector_Min) PT_Sector_Min = PT_Sector[iAngle];
      }

      /*--- Set the value of the distortion, it only works for one surface ---*/

      Mach_Inf           = config->GetMach();
      Gamma              = config->GetGamma();
      TotalPressure_Inf  = config->GetPressure_FreeStreamND() * pow( 1.0 + Mach_Inf * Mach_Inf *
                                                                    0.5 * (Gamma - 1.0), Gamma    / (Gamma - 1.0));

      if (q_Mean != 0.0) DC60 = ((PT_Mean - PT_Sector_Min)*TotalPressure_Inf)/q_Mean;
      else DC60 = 0.0;

      config->SetSurface_DC60(iMarker_Analyze, DC60);

      solver->SetTotal_DC60(DC60);

      /*--- Deallocate the memory ---*/

      delete[] r;

      delete [] PT_Sector;

      for (iAngle = 0; iAngle < nAngle; iAngle++) {
        for (iStation = 0; iStation < nStation; iStation++) {
          delete[] ProbeArray[iAngle][iStation];
        }
      }
      delete[] ProbeArray;


      /*--- Evaluate the IDC, and IDR parameters ---*/

      nStation = SU2_TYPE::Int(config->GetDistortionRack()[0]);
      Theta = SU2_TYPE::Int(config->GetDistortionRack()[1]);
      nAngle = SU2_TYPE::Int(360/float(Theta));

      /*--- Allocate memory ---*/

      r = new su2double [nStation+1];
      ProbeArray = new su2double ** [nAngle];
      for (iAngle = 0; iAngle < nAngle; iAngle++) {
        ProbeArray[iAngle] = new su2double * [nStation];
        for (iStation = 0; iStation < nStation; iStation++) {
          ProbeArray[iAngle][iStation] = new su2double [4];
        }
      }

      /*--- Define the radius for each probe ---*/

      r[0] = HubRadius; r[nStation] = TipRadius;
      for (iStation = 1; iStation < nStation; iStation++) {
        r[iStation] = sqrt(  r[iStation-1]*r[iStation-1] + (r[nStation]*r[nStation] - r[0]*r[0])/float(nStation) );
      }

      /*--- Define the probe rack ---*/

      UpVector[0] = 0.0; UpVector[1] = 0.0; UpVector[2] = 1.0;

      for (iAngle = 0; iAngle < nAngle; iAngle++) {

        radians = -iAngle*Theta*2.0*PI_NUMBER/360;
        RotatedVector[0] =  UpVector[0];
        RotatedVector[1] =  UpVector[1] * cos(radians) - UpVector[2] * sin(radians);
        RotatedVector[2] =  UpVector[1] * sin(radians) + UpVector[2] * cos(radians);

        for (iStation = 1; iStation <= nStation; iStation++) {
          ProbeArray[iAngle][iStation-1][0] = xCoord_CG+RotatedVector[0]*sqrt(0.5*(r[iStation]*r[iStation]+r[iStation-1]*r[iStation-1]));
          ProbeArray[iAngle][iStation-1][1] = yCoord_CG+RotatedVector[1]*sqrt(0.5*(r[iStation]*r[iStation]+r[iStation-1]*r[iStation-1]));
          ProbeArray[iAngle][iStation-1][2] = zCoord_CG+RotatedVector[2]*sqrt(0.5*(r[iStation]*r[iStation]+r[iStation-1]*r[iStation-1]));
        }

      }

      /*--- Compute the Total pressure at each probe, closes grid point to the location ---*/

      for (iAngle = 0; iAngle < nAngle; iAngle++) {
        for (iStation = 0; iStation < nStation; iStation++) {
          xCoord_ = ProbeArray[iAngle][iStation][0];
          yCoord_ = ProbeArray[iAngle][iStation][1];
          zCoord_ = ProbeArray[iAngle][iStation][2];

          MinDistance = 1E6;

          for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
            for (iVertex = 0; iVertex < Buffer_Recv_nVertex[iProcessor]; iVertex++) {

              Total_Index = iProcessor*MaxLocalVertex_Surface+iVertex;
              xCoord = Buffer_Recv_Coord_x[Total_Index];
              yCoord = Buffer_Recv_Coord_y[Total_Index];
              if (nDim == 3) zCoord = Buffer_Recv_Coord_z[Total_Index];

              dx = (xCoord_ - xCoord); dy = (yCoord_ - yCoord);
              if (nDim == 3) dz = (zCoord_ - zCoord);

              Distance = dx*dx + dy*dy; if (nDim == 3) Distance += dz*dz; Distance = sqrt(Distance);

              if (Engine_HalfModel) {

                yCoord = -yCoord;

                dx_ = (xCoord_ - xCoord); dy_ = (yCoord_ - yCoord);
                if (nDim == 3) dz_ = (zCoord_ - zCoord);

                Distance_Mirror = dx_*dx_ + dy_*dy_;
                if (nDim == 3) Distance_Mirror += dz_*dz_;
                Distance_Mirror = sqrt(Distance_Mirror);

                if (Distance_Mirror < Distance) {
                  Distance = Distance_Mirror;
                  dx = dx_; dy = dy_;
                  if (nDim == 3) dz = dz_;
                }

              }

              if (Distance <= MinDistance) {
                MinDistance = Distance;
                ProbeArray[iAngle][iStation][3] = Buffer_Recv_PT[Total_Index];
              }

            }
          }

        }

      }

      /*--- Evaluate the average and min. pressure at each station/radius and fan  ---*/

      PT_Station = new su2double [nStation];
      PT_Station_Min = new su2double [nStation];

      PT_Mean = 0.0;
      for (iStation = 0; iStation < nStation; iStation++) {
        PT_Station[iStation] = 0.0;
        PT_Station_Min[iStation] = ProbeArray[0][iStation][3];
        for (iAngle = 0; iAngle < nAngle; iAngle++) {
          PT = ProbeArray[iAngle][iStation][3];
          PT_Station[iStation] += PT / float(nAngle);
          if (PT <= PT_Station_Min[iStation] ) PT_Station_Min[iStation] = PT;
          PT_Mean += ProbeArray[iAngle][iStation][3]/float(nStation*nAngle);
        }
      }

      /*--- Set the value of the distortion, it only works for one surface ---*/

      IDC = 0.0;
      for (iStation = 0; iStation < nStation-1; iStation++) {
        IDC = max (IDC, 0.5*((PT_Station[iStation] - PT_Station_Min[iStation])/PT_Mean
                             + (PT_Station[iStation+1] - PT_Station_Min[iStation+1])/PT_Mean) );

      }

      config->SetSurface_IDC(iMarker_Analyze, IDC);
      solver->SetTotal_IDC(IDC);

      IDR = 0.0;
      for (iStation = 0; iStation < nStation; iStation++) {
        IDR = max (IDR, (PT_Mean-PT_Station[iStation])/PT_Mean);
      }

      config->SetSurface_IDR(iMarker_Analyze, IDR);

      solver->SetTotal_IDR(IDR);

      /*--- Release IDX parameters ---*/

      delete [] PT_Station_Min;
      delete [] PT_Station;

      /*--- Evaluate the IDC Mach parameter ---*/

      /*--- Compute the Mach number at each probe, closes grid point to the location ---*/

      for (iAngle = 0; iAngle < nAngle; iAngle++) {
        for (iStation = 0; iStation < nStation; iStation++) {
          xCoord_ = ProbeArray[iAngle][iStation][0];
          yCoord_ = ProbeArray[iAngle][iStation][1];
          zCoord_ = ProbeArray[iAngle][iStation][2];

          MinDistance = 1E6;

          for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
            for (iVertex = 0; iVertex < Buffer_Recv_nVertex[iProcessor]; iVertex++) {

              Total_Index = iProcessor*MaxLocalVertex_Surface+iVertex;

              xCoord = Buffer_Recv_Coord_x[Total_Index];
              yCoord = Buffer_Recv_Coord_y[Total_Index];
              if (nDim == 3) zCoord = Buffer_Recv_Coord_z[Total_Index];

              dx = (xCoord_ - xCoord); dy = (yCoord_ - yCoord);
              if (nDim == 3) dz = (zCoord_ - zCoord);

              Distance = dx*dx + dy*dy;
              if (nDim == 3) Distance += dz*dz;
              Distance = sqrt(Distance);

              if (Engine_HalfModel) {

                yCoord = -yCoord;

                dx_ = (xCoord_ - xCoord); dy_ = (yCoord_ - yCoord);
                if (nDim == 3) dz_ = (zCoord_ - zCoord);

                Distance_Mirror = dx_*dx_ + dy_*dy_;
                if (nDim == 3) Distance_Mirror += dz_*dz_;
                Distance_Mirror = sqrt(Distance_Mirror);

                if (Distance_Mirror < Distance) {
                  Distance = Distance_Mirror;
                  dx = dx_; dy = dy_;
                  if (nDim == 3) dz = dz_;
                }

              }

              if (Distance <= MinDistance) {
                MinDistance = Distance;
                ProbeArray[iAngle][iStation][3] = Buffer_Recv_Mach[Total_Index];
              }

            }
          }

        }

      }

      /*--- Evaluate the average and min. pressure at each station/radius and fan face ---*/

      Mach_Station = new su2double [nStation];
      Mach_Station_Min = new su2double [nStation];

      Mach_Mean = 0.0;
      for (iStation = 0; iStation < nStation; iStation++) {
        Mach_Station[iStation] = 0.0;
        Mach_Station_Min[iStation] = ProbeArray[0][iStation][3];
        for (iAngle = 0; iAngle < nAngle; iAngle++) {
          Mach = ProbeArray[iAngle][iStation][3];
          Mach_Station[iStation] += Mach / float(nAngle);
          if (Mach <= Mach_Station_Min[iStation] ) Mach_Station_Min[iStation] = Mach;
          Mach_Mean += ProbeArray[iAngle][iStation][3]/float(nStation*nAngle);
        }
      }

      /*--- Set the value of the distortion, it only works for one surface ---*/

      IDC_Mach = 0.0;
      for (iStation = 0; iStation < nStation-1; iStation++) {
        if (Mach_Mean != 0)
          IDC_Mach = max (IDC_Mach, 0.5*((Mach_Station[iStation] - Mach_Station_Min[iStation])/Mach_Mean
                                         + (Mach_Station[iStation+1] - Mach_Station_Min[iStation+1])/Mach_Mean)   );

      }

      config->SetSurface_IDC_Mach(iMarker_Analyze, IDC_Mach);

      solver->SetTotal_IDC_Mach(IDC_Mach);

      delete [] Mach_Station_Min;
      delete [] Mach_Station;

      /*--- Release distortion parameters ---*/

      delete[] r;
      for (iAngle = 0; iAngle < nAngle; iAngle++) {
        for (iStation = 0; iStation < nStation; iStation++) {
          delete[] ProbeArray[iAngle][iStation];
        }
      }
      delete[] ProbeArray;

      /*--- Create the distortion plot ---*/

      Theta = 10; nStation = 20;

      nAngle = SU2_TYPE::Int(360/float(Theta));
      r = new su2double [nStation+1];

      /*--- Allocate memory ---*/

      ProbeArray = new su2double ** [nAngle];
      for (iAngle = 0; iAngle < nAngle; iAngle++) {
        ProbeArray[iAngle] = new su2double * [nStation];
        for (iStation = 0; iStation < nStation; iStation++) {
          ProbeArray[iAngle][iStation] = new su2double [11];
        }
      }

      /*--- Define the radius for each probe ---*/

      r[0] = HubRadius;
      r[nStation] = TipRadius;

      for (iStation = 1; iStation < nStation; iStation++) {
        r[iStation] = sqrt(  r[iStation-1]*r[iStation-1] + (r[nStation]*r[nStation] - r[0]*r[0])/float(nStation) );
      }

      /*--- Define the probe rack ---*/

      UpVector[0] = 0.0; UpVector[1] = 0.0; UpVector[2] = 1.0;

      for (iAngle = 0; iAngle < nAngle; iAngle++) {

        radians = -iAngle*Theta*2.0*PI_NUMBER/360;
        RotatedVector[0] =  UpVector[0];
        RotatedVector[1] =  UpVector[1] * cos(radians) - UpVector[2] * sin(radians);
        RotatedVector[2] =  UpVector[1] * sin(radians) + UpVector[2] * cos(radians);

        for (iStation = 1; iStation <= nStation; iStation++) {
          ProbeArray[iAngle][iStation-1][0] = xCoord_CG+RotatedVector[0]*sqrt(0.5*(r[iStation]*r[iStation]+r[iStation-1]*r[iStation-1]));
          ProbeArray[iAngle][iStation-1][1] = yCoord_CG+RotatedVector[1]*sqrt(0.5*(r[iStation]*r[iStation]+r[iStation-1]*r[iStation-1]));
          ProbeArray[iAngle][iStation-1][2] = zCoord_CG+RotatedVector[2]*sqrt(0.5*(r[iStation]*r[iStation]+r[iStation-1]*r[iStation-1]));
        }

      }

      /*--- Compute the primitieve variables, closest grid point to the location + gradient ---*/

      for (iAngle = 0; iAngle < nAngle; iAngle++) {
        for (iStation = 0; iStation < nStation; iStation++) {
          xCoord_ = ProbeArray[iAngle][iStation][0];
          yCoord_ = ProbeArray[iAngle][iStation][1];
          zCoord_ = ProbeArray[iAngle][iStation][2];

          MinDistance = 1E6;

          for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
            for (iVertex = 0; iVertex < Buffer_Recv_nVertex[iProcessor]; iVertex++) {

              Total_Index = iProcessor*MaxLocalVertex_Surface+iVertex;
              xCoord = Buffer_Recv_Coord_x[Total_Index];
              yCoord = Buffer_Recv_Coord_y[Total_Index];
              if (nDim == 3) zCoord = Buffer_Recv_Coord_z[Total_Index];

              dx = (xCoord_ - xCoord); dy = (yCoord_ - yCoord);
              if (nDim == 3) dz = (zCoord_ - zCoord);

              Distance = dx*dx + dy*dy; if (nDim == 3) Distance += dz*dz; Distance = sqrt(Distance);

              SignFlip = 1.0;

              if (Engine_HalfModel) {

                yCoord = -yCoord;

                dx_ = (xCoord_ - xCoord);
                dy_ = (yCoord_ - yCoord);
                if (nDim == 3) dz_ = (zCoord_ - zCoord);

                Distance_Mirror = dx_*dx_ + dy_*dy_;
                if (nDim == 3) Distance_Mirror += dz_*dz_;
                Distance_Mirror = sqrt(Distance_Mirror);

                if (Distance_Mirror < Distance) {
                  SignFlip = -1.0;
                  Distance = Distance_Mirror;
                  dx = dx_; dy = dy_;
                  if (nDim == 3) dz = dz_;
                }

              }


              if (Distance <= MinDistance) {
                MinDistance = Distance;
                ProbeArray[iAngle][iStation][3] = Buffer_Recv_PT[Total_Index];
                ProbeArray[iAngle][iStation][4] = Buffer_Recv_TT[Total_Index];
                ProbeArray[iAngle][iStation][5] = Buffer_Recv_P[Total_Index];
                ProbeArray[iAngle][iStation][6] = Buffer_Recv_T[Total_Index];
                ProbeArray[iAngle][iStation][7] = Buffer_Recv_Mach[Total_Index];
                ProbeArray[iAngle][iStation][8] = Buffer_Recv_Vel_x[Total_Index];
                ProbeArray[iAngle][iStation][9] =  SignFlip * Buffer_Recv_Vel_y[Total_Index];
                if (nDim == 3) ProbeArray[iAngle][iStation][10] = Buffer_Recv_Vel_z[Total_Index];
              }

            }
          }

        }

      }

      /*--- Reverse in the Y direction to move the solution from 3D to 2D ---*/

      yCoord_CG = -yCoord_CG;
      for (iAngle = 0; iAngle < nAngle; iAngle++) {
        for (iStation = 0; iStation < nStation; iStation++) {
          ProbeArray[iAngle][iStation][9] = -ProbeArray[iAngle][iStation][9];
          ProbeArray[iAngle][iStation][1] = -ProbeArray[iAngle][iStation][1];
        }
      }

      if (output) {

        if (config->GetTabular_FileFormat() == TAB_OUTPUT::TAB_CSV) {

          SurfFlow_file << "\nDATASET UNSTRUCTURED_GRID" << endl;
          SurfFlow_file <<"POINTS " << nAngle*nStation << " float" << endl;
          for (iAngle = 0; iAngle < nAngle; iAngle++) {
            for (iStation = 0; iStation < nStation; iStation++) {
              SurfFlow_file << ProbeArray[iAngle][iStation][1]-yCoord_CG << " " << ProbeArray[iAngle][iStation][2]-zCoord_CG << " 0.0 " <<" ";
            }
          }

          SurfFlow_file <<"\nCELLS " << nAngle*(nStation-1) <<" "<< nAngle*(nStation-1)*5 << endl;
          for (iAngle = 0; iAngle < nAngle; iAngle++) {
            for (iStation = 0; iStation < nStation-1; iStation++) {
              a = iAngle*nStation+iStation; b = a + nStation; c = b+1; d = a +1;
              if (iAngle == nAngle-1) { b = iStation; c = b+1;   }
              SurfFlow_file << "4 " << a  <<" "<< b <<" "<< c <<" "<< d <<" ";
            }
          }

          SurfFlow_file <<"\nCELL_TYPES " << nAngle*(nStation-1) << endl;
          for (iAngle = 0; iAngle < nAngle; iAngle++) {
            for (iStation = 0; iStation < nStation-1; iStation++) {
              SurfFlow_file << "9 " ;
            }
          }

          SurfFlow_file <<"\nPOINT_DATA " << nAngle*nStation << endl;
          SurfFlow_file <<"SCALARS PT/PT_inf float" << endl;
          SurfFlow_file <<"LOOKUP_TABLE default" << endl;

          for (iAngle = 0; iAngle < nAngle; iAngle++) {
            for (iStation = 0; iStation < nStation; iStation++) {
              SurfFlow_file << ProbeArray[iAngle][iStation][3] << " ";
            }
          }

          SurfFlow_file <<"SCALARS TT/TT_inf float" << endl;
          SurfFlow_file <<"LOOKUP_TABLE default" << endl;

          for (iAngle = 0; iAngle < nAngle; iAngle++) {
            for (iStation = 0; iStation < nStation; iStation++) {
              SurfFlow_file << ProbeArray[iAngle][iStation][4] << " ";
            }
          }

          SurfFlow_file <<"SCALARS Alpha float" << endl;
          SurfFlow_file <<"LOOKUP_TABLE default" << endl;

          for (iAngle = 0; iAngle < nAngle; iAngle++) {
            for (iStation = 0; iStation < nStation; iStation++) {
              Alpha = atan(ProbeArray[iAngle][iStation][10]/ProbeArray[iAngle][iStation][8])*360.0/(2.0*PI_NUMBER);
              SurfFlow_file << Alpha << " ";
            }
          }

          SurfFlow_file <<"SCALARS Beta float" << endl;
          SurfFlow_file <<"LOOKUP_TABLE default" << endl;

          for (iAngle = 0; iAngle < nAngle; iAngle++) {
            for (iStation = 0; iStation < nStation; iStation++) {
              Beta = atan(ProbeArray[iAngle][iStation][9]/ProbeArray[iAngle][iStation][8])*360.0/(2.0*PI_NUMBER);
              SurfFlow_file << Beta << " ";
            }
          }

          SurfFlow_file <<"SCALARS Mach float" << endl;
          SurfFlow_file <<"LOOKUP_TABLE default" << endl;

          for (iAngle = 0; iAngle < nAngle; iAngle++) {
            for (iStation = 0; iStation < nStation; iStation++) {
              SurfFlow_file << ProbeArray[iAngle][iStation][7] << " ";
            }
          }

          SurfFlow_file <<"VECTORS Velocity float" << endl;

          for (iAngle = 0; iAngle < nAngle; iAngle++) {
            for (iStation = 0; iStation < nStation; iStation++) {
              SurfFlow_file << ProbeArray[iAngle][iStation][8] << " " << ProbeArray[iAngle][iStation][9] << " " << ProbeArray[iAngle][iStation][10] << " ";
            }
          }

        }
        else {

          SurfFlow_file <<"ZONE T= \"" << Analyze_TagBound <<"\", NODES=" << nAngle*nStation << " , ELEMENTS= " << nAngle*(nStation-1) <<", DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL" << endl;

          for (iAngle = 0; iAngle < nAngle; iAngle++) {
            for (iStation = 0; iStation < nStation; iStation++) {

              Alpha = atan(ProbeArray[iAngle][iStation][10]/ProbeArray[iAngle][iStation][8])*360.0/(2.0*PI_NUMBER);
              Beta = atan(ProbeArray[iAngle][iStation][9]/ProbeArray[iAngle][iStation][8])*360.0/(2.0*PI_NUMBER);

              Mach_ij = ProbeArray[iAngle][iStation][7];
              if (iAngle+1 != nAngle) Mach_ip1j = ProbeArray[iAngle+1][iStation][7];
              else Mach_ip1j = ProbeArray[0][iStation][7];
              if (iAngle-1 != -1) Mach_im1j = ProbeArray[iAngle-1][iStation][7];
              else Mach_im1j = ProbeArray[nAngle-1][iStation][7];
              if (iStation+1 != nStation) Mach_ijp1 = ProbeArray[iAngle][iStation+1][7];
              else Mach_ijp1 = ProbeArray[iAngle][0][7];
              if (iStation-1 != -1) Mach_ijm1 = ProbeArray[iAngle][iStation-1][7];
              else Mach_ijm1 = ProbeArray[iAngle][nStation-1][7];
              Filtered_Mach = (4.0*Mach_ij+Mach_ip1j+Mach_im1j+Mach_ijp1+Mach_ijm1)/8.0;

              Alpha_ij = atan(ProbeArray[iAngle][iStation][10]/ProbeArray[iAngle][iStation][8])*360.0/(2.0*PI_NUMBER);
              if (iAngle+1 != nAngle) Alpha_ip1j = atan(ProbeArray[iAngle+1][iStation][10]/ProbeArray[iAngle+1][iStation][8])*360.0/(2.0*PI_NUMBER);
              else Alpha_ip1j = atan(ProbeArray[0][iStation][10]/ProbeArray[0][iStation][8])*360.0/(2.0*PI_NUMBER);
              if (iAngle-1 != -1) Alpha_im1j = atan(ProbeArray[iAngle-1][iStation][10]/ProbeArray[iAngle-1][iStation][8])*360.0/(2.0*PI_NUMBER);
              else Alpha_im1j = atan(ProbeArray[nAngle-1][iStation][10]/ProbeArray[nAngle-1][iStation][8])*360.0/(2.0*PI_NUMBER);
              if (iStation+1 != nStation) Alpha_ijp1 = atan(ProbeArray[iAngle][iStation+1][10]/ProbeArray[iAngle][iStation+1][8])*360.0/(2.0*PI_NUMBER);
              else Alpha_ijp1 = atan(ProbeArray[iAngle][0][10]/ProbeArray[iAngle][0][8])*360.0/(2.0*PI_NUMBER);
              if (iStation-1 != -1) Alpha_ijm1 = atan(ProbeArray[iAngle][iStation-1][10]/ProbeArray[iAngle][iStation-1][8])*360.0/(2.0*PI_NUMBER);
              else Alpha_ijm1 = atan(ProbeArray[iAngle][nStation-1][10]/ProbeArray[iAngle][nStation-1][8])*360.0/(2.0*PI_NUMBER);
              Filtered_Alpha = (4.0*Alpha_ij+Alpha_ip1j+Alpha_im1j+Alpha_ijp1+Alpha_ijm1)/8.0;

              Beta_ij = atan(ProbeArray[iAngle][iStation][9]/ProbeArray[iAngle][iStation][8])*360.0/(2.0*PI_NUMBER);
              if (iAngle+1 != nAngle) Beta_ip1j = atan(ProbeArray[iAngle+1][iStation][9]/ProbeArray[iAngle+1][iStation][8])*360.0/(2.0*PI_NUMBER);
              else Beta_ip1j = atan(ProbeArray[0][iStation][9]/ProbeArray[0][iStation][8])*360.0/(2.0*PI_NUMBER);
              if (iAngle-1 != -1) Beta_im1j = atan(ProbeArray[iAngle-1][iStation][9]/ProbeArray[iAngle-1][iStation][8])*360.0/(2.0*PI_NUMBER);
              else Beta_im1j = atan(ProbeArray[nAngle-1][iStation][9]/ProbeArray[nAngle-1][iStation][8])*360.0/(2.0*PI_NUMBER);
              if (iStation+1 != nStation) Beta_ijp1 = atan(ProbeArray[iAngle][iStation+1][9]/ProbeArray[iAngle][iStation+1][8])*360.0/(2.0*PI_NUMBER);
              else Beta_ijp1 = atan(ProbeArray[iAngle][0][9]/ProbeArray[iAngle][0][8])*360.0/(2.0*PI_NUMBER);
              if (iStation-1 != -1) Beta_ijm1 = atan(ProbeArray[iAngle][iStation-1][9]/ProbeArray[iAngle][iStation-1][8])*360.0/(2.0*PI_NUMBER);
              else Beta_ijm1 = atan(ProbeArray[iAngle][nStation-1][9]/ProbeArray[iAngle][nStation-1][8])*360.0/(2.0*PI_NUMBER);
              Filtered_Beta = (4.0*Beta_ij+Beta_ip1j+Beta_im1j+Beta_ijp1+Beta_ijm1)/8.0;


              SurfFlow_file
              << " "  << ProbeArray[iAngle][iStation][1]-yCoord_CG
              <<" " << ProbeArray[iAngle][iStation][2]-zCoord_CG
              <<" " << ProbeArray[iAngle][iStation][3] <<" " << ProbeArray[iAngle][iStation][4]
              <<" " << ProbeArray[iAngle][iStation][5] <<" " << ProbeArray[iAngle][iStation][6]
              <<" " << ProbeArray[iAngle][iStation][8] <<" " << ProbeArray[iAngle][iStation][9]
              <<" " << ProbeArray[iAngle][iStation][10]
              <<" " << Alpha <<" " << Beta << " " << ProbeArray[iAngle][iStation][7]
              <<" " << Filtered_Alpha <<" " << Filtered_Beta << " " << Filtered_Mach << endl;

            }
          }

          for (iAngle = 0; iAngle < nAngle; iAngle++) {
            for (iStation = 0; iStation < nStation-1; iStation++) {
              a = iAngle*nStation+iStation; b = a + nStation; c = b+1; d = a +1;
              if (iAngle == nAngle-1) { b = iStation; c = b+1;   }
              SurfFlow_file << a+1  <<" "<< b+1  <<" "<< c+1 <<" "<< d+1 << endl;
            }
          }

          /*--- Add extra info ---*/

          SurfFlow_file << "TEXT X=14, Y=86, F=HELV-BOLD, C=BLUE, H=2.0, ";
          unsigned short RackProbes = SU2_TYPE::Int(config->GetDistortionRack()[0]);
          unsigned short RackAngle = SU2_TYPE::Int(config->GetDistortionRack()[1]);
          SurfFlow_file << "T=\"Rack Size: " << RackProbes << " probes at "<< RackAngle << "deg." << "\\" << "\\n";
          SurfFlow_file << "Mach " << config->GetMach() << ", Reynolds " << config->GetReynolds() << ", <greek>a</greek> "
          << config->GetAoA() << "deg, <greek>b</greek> " << config->GetAoS() << "deg." << "\\" << "\\n";
          SurfFlow_file.precision(1);
          SurfFlow_file << fixed << "Net Thrust " << solver->GetTotal_NetThrust() << "lbs, Power " << solver->GetTotal_Power() <<  "HP";
          SurfFlow_file.precision(4);
          SurfFlow_file << ", MassFlow " << config->GetSurface_MassFlow(iMarker_Analyze) << ",\\" << "\\n";
          SurfFlow_file << "IDC " << config->GetSurface_IDC(iMarker_Analyze)*100 << "%, IDCM " << config->GetSurface_IDC_Mach(iMarker_Analyze)*100 << "%, IDR " << config->GetSurface_IDR(iMarker_Analyze)*100 << "%,\\" << "\\n";
          SurfFlow_file << "DC60 " << config->GetSurface_DC60(iMarker_Analyze) << ".\"" << endl;

        }

      }

      /*--- Release the recv buffers on the master node ---*/

      delete [] Buffer_Recv_Coord_x;
      delete [] Buffer_Recv_Coord_y;
      if (nDim == 3) delete [] Buffer_Recv_Coord_z;

      delete [] Buffer_Recv_PT;
      delete [] Buffer_Recv_TT;
      delete [] Buffer_Recv_P;
      delete [] Buffer_Recv_T;
      delete [] Buffer_Recv_Mach;
      delete [] Buffer_Recv_Vel_x;
      delete [] Buffer_Recv_Vel_y;
      if (nDim == 3) delete [] Buffer_Recv_Vel_z;
      delete [] Buffer_Recv_q;

      delete [] Buffer_Recv_Area;

      delete [] Buffer_Recv_nVertex;

      delete[] r;
      for (iAngle = 0; iAngle < nAngle; iAngle++) {
        for (iStation = 0; iStation < nStation; iStation++) {
          delete[] ProbeArray[iAngle][iStation];
        }
      }
      delete[] ProbeArray;

    }

//    if ((rank == MASTER_NODE) && !config->GetDiscrete_Adjoint()) {
//
//      cout << "Surface ("<< Analyze_TagBound << "): ";
//      cout.precision(4);
//      cout.setf(ios::fixed, ios::floatfield);
//      cout << setprecision(1) << "IDC " << 100*config->GetSurface_IDC(iMarker_Analyze)
//      << "%. IDC Mach " << 100*config->GetSurface_IDC_Mach(iMarker_Analyze)
//      << "%. IDR " << 100*config->GetSurface_IDR(iMarker_Analyze)
//      << "%. DC60 " << config->GetSurface_DC60(iMarker_Analyze) << "." << endl;
//
//    }

    /*--- Release the memory for the remaining buffers and exit ---*/

    delete [] Buffer_Send_Coord_x;
    delete [] Buffer_Send_Coord_y;
    if (nDim == 3) delete [] Buffer_Send_Coord_z;

    delete [] Buffer_Send_PT;
    delete [] Buffer_Send_TT;
    delete [] Buffer_Send_P;
    delete [] Buffer_Send_T;
    delete [] Buffer_Send_Mach;
    delete [] Buffer_Send_Vel_x;
    delete [] Buffer_Send_Vel_y;
    if (nDim == 3) delete [] Buffer_Send_Vel_z;
    delete [] Buffer_Send_q;

    delete [] Buffer_Send_Area;

  }

  /*--- Close the tecplot  file ---*/

  if (output) {
    SurfFlow_file.close();
  }

}

void COutputLegacy::WriteTurboPerfConvHistory(CConfig *config){

  unsigned short iMarker_Monitoring;
  string inMarker_Tag, outMarker_Tag, inMarkerTag_Mix;
  unsigned short nZone       = config->GetnZone();
  bool turbulent = ((config->GetKind_Solver() == MAIN_SOLVER::RANS) || (config->GetKind_Solver() == MAIN_SOLVER::DISC_ADJ_RANS));
  bool menter_sst = (config->GetKind_Turb_Model() == TURB_MODEL::SST) || (config->GetKind_Turb_Model() == TURB_MODEL::SST_SUST);

  unsigned short nBladesRow, nStages;
  unsigned short iStage;
  nBladesRow = config->GetnMarker_Turbomachinery();
  nStages    = SU2_TYPE::Int(nBladesRow/2);

  cout << endl << "------------------------- Turbomachinery Summary ------------------------" << endl;
  cout << endl;
  for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Turbomachinery(); iMarker_Monitoring++){
    cout << endl << "----------------------------- Blade " << iMarker_Monitoring + 1 << " -----------------------------------" << endl;
    inMarker_Tag = config->GetMarker_TurboPerf_BoundIn(iMarker_Monitoring);
    outMarker_Tag = config->GetMarker_TurboPerf_BoundOut(iMarker_Monitoring);
    if(iMarker_Monitoring == 0){
      cout << "BC Inlet convergence monitoring marker " << inMarker_Tag << " : "<<endl;
      cout << endl;
      cout << "     Inlet Total Enthalpy" << "     Inlet Total Enthalpy BC" << "     err(%)" <<  endl;
      cout.width(25); cout << TotalEnthalpyIn[iMarker_Monitoring][config->GetnSpan_iZones(iMarker_Monitoring)]*config->GetEnergy_Ref();
      cout.width(25); cout << TotalEnthalpyIn_BC[iMarker_Monitoring][config->GetnSpan_iZones(iMarker_Monitoring)]*config->GetEnergy_Ref();
      cout.width(25); cout << abs((TotalEnthalpyIn[iMarker_Monitoring][config->GetnSpan_iZones(iMarker_Monitoring)] - TotalEnthalpyIn_BC[iMarker_Monitoring][config->GetnSpan_iZones(iMarker_Monitoring)])/TotalEnthalpyIn_BC[iMarker_Monitoring][config->GetnSpan_iZones(iMarker_Monitoring)])*100.0;
      cout << endl;
      cout << endl;
      cout << "     Inlet Entropy" << "            Inlet Entropy BC" << "            err(%)" <<  endl;
      cout.width(25); cout << EntropyIn[iMarker_Monitoring][config->GetnSpan_iZones(iMarker_Monitoring)]*config->GetEnergy_Ref()/config->GetTemperature_Ref();
      cout.width(25); cout << EntropyIn_BC[iMarker_Monitoring][config->GetnSpan_iZones(iMarker_Monitoring)]*config->GetEnergy_Ref()/config->GetTemperature_Ref();
      cout.width(25); cout << abs((EntropyIn[iMarker_Monitoring][config->GetnSpan_iZones(iMarker_Monitoring)] - EntropyIn_BC[iMarker_Monitoring][config->GetnSpan_iZones(iMarker_Monitoring)])/EntropyIn_BC[iMarker_Monitoring][config->GetnSpan_iZones(iMarker_Monitoring)])*100.0;
      cout << endl;
      cout << endl;
      cout << "     Inlet Absolute Angle" << "     Inlet Absolute Angle BC" << "     err(%)" <<  endl;
      cout.width(25); cout << 180.0/PI_NUMBER*AbsFlowAngleIn[iMarker_Monitoring][config->GetnSpan_iZones(iMarker_Monitoring)];
      cout.width(25); cout << 180.0/PI_NUMBER*FlowAngleIn_BC[iMarker_Monitoring][config->GetnSpan_iZones(iMarker_Monitoring)];
      cout.width(25); cout << abs((AbsFlowAngleIn[iMarker_Monitoring][config->GetnSpan_iZones(iMarker_Monitoring)] - FlowAngleIn_BC[iMarker_Monitoring][config->GetnSpan_iZones(iMarker_Monitoring)])/FlowAngleIn_BC[iMarker_Monitoring][config->GetnSpan_iZones(iMarker_Monitoring)])*100.0;
      cout << endl;
      cout << endl;
      if(turbulent){
        if(menter_sst){
          cout << "     Inlet TurbIntensity" << "      Inlet TurbIntensity BC" << "      err(%)" <<  endl;
          cout.width(25); cout << TurbIntensityIn[iMarker_Monitoring][config->GetnSpan_iZones(iMarker_Monitoring)];
          cout.width(25); cout << config->GetTurbulenceIntensity_FreeStream();
          cout.width(25); cout << abs((TurbIntensityIn[iMarker_Monitoring][config->GetnSpan_iZones(iMarker_Monitoring)] - config->GetTurbulenceIntensity_FreeStream())/config->GetTurbulenceIntensity_FreeStream())*100.0;
          cout << endl;
          cout << endl;
          cout << "     Inlet Turb2LamRatio" << "      Inlet Turb2LamRatio BC" << "      err(%)" <<  endl;
          cout.width(25); cout << Turb2LamViscRatioIn[iMarker_Monitoring][config->GetnSpan_iZones(iMarker_Monitoring)];
          cout.width(25); cout << config->GetTurb2LamViscRatio_FreeStream();
          cout.width(25); cout << abs((Turb2LamViscRatioIn[iMarker_Monitoring][config->GetnSpan_iZones(iMarker_Monitoring)] - config->GetTurb2LamViscRatio_FreeStream())/config->GetTurb2LamViscRatio_FreeStream())*100.0;
          cout << endl;
          cout << endl;
        }
        else{
          cout << "     Inlet Nu Factor" << "          Inlet Nu Factor BC" << "          err(%)" <<  endl;
          cout.width(25); cout << NuFactorIn[iMarker_Monitoring][config->GetnSpan_iZones(iMarker_Monitoring)];
          cout.width(25); cout << config->GetNuFactor_FreeStream();
          cout.width(25); cout << abs((NuFactorIn[iMarker_Monitoring][config->GetnSpan_iZones(iMarker_Monitoring)] - config->GetNuFactor_FreeStream())/config->GetNuFactor_FreeStream())*100.0;
          cout << endl;
          cout << endl;
        }
      }
    }
    if(iMarker_Monitoring == config->GetnMarker_Turbomachinery() -1 ){
      // if BC outlet
      cout << "BC outlet convergence monitoring  marker " << outMarker_Tag << " : "<<endl;
      cout << endl;
      cout << "     Outlet Pressure" << "          Outlet Pressure BC" << "          err(%)" <<  endl;
      cout.width(25); cout << PressureOut[iMarker_Monitoring][config->GetnSpan_iZones(iMarker_Monitoring)]*config->GetPressure_Ref();
      cout.width(25); cout << PressureOut_BC[iMarker_Monitoring][config->GetnSpan_iZones(iMarker_Monitoring)]*config->GetPressure_Ref();
      cout.width(25); cout << abs((PressureOut[iMarker_Monitoring][config->GetnSpan_iZones(iMarker_Monitoring)] - PressureOut_BC[iMarker_Monitoring][config->GetnSpan_iZones(iMarker_Monitoring)])/PressureOut_BC[iMarker_Monitoring][config->GetnSpan_iZones(iMarker_Monitoring)])*100.0;
      cout << endl;
      cout << endl;
    }

    cout << "Convergence monitoring for integral quantities between markers " << inMarker_Tag << " and "<< outMarker_Tag << " : "<<endl;
    cout << endl;
    cout << "     Inlet Mass Flow " << "         Outlet Mass Flow" << "            err(%)" <<  endl;
    cout.width(25); cout << MassFlowIn[iMarker_Monitoring][config->GetnSpan_iZones(iMarker_Monitoring)]*config->GetVelocity_Ref()*config->GetDensity_Ref();
    cout.width(25); cout << MassFlowOut[iMarker_Monitoring][config->GetnSpan_iZones(iMarker_Monitoring)]*config->GetVelocity_Ref()*config->GetDensity_Ref();
    cout.width(25); cout << abs((MassFlowIn[iMarker_Monitoring][config->GetnSpan_iZones(iMarker_Monitoring)] - MassFlowOut[iMarker_Monitoring][config->GetnSpan_iZones(iMarker_Monitoring)])/MassFlowIn[iMarker_Monitoring][config->GetnSpan_iZones(iMarker_Monitoring)])*100.0;
    cout << endl;
    cout << endl;
    //if(stator)
    //cout << "     Inlet Total Enthalpy " << "    Outlet Total Enthalpy" << "     err(%)" <<  endl;
    //else
    cout << "     Inlet Total Rothalpy " << "    Outlet Total Rothalpy" << "       err(%)" <<  endl;
    cout.width(25); cout << RothalpyIn[iMarker_Monitoring][config->GetnSpan_iZones(iMarker_Monitoring)]*config->GetEnergy_Ref();
    cout.width(25); cout << RothalpyOut[iMarker_Monitoring][config->GetnSpan_iZones(iMarker_Monitoring)]*config->GetEnergy_Ref();
    cout.width(25); cout << abs((RothalpyIn[iMarker_Monitoring][config->GetnSpan_iZones(iMarker_Monitoring)] - RothalpyOut[iMarker_Monitoring][config->GetnSpan_iZones(iMarker_Monitoring)])/RothalpyIn[iMarker_Monitoring][config->GetnSpan_iZones(iMarker_Monitoring)])*100.0;
    cout << endl;
    cout << endl;
    cout << "Blade performance between boundaries " << inMarker_Tag << " and "<< outMarker_Tag << " : "<<endl;
    cout << endl;
    cout << "     Total Pressure Loss(%)" << "   Kinetic Energy Loss(%)" << "      Entropy Generation(%)" << endl;
    cout.width(25); cout << TotalPressureLoss[iMarker_Monitoring][config->GetnSpan_iZones(iMarker_Monitoring)]*100.0;
    cout.width(25); cout << KineticEnergyLoss[iMarker_Monitoring][config->GetnSpan_iZones(iMarker_Monitoring)]*100.0;
    cout.width(25); cout << EntropyGen[iMarker_Monitoring][config->GetnSpan_iZones(iMarker_Monitoring)]*100.0;
    cout << endl;
    cout << endl;
    cout << "     Total Inlet Enthalpy" << "     Eulerian Work" << "               Pressure Ratio" <<  endl;
    cout.width(25); cout << TotalEnthalpyIn[iMarker_Monitoring][config->GetnSpan_iZones(iMarker_Monitoring)]*config->GetEnergy_Ref();
    cout.width(25); cout << EulerianWork[iMarker_Monitoring][config->GetnSpan_iZones(iMarker_Monitoring)]*config->GetEnergy_Ref();
    cout.width(25); cout << PressureRatio[iMarker_Monitoring][config->GetnSpan_iZones(iMarker_Monitoring)];
    cout << endl;
    cout << endl;
    cout << "     Inlet Entropy" << "            Outlet Entropy" << "             Outlet Is. Enthalpy" <<  endl;
    cout.width(25); cout << EntropyIn[iMarker_Monitoring][config->GetnSpan_iZones(iMarker_Monitoring)]*config->GetEnergy_Ref()/config->GetTemperature_Ref();
    cout.width(25); cout << EntropyOut[iMarker_Monitoring][config->GetnSpan_iZones(iMarker_Monitoring)]*config->GetEnergy_Ref()/config->GetTemperature_Ref();
    cout.width(25); cout << EnthalpyOutIs[iMarker_Monitoring][config->GetnSpan_iZones(iMarker_Monitoring)]*config->GetEnergy_Ref();
    cout << endl;
    cout << endl;
    cout << "Cinematic quantities between boundaries " << inMarker_Tag << " and "<< outMarker_Tag << " : "<<endl;
    cout << endl;
    cout << "     Inlet Mach"<< "               Inlet Normal Mach" << "            Inlet Tang. Mach" << endl;
    cout.width(25); cout << sqrt(MachIn[iMarker_Monitoring][config->GetnSpan_iZones(iMarker_Monitoring)][0]*MachIn[iMarker_Monitoring][config->GetnSpan_iZones(iMarker_Monitoring)][0] +MachIn[iMarker_Monitoring][config->GetnSpan_iZones(iMarker_Monitoring)][1]*MachIn[iMarker_Monitoring][config->GetnSpan_iZones(iMarker_Monitoring)][1]);
    cout.width(25); cout << MachIn[iMarker_Monitoring][config->GetnSpan_iZones(iMarker_Monitoring)][0];
    cout.width(25); cout << MachIn[iMarker_Monitoring][config->GetnSpan_iZones(iMarker_Monitoring)][1];
    cout << endl;
    cout << endl;
    cout << "     Outlet Mach"<< "              Outlet Normal Mach" << "           Outlet Tang. Mach" << endl;
    cout.width(25); cout << sqrt(MachOut[iMarker_Monitoring][config->GetnSpan_iZones(iMarker_Monitoring)][0]*MachOut[iMarker_Monitoring][config->GetnSpan_iZones(iMarker_Monitoring)][0] +MachOut[iMarker_Monitoring][config->GetnSpan_iZones(iMarker_Monitoring)][1]*MachOut[iMarker_Monitoring][config->GetnSpan_iZones(iMarker_Monitoring)][1]);
    cout.width(25); cout << MachOut[iMarker_Monitoring][config->GetnSpan_iZones(iMarker_Monitoring)][0];
    cout.width(25); cout << MachOut[iMarker_Monitoring][config->GetnSpan_iZones(iMarker_Monitoring)][1];cout << endl;
    cout << endl;
    cout << "     Inlet Flow Angle" << "         Outlet flow Angle  " << endl;
    cout.width(25); cout << 180.0/PI_NUMBER*FlowAngleIn[iMarker_Monitoring][config->GetnSpan_iZones(iMarker_Monitoring)];
    cout.width(25); cout << 180.0/PI_NUMBER*FlowAngleOut[iMarker_Monitoring][config->GetnSpan_iZones(iMarker_Monitoring)];
    cout << endl;
    cout << endl;
    // if gridmov
    cout << "     Inlet Abs Flow Angle" << "     Outlet Abs Flow Angle  " << endl;
    cout.width(25); cout << 180.0/PI_NUMBER*AbsFlowAngleIn[iMarker_Monitoring][config->GetnSpan_iZones(iMarker_Monitoring)];
    cout.width(25); cout << 180.0/PI_NUMBER*AbsFlowAngleOut[iMarker_Monitoring][config->GetnSpan_iZones(iMarker_Monitoring)];
    cout << endl;
    cout << endl << "-------------------------------------------------------------------------" << endl;
    cout << endl;
    if(nZone > 0 && iMarker_Monitoring < config->GetnMarker_Turbomachinery() -1){
      cout << endl << "---------- Mixing-Plane Interface between Blade " << iMarker_Monitoring + 1 << " and Blade " << iMarker_Monitoring + 2 << " -----------" << endl;
      cout << endl;
      inMarkerTag_Mix = config->GetMarker_TurboPerf_BoundIn(iMarker_Monitoring + 1);
      cout << "Convergence monitoring for the outlet  " << outMarker_Tag << " and the inlet  "<< inMarkerTag_Mix << " : "<<endl;
      cout << endl;
      cout << "     Outlet Density " << "          Inlet Density" << "               err(%)" <<  endl;
      cout.width(25); cout << DensityOut[iMarker_Monitoring][config->GetnSpan_iZones(iMarker_Monitoring)]*config->GetDensity_Ref();
      cout.width(25); cout << DensityIn[iMarker_Monitoring + 1][config->GetnSpan_iZones(iMarker_Monitoring + 1)]*config->GetDensity_Ref();
      cout.width(25); cout << abs((DensityIn[iMarker_Monitoring + 1][config->GetnSpan_iZones(iMarker_Monitoring +1)] - DensityOut[iMarker_Monitoring][config->GetnSpan_iZones(iMarker_Monitoring)])/DensityIn[iMarker_Monitoring + 1][config->GetnSpan_iZones(iMarker_Monitoring +1)])*100.0;
      cout << endl;
      cout << endl;
      cout << "     Outlet Pressure " << "         Inlet Pressure" << "              err(%)" <<  endl;
      cout.width(25); cout << PressureOut[iMarker_Monitoring][config->GetnSpan_iZones(iMarker_Monitoring)]*config->GetPressure_Ref();
      cout.width(25); cout << PressureIn[iMarker_Monitoring + 1][config->GetnSpan_iZones(iMarker_Monitoring +1)]*config->GetPressure_Ref();
      cout.width(25); cout << abs((PressureIn[iMarker_Monitoring + 1][config->GetnSpan_iZones(iMarker_Monitoring +1)] - PressureOut[iMarker_Monitoring][config->GetnSpan_iZones(iMarker_Monitoring)])/PressureIn[iMarker_Monitoring + 1][config->GetnSpan_iZones(iMarker_Monitoring +1)])*100.0;
      cout << endl;
      cout << endl;
      cout << "     Outlet Normal Velocity " << "  Inlet Normal Velocity" << "       err(%)" <<  endl;
      cout.width(25); cout << TurboVelocityOut[iMarker_Monitoring][config->GetnSpan_iZones(iMarker_Monitoring)][0]*config->GetVelocity_Ref();
      cout.width(25); cout << TurboVelocityIn[iMarker_Monitoring + 1][config->GetnSpan_iZones(iMarker_Monitoring +1)][0]*config->GetVelocity_Ref();
      cout.width(25); cout << abs((TurboVelocityIn[iMarker_Monitoring + 1][config->GetnSpan_iZones(iMarker_Monitoring +1)][0] - TurboVelocityOut[iMarker_Monitoring][config->GetnSpan_iZones(iMarker_Monitoring)][0])/TurboVelocityIn[iMarker_Monitoring+1][config->GetnSpan_iZones(iMarker_Monitoring +1)][0])*100.0;
      cout << endl;
      cout << endl;
      cout << "     Outlet Tang. Velocity " << "   Inlet Tang. Velocity" << "        err(%)" <<  endl;
      cout.width(25); cout << TurboVelocityOut[iMarker_Monitoring][config->GetnSpan_iZones(iMarker_Monitoring)][1]*config->GetVelocity_Ref();
      cout.width(25); cout << TurboVelocityIn[iMarker_Monitoring + 1][config->GetnSpan_iZones(iMarker_Monitoring +1)][1]*config->GetVelocity_Ref();
      cout.width(25); cout << abs((TurboVelocityIn[iMarker_Monitoring + 1][config->GetnSpan_iZones(iMarker_Monitoring +1)][1] - TurboVelocityOut[iMarker_Monitoring][config->GetnSpan_iZones(iMarker_Monitoring)][1])/TurboVelocityIn[iMarker_Monitoring+1][config->GetnSpan_iZones(iMarker_Monitoring +1)][1])*100.0;
      cout << endl;
      cout << endl;
      cout << "     Outlet Entropy " << "          Inlet Entropy" << "               err(%)" <<  endl;
      cout.width(25); cout << EntropyOut[iMarker_Monitoring][config->GetnSpan_iZones(iMarker_Monitoring)]*config->GetEnergy_Ref()/config->GetTemperature_Ref();
      cout.width(25); cout << EntropyIn[iMarker_Monitoring + 1][config->GetnSpan_iZones(iMarker_Monitoring +1)]*config->GetEnergy_Ref()/config->GetTemperature_Ref();
      cout.width(25); cout << abs((EntropyIn[iMarker_Monitoring + 1][config->GetnSpan_iZones(iMarker_Monitoring+1)] - EntropyOut[iMarker_Monitoring][config->GetnSpan_iZones(iMarker_Monitoring)])/EntropyIn[iMarker_Monitoring + 1][config->GetnSpan_iZones(iMarker_Monitoring+1)])*100.0;
      if(turbulent){
        cout << endl;
        cout << endl;
        if(menter_sst){
          cout << "     Outlet TurbIntensity " << "    Inlet TurbIntensity" << "         err(%)" <<  endl;
          cout.width(25); cout << TurbIntensityOut[iMarker_Monitoring][config->GetnSpan_iZones(iMarker_Monitoring)];
          cout.width(25); cout << TurbIntensityIn[iMarker_Monitoring + 1][config->GetnSpan_iZones(iMarker_Monitoring +1)];
          cout.width(25); cout << abs((TurbIntensityIn[iMarker_Monitoring + 1][config->GetnSpan_iZones(iMarker_Monitoring+1)] - TurbIntensityOut[iMarker_Monitoring][config->GetnSpan_iZones(iMarker_Monitoring)])/TurbIntensityIn[iMarker_Monitoring + 1][config->GetnSpan_iZones(iMarker_Monitoring+1)])*100.0;
          cout << endl;
          cout << endl;
          cout << "     Outlet Turb2LamRatio " << "    Inlet Turb2LamRatio" << "         err(%)" <<  endl;
          cout.width(25); cout << Turb2LamViscRatioOut[iMarker_Monitoring][config->GetnSpan_iZones(iMarker_Monitoring)];
          cout.width(25); cout << Turb2LamViscRatioIn[iMarker_Monitoring + 1][config->GetnSpan_iZones(iMarker_Monitoring +1)];
          cout.width(25); cout << abs((Turb2LamViscRatioIn[iMarker_Monitoring + 1][config->GetnSpan_iZones(iMarker_Monitoring+1)] - Turb2LamViscRatioOut[iMarker_Monitoring][config->GetnSpan_iZones(iMarker_Monitoring)])/Turb2LamViscRatioIn[iMarker_Monitoring + 1][config->GetnSpan_iZones(iMarker_Monitoring+1)])*100.0;
        }
        else{
          cout << "     Outlet Nu Factor " << "        Inlet Nu Factor" << "             err(%)" <<  endl;
          cout.width(25); cout << NuFactorOut[iMarker_Monitoring][config->GetnSpan_iZones(iMarker_Monitoring)];
          cout.width(25); cout << NuFactorIn[iMarker_Monitoring + 1][config->GetnSpan_iZones(iMarker_Monitoring +1)];
          cout.width(25); cout << abs((NuFactorIn[iMarker_Monitoring + 1][config->GetnSpan_iZones(iMarker_Monitoring+1)] - NuFactorOut[iMarker_Monitoring][config->GetnSpan_iZones(iMarker_Monitoring)])/NuFactorIn[iMarker_Monitoring + 1][config->GetnSpan_iZones(iMarker_Monitoring+1)])*100.0;
        }
      }
      cout << endl;
      cout << endl << "-------------------------------------------------------------------------" << endl;
      cout << endl;
    }

  }
  if(nZone > 1){
    /*--- Stage Performance ---*/
    for(iStage = 0; iStage < nStages; iStage++ ){
      cout << endl << "----------------------------- Stage " << iStage + 1 << " -----------------------------------" << endl;
      inMarker_Tag = config->GetMarker_TurboPerf_BoundIn(iStage*2);
      outMarker_Tag = config->GetMarker_TurboPerf_BoundOut(iStage*2+1);
      cout << "Stage performance between boundaries " << inMarker_Tag << " and "<< outMarker_Tag << " : "<<endl;
      cout << endl;
      cout << "     Total-Total Eff.(%)" << "      Total-Static Eff.(%)" << "      Entropy Generation(%)" << endl;
      cout.width(25); cout << TotalTotalEfficiency[nBladesRow + iStage][nSpanWiseSections]*100.0;
      cout.width(25); cout << TotalStaticEfficiency[nBladesRow + iStage][nSpanWiseSections]*100.0;
      cout.width(25); cout << EntropyGen[nBladesRow + iStage][nSpanWiseSections]*100.0;
      cout << endl;
      cout << endl;
      cout << "     Pressure Ratio " << "          Outlet Is. Enthalpy" << "       In-Out MassFlow Diff (%)" <<  endl;
      cout.width(25); cout << PressureRatio[nBladesRow + iStage][nSpanWiseSections];
      cout.width(25); cout << EnthalpyOutIs[nBladesRow + iStage][nSpanWiseSections]*config->GetEnergy_Ref();
      cout.width(25); cout << abs((MassFlowIn[nBladesRow + iStage][nSpanWiseSections] - MassFlowOut[nBladesRow + iStage][nSpanWiseSections])/MassFlowIn[nBladesRow + iStage][nSpanWiseSections])*100.0;
    }
    cout << endl;
    cout << endl << "-------------------------------------------------------------------------" << endl;
    cout << endl;

    /*--- Full Machine Performance ---*/
    // if(turbine)
    cout << endl << "---------------------------- Turbine ------------------------------------" << endl;
    inMarker_Tag = config->GetMarker_TurboPerf_BoundIn(0);
    outMarker_Tag = config->GetMarker_TurboPerf_BoundOut(nBladesRow-1);
    cout << "Turbine performance between boundaries " << inMarker_Tag << " and "<< outMarker_Tag << " : "<<endl;
    cout << endl;
    cout << "     Total-Total Eff.(%)" << "      Total-Static Eff.(%)" << "      Entropy Generation(%)" << endl;
    cout.width(25); cout << TotalTotalEfficiency[nBladesRow + nStages][nSpanWiseSections]*100.0;
    cout.width(25); cout << TotalStaticEfficiency[nBladesRow + nStages][nSpanWiseSections]*100.0;
    cout.width(25); cout << EntropyGen[nBladesRow + nStages][nSpanWiseSections]*100.0;
    cout << endl;
    cout << endl;
    cout << "     Pressure Ratio " << "          Outlet Is. Enthalpy" << "       In-Out MassFlow Diff (%)" <<  endl;
    cout.width(25); cout << PressureRatio[nBladesRow + nStages][nSpanWiseSections];
    cout.width(25); cout << EnthalpyOutIs[nBladesRow + nStages][nSpanWiseSections]*config->GetEnergy_Ref();;
    cout.width(25); cout << abs((MassFlowIn[nBladesRow + nStages][nSpanWiseSections] - MassFlowOut[nBladesRow + nStages][nSpanWiseSections])/MassFlowIn[nBladesRow + nStages][nSpanWiseSections])*100.0;
    cout << endl;
    cout << endl << "-------------------------------------------------------------------------" << endl;
    cout << endl;
  }

}

void COutputLegacy::SpecialOutput_Turbo(CSolver *****solver, CGeometry ****geometry, CConfig **config,
                                       unsigned short val_iZone, bool output) {

  string inMarker_Tag, outMarker_Tag, inMarkerTag_Mix;
  unsigned short nZone       = config[val_iZone]->GetnZone();

  unsigned short iDim, iSpan;

  unsigned long iExtIter = config[val_iZone]->GetInnerIter();
  const su2double* SpanWiseValuesIn, *SpanWiseValuesOut;
  ofstream myfile;
  string spanwise_performance_filename;


  /*--- Start of write file turboperformance spanwise ---*/
  if (rank == MASTER_NODE){
    SpanWiseValuesIn = geometry[val_iZone][INST_0][MESH_0]->GetSpanWiseValue(1);
    SpanWiseValuesOut = geometry[val_iZone][INST_0][MESH_0]->GetSpanWiseValue(2);



    /*--- Writing Span wise inflow thermodynamic quantities. ---*/
    spanwise_performance_filename = "TURBOMACHINERY/inflow_spanwise_thermodynamic_values.dat";
    char buffer[50];
    if (nZone > 1){
      unsigned short lastindex      =  spanwise_performance_filename.find_last_of(".");
      spanwise_performance_filename =  spanwise_performance_filename.substr(0, lastindex);
      SPRINTF (buffer, "_%d.dat", SU2_TYPE::Int(val_iZone));
      spanwise_performance_filename.append(string(buffer));
    }


    myfile.open (spanwise_performance_filename.data(), ios::out | ios::trunc);
    myfile.setf(ios::scientific);
    myfile.precision(12);

    myfile << "TITLE = \"Inflow Spanwise Thermodynamic Values. iExtIter = " << iExtIter << " \"" << endl;
    myfile << "VARIABLES =" << endl;

    myfile.width(30); myfile << "\"SpanWise Value[m]\"";
    myfile.width(15); myfile << "\"iSpan\"";
    myfile.width(30); myfile << "\"Pressure[Pa]\"";
    myfile.width(30); myfile << "\"TotalPressure[Pa]\"";
    myfile.width(30); myfile << "\"Temperature[K]\"";
    myfile.width(30); myfile << "\"TotalTemperature[K]\"";
    myfile.width(30); myfile << "\"Enthalpy[J]\"";
    myfile.width(30); myfile << "\"TotalEnthalpy[J]\"";
    myfile.width(30); myfile << "\"Density[kg/m3]\"";
    myfile.width(30); myfile << "\"Entropy[J/K]\"";
    myfile.width(30); myfile << "\"TurbIntensity[-]\"";
    myfile.width(30); myfile << "\"Turb2LamViscRatio[-]\"";
    myfile.width(30); myfile << "\"NuFactor[-]\"";
    myfile << endl;

    for(iSpan = 0; iSpan < config[ZONE_0]->GetnSpan_iZones(val_iZone); iSpan++){

      myfile.width(30); myfile << SpanWiseValuesIn[iSpan];
      myfile.width(15); myfile << iSpan;
      myfile.width(30); myfile << PressureIn           [val_iZone][iSpan]*config[ZONE_0]->GetPressure_Ref();
      myfile.width(30); myfile << TotalPressureIn      [val_iZone][iSpan]*config[ZONE_0]->GetPressure_Ref();
      myfile.width(30); myfile << TemperatureIn        [val_iZone][iSpan]*config[ZONE_0]->GetTemperature_Ref();
      myfile.width(30); myfile << TotalTemperatureIn   [val_iZone][iSpan]*config[ZONE_0]->GetTemperature_Ref();
      myfile.width(30); myfile << EnthalpyIn           [val_iZone][iSpan]*config[ZONE_0]->GetEnergy_Ref();
      myfile.width(30); myfile << TotalEnthalpyIn      [val_iZone][iSpan]*config[ZONE_0]->GetEnergy_Ref();
      myfile.width(30); myfile << DensityIn            [val_iZone][iSpan]*config[ZONE_0]->GetDensity_Ref();
      myfile.width(30); myfile << EntropyIn            [val_iZone][iSpan]*config[ZONE_0]->GetEnergy_Ref()/config[ZONE_0]->GetTemperature_Ref();
      if(TurbIntensityIn[val_iZone][iSpan] > 1.0){
        myfile.width(30); myfile << TurbIntensityIn      [val_iZone][config[ZONE_0]->GetnSpan_iZones(val_iZone)/2];
      }else{
        myfile.width(30); myfile << TurbIntensityIn      [val_iZone][iSpan];
      }
      myfile.width(30); myfile << Turb2LamViscRatioIn  [val_iZone][iSpan];
      myfile.width(30); myfile << NuFactorIn           [val_iZone][iSpan];
      myfile << endl;
    }

    myfile.close();

    /*--- Writing Span wise outflow thermodynamic quantities. ---*/
    spanwise_performance_filename = "TURBOMACHINERY/outflow_spanwise_thermodynamic_values.dat";
    if (nZone > 1){
      unsigned short lastindex      =  spanwise_performance_filename.find_last_of(".");
      spanwise_performance_filename =  spanwise_performance_filename.substr(0, lastindex);
      SPRINTF (buffer, "_%d.dat", SU2_TYPE::Int(val_iZone));
      spanwise_performance_filename.append(string(buffer));
    }

    myfile.open (spanwise_performance_filename.data(), ios::out | ios::trunc);
    myfile.setf(ios::scientific);
    myfile.precision(12);

    myfile << "TITLE = \"Outflow Span-wise Thermodynamic Values. iExtIter = " << iExtIter << " \"" << endl;
    myfile << "VARIABLES =" << endl;

    myfile.width(30); myfile << "\"SpanWise Value[m]\"";
    myfile.width(15); myfile << "\"iSpan\"";
    myfile.width(30); myfile << "\"Pressure[Pa]\"";
    myfile.width(30); myfile << "\"TotalPressure[Pa]\"";
    myfile.width(30); myfile << "\"Temperature[K]\"";
    myfile.width(30); myfile << "\"TotalTemperature[K]\"";
    myfile.width(30); myfile << "\"Enthalpy[J]\"";
    myfile.width(30); myfile << "\"TotalEnthalpy[J]\"";
    myfile.width(30); myfile << "\"Density[kg/m3]\"";
    myfile.width(30); myfile << "\"Entropy[J/K]\"";
    myfile.width(30); myfile << "\"TurbIntensity[-]\"";
    myfile.width(30); myfile << "\"Turb2LamViscRatio[-]\"";
    myfile.width(30); myfile << "\"NuFactor[-]\"";
    myfile << endl;


    for(iSpan = 0; iSpan < config[ZONE_0]->GetnSpan_iZones(val_iZone); iSpan++){

      myfile.width(30); myfile << SpanWiseValuesOut[iSpan];
      myfile.width(15); myfile << iSpan;
      myfile.width(30); myfile << PressureOut           [val_iZone][iSpan]*config[ZONE_0]->GetPressure_Ref();
      myfile.width(30); myfile << TotalPressureOut      [val_iZone][iSpan]*config[ZONE_0]->GetPressure_Ref();
      myfile.width(30); myfile << TemperatureOut        [val_iZone][iSpan]*config[ZONE_0]->GetTemperature_Ref();
      myfile.width(30); myfile << TotalTemperatureOut   [val_iZone][iSpan]*config[ZONE_0]->GetTemperature_Ref();
      myfile.width(30); myfile << EnthalpyOut           [val_iZone][iSpan]*config[ZONE_0]->GetEnergy_Ref();
      myfile.width(30); myfile << TotalEnthalpyOut      [val_iZone][iSpan]*config[ZONE_0]->GetEnergy_Ref();
      myfile.width(30); myfile << DensityOut            [val_iZone][iSpan]*config[ZONE_0]->GetDensity_Ref();
      myfile.width(30); myfile << EntropyOut            [val_iZone][iSpan]*config[ZONE_0]->GetEnergy_Ref()/config[ZONE_0]->GetTemperature_Ref();
      if(TurbIntensityOut[val_iZone][iSpan] > 1.0){
        myfile.width(30); myfile << TurbIntensityOut      [val_iZone][config[ZONE_0]->GetnSpan_iZones(val_iZone)/2];
      }else{
        myfile.width(30); myfile << TurbIntensityOut      [val_iZone][iSpan];
      }
      myfile.width(30); myfile << Turb2LamViscRatioOut  [val_iZone][iSpan];
      myfile.width(30); myfile << NuFactorOut           [val_iZone][iSpan];
      myfile << endl;
    }

    myfile.close();

    /*--- Writing Span wise inflow kinematic quantities. ---*/
    spanwise_performance_filename = "TURBOMACHINERY/inflow_spanwise_kinematic_values.dat";
    if (nZone > 1){
      unsigned short lastindex      =  spanwise_performance_filename.find_last_of(".");
      spanwise_performance_filename =  spanwise_performance_filename.substr(0, lastindex);
      SPRINTF (buffer, "_%d.dat", SU2_TYPE::Int(val_iZone));
      spanwise_performance_filename.append(string(buffer));
    }

    myfile.open (spanwise_performance_filename.data(), ios::out | ios::trunc);
    myfile.setf(ios::scientific);
    myfile.precision(12);

    myfile << "TITLE = \"Inflow Span-wise Kinematic Values. iExtIter = " << iExtIter << " \"" << endl;
    myfile << "VARIABLES =" << endl;

    myfile.width(30); myfile << "\"SpanWise Value[m]\"";
    myfile.width(15); myfile << "\"iSpan\"";
    myfile.width(30); myfile << "\"Normal Mach[-]\"";
    myfile.width(30); myfile << "\"Tangential Mach[-]\"";
    myfile.width(30); myfile << "\"3rd Component Mach[-]\"";
    myfile.width(30); myfile << "\"Mach Module[-]\"";
    myfile.width(30); myfile << "\"Normal Velocity[m/s]\"";
    myfile.width(30); myfile << "\"Tangential Velocity[m/s]\"";
    myfile.width(30); myfile << "\"3rd Component Velocity[m/s]\"";
    myfile.width(30); myfile << "\"Velocity Module[m/s]\"";
    myfile.width(30); myfile << "\"Absolute Flow Angle[deg]\"";
    myfile.width(30); myfile << "\"Relative Flow Angle[deg]\"";
    myfile << endl;


    for(iSpan = 0; iSpan < config[ZONE_0]->GetnSpan_iZones(val_iZone); iSpan++){

      myfile.width(30); myfile << SpanWiseValuesIn[iSpan];
      myfile.width(15); myfile << iSpan;
      for (iDim = 0; iDim < 4; iDim++){
        myfile.width(30); myfile << MachIn              [val_iZone][iSpan][iDim];
      }
      for (iDim = 0; iDim < 4; iDim++){
        myfile.width(30); myfile << TurboVelocityIn     [val_iZone][iSpan][iDim]*config[ZONE_0]->GetVelocity_Ref();
      }
      if(AbsFlowAngleIn[val_iZone][iSpan] != AbsFlowAngleIn[val_iZone][iSpan]){
        myfile.width(30); myfile << "0.0000";
      }
      else{
        myfile.width(30); myfile << AbsFlowAngleIn     [val_iZone][iSpan]*180.0/PI_NUMBER;
      }
      if(FlowAngleIn[val_iZone][iSpan] != FlowAngleIn[val_iZone][iSpan]){
        myfile.width(30); myfile << "0.0000";
      }
      else{
        myfile.width(30); myfile << FlowAngleIn      [val_iZone][iSpan]*180.0/PI_NUMBER;
      }
      myfile << endl;
    }

    myfile.close();

    /*--- Writing Span wise outflow thermodynamic quantities. ---*/
    spanwise_performance_filename = "TURBOMACHINERY/outflow_spanwise_kinematic_values.dat";
    if (nZone > 1){
      unsigned short lastindex      =  spanwise_performance_filename.find_last_of(".");
      spanwise_performance_filename =  spanwise_performance_filename.substr(0, lastindex);
      SPRINTF (buffer, "_%d.dat", SU2_TYPE::Int(val_iZone));
      spanwise_performance_filename.append(string(buffer));
    }

    myfile.open (spanwise_performance_filename.data(), ios::out | ios::trunc);
    myfile.setf(ios::scientific);
    myfile.precision(12);

    myfile << "TITLE = \"Outflow Span-wise Kinematic Values. iExtIter = " << iExtIter << " \"" << endl;
    myfile << "VARIABLES =" << endl;

    myfile.width(30); myfile << "\"SpanWise Value[m]\"";
    myfile.width(15); myfile << "\"iSpan\"";
    myfile.width(30); myfile << "\"Normal Mach[-]\"";
    myfile.width(30); myfile << "\"Tangential Mach[-]\"";
    myfile.width(30); myfile << "\"3rd Component Mach[-]\"";
    myfile.width(30); myfile << "\"Mach Module[-]\"";
    myfile.width(30); myfile << "\"Normal Velocity[m/s]\"";
    myfile.width(30); myfile << "\"Tangential Velocity[m/s]\"";
    myfile.width(30); myfile << "\"3rd Component Velocity[m/s]\"";
    myfile.width(30); myfile << "\"Velocity Module[m/s]\"";
    myfile.width(30); myfile << "\"Absolute Flow Angle[deg]\"";
    myfile.width(30); myfile << "\"Relative Flow Angle[deg]\"";
    myfile << endl;


    for(iSpan = 0; iSpan < config[ZONE_0]->GetnSpan_iZones(val_iZone); iSpan++){

      myfile.width(30); myfile << SpanWiseValuesOut[iSpan];
      myfile.width(15); myfile << iSpan;
      for (iDim = 0; iDim < 4; iDim++){
        myfile.width(30); myfile << MachOut              [val_iZone][iSpan][iDim];
      }
      for (iDim = 0; iDim < 4; iDim++){
        myfile.width(30); myfile << TurboVelocityOut     [val_iZone][iSpan][iDim]*config[ZONE_0]->GetVelocity_Ref();
      }
      if(AbsFlowAngleOut[val_iZone][iSpan] != AbsFlowAngleOut[val_iZone][iSpan]){
        myfile.width(30); myfile << "0.0000";
      }
      else{
        myfile.width(30); myfile << AbsFlowAngleOut      [val_iZone][iSpan]*180.0/PI_NUMBER;
      }
      if(FlowAngleOut[val_iZone][iSpan] != FlowAngleOut[val_iZone][iSpan]){
        myfile.width(30); myfile << "0.0000";
      }
      else{
        myfile.width(30); myfile << FlowAngleOut      [val_iZone][iSpan]*180.0/PI_NUMBER;
      }
      myfile << endl;
    }

    myfile.close();

  }
}

void COutputLegacy::SpecialOutput_HarmonicBalance(CSolver *****solver, CGeometry ****geometry, CConfig **config, unsigned short iInst, unsigned short val_nInst, bool output) const {

  /*--- Write file with flow quantities for harmonic balance HB ---*/
  ofstream HB_output_file;
  ofstream mean_HB_file;

  /*--- MPI Send/Recv buffers ---*/
  su2double *sbuf_var = nullptr,  *rbuf_var = nullptr;

  /*--- Other variables ---*/
  unsigned short iVar, kInst;
  unsigned short nVar_output = 5;
  unsigned long current_iter = config[ZONE_0]->GetInnerIter();

  /*--- Allocate memory for send buffer ---*/
  sbuf_var = new su2double[nVar_output];

  su2double *averages = new su2double[nVar_output];
  for (iVar = 0; iVar < nVar_output; iVar++)
    averages[iVar] = 0;

  /*--- Allocate memory for receive buffer ---*/
  if (rank == MASTER_NODE) {
    rbuf_var = new su2double[nVar_output];

    HB_output_file.precision(15);
    HB_output_file.open("HB_output.csv", ios::out);
    HB_output_file <<  "\"time_instance\",\"CL\",\"CD\",\"CMx\",\"CMy\",\"CMz\"" << endl;

    mean_HB_file.precision(15);
    if (current_iter == 0 && iInst == 1) {
      mean_HB_file.open("history_HB.plt", ios::trunc);
      mean_HB_file << "TITLE = \"SU2 HARMONIC BALANCE SIMULATION\"" << endl;
      mean_HB_file <<  "VARIABLES = \"Iteration\",\"CL\",\"CD\",\"CMx\",\"CMy\",\"CMz\",\"CT\",\"CQ\",\"CMerit\"" << endl;
      mean_HB_file << "ZONE T= \"Average Convergence History\"" << endl;
    }
    else
      mean_HB_file.open("history_HB.plt", ios::out | ios::app);
  }

  if (rank == MASTER_NODE) {

    /*--- Run through the zones, collecting the output variables
       N.B. Summing across processors within a given zone is being done
       elsewhere. ---*/
    for (kInst = 0; kInst < val_nInst; kInst++) {

      /*--- Flow solution coefficients (parallel) ---*/
      sbuf_var[0] = solver[ZONE_0][kInst][MESH_0][FLOW_SOL]->GetTotal_CL();
      sbuf_var[1] = solver[ZONE_0][kInst][INST_0][FLOW_SOL]->GetTotal_CD();
      sbuf_var[2] = solver[ZONE_0][kInst][INST_0][FLOW_SOL]->GetTotal_CMx();
      sbuf_var[3] = solver[ZONE_0][kInst][INST_0][FLOW_SOL]->GetTotal_CMy();
      sbuf_var[4] = solver[ZONE_0][kInst][INST_0][FLOW_SOL]->GetTotal_CMz();

      for (iVar = 0; iVar < nVar_output; iVar++) {
        rbuf_var[iVar] = sbuf_var[iVar];
      }

      HB_output_file << kInst << ", ";
      for (iVar = 0; iVar < nVar_output; iVar++)
        HB_output_file << rbuf_var[iVar] << ", ";
      HB_output_file << endl;

      /*--- Increment the total contributions from each zone, dividing by nZone as you go ---*/
      for (iVar = 0; iVar < nVar_output; iVar++) {
        averages[iVar] += (1.0/su2double(val_nInst))*rbuf_var[iVar];
      }
    }
  }

  if (rank == MASTER_NODE && iInst == INST_0) {

    mean_HB_file << current_iter << ", ";
    for (iVar = 0; iVar < nVar_output; iVar++) {
      mean_HB_file << averages[iVar];
      if (iVar < nVar_output-1)
        mean_HB_file << ", ";
    }
    mean_HB_file << endl;
  }

  if (rank == MASTER_NODE) {
    HB_output_file.close();
    mean_HB_file.close();
    delete [] rbuf_var;
  }

  delete [] sbuf_var;
  delete [] averages;
}

void COutputLegacy::SetSpecial_Output(CSolver *****solver_container,
                                       CGeometry ****geometry,
                                       CConfig **config,
                                       unsigned long iExtIter,
                                       unsigned short val_nZone) {

  bool special_output = false;
  unsigned short iZone;

  for (iZone = 0; iZone < val_nZone; iZone++){

    special_output = config[iZone]->GetSpecial_Output();

    /*--- Output a file with the forces breakdown. ---*/
    if (config[iZone]->GetWrt_ForcesBreakdown())
      SpecialOutput_ForcesBreakdown(solver_container, geometry, config, iZone, special_output);

  }

}

void COutputLegacy::SpecialOutput_AnalyzeSurface(CSolver *solver, CGeometry *geometry, CConfig *config, bool output) const {

  unsigned short iDim, iMarker, iMarker_Analyze;
  unsigned long iVertex, iPoint;
  su2double Mach = 0.0, Pressure, Temperature = 0.0, TotalPressure = 0.0, TotalTemperature = 0.0,
  Enthalpy, Velocity[3]= {0.0}, TangVel[3], Velocity2, MassFlow, Density, Area,
  AxiFactor = 1.0, SoundSpeed, Vn, Vn2, Vtang2, Weight = 1.0;

  su2double Gas_Constant      = config->GetGas_ConstantND();
  su2double Gamma             = config->GetGamma();
  unsigned short nMarker      = config->GetnMarker_All();
  unsigned short nDim         = geometry->GetnDim();
  unsigned short Kind_Average = config->GetKind_Average();

  bool compressible   = config->GetKind_Regime() == ENUM_REGIME::COMPRESSIBLE;
  bool incompressible = config->GetKind_Regime() == ENUM_REGIME::INCOMPRESSIBLE;
  bool energy         = config->GetEnergy_Equation();


  bool axisymmetric               = config->GetAxisymmetric();
  unsigned short nMarker_Analyze  = config->GetnMarker_Analyze();

  su2double  *Vector                    = new su2double[nDim];
  su2double  *Surface_MassFlow          = new su2double[nMarker];
  su2double  *Surface_Mach              = new su2double[nMarker];
  su2double  *Surface_Temperature       = new su2double[nMarker];
  su2double  *Surface_Density           = new su2double[nMarker];
  su2double  *Surface_Enthalpy          = new su2double[nMarker];
  su2double  *Surface_NormalVelocity    = new su2double[nMarker];
  su2double  *Surface_StreamVelocity2   = new su2double[nMarker];
  su2double  *Surface_TransvVelocity2   = new su2double[nMarker];
  su2double  *Surface_Pressure          = new su2double[nMarker];
  su2double  *Surface_TotalTemperature  = new su2double[nMarker];
  su2double  *Surface_TotalPressure     = new su2double[nMarker];
  su2double  *Surface_VelocityIdeal     = new su2double[nMarker];
  su2double  *Surface_Area              = new su2double[nMarker];
  su2double  *Surface_MassFlow_Abs      = new su2double[nMarker];

  /*--- Compute the numerical fan face Mach number, and the total area of the inflow ---*/

  for (iMarker = 0; iMarker < nMarker; iMarker++) {

    Surface_MassFlow[iMarker]          = 0.0;
    Surface_Mach[iMarker]              = 0.0;
    Surface_Temperature[iMarker]       = 0.0;
    Surface_Density[iMarker]           = 0.0;
    Surface_Enthalpy[iMarker]          = 0.0;
    Surface_NormalVelocity[iMarker]    = 0.0;
    Surface_StreamVelocity2[iMarker]   = 0.0;
    Surface_TransvVelocity2[iMarker]   = 0.0;
    Surface_Pressure[iMarker]          = 0.0;
    Surface_TotalTemperature[iMarker]  = 0.0;
    Surface_TotalPressure[iMarker]     = 0.0;
    Surface_VelocityIdeal[iMarker]     = 0.0;
    Surface_Area[iMarker]              = 0.0;
    Surface_MassFlow_Abs[iMarker]      = 0.0;

    if (config->GetMarker_All_Analyze(iMarker) == YES) {

      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

        if (geometry->nodes->GetDomain(iPoint)) {

          geometry->vertex[iMarker][iVertex]->GetNormal(Vector);

          if (axisymmetric) {
            if (geometry->nodes->GetCoord(iPoint, 1) != 0.0)
              AxiFactor = 2.0*PI_NUMBER*geometry->nodes->GetCoord(iPoint, 1);
            else
              AxiFactor = 1.0;
          } else {
            AxiFactor = 1.0;
          }

          Density = solver->GetNodes()->GetDensity(iPoint);
          Velocity2 = 0.0; Area = 0.0; MassFlow = 0.0; Vn = 0.0; Vtang2 = 0.0;

          for (iDim = 0; iDim < nDim; iDim++) {
            Area += (Vector[iDim] * AxiFactor) * (Vector[iDim] * AxiFactor);
            Velocity[iDim] = solver->GetNodes()->GetVelocity(iPoint,iDim);
            Velocity2 += Velocity[iDim] * Velocity[iDim];
            Vn += Velocity[iDim] * Vector[iDim] * AxiFactor;
            MassFlow += Vector[iDim] * AxiFactor * Density * Velocity[iDim];
          }

          Area       = sqrt (Area);
          if (AxiFactor == 0.0) Vn = 0.0; else Vn /= Area;
          Vn2        = Vn * Vn;
          Pressure   = solver->GetNodes()->GetPressure(iPoint);
          SoundSpeed = solver->GetNodes()->GetSoundSpeed(iPoint);

          for (iDim = 0; iDim < nDim; iDim++) {
            TangVel[iDim] = Velocity[iDim] - Vn*Vector[iDim]*AxiFactor/Area;
            Vtang2       += TangVel[iDim]*TangVel[iDim];
          }

          if (incompressible){
            if (config->GetKind_DensityModel() == INC_DENSITYMODEL::VARIABLE) {
              Mach = sqrt(solver->GetNodes()->GetVelocity2(iPoint))/
              sqrt(solver->GetNodes()->GetSpecificHeatCp(iPoint)*config->GetPressure_ThermodynamicND()/(solver->GetNodes()->GetSpecificHeatCv(iPoint)*solver->GetNodes()->GetDensity(iPoint)));
            } else {
              Mach = sqrt(solver->GetNodes()->GetVelocity2(iPoint))/
              sqrt(config->GetBulk_Modulus()/(solver->GetNodes()->GetDensity(iPoint)));
            }
            Temperature       = solver->GetNodes()->GetTemperature(iPoint);
            Enthalpy          = solver->GetNodes()->GetSpecificHeatCp(iPoint)*Temperature;
            TotalTemperature  = Temperature + 0.5*Velocity2/solver->GetNodes()->GetSpecificHeatCp(iPoint);
            TotalPressure     = Pressure + 0.5*Density*Velocity2;
          }
          else{
            Mach              = sqrt(Velocity2)/SoundSpeed;
            Temperature       = Pressure / (Gas_Constant * Density);
            Enthalpy          = solver->GetNodes()->GetEnthalpy(iPoint);
            TotalTemperature  = Temperature * (1.0 + Mach * Mach * 0.5 * (Gamma - 1.0));
            TotalPressure     = Pressure * pow( 1.0 + Mach * Mach * 0.5 * (Gamma - 1.0), Gamma / (Gamma - 1.0));
          }

          /*--- Compute the mass Surface_MassFlow ---*/

          Surface_Area[iMarker]             += Area;
          Surface_MassFlow[iMarker]         += MassFlow;
          Surface_MassFlow_Abs[iMarker]     += abs(MassFlow);

          if (Kind_Average == AVERAGE_MASSFLUX) Weight = abs(MassFlow);
          else if (Kind_Average == AVERAGE_AREA) Weight = abs(Area);
          else Weight = 1.0;

          Surface_Mach[iMarker]             += Mach*Weight;
          Surface_Temperature[iMarker]      += Temperature*Weight;
          Surface_Density[iMarker]          += Density*Weight;
          Surface_Enthalpy[iMarker]         += Enthalpy*Weight;
          Surface_NormalVelocity[iMarker]   += Vn*Weight;
          Surface_Pressure[iMarker]         += Pressure*Weight;
          Surface_TotalTemperature[iMarker] += TotalTemperature*Weight;
          Surface_TotalPressure[iMarker]    += TotalPressure*Weight;

          /*--- For now, always used the area to weight the uniformities. ---*/

          Weight = abs(Area);

          Surface_StreamVelocity2[iMarker]   += Vn2*Weight;
          Surface_TransvVelocity2[iMarker]   += Vtang2*Weight;

        }
      }

    }

  }

  /*--- Copy to the appropriate structure ---*/

  su2double *Surface_MassFlow_Local          = new su2double [nMarker_Analyze];
  su2double *Surface_Mach_Local              = new su2double [nMarker_Analyze];
  su2double *Surface_Temperature_Local       = new su2double [nMarker_Analyze];
  su2double *Surface_Density_Local           = new su2double [nMarker_Analyze];
  su2double *Surface_Enthalpy_Local          = new su2double [nMarker_Analyze];
  su2double *Surface_NormalVelocity_Local    = new su2double [nMarker_Analyze];
  su2double *Surface_StreamVelocity2_Local   = new su2double [nMarker_Analyze];
  su2double *Surface_TransvVelocity2_Local   = new su2double [nMarker_Analyze];
  su2double *Surface_Pressure_Local          = new su2double [nMarker_Analyze];
  su2double *Surface_TotalTemperature_Local  = new su2double [nMarker_Analyze];
  su2double *Surface_TotalPressure_Local     = new su2double [nMarker_Analyze];
  su2double *Surface_Area_Local              = new su2double [nMarker_Analyze];
  su2double *Surface_MassFlow_Abs_Local      = new su2double [nMarker_Analyze];

  su2double *Surface_MassFlow_Total          = new su2double [nMarker_Analyze];
  su2double *Surface_Mach_Total              = new su2double [nMarker_Analyze];
  su2double *Surface_Temperature_Total       = new su2double [nMarker_Analyze];
  su2double *Surface_Density_Total           = new su2double [nMarker_Analyze];
  su2double *Surface_Enthalpy_Total          = new su2double [nMarker_Analyze];
  su2double *Surface_NormalVelocity_Total    = new su2double [nMarker_Analyze];
  su2double *Surface_StreamVelocity2_Total   = new su2double [nMarker_Analyze];
  su2double *Surface_TransvVelocity2_Total   = new su2double [nMarker_Analyze];
  su2double *Surface_Pressure_Total          = new su2double [nMarker_Analyze];
  su2double *Surface_TotalTemperature_Total  = new su2double [nMarker_Analyze];
  su2double *Surface_TotalPressure_Total     = new su2double [nMarker_Analyze];
  su2double *Surface_Area_Total              = new su2double [nMarker_Analyze];
  su2double *Surface_MassFlow_Abs_Total      = new su2double [nMarker_Analyze];

  su2double *Surface_MomentumDistortion_Total = new su2double [nMarker_Analyze];

  for (iMarker_Analyze = 0; iMarker_Analyze < nMarker_Analyze; iMarker_Analyze++) {
    Surface_MassFlow_Local[iMarker_Analyze]          = 0.0;
    Surface_Mach_Local[iMarker_Analyze]              = 0.0;
    Surface_Temperature_Local[iMarker_Analyze]       = 0.0;
    Surface_Density_Local[iMarker_Analyze]           = 0.0;
    Surface_Enthalpy_Local[iMarker_Analyze]          = 0.0;
    Surface_NormalVelocity_Local[iMarker_Analyze]    = 0.0;
    Surface_StreamVelocity2_Local[iMarker_Analyze]   = 0.0;
    Surface_TransvVelocity2_Local[iMarker_Analyze]   = 0.0;
    Surface_Pressure_Local[iMarker_Analyze]          = 0.0;
    Surface_TotalTemperature_Local[iMarker_Analyze]  = 0.0;
    Surface_TotalPressure_Local[iMarker_Analyze]     = 0.0;
    Surface_Area_Local[iMarker_Analyze]              = 0.0;
    Surface_MassFlow_Abs_Local[iMarker_Analyze]      = 0.0;

    Surface_MassFlow_Total[iMarker_Analyze]          = 0.0;
    Surface_Mach_Total[iMarker_Analyze]              = 0.0;
    Surface_Temperature_Total[iMarker_Analyze]       = 0.0;
    Surface_Density_Total[iMarker_Analyze]           = 0.0;
    Surface_Enthalpy_Total[iMarker_Analyze]          = 0.0;
    Surface_NormalVelocity_Total[iMarker_Analyze]    = 0.0;
    Surface_StreamVelocity2_Total[iMarker_Analyze]   = 0.0;
    Surface_TransvVelocity2_Total[iMarker_Analyze]   = 0.0;
    Surface_Pressure_Total[iMarker_Analyze]          = 0.0;
    Surface_TotalTemperature_Total[iMarker_Analyze]  = 0.0;
    Surface_TotalPressure_Total[iMarker_Analyze]     = 0.0;
    Surface_Area_Total[iMarker_Analyze]              = 0.0;
    Surface_MassFlow_Abs_Total[iMarker_Analyze]      = 0.0;

    Surface_MomentumDistortion_Total[iMarker_Analyze] = 0.0;

  }

  /*--- Compute the numerical fan face Mach number, mach number, temperature and the total area ---*/

  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {

    if (config->GetMarker_All_Analyze(iMarker) == YES)  {

      for (iMarker_Analyze= 0; iMarker_Analyze < nMarker_Analyze; iMarker_Analyze++) {

        /*--- Add the Surface_MassFlow, and Surface_Area to the particular boundary ---*/

        if (config->GetMarker_All_TagBound(iMarker) == config->GetMarker_Analyze_TagBound(iMarker_Analyze)) {
          Surface_MassFlow_Local[iMarker_Analyze]          += Surface_MassFlow[iMarker];
          Surface_Mach_Local[iMarker_Analyze]              += Surface_Mach[iMarker];
          Surface_Temperature_Local[iMarker_Analyze]       += Surface_Temperature[iMarker];
          Surface_Density_Local[iMarker_Analyze]           += Surface_Density[iMarker];
          Surface_Enthalpy_Local[iMarker_Analyze]          += Surface_Enthalpy[iMarker];
          Surface_NormalVelocity_Local[iMarker_Analyze]    += Surface_NormalVelocity[iMarker];
          Surface_StreamVelocity2_Local[iMarker_Analyze]   += Surface_StreamVelocity2[iMarker];
          Surface_TransvVelocity2_Local[iMarker_Analyze]   += Surface_TransvVelocity2[iMarker];
          Surface_Pressure_Local[iMarker_Analyze]          += Surface_Pressure[iMarker];
          Surface_TotalTemperature_Local[iMarker_Analyze]  += Surface_TotalTemperature[iMarker];
          Surface_TotalPressure_Local[iMarker_Analyze]     += Surface_TotalPressure[iMarker];
          Surface_Area_Local[iMarker_Analyze]              += Surface_Area[iMarker];
          Surface_MassFlow_Abs_Local[iMarker_Analyze]      += Surface_MassFlow_Abs[iMarker];
        }

      }

    }

  }

#ifdef HAVE_MPI
  if (config->GetComm_Level() == COMM_FULL) {
    SU2_MPI::Allreduce(Surface_MassFlow_Local, Surface_MassFlow_Total, nMarker_Analyze, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());
    SU2_MPI::Allreduce(Surface_Mach_Local, Surface_Mach_Total, nMarker_Analyze, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());
    SU2_MPI::Allreduce(Surface_Temperature_Local, Surface_Temperature_Total, nMarker_Analyze, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());
    SU2_MPI::Allreduce(Surface_Density_Local, Surface_Density_Total, nMarker_Analyze, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());
    SU2_MPI::Allreduce(Surface_Enthalpy_Local, Surface_Enthalpy_Total, nMarker_Analyze, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());
    SU2_MPI::Allreduce(Surface_NormalVelocity_Local, Surface_NormalVelocity_Total, nMarker_Analyze, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());
    SU2_MPI::Allreduce(Surface_StreamVelocity2_Local, Surface_StreamVelocity2_Total, nMarker_Analyze, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());
    SU2_MPI::Allreduce(Surface_TransvVelocity2_Local, Surface_TransvVelocity2_Total, nMarker_Analyze, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());
    SU2_MPI::Allreduce(Surface_Pressure_Local, Surface_Pressure_Total, nMarker_Analyze, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());
    SU2_MPI::Allreduce(Surface_TotalTemperature_Local, Surface_TotalTemperature_Total, nMarker_Analyze, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());
    SU2_MPI::Allreduce(Surface_TotalPressure_Local, Surface_TotalPressure_Total, nMarker_Analyze, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());
    SU2_MPI::Allreduce(Surface_Area_Local, Surface_Area_Total, nMarker_Analyze, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());
    SU2_MPI::Allreduce(Surface_MassFlow_Abs_Local, Surface_MassFlow_Abs_Total, nMarker_Analyze, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());
  }
#else

  for (iMarker_Analyze = 0; iMarker_Analyze < nMarker_Analyze; iMarker_Analyze++) {
    Surface_MassFlow_Total[iMarker_Analyze]          = Surface_MassFlow_Local[iMarker_Analyze];
    Surface_Mach_Total[iMarker_Analyze]              = Surface_Mach_Local[iMarker_Analyze];
    Surface_Temperature_Total[iMarker_Analyze]       = Surface_Temperature_Local[iMarker_Analyze];
    Surface_Density_Total[iMarker_Analyze]           = Surface_Density_Local[iMarker_Analyze];
    Surface_Enthalpy_Total[iMarker_Analyze]          = Surface_Enthalpy_Local[iMarker_Analyze];
    Surface_NormalVelocity_Total[iMarker_Analyze]    = Surface_NormalVelocity_Local[iMarker_Analyze];
    Surface_StreamVelocity2_Total[iMarker_Analyze]   = Surface_StreamVelocity2_Local[iMarker_Analyze];
    Surface_TransvVelocity2_Total[iMarker_Analyze]   = Surface_TransvVelocity2_Local[iMarker_Analyze];
    Surface_Pressure_Total[iMarker_Analyze]          = Surface_Pressure_Local[iMarker_Analyze];
    Surface_TotalTemperature_Total[iMarker_Analyze]  = Surface_TotalTemperature_Local[iMarker_Analyze];
    Surface_TotalPressure_Total[iMarker_Analyze]     = Surface_TotalPressure_Local[iMarker_Analyze];
    Surface_Area_Total[iMarker_Analyze]              = Surface_Area_Local[iMarker_Analyze];
    Surface_MassFlow_Abs_Total[iMarker_Analyze]      = Surface_MassFlow_Abs_Local[iMarker_Analyze];
  }

#endif

  /*--- Compute the value of Surface_Area_Total, and Surface_Pressure_Total, and
   set the value in the config structure for future use ---*/

  for (iMarker_Analyze = 0; iMarker_Analyze < nMarker_Analyze; iMarker_Analyze++) {

    if (Kind_Average == AVERAGE_MASSFLUX) Weight = Surface_MassFlow_Abs_Total[iMarker_Analyze];
    else if (Kind_Average == AVERAGE_AREA) Weight = abs(Surface_Area_Total[iMarker_Analyze]);
    else Weight = 1.0;

    if (Weight != 0.0) {
      Surface_Mach_Total[iMarker_Analyze]             /= Weight;
      Surface_Temperature_Total[iMarker_Analyze]      /= Weight;
      Surface_Density_Total[iMarker_Analyze]          /= Weight;
      Surface_Enthalpy_Total[iMarker_Analyze]         /= Weight;
      Surface_NormalVelocity_Total[iMarker_Analyze]   /= Weight;
      Surface_Pressure_Total[iMarker_Analyze]         /= Weight;
      Surface_TotalTemperature_Total[iMarker_Analyze] /= Weight;
      Surface_TotalPressure_Total[iMarker_Analyze]    /= Weight;
    }
    else {
      Surface_Mach_Total[iMarker_Analyze]             = 0.0;
      Surface_Temperature_Total[iMarker_Analyze]      = 0.0;
      Surface_Density_Total[iMarker_Analyze]          = 0.0;
      Surface_Enthalpy_Total[iMarker_Analyze]         = 0.0;
      Surface_NormalVelocity_Total[iMarker_Analyze]   = 0.0;
      Surface_Pressure_Total[iMarker_Analyze]         = 0.0;
      Surface_TotalTemperature_Total[iMarker_Analyze] = 0.0;
      Surface_TotalPressure_Total[iMarker_Analyze]    = 0.0;
    }

    /*--- Compute flow uniformity parameters separately (always area for now). ---*/

    Area = fabs(Surface_Area_Total[iMarker_Analyze]);

    if (Area != 0.0) {
      Surface_MomentumDistortion_Total[iMarker_Analyze] = Surface_StreamVelocity2_Total[iMarker_Analyze]/(Surface_NormalVelocity_Total[iMarker_Analyze]*Surface_NormalVelocity_Total[iMarker_Analyze]*Area) - 1.0;
      Surface_StreamVelocity2_Total[iMarker_Analyze] /= Area;
      Surface_TransvVelocity2_Total[iMarker_Analyze] /= Area;
    }
    else {
      Surface_MomentumDistortion_Total[iMarker_Analyze] = 0.0;
      Surface_StreamVelocity2_Total[iMarker_Analyze]    = 0.0;
      Surface_TransvVelocity2_Total[iMarker_Analyze]    = 0.0;
    }

  }

  for (iMarker_Analyze = 0; iMarker_Analyze < nMarker_Analyze; iMarker_Analyze++) {

    su2double MassFlow = Surface_MassFlow_Total[iMarker_Analyze] * config->GetDensity_Ref() * config->GetVelocity_Ref();
    if (config->GetSystemMeasurements() == US) MassFlow *= 32.174;
    config->SetSurface_MassFlow(iMarker_Analyze, MassFlow);

    su2double Mach = Surface_Mach_Total[iMarker_Analyze];
    config->SetSurface_Mach(iMarker_Analyze, Mach);

    su2double Temperature = Surface_Temperature_Total[iMarker_Analyze] * config->GetTemperature_Ref();
    config->SetSurface_Temperature(iMarker_Analyze, Temperature);

    su2double Pressure = Surface_Pressure_Total[iMarker_Analyze] * config->GetPressure_Ref();
    config->SetSurface_Pressure(iMarker_Analyze, Pressure);

    su2double Density = Surface_Density_Total[iMarker_Analyze] * config->GetDensity_Ref();
    config->SetSurface_Density(iMarker_Analyze, Density);

    su2double Enthalpy = Surface_Enthalpy_Total[iMarker_Analyze];
    config->SetSurface_Enthalpy(iMarker_Analyze, Enthalpy);

    su2double NormalVelocity = Surface_NormalVelocity_Total[iMarker_Analyze] * config->GetVelocity_Ref();
    config->SetSurface_NormalVelocity(iMarker_Analyze, NormalVelocity);

    su2double Uniformity = sqrt(Surface_StreamVelocity2_Total[iMarker_Analyze]) * config->GetVelocity_Ref();
    config->SetSurface_Uniformity(iMarker_Analyze, Uniformity);

    su2double SecondaryStrength = sqrt(Surface_TransvVelocity2_Total[iMarker_Analyze]) * config->GetVelocity_Ref();
    config->SetSurface_SecondaryStrength(iMarker_Analyze, SecondaryStrength);

    su2double MomentumDistortion = Surface_MomentumDistortion_Total[iMarker_Analyze];
    config->SetSurface_MomentumDistortion(iMarker_Analyze, MomentumDistortion);

    su2double SecondOverUniform = SecondaryStrength/Uniformity;
    config->SetSurface_SecondOverUniform(iMarker_Analyze, SecondOverUniform);

    su2double TotalTemperature = Surface_TotalTemperature_Total[iMarker_Analyze] * config->GetTemperature_Ref();
    config->SetSurface_TotalTemperature(iMarker_Analyze, TotalTemperature);

    su2double TotalPressure = Surface_TotalPressure_Total[iMarker_Analyze] * config->GetPressure_Ref();
    config->SetSurface_TotalPressure(iMarker_Analyze, TotalPressure);

  }

  /*--- Compute the average static pressure drop between two surfaces. Note
   that this assumes we have two surfaces being analyzed and that the outlet
   is first followed by the inlet. This is because we may also want to choose
   outlet values (temperature, uniformity, etc.) for our design problems,
   which require the outlet to be listed first. This is a simple first version
   that could be generalized to a different orders/lists/etc. ---*/

  for (iMarker_Analyze = 0; iMarker_Analyze < nMarker_Analyze; iMarker_Analyze++) {
    if (nMarker_Analyze == 2) {
      su2double Pressure_Drop = (Surface_Pressure_Total[1]-Surface_Pressure_Total[0]) * config->GetPressure_Ref();
      config->SetSurface_PressureDrop(iMarker_Analyze, Pressure_Drop);
    } else {
      config->SetSurface_PressureDrop(iMarker_Analyze, 0.0);
    }
  }

  if ((rank == MASTER_NODE) && !config->GetDiscrete_Adjoint() && output) {

    cout.precision(6);
    cout.setf(ios::scientific, ios::floatfield);
    cout << endl << "Computing surface mean values." << endl << endl;

    for (iMarker_Analyze = 0; iMarker_Analyze < nMarker_Analyze; iMarker_Analyze++) {
      cout << "Surface "<< config->GetMarker_Analyze_TagBound(iMarker_Analyze) << ":" << endl;

      if (nDim == 3) { if (config->GetSystemMeasurements() == SI) cout << setw(20) << "Area (m^2): "; else cout << setw(20) << "Area (ft^2): "; }
      else { if (config->GetSystemMeasurements() == SI) cout << setw(20) << "Area (m): "; else cout << setw(20) << "Area (ft): "; }

      if (config->GetSystemMeasurements() == SI)      cout << setw(15) << fabs(Surface_Area_Total[iMarker_Analyze]);
      else if (config->GetSystemMeasurements() == US) cout << setw(15) << fabs(Surface_Area_Total[iMarker_Analyze])*12.0*12.0;

      cout << endl;

      su2double MassFlow = config->GetSurface_MassFlow(iMarker_Analyze);
      if (config->GetSystemMeasurements() == SI)      cout << setw(20) << "Mf (kg/s): " << setw(15) << MassFlow;
      else if (config->GetSystemMeasurements() == US) cout << setw(20) << "Mf (lbs/s): " << setw(15) << MassFlow;

      su2double NormalVelocity = config->GetSurface_NormalVelocity(iMarker_Analyze);
      if (config->GetSystemMeasurements() == SI)      cout << setw(20) << "Vn (m/s): " << setw(15) << NormalVelocity;
      else if (config->GetSystemMeasurements() == US) cout << setw(20) << "Vn (ft/s): " << setw(15) << NormalVelocity;

      cout << endl;

      su2double Uniformity = config->GetSurface_Uniformity(iMarker_Analyze);
      if (config->GetSystemMeasurements() == SI)      cout << setw(20) << "Uniformity (m/s): " << setw(15) << Uniformity;
      else if (config->GetSystemMeasurements() == US) cout << setw(20) << "Uniformity (ft/s): " << setw(15) << Uniformity;

      su2double SecondaryStrength = config->GetSurface_SecondaryStrength(iMarker_Analyze);
      if (config->GetSystemMeasurements() == SI)      cout << setw(20) << "Secondary (m/s): " << setw(15) << SecondaryStrength;
      else if (config->GetSystemMeasurements() == US) cout << setw(20) << "Secondary (ft/s): " << setw(15) << SecondaryStrength;

      cout << endl;

      su2double MomentumDistortion = config->GetSurface_MomentumDistortion(iMarker_Analyze);
      cout << setw(20) << "Mom. Distortion: " << setw(15) << MomentumDistortion;

      su2double SecondOverUniform = config->GetSurface_SecondOverUniform(iMarker_Analyze);
      cout << setw(20) << "Second/Uniform: " << setw(15) << SecondOverUniform;

      cout << endl;

      su2double Pressure = config->GetSurface_Pressure(iMarker_Analyze);
      if (config->GetSystemMeasurements() == SI)      cout << setw(20) << "P (Pa): " << setw(15) << Pressure;
      else if (config->GetSystemMeasurements() == US) cout << setw(20) << "P (psf): " << setw(15) << Pressure;

      su2double TotalPressure = config->GetSurface_TotalPressure(iMarker_Analyze);
      if (config->GetSystemMeasurements() == SI)      cout << setw(20) << "PT (Pa): " << setw(15) <<TotalPressure;
      else if (config->GetSystemMeasurements() == US) cout << setw(20) << "PT (psf): " << setw(15) <<TotalPressure;

      cout << endl;

      su2double Mach = config->GetSurface_Mach(iMarker_Analyze);
      cout << setw(20) << "Mach: " << setw(15) << Mach;

      su2double Density = config->GetSurface_Density(iMarker_Analyze);
      if (config->GetSystemMeasurements() == SI)      cout << setw(20) << "Rho (kg/m^3): " << setw(15) << Density;
      else if (config->GetSystemMeasurements() == US) cout << setw(20) << "Rho (lb/ft^3): " << setw(15) << Density*32.174;

      cout << endl;

      if (compressible || energy) {
        su2double Temperature = config->GetSurface_Temperature(iMarker_Analyze);
        if (config->GetSystemMeasurements() == SI)      cout << setw(20) << "T (K): " << setw(15) << Temperature;
        else if (config->GetSystemMeasurements() == US) cout << setw(20) << "T (R): " << setw(15) << Temperature;

        su2double TotalTemperature = config->GetSurface_TotalTemperature(iMarker_Analyze);
        if (config->GetSystemMeasurements() == SI)      cout << setw(20) << "TT (K): " << setw(15) << TotalTemperature;
        else if (config->GetSystemMeasurements() == US) cout << setw(20) << "TT (R): " << setw(15) << TotalTemperature;

        cout << endl;
      }

    }
    cout.unsetf(ios_base::floatfield);

  }

  delete [] Surface_MassFlow_Local;
  delete [] Surface_Mach_Local;
  delete [] Surface_Temperature_Local;
  delete [] Surface_Density_Local;
  delete [] Surface_Enthalpy_Local;
  delete [] Surface_NormalVelocity_Local;
  delete [] Surface_Pressure_Local;
  delete [] Surface_TotalTemperature_Local;
  delete [] Surface_TotalPressure_Local;
  delete [] Surface_Area_Local;
  delete [] Surface_MassFlow_Abs_Local;

  delete [] Surface_MassFlow_Total;
  delete [] Surface_Mach_Total;
  delete [] Surface_Temperature_Total;
  delete [] Surface_Density_Total;
  delete [] Surface_Enthalpy_Total;
  delete [] Surface_NormalVelocity_Total;
  delete [] Surface_Pressure_Total;
  delete [] Surface_TotalTemperature_Total;
  delete [] Surface_TotalPressure_Total;
  delete [] Surface_Area_Total;
  delete [] Surface_MassFlow_Abs_Total;

  delete [] Surface_MassFlow;
  delete [] Surface_Mach;
  delete [] Surface_Temperature;
  delete [] Surface_Density;
  delete [] Surface_Enthalpy;
  delete [] Surface_NormalVelocity;
  delete [] Surface_Pressure;
  delete [] Surface_TotalTemperature;
  delete [] Surface_TotalPressure;
  delete [] Surface_Area;
  delete [] Vector;
  delete [] Surface_VelocityIdeal;
  delete [] Surface_MassFlow_Abs;

}
