/*!
 * \file output_structure.cpp
 * \brief Main subroutines for output solver information
 * \author F. Palacios, T. Economon
 * \version 7.0.4 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2020, SU2 Contributors (cf. AUTHORS.md)
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

  unsigned short iDim, iZone, iSpan, iMarker;

  /*--- Initialize point and connectivity counters to zero. ---*/

  nGlobal_Poin      = 0;
  nSurf_Poin        = 0;
  nGlobal_Elem      = 0;
  nSurf_Elem        = 0;
  nGlobal_Tria      = 0;
  nGlobal_Quad      = 0;
  nGlobal_Tetr      = 0;
  nGlobal_Hexa      = 0;
  nGlobal_Pris      = 0;
  nGlobal_Pyra      = 0;
  nGlobal_Line      = 0;
  nGlobal_BoundTria = 0;
  nGlobal_BoundQuad = 0;

  /*--- Initialize pointers to NULL ---*/

  Coords = NULL;
  Conn_Line = NULL;     Conn_BoundTria = NULL;  Conn_BoundQuad = NULL;
  Conn_Tria = NULL;     Conn_Quad = NULL;       Conn_Tetr = NULL;
  Conn_Hexa = NULL;     Conn_Pris = NULL;       Conn_Pyra = NULL;
  Data = NULL;

  /*--- Initialize parallel pointers to NULL ---*/

  nGlobal_Poin_Par    = 0;
  nGlobal_Elem_Par    = 0;
  nGlobal_Surf_Poin   = 0;
  nParallel_Poin      = 0;
  nSurf_Poin_Par      = 0;
  nSurf_Elem_Par      = 0;
  nParallel_Tria      = 0;
  nParallel_Quad      = 0;
  nParallel_Tetr      = 0;
  nParallel_Hexa      = 0;
  nParallel_Pris      = 0;
  nParallel_Pyra      = 0;
  nParallel_Line      = 0;
  nParallel_BoundTria = 0;
  nParallel_BoundQuad = 0;

  /*--- Initialize pointers to NULL ---*/

  Conn_BoundLine_Par = NULL;  Conn_BoundTria_Par = NULL;  Conn_BoundQuad_Par = NULL;
  Conn_Tria_Par = NULL;  Conn_Quad_Par = NULL;       Conn_Tetr_Par = NULL;
  Conn_Hexa_Par = NULL;  Conn_Pris_Par = NULL;       Conn_Pyra_Par = NULL;

  Local_Data         = NULL;
  Local_Data_Copy    = NULL;
  Parallel_Data      = NULL;
  Parallel_Surf_Data = NULL;

  /*--- Initialize structures for storing linear partitioning offsets ---*/

  nGlobalPoint_Sort = 0;
  nLocalPoint_Sort  = 0;
  nPoint_Restart    = 0;

  Local_Halo_Sort = NULL;

  beg_node = NULL;
  end_node = NULL;

  nPointLinear     = NULL;
  nPointCumulative = NULL;

  /*--- Initialize CGNS write flag ---*/

  wrote_base_file = false;

  /*--- Initialize CGNS write flag ---*/

  wrote_CGNS_base = false;

  /*--- Initialize Tecplot surface flag ---*/

  wrote_surf_file = false;

  /*--- Initialize Paraview write flag ---*/

  wrote_Paraview_base = false;

  /*--- Initialize residual ---*/

  RhoRes_New = EPS;
  RhoRes_Old = new su2double[config->GetnZone()];
  for (iZone = 0; iZone < config->GetnZone(); iZone++) RhoRes_Old[iZone] = EPS;

  wrote_Paraview_base = false;

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

  if (RhoRes_Old != NULL) delete [] RhoRes_Old;

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

void COutputLegacy::SetSurfaceCSV_Flow(CConfig *config, CGeometry *geometry,
                                 CSolver *FlowSolver, unsigned long iExtIter,
                                 unsigned short val_iZone, unsigned short val_iInst) {

  unsigned short iMarker;
  unsigned long iPoint, iVertex, Global_Index;
  su2double PressCoeff = 0.0, SkinFrictionCoeff[3], HeatFlux;
  su2double xCoord = 0.0, yCoord = 0.0, zCoord = 0.0, Mach, Pressure;
  char cstr[200];

  unsigned short solver = config->GetKind_Solver();
  unsigned short nDim = geometry->GetnDim();

#ifndef HAVE_MPI

  unsigned short iDim;
  char buffer [50];
  ofstream SurfFlow_file;

  /*--- Write file name with extension if unsteady ---*/
  strcpy (cstr, config->GetSurfCoeff_FileName().c_str());

  if (config->GetTime_Marching() == HARMONIC_BALANCE) {
    SPRINTF (buffer, "_%d.csv", SU2_TYPE::Int(val_iInst));

  }else if (config->GetTime_Marching() && config->GetTime_Domain()) {
    if ((SU2_TYPE::Int(iExtIter) >= 0)    && (SU2_TYPE::Int(iExtIter) < 10))    SPRINTF (buffer, "_0000%d.csv", SU2_TYPE::Int(iExtIter));
    if ((SU2_TYPE::Int(iExtIter) >= 10)   && (SU2_TYPE::Int(iExtIter) < 100))   SPRINTF (buffer, "_000%d.csv",  SU2_TYPE::Int(iExtIter));
    if ((SU2_TYPE::Int(iExtIter) >= 100)  && (SU2_TYPE::Int(iExtIter) < 1000))  SPRINTF (buffer, "_00%d.csv",   SU2_TYPE::Int(iExtIter));
    if ((SU2_TYPE::Int(iExtIter) >= 1000) && (SU2_TYPE::Int(iExtIter) < 10000)) SPRINTF (buffer, "_0%d.csv",    SU2_TYPE::Int(iExtIter));
    if (SU2_TYPE::Int(iExtIter) >= 10000) SPRINTF (buffer, "_%d.csv", SU2_TYPE::Int(iExtIter));
  }
  else
    SPRINTF (buffer, ".csv");

  strcat (cstr, buffer);
  SurfFlow_file.precision(15);
  SurfFlow_file.open(cstr, ios::out);

  SurfFlow_file << "\"Global_Index\", \"x_coord\", \"y_coord\", ";
  if (nDim == 3) SurfFlow_file << "\"z_coord\", ";
  SurfFlow_file << "\"Pressure\", \"Pressure_Coefficient\", ";

  switch (solver) {
    case EULER : case INC_EULER: SurfFlow_file <<  "\"Mach_Number\"" << "\n"; break;
    case NAVIER_STOKES: case RANS: case INC_NAVIER_STOKES: case INC_RANS:
      if (nDim == 2) SurfFlow_file <<  "\"Skin_Friction_Coefficient_X\", \"Skin_Friction_Coefficient_Y\", \"Heat_Flux\"" << "\n";
      if (nDim == 3) SurfFlow_file <<  "\"Skin_Friction_Coefficient_X\", \"Skin_Friction_Coefficient_Y\", \"Skin_Friction_Coefficient_Z\", \"Heat_Flux\"" << "\n";
      break;
  }

  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_Plotting(iMarker) == YES) {
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        Global_Index = geometry->node[iPoint]->GetGlobalIndex();
        xCoord = geometry->node[iPoint]->GetCoord(0);
        yCoord = geometry->node[iPoint]->GetCoord(1);
        if (nDim == 3) zCoord = geometry->node[iPoint]->GetCoord(2);

        /*--- The output should be in inches ---*/

        if (config->GetSystemMeasurements() == US) {
          xCoord *= 12.0; yCoord *= 12.0;
          if (nDim == 3) zCoord *= 12.0;
        }

        Pressure = FlowSolver->GetNodes()->GetPressure(iPoint);
        PressCoeff = FlowSolver->GetCPressure(iMarker, iVertex);
        SurfFlow_file << scientific << Global_Index << ", " << xCoord << ", " << yCoord << ", ";
        if (nDim == 3) SurfFlow_file << scientific << zCoord << ", ";
        SurfFlow_file << scientific << Pressure << ", " << PressCoeff << ", ";
        switch (solver) {
          case EULER : case FEM_EULER: case INC_EULER:
            Mach = sqrt(FlowSolver->GetNodes()->GetVelocity2(iPoint)) / FlowSolver->GetNodes()->GetSoundSpeed(iPoint);
            SurfFlow_file << scientific << Mach << "\n";
            break;
          case NAVIER_STOKES: case RANS: case FEM_NAVIER_STOKES: case FEM_RANS: case FEM_LES:
          case INC_NAVIER_STOKES: case INC_RANS:
            for (iDim = 0; iDim < nDim; iDim++)
              SkinFrictionCoeff[iDim] = FlowSolver->GetCSkinFriction(iMarker, iVertex, iDim);
            HeatFlux = FlowSolver->GetHeatFlux(iMarker, iVertex);

            if (nDim == 2) SurfFlow_file << scientific << SkinFrictionCoeff[0] << ", " << SkinFrictionCoeff[1] << ", " << HeatFlux << "\n";
            if (nDim == 3) SurfFlow_file << scientific << SkinFrictionCoeff[0] << ", " << SkinFrictionCoeff[1] << ", " << SkinFrictionCoeff[2] << ", " << HeatFlux << "\n";

            break;
        }
      }
    }
  }

  SurfFlow_file.close();

#else

  int iProcessor, nProcessor = size;

  unsigned long Buffer_Send_nVertex[1], *Buffer_Recv_nVertex = NULL;
  unsigned long nVertex_Surface = 0, nLocalVertex_Surface = 0;
  unsigned long MaxLocalVertex_Surface = 0;

  /*--- Find the max number of surface vertices among all
   partitions and set up buffers. The master node will handle the
   writing of the CSV file after gathering all of the data. ---*/

  nLocalVertex_Surface = 0;
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
    if (config->GetMarker_All_Plotting(iMarker) == YES)
      for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        if (geometry->node[iPoint]->GetDomain()) nLocalVertex_Surface++;
      }

  /*--- Communicate the number of local vertices on each partition
   to the master node ---*/

  Buffer_Send_nVertex[0] = nLocalVertex_Surface;
  if (rank == MASTER_NODE) Buffer_Recv_nVertex = new unsigned long [nProcessor];

  SU2_MPI::Allreduce(&nLocalVertex_Surface, &MaxLocalVertex_Surface, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
  SU2_MPI::Gather(&Buffer_Send_nVertex, 1, MPI_UNSIGNED_LONG, Buffer_Recv_nVertex, 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);

  /*--- Send and Recv buffers ---*/

  su2double *Buffer_Send_Coord_x = new su2double [MaxLocalVertex_Surface];
  su2double *Buffer_Recv_Coord_x = NULL;

  su2double *Buffer_Send_Coord_y = new su2double [MaxLocalVertex_Surface];
  su2double *Buffer_Recv_Coord_y = NULL;

  su2double *Buffer_Send_Coord_z = new su2double [MaxLocalVertex_Surface];
  su2double *Buffer_Recv_Coord_z = NULL;

  su2double *Buffer_Send_Press = new su2double [MaxLocalVertex_Surface];
  su2double *Buffer_Recv_Press = NULL;

  su2double *Buffer_Send_CPress = new su2double [MaxLocalVertex_Surface];
  su2double *Buffer_Recv_CPress = NULL;

  su2double *Buffer_Send_Mach = new su2double [MaxLocalVertex_Surface];
  su2double *Buffer_Recv_Mach = NULL;

  su2double *Buffer_Send_SkinFriction_x = new su2double [MaxLocalVertex_Surface];
  su2double *Buffer_Recv_SkinFriction_x = NULL;

  su2double *Buffer_Send_SkinFriction_y = new su2double [MaxLocalVertex_Surface];
  su2double *Buffer_Recv_SkinFriction_y = NULL;

  su2double *Buffer_Send_SkinFriction_z = new su2double [MaxLocalVertex_Surface];
  su2double *Buffer_Recv_SkinFriction_z = NULL;

  su2double *Buffer_Send_HeatTransfer = new su2double [MaxLocalVertex_Surface];
  su2double *Buffer_Recv_HeatTransfer = NULL;

  unsigned long *Buffer_Send_GlobalIndex = new unsigned long [MaxLocalVertex_Surface];
  unsigned long *Buffer_Recv_GlobalIndex = NULL;

  /*--- Prepare the receive buffers on the master node only. ---*/

  if (rank == MASTER_NODE) {
    Buffer_Recv_Coord_x = new su2double [nProcessor*MaxLocalVertex_Surface];
    Buffer_Recv_Coord_y = new su2double [nProcessor*MaxLocalVertex_Surface];
    if (nDim == 3) Buffer_Recv_Coord_z = new su2double [nProcessor*MaxLocalVertex_Surface];
    Buffer_Recv_Press   = new su2double [nProcessor*MaxLocalVertex_Surface];
    Buffer_Recv_CPress  = new su2double [nProcessor*MaxLocalVertex_Surface];
    Buffer_Recv_Mach    = new su2double [nProcessor*MaxLocalVertex_Surface];
    Buffer_Recv_SkinFriction_x = new su2double [nProcessor*MaxLocalVertex_Surface];
    Buffer_Recv_SkinFriction_y = new su2double [nProcessor*MaxLocalVertex_Surface];
    if (nDim == 3) Buffer_Recv_SkinFriction_z = new su2double [nProcessor*MaxLocalVertex_Surface];
    Buffer_Recv_HeatTransfer = new su2double [nProcessor*MaxLocalVertex_Surface];
    Buffer_Recv_GlobalIndex  = new unsigned long [nProcessor*MaxLocalVertex_Surface];
  }

  /*--- Loop over all vertices in this partition and load the
   data of the specified type into the buffer to be sent to
   the master node. ---*/

  nVertex_Surface = 0;
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
    if (config->GetMarker_All_Plotting(iMarker) == YES)
      for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        if (geometry->node[iPoint]->GetDomain()) {
          Buffer_Send_Press[nVertex_Surface] = FlowSolver->GetNodes()->GetPressure(iPoint);
          Buffer_Send_CPress[nVertex_Surface] = FlowSolver->GetCPressure(iMarker, iVertex);
          Buffer_Send_Coord_x[nVertex_Surface] = geometry->node[iPoint]->GetCoord(0);
          Buffer_Send_Coord_y[nVertex_Surface] = geometry->node[iPoint]->GetCoord(1);
          if (nDim == 3) { Buffer_Send_Coord_z[nVertex_Surface] = geometry->node[iPoint]->GetCoord(2); }

          /*--- If US system, the output should be in inches ---*/

          if (config->GetSystemMeasurements() == US) {
            Buffer_Send_Coord_x[nVertex_Surface] *= 12.0;
            Buffer_Send_Coord_y[nVertex_Surface] *= 12.0;
            if (nDim == 3) Buffer_Send_Coord_z[nVertex_Surface] *= 12.0;
          }

          Buffer_Send_GlobalIndex[nVertex_Surface] = geometry->node[iPoint]->GetGlobalIndex();

          if (solver == EULER || solver == FEM_EULER || solver == INC_EULER)
            Buffer_Send_Mach[nVertex_Surface] = sqrt(FlowSolver->GetNodes()->GetVelocity2(iPoint)) / FlowSolver->GetNodes()->GetSoundSpeed(iPoint);

          if (solver == NAVIER_STOKES     || solver == RANS     || solver == INC_NAVIER_STOKES || solver == INC_RANS ||
              solver == FEM_NAVIER_STOKES || solver == FEM_RANS || solver == FEM_LES) {
            Buffer_Send_SkinFriction_x[nVertex_Surface] = FlowSolver->GetCSkinFriction(iMarker, iVertex, 0);
            Buffer_Send_SkinFriction_y[nVertex_Surface] = FlowSolver->GetCSkinFriction(iMarker, iVertex, 1);
            if (nDim == 3) Buffer_Send_SkinFriction_z[nVertex_Surface] = FlowSolver->GetCSkinFriction(iMarker, iVertex, 2);
            Buffer_Send_HeatTransfer[nVertex_Surface] = FlowSolver->GetHeatFlux(iMarker, iVertex);
          }
          nVertex_Surface++;
        }
      }

  /*--- Send the information to the master node ---*/

  SU2_MPI::Gather(Buffer_Send_Coord_x, MaxLocalVertex_Surface, MPI_DOUBLE, Buffer_Recv_Coord_x, MaxLocalVertex_Surface, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  SU2_MPI::Gather(Buffer_Send_Coord_y, MaxLocalVertex_Surface, MPI_DOUBLE, Buffer_Recv_Coord_y, MaxLocalVertex_Surface, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  if (nDim == 3) SU2_MPI::Gather(Buffer_Send_Coord_z, MaxLocalVertex_Surface, MPI_DOUBLE, Buffer_Recv_Coord_z, MaxLocalVertex_Surface, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  SU2_MPI::Gather(Buffer_Send_Press, MaxLocalVertex_Surface, MPI_DOUBLE, Buffer_Recv_Press, MaxLocalVertex_Surface, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  SU2_MPI::Gather(Buffer_Send_CPress, MaxLocalVertex_Surface, MPI_DOUBLE, Buffer_Recv_CPress, MaxLocalVertex_Surface, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  if (solver == EULER || solver == FEM_EULER || solver == INC_EULER) SU2_MPI::Gather(Buffer_Send_Mach, MaxLocalVertex_Surface, MPI_DOUBLE, Buffer_Recv_Mach, MaxLocalVertex_Surface, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  if ((solver == NAVIER_STOKES) || (solver == RANS) || (solver == FEM_NAVIER_STOKES) || (solver == FEM_RANS) || (solver == FEM_LES) || (solver == INC_RANS) || (solver == INC_NAVIER_STOKES)) {
    SU2_MPI::Gather(Buffer_Send_SkinFriction_x, MaxLocalVertex_Surface, MPI_DOUBLE, Buffer_Recv_SkinFriction_x, MaxLocalVertex_Surface, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    SU2_MPI::Gather(Buffer_Send_SkinFriction_y, MaxLocalVertex_Surface, MPI_DOUBLE, Buffer_Recv_SkinFriction_y, MaxLocalVertex_Surface, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    if (nDim == 3) SU2_MPI::Gather(Buffer_Send_SkinFriction_z, MaxLocalVertex_Surface, MPI_DOUBLE, Buffer_Recv_SkinFriction_z, MaxLocalVertex_Surface, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    SU2_MPI::Gather(Buffer_Send_HeatTransfer, MaxLocalVertex_Surface, MPI_DOUBLE, Buffer_Recv_HeatTransfer, MaxLocalVertex_Surface, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  }
  SU2_MPI::Gather(Buffer_Send_GlobalIndex, MaxLocalVertex_Surface, MPI_UNSIGNED_LONG, Buffer_Recv_GlobalIndex, MaxLocalVertex_Surface, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);

  /*--- The master node unpacks the data and writes the surface CSV file ---*/

  if (rank == MASTER_NODE) {

    /*--- Write file name with extension if unsteady ---*/
    char buffer[50];
    string filename = config->GetSurfCoeff_FileName();
    ofstream SurfFlow_file;

    /*--- Write file name with extension if unsteady ---*/
    strcpy (cstr, filename.c_str());
    if (config->GetTime_Marching() == HARMONIC_BALANCE) {
      SPRINTF (buffer, "_%d.csv", SU2_TYPE::Int(val_iInst));

    } else if (config->GetTime_Marching() && config->GetTime_Domain()) {
      if ((SU2_TYPE::Int(iExtIter) >= 0)    && (SU2_TYPE::Int(iExtIter) < 10))    SPRINTF (buffer, "_0000%d.csv", SU2_TYPE::Int(iExtIter));
      if ((SU2_TYPE::Int(iExtIter) >= 10)   && (SU2_TYPE::Int(iExtIter) < 100))   SPRINTF (buffer, "_000%d.csv",  SU2_TYPE::Int(iExtIter));
      if ((SU2_TYPE::Int(iExtIter) >= 100)  && (SU2_TYPE::Int(iExtIter) < 1000))  SPRINTF (buffer, "_00%d.csv",   SU2_TYPE::Int(iExtIter));
      if ((SU2_TYPE::Int(iExtIter) >= 1000) && (SU2_TYPE::Int(iExtIter) < 10000)) SPRINTF (buffer, "_0%d.csv",    SU2_TYPE::Int(iExtIter));
      if (SU2_TYPE::Int(iExtIter) >= 10000) SPRINTF (buffer, "_%d.csv", SU2_TYPE::Int(iExtIter));
    }
    else
      SPRINTF (buffer, ".csv");

    strcat (cstr, buffer);
    SurfFlow_file.precision(15);
    SurfFlow_file.open(cstr, ios::out);

    SurfFlow_file << "\"Global_Index\", \"x_coord\", \"y_coord\", ";
    if (nDim == 3) SurfFlow_file << "\"z_coord\", ";
    SurfFlow_file << "\"Pressure\", \"Pressure_Coefficient\", ";

    switch (solver) {
      case EULER : case FEM_EULER: case INC_EULER: SurfFlow_file <<  "\"Mach_Number\"" << "\n"; break;
      case NAVIER_STOKES: case RANS: case FEM_NAVIER_STOKES: case FEM_RANS: case FEM_LES: case INC_NAVIER_STOKES: case INC_RANS:
        if (nDim == 2) SurfFlow_file << "\"Skin_Friction_Coefficient_X\", \"Skin_Friction_Coefficient_Y\", \"Heat_Flux\"" << "\n";
        if (nDim == 3) SurfFlow_file << "\"Skin_Friction_Coefficient_X\", \"Skin_Friction_Coefficient_Y\", \"Skin_Friction_Coefficient_Z\", \"Heat_Flux\"" << "\n";
        break;
    }

    /*--- Loop through all of the collected data and write each node's values ---*/

    unsigned long Total_Index;
    for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
      for (iVertex = 0; iVertex < Buffer_Recv_nVertex[iProcessor]; iVertex++) {

        /*--- Current index position and global index ---*/
        Total_Index  = iProcessor*MaxLocalVertex_Surface+iVertex;
        Global_Index = Buffer_Recv_GlobalIndex[Total_Index];

        /*--- Retrieve the merged data for this node ---*/
        xCoord = Buffer_Recv_Coord_x[Total_Index];
        yCoord = Buffer_Recv_Coord_y[Total_Index];
        if (nDim == 3) zCoord = Buffer_Recv_Coord_z[Total_Index];
        Pressure   = Buffer_Recv_Press[Total_Index];
        PressCoeff = Buffer_Recv_CPress[Total_Index];

        /*--- Write the first part of the data ---*/
        SurfFlow_file << scientific << Global_Index << ", " << xCoord << ", " << yCoord << ", ";
        if (nDim == 3) SurfFlow_file << scientific << zCoord << ", ";
        SurfFlow_file << scientific << Pressure << ", " << PressCoeff << ", ";

        /*--- Write the solver-dependent part of the data ---*/
        switch (solver) {
          case EULER :
          case FEM_EULER:
          case INC_EULER:
            Mach = Buffer_Recv_Mach[Total_Index];
            SurfFlow_file << scientific << Mach << "\n";
            break;
          case NAVIER_STOKES: case RANS:
          case FEM_NAVIER_STOKES: case FEM_RANS: case FEM_LES:
          case INC_NAVIER_STOKES: case INC_RANS:
            SkinFrictionCoeff[0] = Buffer_Recv_SkinFriction_x[Total_Index];
            SkinFrictionCoeff[1] = Buffer_Recv_SkinFriction_y[Total_Index];
            if (nDim == 3) SkinFrictionCoeff[2] = Buffer_Recv_SkinFriction_z[Total_Index];
            HeatFlux = Buffer_Recv_HeatTransfer[Total_Index];
            if (nDim == 2) SurfFlow_file << scientific << SkinFrictionCoeff[0] << ", " << SkinFrictionCoeff[1] << ", " << HeatFlux << "\n";
            if (nDim == 3) SurfFlow_file << scientific << SkinFrictionCoeff[0] << ", " << SkinFrictionCoeff[1] << ", " << SkinFrictionCoeff[2] << ", " << HeatFlux << "\n";
            break;
        }
      }
    }

    /*--- Close the CSV file ---*/
    SurfFlow_file.close();

    /*--- Release the recv buffers on the master node ---*/

    delete [] Buffer_Recv_Coord_x;
    delete [] Buffer_Recv_Coord_y;
    if (nDim == 3) delete [] Buffer_Recv_Coord_z;
    delete [] Buffer_Recv_Press;
    delete [] Buffer_Recv_CPress;
    delete [] Buffer_Recv_Mach;
    delete [] Buffer_Recv_SkinFriction_x;
    delete [] Buffer_Recv_SkinFriction_y;
    if (nDim == 3) delete [] Buffer_Recv_SkinFriction_z;
    delete [] Buffer_Recv_HeatTransfer;
    delete [] Buffer_Recv_GlobalIndex;

    delete [] Buffer_Recv_nVertex;

  }

  /*--- Release the memory for the remaining buffers and exit ---*/

  delete [] Buffer_Send_Coord_x;
  delete [] Buffer_Send_Coord_y;
  delete [] Buffer_Send_Coord_z;
  delete [] Buffer_Send_Press;
  delete [] Buffer_Send_CPress;
  delete [] Buffer_Send_Mach;
  delete [] Buffer_Send_SkinFriction_x;
  delete [] Buffer_Send_SkinFriction_y;
  delete [] Buffer_Send_SkinFriction_z;
  delete [] Buffer_Send_HeatTransfer;
  delete [] Buffer_Send_GlobalIndex;

#endif

}

void COutputLegacy::SetSurfaceCSV_Adjoint(CConfig *config, CGeometry *geometry, CSolver *AdjSolver, CSolver *FlowSolution, unsigned long iExtIter, unsigned short val_iZone, unsigned short val_iInst) {

#ifndef HAVE_MPI

  unsigned long iPoint, iVertex, Global_Index;
  su2double *Solution, xCoord, yCoord, zCoord;
  unsigned short iMarker;
  char cstr[200], buffer[50];
  ofstream SurfAdj_file;

  /*--- Write file name with extension if unsteady ---*/

  strcpy (cstr, config->GetSurfAdjCoeff_FileName().c_str());

  if (config->GetTime_Marching() == HARMONIC_BALANCE) {
    SPRINTF (buffer, "_%d.csv", SU2_TYPE::Int(val_iInst));

  } else if (config->GetTime_Marching() && config->GetTime_Domain()) {
    if ((SU2_TYPE::Int(iExtIter) >= 0)    && (SU2_TYPE::Int(iExtIter) < 10))    SPRINTF (buffer, "_0000%d.csv", SU2_TYPE::Int(iExtIter));
    if ((SU2_TYPE::Int(iExtIter) >= 10)   && (SU2_TYPE::Int(iExtIter) < 100))   SPRINTF (buffer, "_000%d.csv",  SU2_TYPE::Int(iExtIter));
    if ((SU2_TYPE::Int(iExtIter) >= 100)  && (SU2_TYPE::Int(iExtIter) < 1000))  SPRINTF (buffer, "_00%d.csv",   SU2_TYPE::Int(iExtIter));
    if ((SU2_TYPE::Int(iExtIter) >= 1000) && (SU2_TYPE::Int(iExtIter) < 10000)) SPRINTF (buffer, "_0%d.csv",    SU2_TYPE::Int(iExtIter));
    if (SU2_TYPE::Int(iExtIter) >= 10000) SPRINTF (buffer, "_%d.csv", SU2_TYPE::Int(iExtIter));
  }
  else
    SPRINTF (buffer, ".csv");

  strcat(cstr, buffer);
  SurfAdj_file.precision(15);
  SurfAdj_file.open(cstr, ios::out);

  SurfAdj_file << "SENS_AOA=" << AdjSolver->GetTotal_Sens_AoA() * PI_NUMBER / 180.0 << endl;

  if (geometry->GetnDim() == 2) {
    if (config ->GetKind_Regime() == COMPRESSIBLE)
      SurfAdj_file <<  "\"Point\",\"Sensitivity\",\"PsiRho\",\"Phi_x\",\"Phi_y\",\"PsiE\",\"x_coord\",\"y_coord\"";
    else if (config ->GetKind_Regime() == INCOMPRESSIBLE)
      SurfAdj_file <<  "\"Point\",\"Sensitivity\",\"PsiRho\",\"Phi_x\",\"Phi_y\",\"x_coord\",\"y_coord\"";

    if (config->GetDiscrete_Adjoint()) {
      SurfAdj_file << ",\"x_Sens\",\"y_Sens\"";
    }
    SurfAdj_file << "\n";

    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if (config->GetMarker_All_Plotting(iMarker) == YES)
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          Global_Index = geometry->node[iPoint]->GetGlobalIndex();
          Solution = AdjSolver->GetNodes()->GetSolution(iPoint);
          xCoord = geometry->node[iPoint]->GetCoord(0);
          yCoord = geometry->node[iPoint]->GetCoord(1);

          /*--- If US system, the output should be in inches ---*/

          if (config->GetSystemMeasurements() == US) {
            xCoord *= 12.0;
            yCoord *= 12.0;
          }
          if (config ->GetKind_Regime() == COMPRESSIBLE)
            SurfAdj_file << scientific << Global_Index << ", " << AdjSolver->GetCSensitivity(iMarker, iVertex) << ", " << Solution[0] << ", "
                         << Solution[1] << ", " << Solution[2] << ", " << Solution[3] <<", " << xCoord <<", "<< yCoord;
          else if (config ->GetKind_Regime() == INCOMPRESSIBLE)
              SurfAdj_file << scientific << Global_Index << ", " << AdjSolver->GetCSensitivity(iMarker, iVertex) << ", " << Solution[0] << ", "
                           << Solution[1] << ", " << Solution[2] <<", " << xCoord <<", "<< yCoord;
          if (config->GetDiscrete_Adjoint()) {
            SurfAdj_file << ", " << AdjSolver->GetNodes()->GetSensitivity(iPoint,0) << ", " << AdjSolver->GetNodes()->GetSensitivity(iPoint,1);
          }
          SurfAdj_file << "\n";
        }
    }
  }

  if (geometry->GetnDim() == 3) {
    if (config ->GetKind_Regime() == COMPRESSIBLE)
      SurfAdj_file <<  "\"Point\",\"Sensitivity\",\"PsiRho\",\"Phi_x\",\"Phi_y\",\"Phi_z\",\"PsiE\",\"x_coord\",\"y_coord\",\"z_coord\"";
    else if (config ->GetKind_Regime() == INCOMPRESSIBLE)
      SurfAdj_file <<  "\"Point\",\"Sensitivity\",\"PsiRho\",\"Phi_x\",\"Phi_y\",\"Phi_z\",\"x_coord\",\"y_coord\",\"z_coord\"";

    if (config->GetDiscrete_Adjoint()) {
      SurfAdj_file << ",\"x_Sens\",\"y_Sens\",\"z_Sens\"";
    }
    SurfAdj_file << "\n";
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if (config->GetMarker_All_Plotting(iMarker) == YES)
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          Global_Index = geometry->node[iPoint]->GetGlobalIndex();
          Solution = AdjSolver->GetNodes()->GetSolution(iPoint);

          xCoord = geometry->node[iPoint]->GetCoord(0);
          yCoord = geometry->node[iPoint]->GetCoord(1);
          zCoord = geometry->node[iPoint]->GetCoord(2);

          /*--- If US system, the output should be in inches ---*/

          if (config->GetSystemMeasurements() == US) {
            xCoord *= 12.0;
            yCoord *= 12.0;
            zCoord *= 12.0;
          }
          if (config ->GetKind_Regime() == COMPRESSIBLE)
            SurfAdj_file << scientific << Global_Index << ", " << AdjSolver->GetCSensitivity(iMarker, iVertex) << ", " << Solution[0] << ", "
                         << Solution[1] << ", " << Solution[2] << ", " << Solution[3] << ", " << Solution[4] << ", "<< xCoord <<", "<< yCoord <<", "<< zCoord;
          else if (config ->GetKind_Regime() == INCOMPRESSIBLE)
            SurfAdj_file << scientific << Global_Index << ", " << AdjSolver->GetCSensitivity(iMarker, iVertex) << ", " << Solution[0] << ", "
                         << Solution[1] << ", " << Solution[2] << ", " << Solution[3] << ", " << xCoord <<", "<< yCoord <<", "<< zCoord;
          if (config->GetDiscrete_Adjoint()) {
            SurfAdj_file << ", " << AdjSolver->GetNodes()->GetSensitivity(iPoint,0) << ", " << AdjSolver->GetNodes()->GetSensitivity(iPoint,1)
            << ", " << AdjSolver->GetNodes()->GetSensitivity(iPoint, 2);
          }
          SurfAdj_file << "\n";
        }
    }
  }

  SurfAdj_file.close();

#else
  int iProcessor, nProcessor = size;

  unsigned short nDim = geometry->GetnDim(), iMarker;
  su2double *Solution, *Coord;
  unsigned long Buffer_Send_nVertex[1], iVertex, iPoint, nVertex_Surface = 0, nLocalVertex_Surface = 0,
  MaxLocalVertex_Surface = 0, nBuffer_Scalar;
  unsigned long *Buffer_Receive_nVertex = NULL;
  ofstream SurfAdj_file;

  /*--- Write the surface .csv file ---*/
  nLocalVertex_Surface = 0;
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
    if (config->GetMarker_All_Plotting(iMarker) == YES)
      for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        if (geometry->node[iPoint]->GetDomain()) nLocalVertex_Surface ++;
      }

  if (rank == MASTER_NODE)
    Buffer_Receive_nVertex = new unsigned long [nProcessor];

  Buffer_Send_nVertex[0] = nLocalVertex_Surface;

  SU2_MPI::Allreduce(&nLocalVertex_Surface, &MaxLocalVertex_Surface, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
  SU2_MPI::Gather(&Buffer_Send_nVertex, 1, MPI_UNSIGNED_LONG, Buffer_Receive_nVertex, 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);

  su2double *Buffer_Send_Coord_x = new su2double[MaxLocalVertex_Surface];
  su2double *Buffer_Send_Coord_y= new su2double[MaxLocalVertex_Surface];
  su2double *Buffer_Send_Coord_z= new su2double[MaxLocalVertex_Surface];
  unsigned long *Buffer_Send_GlobalPoint= new unsigned long[MaxLocalVertex_Surface];
  su2double *Buffer_Send_Sensitivity= new su2double[MaxLocalVertex_Surface];
  su2double *Buffer_Send_PsiRho= new su2double[MaxLocalVertex_Surface];
  su2double *Buffer_Send_Phi_x= new su2double[MaxLocalVertex_Surface];
  su2double *Buffer_Send_Phi_y= new su2double[MaxLocalVertex_Surface];
  su2double *Buffer_Send_Phi_z= new su2double[MaxLocalVertex_Surface];
  su2double *Buffer_Send_PsiE = NULL;

  if (config ->GetKind_Regime() == COMPRESSIBLE)
    Buffer_Send_PsiE =  new su2double[MaxLocalVertex_Surface];

  su2double *Buffer_Send_Sens_x = NULL, *Buffer_Send_Sens_y = NULL, *Buffer_Send_Sens_z = NULL;

  if (config->GetDiscrete_Adjoint()) {
    Buffer_Send_Sens_x = new su2double[MaxLocalVertex_Surface];
    Buffer_Send_Sens_y = new su2double[MaxLocalVertex_Surface];
    if (nDim == 3) {
      Buffer_Send_Sens_z = new su2double[MaxLocalVertex_Surface];
    }
  }

  nVertex_Surface = 0;
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
    if (config->GetMarker_All_Plotting(iMarker) == YES)
      for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        if (geometry->node[iPoint]->GetDomain()) {
          Solution = AdjSolver->GetNodes()->GetSolution(iPoint);
          //Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
          Coord = geometry->node[iPoint]->GetCoord();
          //d = AdjSolver->GetNodes()->GetForceProj_Vector(iPoint);
          Buffer_Send_GlobalPoint[nVertex_Surface] = geometry->node[iPoint]->GetGlobalIndex();
          Buffer_Send_Coord_x[nVertex_Surface] = Coord[0];
          Buffer_Send_Coord_y[nVertex_Surface] = Coord[1];
          Buffer_Send_Sensitivity[nVertex_Surface] =  AdjSolver->GetCSensitivity(iMarker, iVertex);
          Buffer_Send_PsiRho[nVertex_Surface] = Solution[0];
          Buffer_Send_Phi_x[nVertex_Surface] = Solution[1];
          Buffer_Send_Phi_y[nVertex_Surface] = Solution[2];
          if ((nDim == 2) && (config->GetKind_Regime() == COMPRESSIBLE)) Buffer_Send_PsiE[nVertex_Surface] = Solution[3];
          if (nDim == 3) {
            Buffer_Send_Coord_z[nVertex_Surface] = Coord[2];
            Buffer_Send_Phi_z[nVertex_Surface] = Solution[3];
            if(config->GetKind_Regime() == COMPRESSIBLE) Buffer_Send_PsiE[nVertex_Surface] = Solution[4];
          }
          if (config->GetDiscrete_Adjoint()) {
            Buffer_Send_Sens_x[nVertex_Surface] = AdjSolver->GetNodes()->GetSensitivity(iPoint, 0);
            Buffer_Send_Sens_y[nVertex_Surface] = AdjSolver->GetNodes()->GetSensitivity(iPoint, 1);
            if (nDim == 3) {
              Buffer_Send_Sens_z[nVertex_Surface] = AdjSolver->GetNodes()->GetSensitivity(iPoint, 2);
            }
          }

          /*--- If US system, the output should be in inches ---*/

          if (config->GetSystemMeasurements() == US) {
            Buffer_Send_Coord_x[nVertex_Surface] *= 12.0;
            Buffer_Send_Coord_y[nVertex_Surface] *= 12.0;
            if (nDim == 3) Buffer_Send_Coord_z[nVertex_Surface] *= 12.0;
          }

          nVertex_Surface++;
        }
      }

  su2double *Buffer_Receive_Coord_x = NULL, *Buffer_Receive_Coord_y = NULL, *Buffer_Receive_Coord_z = NULL, *Buffer_Receive_Sensitivity = NULL,
  *Buffer_Receive_PsiRho = NULL, *Buffer_Receive_Phi_x = NULL, *Buffer_Receive_Phi_y = NULL, *Buffer_Receive_Phi_z = NULL,
  *Buffer_Receive_PsiE = NULL, *Buffer_Receive_Sens_x = NULL, *Buffer_Receive_Sens_y = NULL, *Buffer_Receive_Sens_z = NULL;
  unsigned long *Buffer_Receive_GlobalPoint = NULL;

  if (rank == MASTER_NODE) {
    Buffer_Receive_Coord_x = new su2double [nProcessor*MaxLocalVertex_Surface];
    Buffer_Receive_Coord_y = new su2double [nProcessor*MaxLocalVertex_Surface];
    if (nDim == 3) Buffer_Receive_Coord_z = new su2double [nProcessor*MaxLocalVertex_Surface];
    Buffer_Receive_GlobalPoint = new unsigned long [nProcessor*MaxLocalVertex_Surface];
    Buffer_Receive_Sensitivity = new su2double [nProcessor*MaxLocalVertex_Surface];
    Buffer_Receive_PsiRho = new su2double [nProcessor*MaxLocalVertex_Surface];
    Buffer_Receive_Phi_x = new su2double [nProcessor*MaxLocalVertex_Surface];
    Buffer_Receive_Phi_y = new su2double [nProcessor*MaxLocalVertex_Surface];
    if (nDim == 3) Buffer_Receive_Phi_z = new su2double [nProcessor*MaxLocalVertex_Surface];
    if  (config->GetKind_Regime() == COMPRESSIBLE)
      Buffer_Receive_PsiE = new su2double [nProcessor*MaxLocalVertex_Surface];
    if (config->GetDiscrete_Adjoint()) {
      Buffer_Receive_Sens_x = new su2double[nProcessor*MaxLocalVertex_Surface];
      Buffer_Receive_Sens_y = new su2double[nProcessor*MaxLocalVertex_Surface];
      if (nDim == 3) {
        Buffer_Receive_Sens_z = new su2double[nProcessor*MaxLocalVertex_Surface];
      }
    }
  }

  nBuffer_Scalar = MaxLocalVertex_Surface;

  /*--- Send the information to the Master node ---*/
  SU2_MPI::Gather(Buffer_Send_Coord_x, nBuffer_Scalar, MPI_DOUBLE, Buffer_Receive_Coord_x, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  SU2_MPI::Gather(Buffer_Send_Coord_y, nBuffer_Scalar, MPI_DOUBLE, Buffer_Receive_Coord_y, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  if (nDim == 3) SU2_MPI::Gather(Buffer_Send_Coord_z, nBuffer_Scalar, MPI_DOUBLE, Buffer_Receive_Coord_z, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  SU2_MPI::Gather(Buffer_Send_GlobalPoint, nBuffer_Scalar, MPI_UNSIGNED_LONG, Buffer_Receive_GlobalPoint, nBuffer_Scalar, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
  SU2_MPI::Gather(Buffer_Send_Sensitivity, nBuffer_Scalar, MPI_DOUBLE, Buffer_Receive_Sensitivity, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  SU2_MPI::Gather(Buffer_Send_PsiRho, nBuffer_Scalar, MPI_DOUBLE, Buffer_Receive_PsiRho, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  SU2_MPI::Gather(Buffer_Send_Phi_x, nBuffer_Scalar, MPI_DOUBLE, Buffer_Receive_Phi_x, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  SU2_MPI::Gather(Buffer_Send_Phi_y, nBuffer_Scalar, MPI_DOUBLE, Buffer_Receive_Phi_y, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  if (nDim == 3) SU2_MPI::Gather(Buffer_Send_Phi_z, nBuffer_Scalar, MPI_DOUBLE, Buffer_Receive_Phi_z, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  if (config->GetKind_Regime() == COMPRESSIBLE)
      SU2_MPI::Gather(Buffer_Send_PsiE, nBuffer_Scalar, MPI_DOUBLE, Buffer_Receive_PsiE, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  if (config->GetDiscrete_Adjoint()) {
    SU2_MPI::Gather(Buffer_Send_Sens_x, nBuffer_Scalar, MPI_DOUBLE, Buffer_Receive_Sens_x, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    SU2_MPI::Gather(Buffer_Send_Sens_y, nBuffer_Scalar, MPI_DOUBLE, Buffer_Receive_Sens_y, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    if (nDim == 3) {
      SU2_MPI::Gather(Buffer_Send_Sens_z, nBuffer_Scalar, MPI_DOUBLE, Buffer_Receive_Sens_z, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    }
  }

  /*--- The master node is the one who writes the surface files ---*/
  if (rank == MASTER_NODE) {
    unsigned long iVertex, GlobalPoint, position;
    char cstr[200], buffer[50];
    ofstream SurfAdj_file;
    string filename = config->GetSurfAdjCoeff_FileName();

    /*--- Write file name with extension if unsteady ---*/
    strcpy (cstr, filename.c_str());

    if (config->GetTime_Marching() == HARMONIC_BALANCE) {
      SPRINTF (buffer, "_%d.csv", SU2_TYPE::Int(val_iInst));

    } else if (config->GetTime_Marching() && config->GetTime_Domain()) {
      if ((SU2_TYPE::Int(iExtIter) >= 0) && (SU2_TYPE::Int(iExtIter) < 10)) SPRINTF (buffer, "_0000%d.csv", SU2_TYPE::Int(iExtIter));
      if ((SU2_TYPE::Int(iExtIter) >= 10) && (SU2_TYPE::Int(iExtIter) < 100)) SPRINTF (buffer, "_000%d.csv", SU2_TYPE::Int(iExtIter));
      if ((SU2_TYPE::Int(iExtIter) >= 100) && (SU2_TYPE::Int(iExtIter) < 1000)) SPRINTF (buffer, "_00%d.csv", SU2_TYPE::Int(iExtIter));
      if ((SU2_TYPE::Int(iExtIter) >= 1000) && (SU2_TYPE::Int(iExtIter) < 10000)) SPRINTF (buffer, "_0%d.csv", SU2_TYPE::Int(iExtIter));
      if (SU2_TYPE::Int(iExtIter) >= 10000) SPRINTF (buffer, "_%d.csv", SU2_TYPE::Int(iExtIter));
    }
    else
      SPRINTF (buffer, ".csv");

    strcat (cstr, buffer);
    SurfAdj_file.open(cstr, ios::out);
    SurfAdj_file.precision(15);

    SurfAdj_file << "SENS_AOA=" << AdjSolver->GetTotal_Sens_AoA() * PI_NUMBER / 180.0 << endl;

    /*--- Write the 2D surface flow coefficient file ---*/
    if (geometry->GetnDim() == 2) {
      if (config->GetKind_Regime() == COMPRESSIBLE)
        SurfAdj_file <<  "\"Point\",\"Sensitivity\",\"PsiRho\",\"Phi_x\",\"Phi_y\",\"PsiE\",\"x_coord\",\"y_coord\"";
      else if (config->GetKind_Regime() == INCOMPRESSIBLE)
        SurfAdj_file <<  "\"Point\",\"Sensitivity\",\"PsiRho\",\"Phi_x\",\"Phi_y\",\"x_coord\",\"y_coord\"";

      if (config->GetDiscrete_Adjoint()) {
        SurfAdj_file << ",\" x_Sens\",\"y_Sens\"";
      }
      SurfAdj_file << "\n";

      for (iProcessor = 0; iProcessor < nProcessor; iProcessor++)
        for (iVertex = 0; iVertex < Buffer_Receive_nVertex[iProcessor]; iVertex++) {

          position = iProcessor*MaxLocalVertex_Surface+iVertex;
          GlobalPoint = Buffer_Receive_GlobalPoint[position];

          if (config->GetKind_Regime() == COMPRESSIBLE)
            SurfAdj_file << scientific << GlobalPoint <<
                            ", " << Buffer_Receive_Sensitivity[position] << ", " << Buffer_Receive_PsiRho[position] <<
                            ", " << Buffer_Receive_Phi_x[position] << ", " << Buffer_Receive_Phi_y[position] <<
                            ", " << Buffer_Receive_PsiE[position] << ", " << Buffer_Receive_Coord_x[position] <<
                            ", "<< Buffer_Receive_Coord_y[position];
          else if (config->GetKind_Regime() == INCOMPRESSIBLE)
            SurfAdj_file << scientific << GlobalPoint <<
                            ", " << Buffer_Receive_Sensitivity[position] << ", " << Buffer_Receive_PsiRho[position] <<
                            ", " << Buffer_Receive_Phi_x[position] << ", " << Buffer_Receive_Phi_y[position] <<
                            ", " << Buffer_Receive_Coord_x[position] <<
                            ", "<< Buffer_Receive_Coord_y[position];
          if (config->GetDiscrete_Adjoint()) {
            SurfAdj_file << ", " << Buffer_Receive_Sens_x[position] << ", " << Buffer_Receive_Sens_y[position];
          }
          SurfAdj_file << "\n";
        }
    }

    /*--- Write the 3D surface flow coefficient file ---*/
    if (geometry->GetnDim() == 3) {
      if (config->GetKind_Regime() == COMPRESSIBLE)
        SurfAdj_file <<  "\"Point\",\"Sensitivity\",\"PsiRho\",\"Phi_x\",\"Phi_y\",\"Phi_z\",\"PsiE\",\"x_coord\",\"y_coord\",\"z_coord\"";
      else if (config->GetKind_Regime() == INCOMPRESSIBLE)
        SurfAdj_file <<  "\"Point\",\"Sensitivity\",\"PsiRho\",\"Phi_x\",\"Phi_y\",\"Phi_z\",\"x_coord\",\"y_coord\",\"z_coord\"";

      if (config->GetDiscrete_Adjoint()) {
        SurfAdj_file << ",\"x_Sens\",\"y_Sens\",\"z_Sens\"";
      }
      SurfAdj_file << "\n";

      for (iProcessor = 0; iProcessor < nProcessor; iProcessor++)
        for (iVertex = 0; iVertex < Buffer_Receive_nVertex[iProcessor]; iVertex++) {
          position = iProcessor*MaxLocalVertex_Surface+iVertex;
          GlobalPoint = Buffer_Receive_GlobalPoint[position];

          if (config->GetKind_Regime() == COMPRESSIBLE)
            SurfAdj_file << scientific << GlobalPoint <<
                            ", " << Buffer_Receive_Sensitivity[position] << ", " << Buffer_Receive_PsiRho[position] <<
                            ", " << Buffer_Receive_Phi_x[position] << ", " << Buffer_Receive_Phi_y[position] << ", " << Buffer_Receive_Phi_z[position] <<
                            ", " << Buffer_Receive_PsiE[position] <<", "<< Buffer_Receive_Coord_x[position] <<
                            ", "<< Buffer_Receive_Coord_y[position] <<", "<< Buffer_Receive_Coord_z[position];
          else if (config->GetKind_Regime() == INCOMPRESSIBLE)
            SurfAdj_file << scientific << GlobalPoint <<
                            ", " << Buffer_Receive_Sensitivity[position] << ", " << Buffer_Receive_PsiRho[position] <<
                            ", " << Buffer_Receive_Phi_x[position] << ", " << Buffer_Receive_Phi_y[position] << ", " << Buffer_Receive_Phi_z[position] <<
                            ", "<< Buffer_Receive_Coord_x[position] <<
                            ", "<< Buffer_Receive_Coord_y[position] <<", "<< Buffer_Receive_Coord_z[position];

          if (config->GetDiscrete_Adjoint()) {
            SurfAdj_file << ", " << Buffer_Receive_Sens_x[position] << ", " << Buffer_Receive_Sens_y[position] << ", " << Buffer_Receive_Sens_z[position];
          }
          SurfAdj_file << "\n";
        }
    }

  }

  if (rank == MASTER_NODE) {
    delete [] Buffer_Receive_nVertex;
    delete [] Buffer_Receive_Coord_x;
    delete [] Buffer_Receive_Coord_y;
    if (nDim == 3) delete [] Buffer_Receive_Coord_z;
    delete [] Buffer_Receive_Sensitivity;
    delete [] Buffer_Receive_PsiRho;
    delete [] Buffer_Receive_Phi_x;
    delete [] Buffer_Receive_Phi_y;
    if (nDim == 3) delete [] Buffer_Receive_Phi_z;
    if (config->GetKind_Regime() == COMPRESSIBLE)
      delete [] Buffer_Receive_PsiE;
    delete [] Buffer_Receive_GlobalPoint;
    if (config->GetDiscrete_Adjoint()) {
      delete [] Buffer_Receive_Sens_x;
      delete [] Buffer_Receive_Sens_y;
      if (nDim == 3) {
        delete [] Buffer_Receive_Sens_z;
      }
    }
  }

  delete [] Buffer_Send_Coord_x;
  delete [] Buffer_Send_Coord_y;
  delete [] Buffer_Send_Coord_z;
  delete [] Buffer_Send_GlobalPoint;
  delete [] Buffer_Send_Sensitivity;
  delete [] Buffer_Send_PsiRho;
  delete [] Buffer_Send_Phi_x;
  delete [] Buffer_Send_Phi_y;
  delete [] Buffer_Send_Phi_z;
  delete [] Buffer_Send_PsiE;
  if (Buffer_Send_Sens_x != NULL) delete [] Buffer_Send_Sens_x;
  if (Buffer_Send_Sens_y != NULL) delete [] Buffer_Send_Sens_y;
  if (Buffer_Send_Sens_z != NULL) delete [] Buffer_Send_Sens_z;

  SurfAdj_file.close();

#endif
}

void COutputLegacy::MergeConnectivity(CConfig *config, CGeometry *geometry, unsigned short val_iZone) {

  /*--- Flags identifying the types of files to be written. ---*/

  bool Wrt_Vol = config->GetWrt_Vol_Sol();
  bool Wrt_Srf = config->GetWrt_Srf_Sol();

  /*--- Merge connectivity for each type of element (excluding halos). Note
   that we only need to merge the connectivity once, as it does not change
   during computation. Check whether the base file has been written. ---*/

  /*--- Merge volumetric grid. ---*/

  if (Wrt_Vol) {

    if ((rank == MASTER_NODE) && (size != SINGLE_NODE) && (nGlobal_Tria != 0))
      cout <<"Merging volumetric triangle grid connectivity." << endl;
    MergeVolumetricConnectivity(config, geometry, TRIANGLE    );

    if ((rank == MASTER_NODE) && (size != SINGLE_NODE) && (nGlobal_Quad != 0))
      cout <<"Merging volumetric quadrilateral grid connectivity." << endl;
    MergeVolumetricConnectivity(config, geometry, QUADRILATERAL   );

    if ((rank == MASTER_NODE) && (size != SINGLE_NODE) && (nGlobal_Tetr != 0))
      cout <<"Merging volumetric tetrahedron grid connectivity." << endl;
    MergeVolumetricConnectivity(config, geometry, TETRAHEDRON );

    if ((rank == MASTER_NODE) && (size != SINGLE_NODE) && (nGlobal_Hexa != 0))
      cout <<"Merging volumetric hexahedron grid connectivity." << endl;
    MergeVolumetricConnectivity(config, geometry, HEXAHEDRON  );

    if ((rank == MASTER_NODE) && (size != SINGLE_NODE) && (nGlobal_Pris != 0))
      cout <<"Merging volumetric prism grid connectivity." << endl;
    MergeVolumetricConnectivity(config, geometry, PRISM       );

    if ((rank == MASTER_NODE) && (size != SINGLE_NODE) && (nGlobal_Pyra != 0))
      cout <<"Merging volumetric pyramid grid connectivity." << endl;
    MergeVolumetricConnectivity(config, geometry, PYRAMID     );

  }

  /*--- Merge surface grid. ---*/

  if (Wrt_Srf) {

    if ((rank == MASTER_NODE) && (size != SINGLE_NODE) && (nGlobal_Line != 0))
      cout <<"Merging surface line grid connectivity." << endl;
    MergeSurfaceConnectivity(config, geometry, LINE);

    if ((rank == MASTER_NODE) && (size != SINGLE_NODE) && (nGlobal_BoundTria != 0))
      cout <<"Merging surface triangle grid connectivity." << endl;
    MergeSurfaceConnectivity(config, geometry, TRIANGLE);

    if ((rank == MASTER_NODE) && (size != SINGLE_NODE) && (nGlobal_BoundQuad != 0))
      cout <<"Merging surface quadrilateral grid connectivity." << endl;
    MergeSurfaceConnectivity(config, geometry, QUADRILATERAL);

  }

  /*--- Update total number of volume elements after merge. ---*/

  nGlobal_Elem = nGlobal_Tria + nGlobal_Quad + nGlobal_Tetr +
  nGlobal_Hexa + nGlobal_Pyra + nGlobal_Pris;

  /*--- Update total number of surface elements after merge. ---*/

  nSurf_Elem = nGlobal_Line + nGlobal_BoundTria + nGlobal_BoundQuad;

}

void COutputLegacy::MergeCoordinates(CConfig *config, CGeometry *geometry) {

  /*--- Local variables needed on all processors ---*/

  unsigned short iDim, nDim = geometry->GetnDim();
  unsigned long iPoint;

  unsigned short kind_SU2 = config->GetKind_SU2();

#ifndef HAVE_MPI

  /*--- In serial, the single process has access to all geometry, so simply
   load the coordinates into the data structure. ---*/

  unsigned short iMarker;
  unsigned long iVertex, nTotalPoints = 0;
  int SendRecv;

  bool isPeriodic;

  /*--- First, create a structure to locate any periodic halo nodes ---*/
  int *Local_Halo = new int[geometry->GetnPoint()];
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
    Local_Halo[iPoint] = !geometry->node[iPoint]->GetDomain();

  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) {
      SendRecv = config->GetMarker_All_SendRecv(iMarker);
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

        /*--- For SU2_CFD and SU2_SOL we want to remove the periodic halo nodes,
         * but for SU2_DEF we want them to be included, therefore the definition of a periodic point
         * is different in each case ---*/

        if (kind_SU2 == SU2_DEF) {
          isPeriodic = ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() > 0));
        }else {
          isPeriodic = ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() > 0) &&
                        (geometry->vertex[iMarker][iVertex]->GetRotation_Type() % 2 == 1));
        }

        if (isPeriodic && (SendRecv < 0)) {
          Local_Halo[iPoint] = false;
        }
      }

    }
  }

  /*--- Total number of points in the mesh (this might include periodic points). ---*/
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
    if (!Local_Halo[iPoint]) nTotalPoints++;

  nGlobal_Poin = nTotalPoints;
  nGlobal_Doma = geometry->GetnPointDomain();

  /*--- Allocate the coordinates data structure. ---*/

  Coords = new su2double*[nDim];
  for (iDim = 0; iDim < nDim; iDim++) {
    Coords[iDim] = new su2double[nGlobal_Poin];
  }

  /*--- Loop over the mesh to collect the coords of the local points ---*/

  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {

    /*--- Check if the node belongs to the domain (i.e, not a halo node).
     Sort by the global index, even in serial there is a renumbering (e.g. RCM). ---*/

    if (!Local_Halo[iPoint]) {

      /*--- Retrieve the current coordinates at this node. ---*/

      unsigned long iGlobal_Index = geometry->node[iPoint]->GetGlobalIndex();

      for (iDim = 0; iDim < nDim; iDim++) {
        Coords[iDim][iGlobal_Index] = geometry->node[iPoint]->GetCoord(iDim);

        /*--- If US system, the output should be in inches ---*/

        if ((config->GetSystemMeasurements() == US) && (config->GetKind_SU2() != SU2_DEF)) {
          Coords[iDim][iGlobal_Index] *= 12.0;
        }

      }

    }
  }


  delete [] Local_Halo;

#else

  /*--- MPI preprocessing ---*/
  int iProcessor, nProcessor = size;
  unsigned long jPoint;

  bool Wrt_Halo = config->GetWrt_Halo(), isPeriodic;

  /*--- Local variables needed for merging the geometry with MPI. ---*/

  unsigned long iVertex, iMarker;
  unsigned long Buffer_Send_nPoin[1], *Buffer_Recv_nPoin = NULL;
  unsigned long nLocalPoint = 0, MaxLocalPoint = 0;
  unsigned long iGlobal_Index = 0, nBuffer_Scalar = 0;

  if (rank == MASTER_NODE) Buffer_Recv_nPoin = new unsigned long[nProcessor];

  int *Local_Halo = new int[geometry->GetnPoint()];
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
    Local_Halo[iPoint] = !geometry->node[iPoint]->GetDomain();

  /*--- Search all send/recv boundaries on this partition for any periodic
   nodes that were part of the original domain. We want to recover these
   for visualization purposes. ---*/

  if (Wrt_Halo) {
    nLocalPoint = geometry->GetnPoint();
  } else {
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) {

        /*--- Checking for less than or equal to the rank, because there may
         be some periodic halo nodes that send info to the same rank. ---*/

        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

          /*--- For SU2_CFD and SU2_SOL we want to remove the periodic halo nodes,
           * but for SU2_DEF we want them to be included, therefore the definition of a periodic point
           * is different in each case ---*/

          if (kind_SU2 == SU2_DEF) {
            isPeriodic = ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() > 0));
          }else {
          isPeriodic = ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() > 0) &&
                        (geometry->vertex[iMarker][iVertex]->GetRotation_Type() % 2 == 1));
          }
          if (isPeriodic) {
            Local_Halo[iPoint] = false;
          }
        }
      }
    }

    /*--- Sum total number of nodes that belong to the domain ---*/

    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
      if (Local_Halo[iPoint] == false)
        nLocalPoint++;
  }
  Buffer_Send_nPoin[0] = nLocalPoint;

  /*--- Communicate the total number of nodes on this domain. ---*/

  SU2_MPI::Gather(&Buffer_Send_nPoin, 1, MPI_UNSIGNED_LONG,
                  Buffer_Recv_nPoin, 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&nLocalPoint, &MaxLocalPoint, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);

  if (rank == MASTER_NODE) {
    nGlobal_Doma = 0;
    for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
      nGlobal_Doma += Buffer_Recv_nPoin[iProcessor];
    }
  }
  nBuffer_Scalar = MaxLocalPoint;

  /*--- Send and Recv buffers. ---*/

  su2double *Buffer_Send_X = new su2double[MaxLocalPoint];
  su2double *Buffer_Recv_X = NULL;

  su2double *Buffer_Send_Y = new su2double[MaxLocalPoint];
  su2double *Buffer_Recv_Y = NULL;

  su2double *Buffer_Send_Z = NULL, *Buffer_Recv_Z = NULL;
  if (nDim == 3) Buffer_Send_Z = new su2double[MaxLocalPoint];

  unsigned long *Buffer_Send_GlobalIndex = new unsigned long[MaxLocalPoint];
  unsigned long *Buffer_Recv_GlobalIndex = NULL;

  /*--- Prepare the receive buffers in the master node only. ---*/

  if (rank == MASTER_NODE) {

    Buffer_Recv_X = new su2double[nProcessor*MaxLocalPoint];
    Buffer_Recv_Y = new su2double[nProcessor*MaxLocalPoint];
    if (nDim == 3) Buffer_Recv_Z = new su2double[nProcessor*MaxLocalPoint];
    Buffer_Recv_GlobalIndex = new unsigned long[nProcessor*MaxLocalPoint];

    /*--- Sum total number of nodes to be written and allocate arrays ---*/
    nGlobal_Poin = 0;
    for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
      nGlobal_Poin += Buffer_Recv_nPoin[iProcessor];
    }
    Coords = new su2double*[nDim];
    for (iDim = 0; iDim < nDim; iDim++) {
      Coords[iDim] = new su2double[nGlobal_Poin];
    }
  }

  /*--- Main communication routine. Loop over each coordinate and perform
   the MPI comm. Temporary 1-D buffers are used to send the coordinates at
   all nodes on each partition to the master node. These are then unpacked
   by the master and sorted by global index in one large n-dim. array. ---*/

  /*--- Loop over this partition to collect the coords of the local points. ---*/
  su2double *Coords_Local; jPoint = 0;
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {

    /*--- Check for halos and write only if requested ---*/
    if (!Local_Halo[iPoint] || Wrt_Halo) {

      /*--- Retrieve local coordinates at this node. ---*/
      Coords_Local = geometry->node[iPoint]->GetCoord();

      /*--- Load local coords into the temporary send buffer. ---*/
      Buffer_Send_X[jPoint] = Coords_Local[0];
      Buffer_Send_Y[jPoint] = Coords_Local[1];
      if (nDim == 3) Buffer_Send_Z[jPoint] = Coords_Local[2];

      /*--- If US system, the output should be in inches ---*/

      if ((config->GetSystemMeasurements() == US) && (config->GetKind_SU2() != SU2_DEF)) {
        Buffer_Send_X[jPoint] *= 12.0;
        Buffer_Send_Y[jPoint] *= 12.0;
        if (nDim == 3) Buffer_Send_Z[jPoint] *= 12.0;
      }

      /*--- Store the global index for this local node. ---*/
      Buffer_Send_GlobalIndex[jPoint] = geometry->node[iPoint]->GetGlobalIndex();

      /*--- Increment jPoint as the counter. We need this because iPoint
       may include halo nodes that we skip over during this loop. ---*/
      jPoint++;
    }
  }

  /*--- Gather the coordinate data on the master node using MPI. ---*/

  SU2_MPI::Gather(Buffer_Send_X, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_X, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  SU2_MPI::Gather(Buffer_Send_Y, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Y, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  if (nDim == 3) {
    SU2_MPI::Gather(Buffer_Send_Z, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Z, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  }
  SU2_MPI::Gather(Buffer_Send_GlobalIndex, nBuffer_Scalar, MPI_UNSIGNED_LONG, Buffer_Recv_GlobalIndex, nBuffer_Scalar, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);

  /*--- The master node unpacks and sorts this variable by global index ---*/

  if (rank == MASTER_NODE) {
    jPoint = 0;
    for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
      for (iPoint = 0; iPoint < Buffer_Recv_nPoin[iProcessor]; iPoint++) {
        /*--- Get global index, then loop over each variable and store ---*/
        iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
        if (iGlobal_Index >= nGlobal_Poin) {
          cout << iGlobal_Index << " " << nGlobal_Poin << endl;
        }
        Coords[0][iGlobal_Index] = Buffer_Recv_X[jPoint];
        Coords[1][iGlobal_Index] = Buffer_Recv_Y[jPoint];
        if (nDim == 3) Coords[2][iGlobal_Index] = Buffer_Recv_Z[jPoint];
        jPoint++;
      }
      /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
      jPoint = (iProcessor+1)*nBuffer_Scalar;
    }
  }

  /*--- Immediately release the temporary data buffers. ---*/

  delete [] Local_Halo;
  delete [] Buffer_Send_X;
  delete [] Buffer_Send_Y;
  if (Buffer_Send_Z != NULL) delete [] Buffer_Send_Z;
  delete [] Buffer_Send_GlobalIndex;
  if (rank == MASTER_NODE) {
    delete [] Buffer_Recv_X;
    delete [] Buffer_Recv_Y;
    if (Buffer_Recv_Z != NULL)  delete [] Buffer_Recv_Z;
    delete [] Buffer_Recv_GlobalIndex;
    delete [] Buffer_Recv_nPoin;
  }

#endif

}

void COutputLegacy::MergeVolumetricConnectivity(CConfig *config, CGeometry *geometry, unsigned short Elem_Type) {

  int iProcessor;
  unsigned short NODES_PER_ELEMENT = 0;
  unsigned long iPoint, iNode, jNode;
  unsigned long iElem = 0;
  unsigned long nLocalElem = 0, nElem_Total = 0;

  unsigned long iVertex, iMarker;
  unsigned long jElem;
  int SendRecv, RecvFrom;

  unsigned long Buffer_Send_nElem[1], *Buffer_Recv_nElem = NULL;
  unsigned long nBuffer_Scalar = 0;
  unsigned long kNode = 0, kElem = 0;
  unsigned long MaxLocalElem = 0, iGlobal_Index, jPoint, kPoint;

  bool Wrt_Halo = config->GetWrt_Halo();
  bool *Write_Elem = NULL, notPeriodic, notHalo, addedPeriodic, isPeriodic;

  unsigned short kind_SU2 = config->GetKind_SU2();

  int *Conn_Elem = NULL;

  /*--- Store the local number of this element type and the number of nodes
   per this element type. In serial, this will be the total number of this
   element type in the entire mesh. In parallel, it is the number on only
   the current partition. ---*/

  switch (Elem_Type) {
    case TRIANGLE:
      nLocalElem = geometry->GetnElemTria();
      NODES_PER_ELEMENT = N_POINTS_TRIANGLE;
      break;
    case QUADRILATERAL:
      nLocalElem = geometry->GetnElemQuad();
      NODES_PER_ELEMENT = N_POINTS_QUADRILATERAL;
      break;
    case TETRAHEDRON:
      nLocalElem = geometry->GetnElemTetr();
      NODES_PER_ELEMENT = N_POINTS_TETRAHEDRON;
      break;
    case HEXAHEDRON:
      nLocalElem = geometry->GetnElemHexa();
      NODES_PER_ELEMENT = N_POINTS_HEXAHEDRON;
      break;
    case PRISM:
      nLocalElem = geometry->GetnElemPris();
      NODES_PER_ELEMENT = N_POINTS_PRISM;
      break;
    case PYRAMID:
      nLocalElem = geometry->GetnElemPyra();
      NODES_PER_ELEMENT = N_POINTS_PYRAMID;
      break;
    default:
      SU2_MPI::Error("Unrecognized element type", CURRENT_FUNCTION);
  }

  /*--- Find the max number of this element type among all
   partitions and set up buffers. ---*/

  Buffer_Send_nElem[0] = nLocalElem;
  if (rank == MASTER_NODE) Buffer_Recv_nElem = new unsigned long[size];

#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&nLocalElem, &MaxLocalElem, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
  SU2_MPI::Gather(&Buffer_Send_nElem, 1, MPI_UNSIGNED_LONG, Buffer_Recv_nElem, 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
#else
  MaxLocalElem = nLocalElem;
  Buffer_Recv_nElem[0] = Buffer_Send_nElem[0];
#endif

  nBuffer_Scalar = MaxLocalElem*NODES_PER_ELEMENT;

  /*--- Send and Recv buffers ---*/

  unsigned long *Buffer_Send_Elem = new unsigned long[nBuffer_Scalar];
  unsigned long *Buffer_Recv_Elem = NULL;

  unsigned short *Buffer_Send_Halo = new unsigned short[MaxLocalElem];
  unsigned short *Buffer_Recv_Halo = NULL;

  /*--- Prepare the receive buffers on the master node only. ---*/

  if (rank == MASTER_NODE) {
    Buffer_Recv_Elem = new unsigned long[size*nBuffer_Scalar];
    Buffer_Recv_Halo = new unsigned short[size*MaxLocalElem];
    if (MaxLocalElem > 0) Conn_Elem = new int[size*MaxLocalElem*NODES_PER_ELEMENT];
  }

  /*--- Force the removal of all added periodic elements (use global index).
   First, we isolate and create a list of all added periodic points, excluding
   those that we part of the original domain (we want these to be in the
   output files). ---*/

  vector<unsigned long> Added_Periodic;
  Added_Periodic.clear();

  if (kind_SU2 != SU2_DEF) {
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) {
        SendRecv = config->GetMarker_All_SendRecv(iMarker);
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

          if ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() > 0) &&
              (geometry->vertex[iMarker][iVertex]->GetRotation_Type() % 2 == 0) &&
              (SendRecv < 0)) {
            Added_Periodic.push_back(geometry->node[iPoint]->GetGlobalIndex());
          }
        }
      }
    }
  }

  /*--- Now we communicate this information to all processors, so that they
   can force the removal of these particular nodes by flagging them as halo
   points. In general, this should be a small percentage of the total mesh,
   so the communication/storage costs here shouldn't be prohibitive. ---*/

  /*--- First communicate the number of points that each rank has found ---*/
  unsigned long nAddedPeriodic = 0, maxAddedPeriodic = 0;
  unsigned long Buffer_Send_nAddedPeriodic[1], *Buffer_Recv_nAddedPeriodic = NULL;
  Buffer_Recv_nAddedPeriodic = new unsigned long[size];

  nAddedPeriodic = Added_Periodic.size();
  Buffer_Send_nAddedPeriodic[0] = nAddedPeriodic;

#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&nAddedPeriodic, &maxAddedPeriodic, 1, MPI_UNSIGNED_LONG,
                     MPI_MAX, MPI_COMM_WORLD);
  SU2_MPI::Allgather(&Buffer_Send_nAddedPeriodic, 1, MPI_UNSIGNED_LONG,
                     Buffer_Recv_nAddedPeriodic,  1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
#else
  maxAddedPeriodic = nAddedPeriodic;
  Buffer_Recv_nAddedPeriodic[0] = Buffer_Send_nAddedPeriodic[0];
#endif

  /*--- Communicate the global index values of all added periodic nodes. ---*/
  unsigned long *Buffer_Send_AddedPeriodic = new unsigned long[maxAddedPeriodic];
  unsigned long *Buffer_Recv_AddedPeriodic = new unsigned long[size*maxAddedPeriodic];

  for (iPoint = 0; iPoint < Added_Periodic.size(); iPoint++) {
    Buffer_Send_AddedPeriodic[iPoint] = Added_Periodic[iPoint];
  }

  /*--- Gather the element connectivity information. All processors will now
   have a copy of the global index values for all added periodic points. ---*/

#ifdef HAVE_MPI
  SU2_MPI::Allgather(Buffer_Send_AddedPeriodic, maxAddedPeriodic, MPI_UNSIGNED_LONG,
                     Buffer_Recv_AddedPeriodic, maxAddedPeriodic, MPI_UNSIGNED_LONG,
                     MPI_COMM_WORLD);
#else
  for (iPoint = 0; iPoint < maxAddedPeriodic; iPoint++) Buffer_Recv_AddedPeriodic[iPoint] = Buffer_Send_AddedPeriodic[iPoint];
#endif

  /*--- Search all send/recv boundaries on this partition for halo cells. In
   particular, consider only the recv conditions (these are the true halo
   nodes). Check the ranks of the processors that are communicating and
   choose to keep only the halo cells from the higher rank processor. Here,
   we are also choosing to keep periodic nodes that were part of the original
   domain. We will check the communicated list of added periodic points. ---*/

  int *Local_Halo = new int[geometry->GetnPoint()];
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
    Local_Halo[iPoint] = !geometry->node[iPoint]->GetDomain();

  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) {
      SendRecv = config->GetMarker_All_SendRecv(iMarker);
      RecvFrom = abs(SendRecv)-1;

      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        iGlobal_Index = geometry->node[iPoint]->GetGlobalIndex();

        /*--- We need to keep one copy of overlapping halo cells. ---*/
        notHalo = ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() == 0) &&
                   (SendRecv < 0) && (rank > RecvFrom));

        /*--- We want to keep the periodic nodes that were part of the original domain.
         For SU2_DEF we want to keep all periodic nodes. ---*/

        if (kind_SU2 == SU2_DEF) {
          isPeriodic = ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() > 0));
        }else {
          isPeriodic = ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() > 0) &&
                        (geometry->vertex[iMarker][iVertex]->GetRotation_Type() % 2 == 1));
        }

        notPeriodic = (isPeriodic && (SendRecv < 0));

        /*--- Lastly, check that this isn't an added periodic point that
         we will forcibly remove. Use the communicated list of these points. ---*/
        addedPeriodic = false; kPoint = 0;
        for (iProcessor = 0; iProcessor < size; iProcessor++) {
          for (jPoint = 0; jPoint < Buffer_Recv_nAddedPeriodic[iProcessor]; jPoint++) {
            if (iGlobal_Index == Buffer_Recv_AddedPeriodic[kPoint+jPoint])
              addedPeriodic = true;
          }
          /*--- Adjust jNode to index of next proc's data in the buffers. ---*/
          kPoint = (iProcessor+1)*maxAddedPeriodic;
        }

        /*--- If we found either of these types of nodes, flag them to be kept. ---*/
        if ((notHalo || notPeriodic) && !addedPeriodic) {
          Local_Halo[iPoint] = false;
        }
      }
    }
  }

  /*--- Loop over all elements in this partition and load the
   elements of the current type into the buffer to be sent to
   the master node. ---*/

  jNode = 0; jElem = 0;
  for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
    if (geometry->elem[iElem]->GetVTK_Type() == Elem_Type) {

      /*--- Loop over all nodes in this element and load the
       connectivity into the send buffer. ---*/

      Buffer_Send_Halo[jElem] = false;
      for (iNode = 0; iNode < NODES_PER_ELEMENT; iNode++) {

        /*--- Store the global index values directly. ---*/

        iPoint = geometry->elem[iElem]->GetNode(iNode);
        Buffer_Send_Elem[jNode] = geometry->node[iPoint]->GetGlobalIndex();

        /*--- Check if this is a halo node. If so, flag this element
         as a halo cell. We will use this later to sort and remove
         any duplicates from the connectivity list. ---*/

        if (Local_Halo[iPoint]) {
          Buffer_Send_Halo[jElem] = true;
        }

        /*--- Increment jNode as the counter. We need this because iElem
         may include other elements that we skip over during this loop. ---*/

        jNode++;
      }
      jElem++;
    }
  }

  /*--- Gather the element connectivity information. ---*/

#ifdef HAVE_MPI
  SU2_MPI::Gather(Buffer_Send_Elem, nBuffer_Scalar, MPI_UNSIGNED_LONG, Buffer_Recv_Elem, nBuffer_Scalar, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
  SU2_MPI::Gather(Buffer_Send_Halo, MaxLocalElem, MPI_UNSIGNED_SHORT, Buffer_Recv_Halo, MaxLocalElem, MPI_UNSIGNED_SHORT, MASTER_NODE, MPI_COMM_WORLD);
#else
  for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Elem[iPoint] = Buffer_Send_Elem[iPoint];
  for (iPoint = 0; iPoint < MaxLocalElem; iPoint++) Buffer_Recv_Halo[iPoint] = Buffer_Send_Halo[iPoint];
#endif

  /*--- The master node unpacks and sorts the connectivity. ---*/

  if (rank == MASTER_NODE) {

    /*---  We need to remove any duplicate elements (halo cells) that
     exist on multiple partitions. Start by initializing all elements
     to the "write" state by using a boolean array. ---*/

    Write_Elem = new bool[size*MaxLocalElem];
    for (iElem = 0; iElem < size*MaxLocalElem; iElem++) {
      Write_Elem[iElem] = true;
    }

    /*--- Remove the rind layer from the solution only if requested ---*/

    if (!Wrt_Halo) {

      /*--- Loop for flagging duplicate elements so that they are not
       included in the final connectivity list. ---*/

      kElem = 0;
      for (iProcessor = 0; iProcessor < size; iProcessor++) {
        for (iElem = 0; iElem < Buffer_Recv_nElem[iProcessor]; iElem++) {

          /*--- Check if this element was marked as a halo. ---*/
          if (Buffer_Recv_Halo[kElem+iElem])
            Write_Elem[kElem+iElem] = false;

        }
        kElem = (iProcessor+1)*MaxLocalElem;
      }
    }

    /*--- Store the unique connectivity list for this element type. ---*/

    jNode = 0; kNode = 0; jElem = 0; nElem_Total = 0;
    for (iProcessor = 0; iProcessor < size; iProcessor++) {
      for (iElem = 0; iElem < Buffer_Recv_nElem[iProcessor]; iElem++) {

        /*--- Only write the elements that were flagged for it. ---*/
        if (Write_Elem[jElem+iElem]) {

          /*--- Increment total count for this element type ---*/
          nElem_Total++;

          /*--- Get global index, then loop over each variable and store.
           Note that we are adding one to the index value because CGNS/Tecplot
           use 1-based indexing.---*/

          for (iNode = 0; iNode < NODES_PER_ELEMENT; iNode++) {
            Conn_Elem[kNode] = (int)Buffer_Recv_Elem[jNode+iElem*NODES_PER_ELEMENT+iNode] + 1;
            kNode++;
          }
        }
      }
      /*--- Adjust jNode to index of next proc's data in the buffers. ---*/
      jElem = (iProcessor+1)*MaxLocalElem;
      jNode = (iProcessor+1)*nBuffer_Scalar;
    }
  }

  /*--- Immediately release the temporary buffers. ---*/
  delete [] Buffer_Send_Elem;
  delete [] Buffer_Send_Halo;
  delete [] Buffer_Recv_nAddedPeriodic;
  delete [] Buffer_Send_AddedPeriodic;
  delete [] Buffer_Recv_AddedPeriodic;
  delete [] Local_Halo;
  if (rank == MASTER_NODE) {
    delete [] Buffer_Recv_nElem;
    delete [] Buffer_Recv_Elem;
    delete [] Buffer_Recv_Halo;
    delete [] Write_Elem;
  }

  /*--- Store the particular global element count in the class data,
   and set the class data pointer to the connectivity array. ---*/

  if (rank == MASTER_NODE) {
    switch (Elem_Type) {
      case TRIANGLE:
        nGlobal_Tria = nElem_Total;
        if (nGlobal_Tria > 0) Conn_Tria = Conn_Elem;
        break;
      case QUADRILATERAL:
        nGlobal_Quad = nElem_Total;
        if (nGlobal_Quad > 0) Conn_Quad = Conn_Elem;
        break;
      case TETRAHEDRON:
        nGlobal_Tetr = nElem_Total;
        if (nGlobal_Tetr > 0) Conn_Tetr = Conn_Elem;
        break;
      case HEXAHEDRON:
        nGlobal_Hexa = nElem_Total;
        if (nGlobal_Hexa > 0) Conn_Hexa = Conn_Elem;
        break;
      case PRISM:
        nGlobal_Pris = nElem_Total;
        if (nGlobal_Pris > 0) Conn_Pris = Conn_Elem;
        break;
      case PYRAMID:
        nGlobal_Pyra = nElem_Total;
        if (nGlobal_Pyra > 0) Conn_Pyra = Conn_Elem;
        break;
      default:
        SU2_MPI::Error("Unrecognized element type", CURRENT_FUNCTION);
        break;
    }
  }

}

void COutputLegacy::MergeSurfaceConnectivity(CConfig *config, CGeometry *geometry, unsigned short Elem_Type) {

  unsigned short NODES_PER_ELEMENT;

  unsigned short iMarker;
  unsigned long iPoint, iNode, jNode;
  unsigned long iElem = 0;
  unsigned long nLocalElem = 0, nElem_Total = 0;

  int iProcessor;
  unsigned long jElem;

  unsigned long iVertex;

  int SendRecv, RecvFrom;

  unsigned long Buffer_Send_nElem[1], *Buffer_Recv_nElem = NULL;
  unsigned long nBuffer_Scalar = 0;
  unsigned long kNode = 0, kElem = 0;
  unsigned long MaxLocalElem = 0, iGlobal_Index, jPoint, kPoint;

  bool Wrt_Halo = config->GetWrt_Halo();
  bool *Write_Elem = NULL, notPeriodic, notHalo, addedPeriodic;


  int *Conn_Elem = NULL;

  /*--- Store the local number of this element type and the number of nodes
   per this element type. In serial, this will be the total number of this
   element type in the entire mesh. In parallel, it is the number on only
   the current partition. ---*/

  nLocalElem = 0;

  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_Plotting(iMarker) == YES) {
      for (iElem = 0; iElem < geometry->GetnElem_Bound(iMarker); iElem++) {
        if (geometry->bound[iMarker][iElem]->GetVTK_Type() == Elem_Type) {
          nLocalElem++;
        }
      }
    }
  }

  switch (Elem_Type) {
    case LINE:
      NODES_PER_ELEMENT = N_POINTS_LINE;
      break;
    case TRIANGLE:
      NODES_PER_ELEMENT = N_POINTS_TRIANGLE;
      break;
    case QUADRILATERAL:
      NODES_PER_ELEMENT = N_POINTS_QUADRILATERAL;
      break;
    default:
      SU2_MPI::Error("Unrecognized element type", CURRENT_FUNCTION);
      NODES_PER_ELEMENT = 0;
      break;
  }

  /*--- Find the max number of this element type among all
   partitions and set up buffers. ---*/

  Buffer_Send_nElem[0] = nLocalElem;
  if (rank == MASTER_NODE) Buffer_Recv_nElem = new unsigned long[size];

#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&nLocalElem, &MaxLocalElem, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
  SU2_MPI::Gather(&Buffer_Send_nElem, 1, MPI_UNSIGNED_LONG, Buffer_Recv_nElem, 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
#else
  MaxLocalElem = nLocalElem;
  Buffer_Recv_nElem[0] = Buffer_Send_nElem[0];
#endif

  nBuffer_Scalar = MaxLocalElem*NODES_PER_ELEMENT;

  /*--- Send and Recv buffers ---*/

  unsigned long *Buffer_Send_Elem = new unsigned long[nBuffer_Scalar];
  unsigned long *Buffer_Recv_Elem = NULL;

  unsigned short *Buffer_Send_Halo = new unsigned short[MaxLocalElem];
  unsigned short *Buffer_Recv_Halo = NULL;

  /*--- Prepare the receive buffers on the master node only. ---*/

  if (rank == MASTER_NODE) {
    Buffer_Recv_Elem = new unsigned long[size*nBuffer_Scalar];
    Buffer_Recv_Halo = new unsigned short[size*MaxLocalElem];
    if (MaxLocalElem > 0) Conn_Elem = new int[size*MaxLocalElem*NODES_PER_ELEMENT];
  }

  /*--- Force the removal of all added periodic elements (use global index).
   First, we isolate and create a list of all added periodic points, excluding
   those that we part of the original domain (we want these to be in the
   output files). ---*/

  vector<unsigned long> Added_Periodic;
  Added_Periodic.clear();
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) {
      SendRecv = config->GetMarker_All_SendRecv(iMarker);
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        if ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() > 0) &&
            (geometry->vertex[iMarker][iVertex]->GetRotation_Type() % 2 == 0) &&
            (SendRecv < 0)) {
          Added_Periodic.push_back(geometry->node[iPoint]->GetGlobalIndex());
        }
      }
    }
  }

  /*--- Now we communicate this information to all processors, so that they
   can force the removal of these particular nodes by flagging them as halo
   points. In general, this should be a small percentage of the total mesh,
   so the communication/storage costs here shouldn't be prohibitive. ---*/

  /*--- First communicate the number of points that each rank has found ---*/
  unsigned long nAddedPeriodic = 0, maxAddedPeriodic = 0;
  unsigned long Buffer_Send_nAddedPeriodic[1], *Buffer_Recv_nAddedPeriodic = NULL;
  Buffer_Recv_nAddedPeriodic = new unsigned long[size];

  nAddedPeriodic = Added_Periodic.size();
  Buffer_Send_nAddedPeriodic[0] = nAddedPeriodic;

#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&nAddedPeriodic, &maxAddedPeriodic, 1, MPI_UNSIGNED_LONG,
                     MPI_MAX, MPI_COMM_WORLD);
  SU2_MPI::Allgather(&Buffer_Send_nAddedPeriodic, 1, MPI_UNSIGNED_LONG,
                     Buffer_Recv_nAddedPeriodic,  1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
#else
  maxAddedPeriodic = nAddedPeriodic;
  Buffer_Recv_nAddedPeriodic[0] = Buffer_Send_nAddedPeriodic[0];
#endif

  /*--- Communicate the global index values of all added periodic nodes. ---*/
  unsigned long *Buffer_Send_AddedPeriodic = new unsigned long[maxAddedPeriodic];
  unsigned long *Buffer_Recv_AddedPeriodic = new unsigned long[size*maxAddedPeriodic];

  for (iPoint = 0; iPoint < Added_Periodic.size(); iPoint++) {
    Buffer_Send_AddedPeriodic[iPoint] = Added_Periodic[iPoint];
  }

  /*--- Gather the element connectivity information. All processors will now
   have a copy of the global index values for all added periodic points. ---*/

#ifdef HAVE_MPI
  SU2_MPI::Allgather(Buffer_Send_AddedPeriodic, maxAddedPeriodic, MPI_UNSIGNED_LONG,
                     Buffer_Recv_AddedPeriodic, maxAddedPeriodic, MPI_UNSIGNED_LONG,
                     MPI_COMM_WORLD);
#else
  for (iPoint = 0; iPoint < maxAddedPeriodic; iPoint++) Buffer_Recv_AddedPeriodic[iPoint] = Buffer_Send_AddedPeriodic[iPoint];
#endif

  /*--- Search all send/recv boundaries on this partition for halo cells. In
   particular, consider only the recv conditions (these are the true halo
   nodes). Check the ranks of the processors that are communicating and
   choose to keep only the halo cells from the higher rank processor. Here,
   we are also choosing to keep periodic nodes that were part of the original
   domain. We will check the communicated list of added periodic points. ---*/

  int *Local_Halo = new int[geometry->GetnPoint()];
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
    Local_Halo[iPoint] = !geometry->node[iPoint]->GetDomain();

  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) {
      SendRecv = config->GetMarker_All_SendRecv(iMarker);
      RecvFrom = abs(SendRecv)-1;

      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        iGlobal_Index = geometry->node[iPoint]->GetGlobalIndex();

        /*--- We need to keep one copy of overlapping halo cells. ---*/
        notHalo = ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() == 0) &&
                   (SendRecv < 0) && (rank > RecvFrom));

        /*--- We want to keep the periodic nodes that were part of the original domain ---*/
        notPeriodic = ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() > 0) &&
                       (geometry->vertex[iMarker][iVertex]->GetRotation_Type() % 2 == 1) &&
                       (SendRecv < 0));

        /*--- Lastly, check that this isn't an added periodic point that
         we will forcibly remove. Use the communicated list of these points. ---*/
        addedPeriodic = false; kPoint = 0;
        for (iProcessor = 0; iProcessor < size; iProcessor++) {
          for (jPoint = 0; jPoint < Buffer_Recv_nAddedPeriodic[iProcessor]; jPoint++) {
            if (iGlobal_Index == Buffer_Recv_AddedPeriodic[kPoint+jPoint])
              addedPeriodic = true;
          }
          /*--- Adjust jNode to index of next proc's data in the buffers. ---*/
          kPoint = (iProcessor+1)*maxAddedPeriodic;
        }

        /*--- If we found either of these types of nodes, flag them to be kept. ---*/
        if ((notHalo || notPeriodic) && !addedPeriodic) {
          Local_Halo[iPoint] = false;
        }
      }
    }
  }

  /*--- Loop over all elements in this partition and load the
   elements of the current type into the buffer to be sent to
   the master node. ---*/
  jNode = 0; jElem = 0;
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
    if (config->GetMarker_All_Plotting(iMarker) == YES)
      for (iElem = 0; iElem < geometry->GetnElem_Bound(iMarker); iElem++) {

        if (geometry->bound[iMarker][iElem]->GetVTK_Type() == Elem_Type) {

          /*--- Loop over all nodes in this element and load the
           connectivity into the send buffer. ---*/

          Buffer_Send_Halo[jElem] = false;
          for (iNode = 0; iNode < NODES_PER_ELEMENT; iNode++) {

            /*--- Store the global index values directly. ---*/

            iPoint = geometry->bound[iMarker][iElem]->GetNode(iNode);
            Buffer_Send_Elem[jNode] = geometry->node[iPoint]->GetGlobalIndex();

            /*--- Check if this is a halo node. If so, flag this element
             as a halo cell. We will use this later to sort and remove
             any duplicates from the connectivity list. ---*/

            if (Local_Halo[iPoint])
              Buffer_Send_Halo[jElem] = true;

            /*--- Increment jNode as the counter. We need this because iElem
             may include other elements that we skip over during this loop. ---*/

            jNode++;
          }
          jElem++;
        }
      }

  /*--- Gather the element connectivity information. ---*/

#ifdef HAVE_MPI
  SU2_MPI::Gather(Buffer_Send_Elem, nBuffer_Scalar, MPI_UNSIGNED_LONG, Buffer_Recv_Elem, nBuffer_Scalar, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
  SU2_MPI::Gather(Buffer_Send_Halo, MaxLocalElem, MPI_UNSIGNED_SHORT, Buffer_Recv_Halo, MaxLocalElem, MPI_UNSIGNED_SHORT, MASTER_NODE, MPI_COMM_WORLD);
#else
  for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Elem[iPoint] = Buffer_Send_Elem[iPoint];
  for (iPoint = 0; iPoint < MaxLocalElem; iPoint++) Buffer_Recv_Halo[iPoint] = Buffer_Send_Halo[iPoint];
#endif

  /*--- The master node unpacks and sorts the connectivity. ---*/

  if (rank == MASTER_NODE) {

    /*---  We need to remove any duplicate elements (halo cells) that
     exist on multiple partitions. Start by initializing all elements
     to the "write" state by using a boolean array. ---*/

    Write_Elem = new bool[size*MaxLocalElem];
    for (iElem = 0; iElem < size*MaxLocalElem; iElem++) {
      Write_Elem[iElem] = true;
    }

    /*--- Remove the rind layer from the solution only if requested ---*/

    if (!Wrt_Halo) {

      /*--- Loop for flagging duplicate elements so that they are not
       included in the final connectivity list. ---*/

      kElem = 0;
      for (iProcessor = 0; iProcessor < size; iProcessor++) {
        for (iElem = 0; iElem < Buffer_Recv_nElem[iProcessor]; iElem++) {

          /*--- Check if this element was marked as a halo. ---*/
          if (Buffer_Recv_Halo[kElem+iElem])
            Write_Elem[kElem+iElem] = false;

        }
        kElem = (iProcessor+1)*MaxLocalElem;
      }
    }

    /*--- Store the unique connectivity list for this element type. ---*/

    jNode = 0; kNode = 0; jElem = 0; nElem_Total = 0;
    for (iProcessor = 0; iProcessor < size; iProcessor++) {
      for (iElem = 0; iElem < Buffer_Recv_nElem[iProcessor]; iElem++) {

        /*--- Only write the elements that were flagged for it. ---*/
        if (Write_Elem[jElem+iElem]) {

          /*--- Increment total count for this element type ---*/
          nElem_Total++;

          /*--- Get global index, then loop over each variable and store.
           Note that we are adding one to the index value because CGNS/Tecplot
           use 1-based indexing.---*/

          for (iNode = 0; iNode < NODES_PER_ELEMENT; iNode++) {
            Conn_Elem[kNode] = (int)Buffer_Recv_Elem[jNode+iElem*NODES_PER_ELEMENT+iNode] + 1;
            kNode++;
          }
        }
      }
      /*--- Adjust jNode to index of next proc's data in the buffers. ---*/
      jElem = (iProcessor+1)*MaxLocalElem;
      jNode = (iProcessor+1)*nBuffer_Scalar;
    }
  }

  /*--- Immediately release the temporary buffers. ---*/
  delete [] Buffer_Send_Elem;
  delete [] Buffer_Send_Halo;
  delete [] Buffer_Recv_nAddedPeriodic;
  delete [] Buffer_Send_AddedPeriodic;
  delete [] Buffer_Recv_AddedPeriodic;
  delete [] Local_Halo;
  if (rank == MASTER_NODE) {
    delete [] Buffer_Recv_nElem;
    delete [] Buffer_Recv_Elem;
    delete [] Buffer_Recv_Halo;
    delete [] Write_Elem;
  }

  /*--- Store the particular global element count in the class data,
   and set the class data pointer to the connectivity array. ---*/

  if (rank == MASTER_NODE) {
    switch (Elem_Type) {
      case LINE:
        nGlobal_Line = nElem_Total;
        if (nGlobal_Line > 0) Conn_Line = Conn_Elem;
        break;
      case TRIANGLE:
        nGlobal_BoundTria = nElem_Total;
        if (nGlobal_BoundTria > 0) Conn_BoundTria = Conn_Elem;
        break;
      case QUADRILATERAL:
        nGlobal_BoundQuad = nElem_Total;
        if (nGlobal_BoundQuad > 0) Conn_BoundQuad = Conn_Elem;
        break;
      default:
        SU2_MPI::Error("Unrecognized element type", CURRENT_FUNCTION);
        break;
    }
  }

}

void COutputLegacy::MergeSolution(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned short val_iZone) {

  unsigned short Kind_Solver  = config->GetKind_Solver();
  unsigned short iVar = 0, jVar = 0, FirstIndex = NONE, SecondIndex = NONE, ThirdIndex = NONE;
  unsigned short nVar_First = 0, nVar_Second = 0, nVar_Third = 0;
  unsigned short iVar_GridVel = 0, iVar_PressCp = 0, iVar_Lam = 0, iVar_MachMean = 0,
  iVar_ViscCoeffs = 0, iVar_HeatCoeffs = 0, iVar_Sens = 0, iVar_Extra = 0, iVar_Eddy = 0, iVar_Sharp = 0,
  iVar_FEA_Vel = 0, iVar_FEA_Accel = 0, iVar_FEA_Stress = 0, iVar_FEA_Stress_3D = 0,
  iVar_FEA_Extra = 0, iVar_SensDim = 0;
  unsigned long iPoint = 0, jPoint = 0, iVertex = 0, iMarker = 0;
  su2double Gas_Constant, Mach2Vel, Mach_Motion, RefDensity, RefPressure = 0.0, factor = 0.0;

  su2double *Aux_Frict_x = NULL, *Aux_Frict_y = NULL, *Aux_Frict_z = NULL, *Aux_Heat = NULL, *Aux_yPlus = NULL, *Aux_Sens = NULL;

  unsigned short CurrentIndex;
  int *Local_Halo;
  unsigned long Buffer_Send_nPoint[1], *Buffer_Recv_nPoint = NULL;
  unsigned long nLocalPoint = 0, MaxLocalPoint = 0;
  unsigned long iGlobal_Index = 0, nBuffer_Scalar = 0;
  bool Wrt_Halo = config->GetWrt_Halo(), isPeriodic;

  int iProcessor;

  bool dynamic_grid = config->GetDynamic_Grid();
  bool compressible   = (config->GetKind_Regime() == COMPRESSIBLE);
  bool incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  bool transition     = (config->GetKind_Trans_Model() == LM);
  bool flow           = (( config->GetKind_Solver() == EULER             ) ||
                         ( config->GetKind_Solver() == NAVIER_STOKES     ) ||
                         ( config->GetKind_Solver() == RANS              ) ||
                         ( config->GetKind_Solver() == INC_EULER         ) ||
                         ( config->GetKind_Solver() == INC_NAVIER_STOKES ) ||
                         ( config->GetKind_Solver() == INC_RANS          ) ||
                         ( config->GetKind_Solver() == FEM_EULER         ) ||
                         ( config->GetKind_Solver() == FEM_NAVIER_STOKES ) ||
                         ( config->GetKind_Solver() == FEM_RANS          ) ||
                         ( config->GetKind_Solver() == FEM_LES           ) ||
                         ( config->GetKind_Solver() == ADJ_EULER         ) ||
                         ( config->GetKind_Solver() == ADJ_NAVIER_STOKES ) ||
                         ( config->GetKind_Solver() == ADJ_RANS          )   );
  bool fem = (config->GetKind_Solver() == FEM_ELASTICITY);

  unsigned short iDim;
  unsigned short nDim = geometry->GetnDim();
  su2double RefArea = config->GetRefArea();
  su2double Gamma = config->GetGamma();
  su2double RefVel2, *Normal, Area;

  /*--- Set the non-dimensionalization ---*/
  if (flow) {
    if (dynamic_grid) {
      Gas_Constant = config->GetGas_ConstantND();
      Mach2Vel = sqrt(Gamma*Gas_Constant*config->GetTemperature_FreeStreamND());
      Mach_Motion = config->GetMach_Motion();
      RefVel2 = (Mach_Motion*Mach2Vel)*(Mach_Motion*Mach2Vel);
    }
    else {
      RefVel2 = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        RefVel2  += solver[FLOW_SOL]->GetVelocity_Inf(iDim)*solver[FLOW_SOL]->GetVelocity_Inf(iDim);
    }
    RefDensity  = solver[FLOW_SOL]->GetDensity_Inf();
    RefPressure = solver[FLOW_SOL]->GetPressure_Inf();
    factor = 1.0 / (0.5*RefDensity*RefArea*RefVel2);
  }

  /*--- Prepare send buffers for the conservative variables. Need to
   find the total number of conservative variables and also the
   index for their particular solution container. ---*/

  switch (Kind_Solver) {
    case EULER : case NAVIER_STOKES: case INC_EULER: case INC_NAVIER_STOKES: FirstIndex = FLOW_SOL; if(config->GetWeakly_Coupled_Heat()) SecondIndex = HEAT_SOL; else SecondIndex = NONE; ThirdIndex = NONE; break;
    case RANS : case INC_RANS: FirstIndex = FLOW_SOL; SecondIndex = TURB_SOL; if (transition) ThirdIndex=TRANS_SOL; else ThirdIndex = NONE;  if(config->GetWeakly_Coupled_Heat()) ThirdIndex = HEAT_SOL; else ThirdIndex = NONE; break;
    case FEM_EULER : case FEM_NAVIER_STOKES: case FEM_LES: FirstIndex = FLOW_SOL; SecondIndex = NONE; ThirdIndex = NONE; break;
    case FEM_RANS : FirstIndex = FLOW_SOL; SecondIndex = TURB_SOL; break;
    case FEM_ELASTICITY: FirstIndex = FEA_SOL; SecondIndex = NONE; ThirdIndex = NONE; break;
    case ADJ_EULER : case ADJ_NAVIER_STOKES : FirstIndex = ADJFLOW_SOL; SecondIndex = NONE; ThirdIndex = NONE; break;
    case ADJ_RANS : FirstIndex = ADJFLOW_SOL; if (config->GetFrozen_Visc_Cont()) SecondIndex = NONE; else SecondIndex = ADJTURB_SOL; ThirdIndex = NONE; break;
    case DISC_ADJ_EULER: case DISC_ADJ_NAVIER_STOKES: case DISC_ADJ_INC_EULER: case DISC_ADJ_INC_NAVIER_STOKES: FirstIndex = ADJFLOW_SOL; SecondIndex = NONE; ThirdIndex = NONE; break;
    case DISC_ADJ_RANS: case DISC_ADJ_INC_RANS: FirstIndex = ADJFLOW_SOL; if (config->GetFrozen_Visc_Disc()) SecondIndex = NONE; else SecondIndex = ADJTURB_SOL; ThirdIndex = NONE; break;
    case DISC_ADJ_FEM: FirstIndex = ADJFEA_SOL; SecondIndex = NONE; ThirdIndex = NONE; break;
    default: SecondIndex = NONE; ThirdIndex = NONE; break;
  }

  nVar_First = solver[FirstIndex]->GetnVar();
  if (SecondIndex != NONE) nVar_Second = solver[SecondIndex]->GetnVar();
  if (ThirdIndex != NONE) nVar_Third = solver[ThirdIndex]->GetnVar();
  nVar_Consv = nVar_First + nVar_Second + nVar_Third;
  nVar_Total = nVar_Consv;


  /*--- Add the limiters ---*/

  if (config->GetWrt_Limiters()) nVar_Total += nVar_Consv;

  /*--- Add the residuals ---*/

  if (config->GetWrt_Residuals()) nVar_Total += nVar_Consv;

  /*--- Add the grid velocity to the restart file for the unsteady adjoint ---*/

  if (dynamic_grid && !fem) {
    iVar_GridVel = nVar_Total;
    if (geometry->GetnDim() == 2) nVar_Total += 2;
    else if (geometry->GetnDim() == 3) nVar_Total += 3;
  }

  /*--- Add Pressure, Temperature, Cp, Mach to the restart file ---*/

  if (Kind_Solver == EULER     || Kind_Solver == NAVIER_STOKES     || Kind_Solver == RANS ||
      Kind_Solver == INC_EULER || Kind_Solver == INC_NAVIER_STOKES || Kind_Solver == INC_RANS ||
      Kind_Solver == FEM_EULER || Kind_Solver == FEM_NAVIER_STOKES || Kind_Solver == FEM_RANS ||
      Kind_Solver == FEM_LES) {
    iVar_PressCp = nVar_Total; nVar_Total += 3;
    iVar_MachMean = nVar_Total; nVar_Total += 1;
  }

  /*--- Add Laminar Viscosity, Skin Friction, Heat Flux, & yPlus to the restart file ---*/

  if (Kind_Solver == NAVIER_STOKES     || Kind_Solver == RANS     ||
      Kind_Solver == INC_NAVIER_STOKES     || Kind_Solver == INC_RANS     ||
      Kind_Solver == FEM_NAVIER_STOKES || Kind_Solver == FEM_RANS || Kind_Solver == FEM_LES) {
    iVar_Lam = nVar_Total;
    nVar_Total += 1;
    iVar_ViscCoeffs = nVar_Total;
    if (geometry->GetnDim() == 2) nVar_Total += 2;
    else if (geometry->GetnDim() == 3) nVar_Total += 3;
    iVar_HeatCoeffs = nVar_Total;
    nVar_Total += 2;
  }

  /*--- Add Eddy Viscosity to the restart file ---*/

  if (Kind_Solver == RANS || Kind_Solver == FEM_RANS || Kind_Solver == FEM_LES || Kind_Solver == INC_RANS) {
    iVar_Eddy = nVar_Total; nVar_Total += 1;
  }

  /*--- Add Sharp edges to the restart file ---*/

  if (config->GetWrt_SharpEdges()) {
    if (((Kind_Solver == EULER) || (Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) ||
        ((Kind_Solver == INC_EULER) || (Kind_Solver == INC_NAVIER_STOKES) || (Kind_Solver == INC_RANS)) ||
        ((Kind_Solver == FEM_EULER) || (Kind_Solver == FEM_NAVIER_STOKES) || (Kind_Solver == FEM_RANS) || (Kind_Solver == FEM_LES)))  {
      iVar_Sharp = nVar_Total; nVar_Total += 1;
    }
  }


  if (( Kind_Solver == ADJ_EULER              ) ||
      ( Kind_Solver == ADJ_NAVIER_STOKES      ) ||
      ( Kind_Solver == ADJ_RANS               )) {
    iVar_Sens   = nVar_Total; nVar_Total += 2;
  }

  if (Kind_Solver == FEM_ELASTICITY)  {
    /*--- If the analysis is dynamic... ---*/
    if (config->GetTime_Domain()) {
      /*--- Velocities ---*/
      iVar_FEA_Vel = nVar_Total;
      if (geometry->GetnDim() == 2) nVar_Total += 2;
      else if (geometry->GetnDim() == 3) nVar_Total += 3;
      /*--- Accelerations ---*/
      iVar_FEA_Accel = nVar_Total;
      if (geometry->GetnDim() == 2) nVar_Total += 2;
      else if (geometry->GetnDim() == 3) nVar_Total += 3;
    }
    iVar_FEA_Stress  = nVar_Total; nVar_Total += 3;
    if (geometry->GetnDim() == 3) {iVar_FEA_Stress_3D = nVar_Total; nVar_Total += 3;}
    iVar_FEA_Extra = nVar_Total; nVar_Total += 1;
  }

  if ((Kind_Solver == DISC_ADJ_EULER)         ||
      (Kind_Solver == DISC_ADJ_NAVIER_STOKES) ||
      (Kind_Solver == DISC_ADJ_RANS) ||
      (Kind_Solver == DISC_ADJ_INC_EULER)         ||
      (Kind_Solver == DISC_ADJ_INC_NAVIER_STOKES) ||
      (Kind_Solver == DISC_ADJ_INC_RANS)) {
    iVar_Sens    = nVar_Total; nVar_Total += 1;
    iVar_SensDim = nVar_Total; nVar_Total += nDim;
  }

  if ((Kind_Solver == DISC_ADJ_FEM) && (config->GetFSI_Simulation())){
    iVar_FEA_Extra = nVar_Total; nVar_Total += 2;
  }

  if (config->GetExtraOutput()) {
    if (Kind_Solver == RANS || Kind_Solver == INC_RANS) {
      iVar_Extra  = nVar_Total; nVar_Extra  = solver[TURB_SOL]->GetnOutputVariables(); nVar_Total += nVar_Extra;
    }
  }



  Local_Halo = new int[geometry->GetnPoint()];
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
    Local_Halo[iPoint] = !geometry->node[iPoint]->GetDomain();

  /*--- Search all send/recv boundaries on this partition for any periodic
   nodes that were part of the original domain. We want to recover these
   for visualization purposes. ---*/

  if (Wrt_Halo) {
    nLocalPoint = geometry->GetnPoint();
  } else {
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) {

        /*--- Checking for less than or equal to the rank, because there may
         be some periodic halo nodes that send info to the same rank. ---*/

        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          isPeriodic = ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() > 0) &&
                        (geometry->vertex[iMarker][iVertex]->GetRotation_Type() % 2 == 1));
          if (isPeriodic) Local_Halo[iPoint] = false;
        }
      }
    }

    /*--- Sum total number of nodes that belong to the domain ---*/

    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
      if (Local_Halo[iPoint] == false)
        nLocalPoint++;

  }
  Buffer_Send_nPoint[0] = nLocalPoint;

  /*--- Each processor sends its local number of nodes to the master. ---*/

  if (rank == MASTER_NODE) Buffer_Recv_nPoint = new unsigned long[size];

#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&nLocalPoint, &MaxLocalPoint, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
  SU2_MPI::Gather(&Buffer_Send_nPoint, 1, MPI_UNSIGNED_LONG, Buffer_Recv_nPoint, 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
#else
  MaxLocalPoint = nLocalPoint;
  Buffer_Recv_nPoint[0] = Buffer_Send_nPoint[0];
#endif

  nBuffer_Scalar = MaxLocalPoint;

  /*--- Send and Recv buffers. ---*/

  su2double *Buffer_Send_Var = new su2double[MaxLocalPoint];
  su2double *Buffer_Recv_Var = NULL;

  su2double *Buffer_Send_Res = new su2double[MaxLocalPoint];
  su2double *Buffer_Recv_Res = NULL;

  su2double *Buffer_Send_Vol = new su2double[MaxLocalPoint];
  su2double *Buffer_Recv_Vol = NULL;

  unsigned long *Buffer_Send_GlobalIndex = new unsigned long[MaxLocalPoint];
  unsigned long *Buffer_Recv_GlobalIndex = NULL;

  /*--- Auxiliary vectors for surface coefficients ---*/

  if (((Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS))  ||
      ((Kind_Solver == INC_NAVIER_STOKES) || (Kind_Solver == INC_RANS)) ||
      ((Kind_Solver == FEM_NAVIER_STOKES) || (Kind_Solver == FEM_RANS) || (Kind_Solver == FEM_LES)))  {
    Aux_Frict_x = new su2double[geometry->GetnPoint()];
    Aux_Frict_y = new su2double[geometry->GetnPoint()];
    Aux_Frict_z = new su2double[geometry->GetnPoint()];
    Aux_Heat  = new su2double[geometry->GetnPoint()];
    Aux_yPlus = new su2double[geometry->GetnPoint()];
  }

  if ((Kind_Solver == ADJ_EULER) ||
      (Kind_Solver == ADJ_NAVIER_STOKES) ||
      (Kind_Solver == ADJ_RANS)  ||
      (Kind_Solver == DISC_ADJ_EULER) ||
      (Kind_Solver == DISC_ADJ_NAVIER_STOKES) ||
      (Kind_Solver == DISC_ADJ_RANS)  ||
      (Kind_Solver == DISC_ADJ_INC_EULER) ||
      (Kind_Solver == DISC_ADJ_INC_NAVIER_STOKES) ||
      (Kind_Solver == DISC_ADJ_INC_RANS)) {
    Aux_Sens = new su2double[geometry->GetnPoint()];
  }

  /*--- Prepare the receive buffers in the master node only. ---*/

  if (rank == MASTER_NODE) {

    Buffer_Recv_Var = new su2double[size*MaxLocalPoint];
    Buffer_Recv_Res = new su2double[size*MaxLocalPoint];
    Buffer_Recv_Vol = new su2double[size*MaxLocalPoint];
    Buffer_Recv_GlobalIndex = new unsigned long[size*MaxLocalPoint];

    /*--- Sum total number of nodes to be written and allocate arrays ---*/
    nGlobal_Poin = 0;
    for (iProcessor = 0; iProcessor < size; iProcessor++) {
      nGlobal_Poin += Buffer_Recv_nPoint[iProcessor];
    }
    Data = new su2double*[nVar_Total];
    for (iVar = 0; iVar < nVar_Total; iVar++) {
      Data[iVar] = new su2double[nGlobal_Poin];
    }
  }

  /*--- Main communication routine. Loop over each variable that has
   been requested by the user and perform the MPI comm. Temporary
   1-D buffers are used to send the solution for each variable at all
   nodes on each partition to the master node. These are then unpacked
   by the master and sorted by global index in one large n-dim. array. ---*/

  for (iVar = 0; iVar < nVar_Consv; iVar++) {

    /*--- Logic for which solution class to draw from. ---*/

    jVar = iVar;
    CurrentIndex = FirstIndex;
    if ((SecondIndex != NONE) && (iVar > nVar_First-1)) {
      jVar = iVar - nVar_First;
      CurrentIndex = SecondIndex;
    }
    if ((SecondIndex != NONE) && (ThirdIndex != NONE) && (iVar > (nVar_First + nVar_Second-1))) {
      jVar = iVar - nVar_First - nVar_Second;
      CurrentIndex = ThirdIndex;
    }

    /*--- Loop over this partition to collect the current variable ---*/

    jPoint = 0;
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {

      /*--- Check for halos & write only if requested ---*/

      if (!Local_Halo[iPoint] || Wrt_Halo) {

        /*--- Get this variable into the temporary send buffer. ---*/

        Buffer_Send_Var[jPoint] = solver[CurrentIndex]->GetNodes()->GetSolution(iPoint, jVar);


        if (config->GetWrt_Limiters()) {
          Buffer_Send_Vol[jPoint] = solver[CurrentIndex]->GetNodes()->GetLimiter_Primitive(iPoint, jVar);
        }

        if (config->GetWrt_Residuals()) {
          if (!config->GetDiscrete_Adjoint()) {
            Buffer_Send_Res[jPoint] = solver[CurrentIndex]->LinSysRes(iPoint, jVar);
          } else {
            Buffer_Send_Res[jPoint] = solver[CurrentIndex]->GetNodes()->GetSolution(iPoint, jVar) -
                                      solver[CurrentIndex]->GetNodes()->GetSolution_Old(iPoint, jVar);
          }
        }


        /*--- Only send/recv the volumes & global indices during the first loop ---*/

        if (iVar == 0) {
          Buffer_Send_GlobalIndex[jPoint] = geometry->node[iPoint]->GetGlobalIndex();
        }

        jPoint++;

      }
    }

    /*--- Gather the data on the master node. ---*/

#ifdef HAVE_MPI
    SU2_MPI::Gather(Buffer_Send_Var, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Var, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
#else
    for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Var[iPoint] = Buffer_Send_Var[iPoint];
#endif

    if (config->GetWrt_Limiters()) {
#ifdef HAVE_MPI
      SU2_MPI::Gather(Buffer_Send_Vol, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Vol, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
#else
      for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Vol[iPoint] = Buffer_Send_Vol[iPoint];
#endif
    }

    if (config->GetWrt_Residuals()) {
#ifdef HAVE_MPI
      SU2_MPI::Gather(Buffer_Send_Res, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Res, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
#else
      for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Res[iPoint] = Buffer_Send_Res[iPoint];
#endif
    }



    if (iVar == 0) {
#ifdef HAVE_MPI
      SU2_MPI::Gather(Buffer_Send_GlobalIndex, nBuffer_Scalar, MPI_UNSIGNED_LONG, Buffer_Recv_GlobalIndex, nBuffer_Scalar, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
#else
      for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_GlobalIndex[iPoint] = Buffer_Send_GlobalIndex[iPoint];
#endif
    }

    /*--- The master node unpacks and sorts this variable by global index ---*/

    if (rank == MASTER_NODE) {
      jPoint = 0;
      for (iProcessor = 0; iProcessor < size; iProcessor++) {
        for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++) {

          /*--- Get global index, then loop over each variable and store ---*/

          iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];

          Data[iVar][iGlobal_Index] = Buffer_Recv_Var[jPoint];


          if (config->GetWrt_Limiters()) {
            Data[iVar+nVar_Consv][iGlobal_Index] = Buffer_Recv_Vol[jPoint];
          }

          if (config->GetWrt_Residuals()) {
            unsigned short ExtraIndex;
            ExtraIndex = nVar_Consv;
            if (config->GetWrt_Limiters()) ExtraIndex = 2*nVar_Consv;
            Data[iVar+ExtraIndex][iGlobal_Index] = Buffer_Recv_Res[jPoint];
          }



          jPoint++;
        }
        /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
        jPoint = (iProcessor+1)*nBuffer_Scalar;
      }
    }

  }

    /*--- Additional communication routine for the grid velocity. Note that
     we are reusing the same temporary buffers from above for efficiency.
     Also, in the future more routines like this could be used to write
     an arbitrary number of additional variables to the file. ---*/

    if (dynamic_grid && !fem) {

      /*--- Loop over this partition to collect the current variable ---*/

      jPoint = 0; su2double *Grid_Vel;
      for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {

        /*--- Check for halos & write only if requested ---*/

        if (!Local_Halo[iPoint] || Wrt_Halo) {

          /*--- Load buffers with the three grid velocity components. ---*/

          Grid_Vel = geometry->node[iPoint]->GetGridVel();
          Buffer_Send_Var[jPoint] = Grid_Vel[0];
          Buffer_Send_Res[jPoint] = Grid_Vel[1];
          if (geometry->GetnDim() == 3) Buffer_Send_Vol[jPoint] = Grid_Vel[2];
          jPoint++;
        }
      }

      /*--- Gather the data on the master node. ---*/

#ifdef HAVE_MPI
      SU2_MPI::Gather(Buffer_Send_Var, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Var, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
      SU2_MPI::Gather(Buffer_Send_Res, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Res, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
      if (geometry->GetnDim() == 3) {
        SU2_MPI::Gather(Buffer_Send_Vol, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Vol, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
      }
#else
      for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Var[iPoint] = Buffer_Send_Var[iPoint];
      for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Res[iPoint] = Buffer_Send_Res[iPoint];
      if (geometry->GetnDim() == 3) {
        for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Vol[iPoint] = Buffer_Send_Vol[iPoint];
      }
#endif

      /*--- The master node unpacks and sorts this variable by global index ---*/

      if (rank == MASTER_NODE) {
        jPoint = 0; iVar = iVar_GridVel;
        for (iProcessor = 0; iProcessor < size; iProcessor++) {
          for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++) {

            /*--- Get global index, then loop over each variable and store ---*/

            iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
            Data[iVar][iGlobal_Index]   = Buffer_Recv_Var[jPoint];
            Data[iVar+1][iGlobal_Index] = Buffer_Recv_Res[jPoint];
            if (geometry->GetnDim() == 3)
              Data[iVar+2][iGlobal_Index] = Buffer_Recv_Vol[jPoint];
            jPoint++;
          }

          /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/

          jPoint = (iProcessor+1)*nBuffer_Scalar;
        }
      }
    }

    /*--- Communicate Pressure, Cp, and Mach ---*/

    if (((Kind_Solver == EULER) || (Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) ||
        ((Kind_Solver == INC_EULER) || (Kind_Solver == INC_NAVIER_STOKES) || (Kind_Solver == INC_RANS)) ||
      ((Kind_Solver == FEM_EULER) || (Kind_Solver == FEM_NAVIER_STOKES) || (Kind_Solver == FEM_RANS) || (Kind_Solver == FEM_LES)))  {

      /*--- First, loop through the mesh in order to find and store the
       value of the coefficient of pressure at any surface nodes. They
       will be placed in an auxiliary vector and then communicated like
       all other volumetric variables. ---*/

      /*--- Loop over this partition to collect the current variable ---*/

      jPoint = 0;
      for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {

        /*--- Check for halos & write only if requested ---*/

        if (!Local_Halo[iPoint] || Wrt_Halo) {

          /*--- Load buffers with the pressure, Cp, and mach variables. ---*/

          Buffer_Send_Var[jPoint] = solver[FLOW_SOL]->GetNodes()->GetPressure(iPoint);
          if (compressible){
            Buffer_Send_Res[jPoint] = solver[FLOW_SOL]->GetNodes()->GetTemperature(iPoint);
          } else{
            Buffer_Send_Res[jPoint] =  0.0;
          }
          Buffer_Send_Vol[jPoint] = (solver[FLOW_SOL]->GetNodes()->GetPressure(iPoint) - RefPressure)*factor*RefArea;

          jPoint++;
        }
      }

      /*--- Gather the data on the master node. ---*/

#ifdef HAVE_MPI
      SU2_MPI::Gather(Buffer_Send_Var, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Var, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
      SU2_MPI::Gather(Buffer_Send_Res, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Res, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
      SU2_MPI::Gather(Buffer_Send_Vol, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Vol, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
#else
      for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Var[iPoint] = Buffer_Send_Var[iPoint];
      for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Res[iPoint] = Buffer_Send_Res[iPoint];
      for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Vol[iPoint] = Buffer_Send_Vol[iPoint];
#endif

      /*--- The master node unpacks and sorts this variable by global index ---*/

      if (rank == MASTER_NODE) {
        jPoint = 0; iVar = iVar_PressCp;
        for (iProcessor = 0; iProcessor < size; iProcessor++) {
          for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++) {

            /*--- Get global index, then loop over each variable and store ---*/

            iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
            Data[iVar][iGlobal_Index]   = Buffer_Recv_Var[jPoint];
            Data[iVar+1][iGlobal_Index] = Buffer_Recv_Res[jPoint];
            Data[iVar+2][iGlobal_Index] = Buffer_Recv_Vol[jPoint];
            jPoint++;
          }

          /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/

          jPoint = (iProcessor+1)*nBuffer_Scalar;
        }
      }
    }

    /*--- Communicate Mach---*/

    if (((Kind_Solver == EULER) || (Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) ||
        ((Kind_Solver == INC_EULER) || (Kind_Solver == INC_NAVIER_STOKES) || (Kind_Solver == INC_RANS)) ||
      ((Kind_Solver == FEM_EULER) || (Kind_Solver == FEM_NAVIER_STOKES) || (Kind_Solver == FEM_RANS) || (Kind_Solver == FEM_LES)))  {

      /*--- Loop over this partition to collect the current variable ---*/

      jPoint = 0;
      for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {

        /*--- Check for halos & write only if requested ---*/

        if (!Local_Halo[iPoint] || Wrt_Halo) {

          /*--- Load buffers with the temperature and laminar viscosity variables. ---*/

          if (compressible) {
            Buffer_Send_Var[jPoint] = sqrt(solver[FLOW_SOL]->GetNodes()->GetVelocity2(iPoint))/
            solver[FLOW_SOL]->GetNodes()->GetSoundSpeed(iPoint);
          }
          if (incompressible) {
            Buffer_Send_Var[jPoint] = sqrt(solver[FLOW_SOL]->GetNodes()->GetVelocity2(iPoint))*config->GetVelocity_Ref()/
            sqrt(config->GetBulk_Modulus()/(solver[FLOW_SOL]->GetNodes()->GetDensity(iPoint)*config->GetDensity_Ref()));
          }
          jPoint++;
        }
      }

      /*--- Gather the data on the master node. ---*/

#ifdef HAVE_MPI
      SU2_MPI::Gather(Buffer_Send_Var, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Var, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
#else
      for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Var[iPoint] = Buffer_Send_Var[iPoint];
#endif

      /*--- The master node unpacks and sorts this variable by global index ---*/

      if (rank == MASTER_NODE) {
        jPoint = 0; iVar = iVar_MachMean;
        for (iProcessor = 0; iProcessor < size; iProcessor++) {
          for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++) {

            /*--- Get global index, then loop over each variable and store ---*/

            iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
            Data[iVar][iGlobal_Index]   = Buffer_Recv_Var[jPoint];
            jPoint++;
          }

          /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/

          jPoint = (iProcessor+1)*nBuffer_Scalar;
        }
      }
    }

    /*--- Laminar Viscosity ---*/

    if (((Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) ||
        ((Kind_Solver == INC_NAVIER_STOKES) || (Kind_Solver == INC_RANS)) ||
        ((Kind_Solver == FEM_NAVIER_STOKES) || (Kind_Solver == FEM_RANS) || (Kind_Solver == FEM_LES)))  {

      /*--- Loop over this partition to collect the current variable ---*/

      jPoint = 0;
      for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {

        /*--- Check for halos & write only if requested ---*/

        if (!Local_Halo[iPoint] || Wrt_Halo) {

          /*--- Load buffers with the temperature and laminar viscosity variables. ---*/

          Buffer_Send_Res[jPoint] = solver[FLOW_SOL]->GetNodes()->GetLaminarViscosity(iPoint);

          jPoint++;
        }
      }

      /*--- Gather the data on the master node. ---*/

#ifdef HAVE_MPI
      SU2_MPI::Gather(Buffer_Send_Res, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Res, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
#else
      for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Res[iPoint] = Buffer_Send_Res[iPoint];
#endif

      /*--- The master node unpacks and sorts this variable by global index ---*/

      if (rank == MASTER_NODE) {
        jPoint = 0; iVar = iVar_Lam;
        for (iProcessor = 0; iProcessor < size; iProcessor++) {
          for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++) {

            /*--- Get global index, then loop over each variable and store ---*/

            iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
            Data[iVar][iGlobal_Index] = Buffer_Recv_Res[jPoint];
            jPoint++;
          }

          /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/

          jPoint = (iProcessor+1)*nBuffer_Scalar;
        }
      }

      /*--- Communicate skin friction ---*/

      /*--- First, loop through the mesh in order to find and store the
       value of the viscous coefficients at any surface nodes. They
       will be placed in an auxiliary vector and then communicated like
       all other volumetric variables. ---*/

      for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
        Aux_Frict_x[iPoint] = 0.0;
        Aux_Frict_y[iPoint] = 0.0;
        Aux_Frict_z[iPoint] = 0.0;
      }
      for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
        if (config->GetMarker_All_Plotting(iMarker) == YES) {
          for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
            iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
            Aux_Frict_x[iPoint] = solver[FLOW_SOL]->GetCSkinFriction(iMarker, iVertex, 0);
            Aux_Frict_y[iPoint] = solver[FLOW_SOL]->GetCSkinFriction(iMarker, iVertex, 1);
            if (geometry->GetnDim() == 3) Aux_Frict_z[iPoint] = solver[FLOW_SOL]->GetCSkinFriction(iMarker, iVertex, 2);
          }
        }

      /*--- Loop over this partition to collect the current variable ---*/

      jPoint = 0;
      for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {

        /*--- Check for halos & write only if requested ---*/

        if (!Local_Halo[iPoint] || Wrt_Halo) {

          /*--- Load buffers with the three grid velocity components. ---*/

          Buffer_Send_Var[jPoint] = Aux_Frict_x[iPoint];
          Buffer_Send_Res[jPoint] = Aux_Frict_y[iPoint];
          if (geometry->GetnDim() == 3)
            Buffer_Send_Vol[jPoint] = Aux_Frict_z[iPoint];
          jPoint++;
        }
      }

      /*--- Gather the data on the master node. ---*/

#ifdef HAVE_MPI
      SU2_MPI::Gather(Buffer_Send_Var, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Var, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
      SU2_MPI::Gather(Buffer_Send_Res, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Res, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
      if (geometry->GetnDim() == 3) {
        SU2_MPI::Gather(Buffer_Send_Vol, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Vol, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
      }
#else
      for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++)
        Buffer_Recv_Var[iPoint] = Buffer_Send_Var[iPoint];
      for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++)
        Buffer_Recv_Res[iPoint] = Buffer_Send_Res[iPoint];
      if (geometry->GetnDim() == 3) {
        for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++)
          Buffer_Recv_Vol[iPoint] = Buffer_Send_Vol[iPoint];
      }
#endif

      /*--- The master node unpacks and sorts this variable by global index ---*/

      if (rank == MASTER_NODE) {
        jPoint = 0;
        iVar = iVar_ViscCoeffs;
        for (iProcessor = 0; iProcessor < size; iProcessor++) {
          for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++) {

            /*--- Get global index, then loop over each variable and store ---*/

            iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
            Data[iVar][iGlobal_Index] = Buffer_Recv_Var[jPoint];
            Data[iVar + 1][iGlobal_Index] = Buffer_Recv_Res[jPoint];
            if (geometry->GetnDim() == 3)
              Data[iVar + 2][iGlobal_Index] = Buffer_Recv_Vol[jPoint];
            jPoint++;
          }

          /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/

          jPoint = (iProcessor + 1) * nBuffer_Scalar;
        }
      }

      /*--- Communicate heat transfer, y+ ---*/

      /*--- First, loop through the mesh in order to find and store the
       value of the viscous coefficients at any surface nodes. They
       will be placed in an auxiliary vector and then communicated like
       all other volumetric variables. ---*/

      for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
        Aux_Heat[iPoint] = 0.0;
        Aux_yPlus[iPoint] = 0.0;
      }
      for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
        if (config->GetMarker_All_Plotting(iMarker) == YES) {
          for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
            iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
            Aux_Heat[iPoint] = solver[FLOW_SOL]->GetHeatFlux(iMarker, iVertex);
            Aux_yPlus[iPoint] = solver[FLOW_SOL]->GetYPlus(iMarker, iVertex);
          }
        }

      /*--- Loop over this partition to collect the current variable ---*/

      jPoint = 0;
      for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {

        /*--- Check for halos & write only if requested ---*/

        if (!Local_Halo[iPoint] || Wrt_Halo) {

          /*--- Load buffers with the skin friction, heat transfer, y+ variables. ---*/

          if (compressible) {
            Buffer_Send_Res[jPoint] = Aux_Heat[iPoint];
            Buffer_Send_Vol[jPoint] = Aux_yPlus[iPoint];
          }
          if (incompressible) {
            Buffer_Send_Res[jPoint] = Aux_Heat[iPoint];
            Buffer_Send_Vol[jPoint] = Aux_yPlus[iPoint];
          }
          jPoint++;
        }
      }

      /*--- Gather the data on the master node. ---*/

#ifdef HAVE_MPI
      SU2_MPI::Gather(Buffer_Send_Res, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Res, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
      SU2_MPI::Gather(Buffer_Send_Vol, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Vol, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
#else
      for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++)
        Buffer_Recv_Res[iPoint] = Buffer_Send_Res[iPoint];
      for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++)
        Buffer_Recv_Vol[iPoint] = Buffer_Send_Vol[iPoint];
#endif

      /*--- The master node unpacks and sorts this variable by global index ---*/

      if (rank == MASTER_NODE) {
        jPoint = 0;
        iVar = iVar_HeatCoeffs;

        for (iProcessor = 0; iProcessor < size; iProcessor++) {
          for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++) {

            /*--- Get global index, then loop over each variable and store ---*/

            iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
            Data[iVar + 0][iGlobal_Index] = Buffer_Recv_Res[jPoint];
            Data[iVar + 1][iGlobal_Index] = Buffer_Recv_Vol[jPoint];
            jPoint++;
          }

          /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/

          jPoint = (iProcessor + 1) * nBuffer_Scalar;
        }
      }
    }


    /*--- Communicate the Eddy Viscosity ---*/

    if (Kind_Solver == RANS || Kind_Solver == FEM_RANS || Kind_Solver == FEM_LES || Kind_Solver == INC_RANS) {

      /*--- Loop over this partition to collect the current variable ---*/

      jPoint = 0;
      for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {

        /*--- Check for halos & write only if requested ---*/

        if (!Local_Halo[iPoint] || Wrt_Halo) {

          /*--- Load buffers with the pressure and mach variables. ---*/

          Buffer_Send_Var[jPoint] = solver[FLOW_SOL]->GetNodes()->GetEddyViscosity(iPoint);

          jPoint++;
        }
      }

      /*--- Gather the data on the master node. ---*/

#ifdef HAVE_MPI
      SU2_MPI::Gather(Buffer_Send_Var, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Var, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
#else
      for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Var[iPoint] = Buffer_Send_Var[iPoint];
#endif

      /*--- The master node unpacks and sorts this variable by global index ---*/

      if (rank == MASTER_NODE) {
        jPoint = 0; iVar = iVar_Eddy;
        for (iProcessor = 0; iProcessor < size; iProcessor++) {
          for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++) {

            /*--- Get global index, then loop over each variable and store ---*/

            iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
            Data[iVar][iGlobal_Index] = Buffer_Recv_Var[jPoint];
            jPoint++;
          }

          /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/

          jPoint = (iProcessor+1)*nBuffer_Scalar;
        }
      }

    }

    /*--- Communicate the Sharp Edges ---*/

    if (config->GetWrt_SharpEdges()) {

      if (((Kind_Solver == EULER) || (Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) ||
          ((Kind_Solver == INC_EULER) || (Kind_Solver == INC_NAVIER_STOKES) || (Kind_Solver == INC_RANS)) ||
        ((Kind_Solver == FEM_EULER) || (Kind_Solver == FEM_NAVIER_STOKES) || (Kind_Solver == FEM_RANS) || (Kind_Solver == FEM_LES)))  {

        /*--- Loop over this partition to collect the current variable ---*/
        jPoint = 0;
        for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {

          /*--- Check for halos & write only if requested ---*/

          if (!Local_Halo[iPoint] || Wrt_Halo) {

            /*--- Load buffers with the pressure and mach variables. ---*/

            Buffer_Send_Var[jPoint] = geometry->node[iPoint]->GetSharpEdge_Distance();
            jPoint++;
          }
        }

        /*--- Gather the data on the master node. ---*/

#ifdef HAVE_MPI
        SU2_MPI::Gather(Buffer_Send_Var, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Var, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
#else
        for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Var[iPoint] = Buffer_Send_Var[iPoint];
#endif

        /*--- The master node unpacks and sorts this variable by global index ---*/

        if (rank == MASTER_NODE) {
          jPoint = 0; iVar = iVar_Sharp;
          for (iProcessor = 0; iProcessor < size; iProcessor++) {
            for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++) {

              /*--- Get global index, then loop over each variable and store ---*/

              iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
              Data[iVar][iGlobal_Index] = Buffer_Recv_Var[jPoint];
              jPoint++;
            }

            /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/

            jPoint = (iProcessor+1)*nBuffer_Scalar;
          }
        }
      }
    }

    /*--- Communicate the surface sensitivity ---*/

    if ((Kind_Solver == ADJ_EULER)         ||
        (Kind_Solver == ADJ_NAVIER_STOKES) ||
        (Kind_Solver == ADJ_RANS)          ||
        (Kind_Solver == DISC_ADJ_EULER)    ||
        (Kind_Solver == DISC_ADJ_NAVIER_STOKES) ||
        (Kind_Solver == DISC_ADJ_RANS)  ||
        (Kind_Solver == DISC_ADJ_INC_EULER)    ||
        (Kind_Solver == DISC_ADJ_INC_NAVIER_STOKES) ||
        (Kind_Solver == DISC_ADJ_INC_RANS)) {

      /*--- First, loop through the mesh in order to find and store the
       value of the surface sensitivity at any surface nodes. They
       will be placed in an auxiliary vector and then communicated like
       all other volumetric variables. ---*/

      for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) Aux_Sens[iPoint] = 0.0;
      for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
        if (config->GetMarker_All_Plotting(iMarker) == YES) {
          for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
            iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
            Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
            Area = 0.0;
            for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim];
            Area = sqrt (Area);
            Aux_Sens[iPoint] = solver[ADJFLOW_SOL]->GetCSensitivity(iMarker, iVertex)/Area;
          }
        }

      /*--- Loop over this partition to collect the current variable ---*/

      jPoint = 0;
      for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {

        /*--- Check for halos & write only if requested ---*/

        if (!Local_Halo[iPoint] || Wrt_Halo) {

          /*--- Load buffers with the skin friction, heat transfer, y+ variables. ---*/

          Buffer_Send_Var[jPoint] = Aux_Sens[iPoint];
          if ((config->GetKind_ConvNumScheme() == SPACE_CENTERED) && (!config->GetDiscrete_Adjoint()))
            Buffer_Send_Res[jPoint] = solver[ADJFLOW_SOL]->GetNodes()->GetSensor(iPoint, iPoint);
          if ((config->GetKind_ConvNumScheme() == SPACE_UPWIND) && (!config->GetDiscrete_Adjoint()))
            Buffer_Send_Res[jPoint] = solver[ADJFLOW_SOL]->GetNodes()->GetLimiter(iPoint, 0);

          jPoint++;
        }
      }

      /*--- Gather the data on the master node. ---*/

#ifdef HAVE_MPI
      SU2_MPI::Gather(Buffer_Send_Var, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Var, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
      if (!config->GetDiscrete_Adjoint())
        SU2_MPI::Gather(Buffer_Send_Res, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Res, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
#else
      for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Var[iPoint] = Buffer_Send_Var[iPoint];
      if (!config->GetDiscrete_Adjoint())
        for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Res[iPoint] = Buffer_Send_Res[iPoint];
#endif

      /*--- The master node unpacks and sorts this variable by global index ---*/

      if (rank == MASTER_NODE) {
        jPoint = 0; iVar = iVar_Sens;
        for (iProcessor = 0; iProcessor < size; iProcessor++) {
          for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++) {

            /*--- Get global index, then loop over each variable and store ---*/

            iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
            Data[iVar+0][iGlobal_Index] = Buffer_Recv_Var[jPoint];
            if (!config->GetDiscrete_Adjoint())
              Data[iVar+1][iGlobal_Index] = Buffer_Recv_Res[jPoint];
            jPoint++;
          }

          /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/

          jPoint = (iProcessor+1)*nBuffer_Scalar;
        }
      }
    }

    if ((Kind_Solver == DISC_ADJ_EULER)    ||
        (Kind_Solver == DISC_ADJ_NAVIER_STOKES) ||
        (Kind_Solver == DISC_ADJ_RANS) ||
        (Kind_Solver == DISC_ADJ_INC_EULER)    ||
        (Kind_Solver == DISC_ADJ_INC_NAVIER_STOKES) ||
        (Kind_Solver == DISC_ADJ_INC_RANS)) {
      /*--- Loop over this partition to collect the current variable ---*/

      jPoint = 0;
      for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {

        /*--- Check for halos & write only if requested ---*/

        if (!Local_Halo[iPoint] || Wrt_Halo) {

          /*--- Load buffers with the skin friction, heat transfer, y+ variables. ---*/

          Buffer_Send_Var[jPoint] = solver[ADJFLOW_SOL]->GetNodes()->GetSensitivity(iPoint, 0);
          Buffer_Send_Res[jPoint] = solver[ADJFLOW_SOL]->GetNodes()->GetSensitivity(iPoint, 1);
          if (nDim == 3)
            Buffer_Send_Vol[jPoint] = solver[ADJFLOW_SOL]->GetNodes()->GetSensitivity(iPoint, 2);
          jPoint++;
        }
      }

      /*--- Gather the data on the master node. ---*/

#ifdef HAVE_MPI
      SU2_MPI::Gather(Buffer_Send_Var, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Var, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
      SU2_MPI::Gather(Buffer_Send_Res, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Res, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
      if (nDim == 3)
        SU2_MPI::Gather(Buffer_Send_Vol, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Vol, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
#else
      for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Var[iPoint] = Buffer_Send_Var[iPoint];
      for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Res[iPoint] = Buffer_Send_Res[iPoint];
      if (nDim == 3)
        for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Vol[iPoint] = Buffer_Send_Vol[iPoint];
#endif

      /*--- The master node unpacks and sorts this variable by global index ---*/

      if (rank == MASTER_NODE) {
        jPoint = 0; iVar = iVar_SensDim;
        for (iProcessor = 0; iProcessor < size; iProcessor++) {
          for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++) {

            /*--- Get global index, then loop over each variable and store ---*/

            iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
            Data[iVar+0][iGlobal_Index] = Buffer_Recv_Var[jPoint];
            Data[iVar+1][iGlobal_Index] = Buffer_Recv_Res[jPoint];
            if (nDim == 3)
              Data[iVar+2][iGlobal_Index] = Buffer_Recv_Vol[jPoint];
            jPoint++;
          }

          /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/

          jPoint = (iProcessor+1)*nBuffer_Scalar;
        }
      }
    }


    /*--- Communicate the Velocities for dynamic FEM problem ---*/

    if ((Kind_Solver == FEM_ELASTICITY) && (config->GetTime_Domain())) {

      /*--- Loop over this partition to collect the current variable ---*/

      jPoint = 0; su2double *Node_Vel;
      for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {

        /*--- Check for halos & write only if requested ---*/

        if (!Local_Halo[iPoint] || Wrt_Halo) {

          /*--- Load buffers with the three grid velocity components. ---*/

          Node_Vel = solver[FEA_SOL]->GetNodes()->GetSolution_Vel(iPoint);
          Buffer_Send_Var[jPoint] = Node_Vel[0];
          Buffer_Send_Res[jPoint] = Node_Vel[1];
          if (geometry->GetnDim() == 3) Buffer_Send_Vol[jPoint] = Node_Vel[2];
          jPoint++;
        }
      }

      /*--- Gather the data on the master node. ---*/

#ifdef HAVE_MPI
      SU2_MPI::Gather(Buffer_Send_Var, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Var, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
      SU2_MPI::Gather(Buffer_Send_Res, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Res, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
      if (geometry->GetnDim() == 3) {
        SU2_MPI::Gather(Buffer_Send_Vol, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Vol, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
      }
#else
      for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Var[iPoint] = Buffer_Send_Var[iPoint];
      for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Res[iPoint] = Buffer_Send_Res[iPoint];
      if (geometry->GetnDim() == 3) {
        for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Vol[iPoint] = Buffer_Send_Vol[iPoint];
      }
#endif

      /*--- The master node unpacks and sorts this variable by global index ---*/

      if (rank == MASTER_NODE) {
        jPoint = 0; iVar = iVar_FEA_Vel;
        for (iProcessor = 0; iProcessor < size; iProcessor++) {
          for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++) {

            /*--- Get global index, then loop over each variable and store ---*/

            iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
            Data[iVar][iGlobal_Index]   = Buffer_Recv_Var[jPoint];
            Data[iVar+1][iGlobal_Index] = Buffer_Recv_Res[jPoint];
            if (geometry->GetnDim() == 3)
              Data[iVar+2][iGlobal_Index] = Buffer_Recv_Vol[jPoint];
            jPoint++;
          }

          /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/

          jPoint = (iProcessor+1)*nBuffer_Scalar;
        }
      }
    }

    /*--- Communicate the Accelerations for dynamic FEM problem ---*/

    if ((Kind_Solver == FEM_ELASTICITY) && (config->GetTime_Domain())) {

      /*--- Loop over this partition to collect the current variable ---*/

      jPoint = 0; su2double *Node_Accel;
      for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {

        /*--- Check for halos & write only if requested ---*/

        if (!Local_Halo[iPoint] || Wrt_Halo) {

          /*--- Load buffers with the three grid velocity components. ---*/

          Node_Accel = solver[FEA_SOL]->GetNodes()->GetSolution_Accel(iPoint);
          Buffer_Send_Var[jPoint] = Node_Accel[0];
          Buffer_Send_Res[jPoint] = Node_Accel[1];
          if (geometry->GetnDim() == 3) Buffer_Send_Vol[jPoint] = Node_Accel[2];
          jPoint++;
        }
      }

      /*--- Gather the data on the master node. ---*/

#ifdef HAVE_MPI
      SU2_MPI::Gather(Buffer_Send_Var, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Var, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
      SU2_MPI::Gather(Buffer_Send_Res, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Res, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
      if (geometry->GetnDim() == 3) {
        SU2_MPI::Gather(Buffer_Send_Vol, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Vol, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
      }
#else
      for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Var[iPoint] = Buffer_Send_Var[iPoint];
      for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Res[iPoint] = Buffer_Send_Res[iPoint];
      if (geometry->GetnDim() == 3) {
        for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Vol[iPoint] = Buffer_Send_Vol[iPoint];
      }
#endif

      /*--- The master node unpacks and sorts this variable by global index ---*/

      if (rank == MASTER_NODE) {
        jPoint = 0; iVar = iVar_FEA_Accel;
        for (iProcessor = 0; iProcessor < size; iProcessor++) {
          for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++) {

            /*--- Get global index, then loop over each variable and store ---*/

            iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
            Data[iVar][iGlobal_Index]   = Buffer_Recv_Var[jPoint];
            Data[iVar+1][iGlobal_Index] = Buffer_Recv_Res[jPoint];
            if (geometry->GetnDim() == 3)
              Data[iVar+2][iGlobal_Index] = Buffer_Recv_Vol[jPoint];
            jPoint++;
          }

          /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/

          jPoint = (iProcessor+1)*nBuffer_Scalar;
        }
      }
    }

    /*--- Communicate the FEM elasticity stresses (2D) - New elasticity solver---*/

    if (Kind_Solver == FEM_ELASTICITY) {

      /*--- Loop over this partition to collect the current variable ---*/

      jPoint = 0; const su2double *Stress;
      for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {

        /*--- Check for halos & write only if requested ---*/

        if (!Local_Halo[iPoint] || Wrt_Halo) {

          /*--- Load buffers with the three grid velocity components. ---*/

          Stress = solver[FEA_SOL]->GetNodes()->GetStress_FEM(iPoint);
          /*--- Sigma xx ---*/
          Buffer_Send_Var[jPoint] = Stress[0];
          /*--- Sigma yy ---*/
          Buffer_Send_Res[jPoint] = Stress[1];
          /*--- Sigma xy ---*/
          Buffer_Send_Vol[jPoint] = Stress[2];
          jPoint++;
        }
      }

      /*--- Gather the data on the master node. ---*/

#ifdef HAVE_MPI
      SU2_MPI::Gather(Buffer_Send_Var, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Var, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
      SU2_MPI::Gather(Buffer_Send_Res, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Res, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
      SU2_MPI::Gather(Buffer_Send_Vol, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Vol, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
#else
      for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Var[iPoint] = Buffer_Send_Var[iPoint];
      for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Res[iPoint] = Buffer_Send_Res[iPoint];
      for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Vol[iPoint] = Buffer_Send_Vol[iPoint];
#endif

      /*--- The master node unpacks and sorts this variable by global index ---*/

      if (rank == MASTER_NODE) {
        jPoint = 0; iVar = iVar_FEA_Stress;
        for (iProcessor = 0; iProcessor < size; iProcessor++) {
          for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++) {

            /*--- Get global index, then loop over each variable and store ---*/

            iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
            Data[iVar][iGlobal_Index]   = Buffer_Recv_Var[jPoint];
            Data[iVar+1][iGlobal_Index] = Buffer_Recv_Res[jPoint];
            Data[iVar+2][iGlobal_Index] = Buffer_Recv_Vol[jPoint];
            jPoint++;
          }

          /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/

          jPoint = (iProcessor+1)*nBuffer_Scalar;
        }
      }
    }

    /*--- Communicate the FEM elasticity stresses (3D) - New elasticity solver---*/

    if ((Kind_Solver == FEM_ELASTICITY) && (geometry->GetnDim() == 3)) {

      /*--- Loop over this partition to collect the current variable ---*/

      jPoint = 0; const su2double *Stress;
      for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {

        /*--- Check for halos & write only if requested ---*/

        if (!Local_Halo[iPoint] || Wrt_Halo) {

          /*--- Load buffers with the three grid velocity components. ---*/

          Stress = solver[FEA_SOL]->GetNodes()->GetStress_FEM(iPoint);
          /*--- Sigma zz ---*/
          Buffer_Send_Var[jPoint] = Stress[3];
          /*--- Sigma xz ---*/
          Buffer_Send_Res[jPoint] = Stress[4];
          /*--- Sigma yz ---*/
          Buffer_Send_Vol[jPoint] = Stress[5];
          jPoint++;
        }
      }

      /*--- Gather the data on the master node. ---*/

#ifdef HAVE_MPI
      SU2_MPI::Gather(Buffer_Send_Var, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Var, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
      SU2_MPI::Gather(Buffer_Send_Res, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Res, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
      SU2_MPI::Gather(Buffer_Send_Vol, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Vol, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
#else
      for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Var[iPoint] = Buffer_Send_Var[iPoint];
      for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Res[iPoint] = Buffer_Send_Res[iPoint];
      for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Vol[iPoint] = Buffer_Send_Vol[iPoint];

#endif

      /*--- The master node unpacks and sorts this variable by global index ---*/

      if (rank == MASTER_NODE) {
        jPoint = 0; iVar = iVar_FEA_Stress_3D;
        for (iProcessor = 0; iProcessor < size; iProcessor++) {
          for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++) {

            /*--- Get global index, then loop over each variable and store ---*/

            iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
            Data[iVar][iGlobal_Index]   = Buffer_Recv_Var[jPoint];
            Data[iVar+1][iGlobal_Index] = Buffer_Recv_Res[jPoint];
            Data[iVar+2][iGlobal_Index] = Buffer_Recv_Vol[jPoint];
            jPoint++;
          }

          /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/

          jPoint = (iProcessor+1)*nBuffer_Scalar;
        }
      }
    }


    /*--- Communicate the Linear elasticity ---*/

    if ( Kind_Solver == FEM_ELASTICITY ) {

      /*--- Loop over this partition to collect the current variable ---*/
      jPoint = 0;
      for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {

        /*--- Check for halos & write only if requested ---*/

        if (!Local_Halo[iPoint] || Wrt_Halo) {

          /*--- Load buffers with the temperature and laminar viscosity variables. ---*/

          Buffer_Send_Var[jPoint] = solver[FEA_SOL]->GetNodes()->GetVonMises_Stress(iPoint);
          jPoint++;
        }
      }

      /*--- Gather the data on the master node. ---*/

#ifdef HAVE_MPI
      SU2_MPI::Gather(Buffer_Send_Var, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Var, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
#else
      for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Var[iPoint] = Buffer_Send_Var[iPoint];
#endif

      /*--- The master node unpacks and sorts this variable by global index ---*/

      if (rank == MASTER_NODE) {
        jPoint = 0; iVar = iVar_FEA_Extra;
        for (iProcessor = 0; iProcessor < size; iProcessor++) {
          for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++) {

            /*--- Get global index, then loop over each variable and store ---*/

            iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
            Data[iVar][iGlobal_Index] = Buffer_Recv_Var[jPoint];
            jPoint++;
          }

          /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/

          jPoint = (iProcessor+1)*nBuffer_Scalar;
        }
      }
    }

    if (config->GetExtraOutput()) {

      for (jVar = 0; jVar < nVar_Extra; jVar++) {

        /*--- Loop over this partition to collect the current variable ---*/

        jPoint = 0;
        for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {

          /*--- Check for halos & write only if requested ---*/

          if (!Local_Halo[iPoint] || Wrt_Halo) {

            /*--- Get this variable into the temporary send buffer. ---*/

            if (Kind_Solver == RANS) {
              Buffer_Send_Var[jPoint] = solver[TURB_SOL]->OutputVariables[iPoint*nVar_Extra+jVar];
            }
            jPoint++;

          }
        }

        /*--- Gather the data on the master node. ---*/

#ifdef HAVE_MPI
        SU2_MPI::Gather(Buffer_Send_Var, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Var, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
#else
        for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Var[iPoint] = Buffer_Send_Var[iPoint];
#endif

        /*--- The master node unpacks and sorts this variable by global index ---*/

        if (rank == MASTER_NODE) {
          jPoint = 0; iVar = iVar_Extra;
          for (iProcessor = 0; iProcessor < size; iProcessor++) {
            for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++) {

              /*--- Get global index, then loop over each variable and store ---*/

              iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
              Data[iVar+jVar][iGlobal_Index] = Buffer_Recv_Var[jPoint];
              jPoint++;
            }

            /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/

            jPoint = (iProcessor+1)*nBuffer_Scalar;
          }
        }
      }
    }



  /*--- Immediately release the temporary buffers. ---*/

  delete [] Buffer_Send_Var;
  delete [] Buffer_Send_Res;
  delete [] Buffer_Send_Vol;
  delete [] Buffer_Send_GlobalIndex;
  if (rank == MASTER_NODE) {
    delete [] Buffer_Recv_nPoint;
    delete [] Buffer_Recv_Var;
    delete [] Buffer_Recv_Res;
    delete [] Buffer_Recv_Vol;
    delete [] Buffer_Recv_GlobalIndex;
  }

  /*--- Release memory needed for surface coefficients ---*/

  delete [] Local_Halo;

  if (((Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) ||
      ((Kind_Solver == INC_NAVIER_STOKES) || (Kind_Solver == INC_RANS)) ||
    ((Kind_Solver == FEM_NAVIER_STOKES) || (Kind_Solver == FEM_RANS) || (Kind_Solver == FEM_LES)))  {
    delete[] Aux_Frict_x; delete[] Aux_Frict_y; delete[] Aux_Frict_z;
    delete [] Aux_Heat; delete [] Aux_yPlus;
  }
  if (( Kind_Solver == ADJ_EULER              ) ||
      ( Kind_Solver == ADJ_NAVIER_STOKES      ) ||
      ( Kind_Solver == ADJ_RANS               ) ||
      ( Kind_Solver == DISC_ADJ_EULER         ) ||
      ( Kind_Solver == DISC_ADJ_NAVIER_STOKES ) ||
      ( Kind_Solver == DISC_ADJ_RANS          ) ||
      ( Kind_Solver == DISC_ADJ_INC_EULER         ) ||
      ( Kind_Solver == DISC_ADJ_INC_NAVIER_STOKES ) ||
      ( Kind_Solver == DISC_ADJ_INC_RANS          )) {
    delete [] Aux_Sens;
  }

}

void COutputLegacy::MergeBaselineSolution(CConfig *config, CGeometry *geometry, CSolver *solver, unsigned short val_iZone) {

  /*--- Local variables needed on all processors ---*/
  unsigned short iVar;
  unsigned long iPoint = 0, jPoint = 0;

  nVar_Total = config->fields.size() - 1;

  /*--- Merge the solution either in serial or parallel. ---*/

#ifndef HAVE_MPI

  /*--- In serial, the single process has access to all solution data,
   so it is simple to retrieve and store inside Solution_Data. ---*/

  unsigned short iMarker;
  unsigned long iVertex, nTotalPoints = 0;
  int SendRecv;

  /*--- First, create a structure to locate any periodic halo nodes ---*/
  int *Local_Halo = new int[geometry->GetnPoint()];
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
    Local_Halo[iPoint] = !geometry->node[iPoint]->GetDomain();

  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) {
      SendRecv = config->GetMarker_All_SendRecv(iMarker);
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        if ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() > 0) &&
            (geometry->vertex[iMarker][iVertex]->GetRotation_Type() % 2 == 1) &&
            (SendRecv < 0)) {
          Local_Halo[iPoint] = false;
        }
      }

    }
  }

  /*--- Total number of points in the mesh (this might include periodic points). ---*/

  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
    if (!Local_Halo[iPoint]) nTotalPoints++;

  nGlobal_Poin = nTotalPoints;
  Data = new su2double*[nVar_Total];
  for (iVar = 0; iVar < nVar_Total; iVar++) {
    Data[iVar] = new su2double[nGlobal_Poin];
  }

  /*--- Loop over all points in the mesh, but only write data
   for nodes in the domain (ignore periodic halo nodes). ---*/

  jPoint = 0;
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
    if (!Local_Halo[iPoint]) {

      /*--- Solution (first, and second system of equations) ---*/

      unsigned short jVar = 0;
      for (iVar = 0; iVar < nVar_Total; iVar++) {
        Data[jVar][jPoint] = solver->GetNodes()->GetSolution(iPoint,iVar);
        jVar++;
      }
    }

    /*--- Increment jPoint as the counter. We need this because iPoint
     may include halo nodes that we skip over during this loop. ---*/

    jPoint++;

  }

#else

  /*--- MPI preprocessing ---*/

  int nProcessor = size, iProcessor;

  /*--- Local variables needed for merging with MPI ---*/

  unsigned long iVertex, iMarker;
  unsigned long Buffer_Send_nPoint[1], *Buffer_Recv_nPoint = NULL;
  unsigned long nLocalPoint = 0, MaxLocalPoint = 0;
  unsigned long iGlobal_Index = 0, nBuffer_Scalar = 0;

  int *Local_Halo = new int[geometry->GetnPoint()];
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
    Local_Halo[iPoint] = !geometry->node[iPoint]->GetDomain();

  bool Wrt_Halo = config->GetWrt_Halo(), isPeriodic;

  /*--- Search all send/recv boundaries on this partition for any periodic
   nodes that were part of the original domain. We want to recover these
   for visualization purposes. ---*/

  if (Wrt_Halo) {
    nLocalPoint = geometry->GetnPoint();
  } else {
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) {

        /*--- Checking for less than or equal to the rank, because there may
         be some periodic halo nodes that send info to the same rank. ---*/

        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          isPeriodic = ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() > 0) &&
                        (geometry->vertex[iMarker][iVertex]->GetRotation_Type() % 2 == 1));
          if (isPeriodic) Local_Halo[iPoint] = false;
        }
      }
    }

    /*--- Sum total number of nodes that belong to the domain ---*/

    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
      if (Local_Halo[iPoint] == false)
        nLocalPoint++;

  }
  Buffer_Send_nPoint[0] = nLocalPoint;

  if (rank == MASTER_NODE) Buffer_Recv_nPoint = new unsigned long[nProcessor];

  SU2_MPI::Allreduce(&nLocalPoint, &MaxLocalPoint, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
  SU2_MPI::Gather(&Buffer_Send_nPoint, 1, MPI_UNSIGNED_LONG, Buffer_Recv_nPoint, 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);

  nBuffer_Scalar = MaxLocalPoint;

  /*--- Send and Recv buffers. ---*/

  su2double *Buffer_Send_Var = new su2double[MaxLocalPoint];
  su2double *Buffer_Recv_Var = NULL;

  unsigned long *Buffer_Send_GlobalIndex = new unsigned long[MaxLocalPoint];
  unsigned long *Buffer_Recv_GlobalIndex = NULL;

  /*--- Prepare the receive buffers in the master node only. ---*/
  if (rank == MASTER_NODE) {

    Buffer_Recv_Var = new su2double[nProcessor*MaxLocalPoint];
    Buffer_Recv_GlobalIndex = new unsigned long[nProcessor*MaxLocalPoint];

    /*--- Sum total number of nodes to be written and allocate arrays ---*/
    nGlobal_Poin = 0;
    for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
      nGlobal_Poin += Buffer_Recv_nPoint[iProcessor];
    }
    Data = new su2double*[nVar_Total];
    for (iVar = 0; iVar < nVar_Total; iVar++) {
      Data[iVar] = new su2double[nGlobal_Poin];
    }

  }

  /*--- Main communication routine. Loop over each variable that has
   been requested by the user and perform the MPI comm. Temporary
   1-D buffers are used to send the solution for each variable at all
   nodes on each partition to the master node. These are then unpacked
   by the master and sorted by global index in one large n-dim. array. ---*/

  for (iVar = 0; iVar < nVar_Total; iVar++) {

    /*--- Loop over this partition to collect the current variable ---*/
    jPoint = 0;
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {

      /*--- Check for halos and write only if requested ---*/
      if (!Local_Halo[iPoint] || Wrt_Halo) {

        /*--- Get this variable into the temporary send buffer. ---*/
        Buffer_Send_Var[jPoint] = solver->GetNodes()->GetSolution(iPoint,iVar);

        /*--- Only send/recv the volumes & global indices during the first loop ---*/
        if (iVar == 0) {
          Buffer_Send_GlobalIndex[jPoint] = geometry->node[iPoint]->GetGlobalIndex();
        }
        jPoint++;
      }
    }

    /*--- Gather the data on the master node. ---*/

    SU2_MPI::Gather(Buffer_Send_Var, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Var, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    if (iVar == 0) {
      SU2_MPI::Gather(Buffer_Send_GlobalIndex, nBuffer_Scalar, MPI_UNSIGNED_LONG, Buffer_Recv_GlobalIndex, nBuffer_Scalar, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
    }

    /*--- The master node unpacks and sorts this variable by global index ---*/
    if (rank == MASTER_NODE) {
      jPoint = 0;
      for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
        for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++) {

          /*--- Get global index, then loop over each variable and store ---*/
          iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
          Data[iVar][iGlobal_Index] = Buffer_Recv_Var[jPoint];
          jPoint++;
        }
        /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
        jPoint = (iProcessor+1)*nBuffer_Scalar;
      }
    }
  }

  /*--- Immediately release the temporary buffers. ---*/

  delete [] Buffer_Send_Var;
  delete [] Buffer_Send_GlobalIndex;
  if (rank == MASTER_NODE) {
    delete [] Buffer_Recv_Var;
    delete [] Buffer_Recv_GlobalIndex;
  }

#endif

  delete [] Local_Halo;

}

void COutputLegacy::SetRestart(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned short val_iZone) {

  /*--- Local variables ---*/

  unsigned short nZone = geometry->GetnZone();
  unsigned short Kind_Solver  = config->GetKind_Solver();
  unsigned short iVar, iDim, nDim = geometry->GetnDim();
  unsigned long iPoint, iExtIter = config->GetInnerIter();
  bool dynamic_grid = config->GetDynamic_Grid();
  bool dynamic_fem = (config->GetTime_Domain());
  bool fem = (config->GetKind_Solver() == FEM_ELASTICITY);
  bool disc_adj_fem = (config->GetKind_Solver() == DISC_ADJ_FEM);
  ofstream restart_file;
  ofstream meta_file;
  string filename, meta_filename;
  bool adjoint = config->GetContinuous_Adjoint() || config->GetDiscrete_Adjoint();

  /*--- Retrieve filename from config ---*/

//  if (((config->GetContinuous_Adjoint()) || (config->GetDiscrete_Adjoint())) && ((config->GetKind_Solver() != DISC_ADJ_FEM)))  {
//    filename = config->GetRestart_AdjFileName();
//    filename = config->GetObjFunc_Extension(filename);
//  } else if (fem) {
//    filename = config->GetRestart_FEMFileName();
//  } else if (disc_adj_fem){
//    filename = config->GetRestart_AdjFEMFileName();
//  } else {
//    filename = config->GetRestart_FileName();
//  }

  /*--- Append the zone number if multizone problems ---*/
  if (nZone > 1)
    filename= config->GetMultizone_FileName(filename, val_iZone, ".dat");

  /*--- Unsteady problems require an iteration number to be appended. ---*/
  if (config->GetTime_Marching() == HARMONIC_BALANCE) {
    filename = config->GetUnsteady_FileName(filename, SU2_TYPE::Int(config->GetiInst()), ".dat");
  } else if (config->GetTime_Domain()) {
    filename = config->GetUnsteady_FileName(filename, SU2_TYPE::Int(iExtIter), ".dat");
  } else if ((fem || disc_adj_fem) && (config->GetWrt_Dynamic())) {
    filename = config->GetUnsteady_FileName(filename, SU2_TYPE::Int(iExtIter), ".dat");
  }

  /*--- Open the restart file and write the solution. ---*/

  restart_file.open(filename.c_str(), ios::out);
  restart_file.precision(15);

  /*--- Write the header line based on the particular solver ----*/

  restart_file << "\"PointID\"";

  /*--- Mesh coordinates are always written to the restart first ---*/

  if (nDim == 2) {
    restart_file << "\t\"x\"\t\"y\"";
  } else {
    restart_file << "\t\"x\"\t\"y\"\t\"z\"";
  }

  for (iVar = 0; iVar < nVar_Consv; iVar++) {
  if (( Kind_Solver == FEM_ELASTICITY ) || ( Kind_Solver == DISC_ADJ_FEM))
      restart_file << "\t\"Displacement_" << iVar+1<<"\"";
    else
      restart_file << "\t\"Conservative_" << iVar+1<<"\"";
  }

  if (config->GetWrt_Limiters()) {
    for (iVar = 0; iVar < nVar_Consv; iVar++) {
      restart_file << "\t\"Limiter_" << iVar+1<<"\"";
    }
  }
  if (config->GetWrt_Residuals()) {
    for (iVar = 0; iVar < nVar_Consv; iVar++) {
      restart_file << "\t\"Residual_" << iVar+1<<"\"";
    }
  }

  /*--- Mesh velocities for dynamic mesh cases ---*/

  if (dynamic_grid && !fem) {
    if (nDim == 2) {
      restart_file << "\t\"Grid_Velx\"\t\"Grid_Vely\"";
    } else {
      restart_file << "\t\"Grid_Velx\"\t\"Grid_Vely\"\t\"Grid_Velz\"";
    }
  }

  if ((Kind_Solver == EULER) || (Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS) ||
      (Kind_Solver == INC_EULER) || (Kind_Solver == INC_NAVIER_STOKES) || (Kind_Solver == INC_RANS)) {
    if (config->GetTabular_FileFormat() == TAB_CSV) {
      restart_file << "\t\"Pressure\"\t\"Temperature\"\t\"Pressure_Coefficient\"\t\"Mach\"";
    } else
      restart_file << "\t\"Pressure\"\t\"Temperature\"\t\"C<sub>p</sub>\"\t\"Mach\"";
  }

  if ((Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS) ||
      (Kind_Solver == INC_NAVIER_STOKES) || (Kind_Solver == INC_RANS)) {
    if (config->GetTabular_FileFormat() == TAB_CSV) {
      if (nDim == 2) restart_file << "\t\"Laminar_Viscosity\"\t\"Skin_Friction_Coefficient_X\"\t\"Skin_Friction_Coefficient_Y\"\t\"Heat_Flux\"\t\"Y_Plus\"";
      if (nDim == 3) restart_file << "\t\"Laminar_Viscosity\"\t\"Skin_Friction_Coefficient_X\"\t\"Skin_Friction_Coefficient_Y\"\t\"Skin_Friction_Coefficient_Z\"\t\"Heat_Flux\"\t\"Y_Plus\"";
    } else {
      if (nDim == 2) restart_file << "\t\"<greek>m</greek>\"\t\"C<sub>f</sub>_x\"\t\"C<sub>f</sub>_y\"\t\"h\"\t\"y<sup>+</sup>\"";
      if (nDim == 3) restart_file << "\t\"<greek>m</greek>\"\t\"C<sub>f</sub>_x\"\t\"C<sub>f</sub>_y\"\t\"C<sub>f</sub>_z\"\t\"h\"\t\"y<sup>+</sup>\"";
    }
  }

  if (Kind_Solver == RANS || Kind_Solver == INC_RANS) {
    if ((config->GetTabular_FileFormat() == TAB_CSV)) {
      restart_file << "\t\"Eddy_Viscosity\"";
    } else
      restart_file << "\t\"<greek>m</greek><sub>t</sub>\"";
  }

  if (config->GetWrt_SharpEdges()) {
    if (((Kind_Solver == EULER) || (Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) ||
        ((Kind_Solver == INC_EULER) || (Kind_Solver == INC_NAVIER_STOKES) || (Kind_Solver == INC_RANS)))  {
      restart_file << "\t\"Sharp_Edge_Dist\"";
    }
  }

  if ((Kind_Solver == ADJ_EULER              ) ||
      (Kind_Solver == ADJ_NAVIER_STOKES      ) ||
      (Kind_Solver == ADJ_RANS               )   ) {
    restart_file << "\t\"Surface_Sensitivity\"\t\"Solution_Sensor\"";
  }
  if (( Kind_Solver == DISC_ADJ_EULER              ) ||
      ( Kind_Solver == DISC_ADJ_NAVIER_STOKES      ) ||
      ( Kind_Solver == DISC_ADJ_RANS               ) ||
      ( Kind_Solver == DISC_ADJ_INC_EULER              ) ||
      ( Kind_Solver == DISC_ADJ_INC_NAVIER_STOKES      ) ||
      ( Kind_Solver == DISC_ADJ_INC_RANS               )) {
    restart_file << "\t\"Surface_Sensitivity\"\t\"Sensitivity_x\"\t\"Sensitivity_y\"";
    if (geometry->GetnDim() == 3) {
      restart_file << "\t\"Sensitivity_z\"";
    }
  }

  if (Kind_Solver == FEM_ELASTICITY) {
    if (!dynamic_fem) {
      if (geometry->GetnDim() == 2)
        restart_file << "\t\"Sxx\"\t\"Syy\"\t\"Sxy\"\t\"Von_Mises_Stress\"";
      if (geometry->GetnDim() == 3)
        restart_file << "\t\"Sxx\"\t\"Syy\"\t\"Sxy\"\t\"Szz\"\t\"Sxz\"\t\"Syz\"\t\"Von_Mises_Stress\"";
    }
    else if (dynamic_fem) {
      if (geometry->GetnDim() == 2) {
        restart_file << "\t\"Velocity_1\"\t\"Velocity_2\"\t\"Acceleration_1\"\t\"Acceleration_2\"";
        restart_file << "\t\"Sxx\"\t\"Syy\"\t\"Sxy\"\t\"Von_Mises_Stress\"";
      }
      if (geometry->GetnDim() == 3) {
        restart_file << "\t\"Velocity_1\"\t\"Velocity_2\"\t\"Velocity_3\"\t\"Acceleration_1\"\t\"Acceleration_2\"\t\"Acceleration_3\"";
        restart_file << "\t\"Sxx\"\t\"Syy\"\t\"Sxy\"\t\"Szz\"\t\"Sxz\"\t\"Syz\"\t\"Von_Mises_Stress\"";
      }
    }
  }

  if ((Kind_Solver == DISC_ADJ_FEM) && (config->GetFSI_Simulation())){
    if (geometry->GetnDim() == 2)
      restart_file << "\t\"CrossTerm_1\"\t\"CrossTerm_2\"";
    if (geometry->GetnDim() == 3)
      restart_file << "\t\"CrossTerm_1\"\t\"CrossTerm_2\"\t\"CrossTerm_3\"";
  }


  if (config->GetExtraOutput()) {
    string *headings = NULL;
    //if (Kind_Solver == RANS) {
    headings = solver[TURB_SOL]->OutputHeadingNames;
    //}

    for (iVar = 0; iVar < nVar_Extra; iVar++) {
      if (headings == NULL) {
        restart_file << "\t\"ExtraOutput_" << iVar+1<<"\"";
      } else {
        restart_file << "\t\""<< headings[iVar] <<"\"";
      }
    }
  }


  restart_file << "\n";

  /*--- Write the restart file ---*/

  /*--- Determine whether or not the FEM solver is used, which decides the
   type of geometry classes that are instantiated. ---*/
  const bool fem_solver = config->GetFEMSolver();

  unsigned long nPointTotal = 0;
  if ( fem_solver ) {
    nPointTotal = solver[FLOW_SOL]->GetnDOFsGlobal();
  } else {
    nPointTotal = geometry->GetGlobal_nPointDomain();
  }

  for (iPoint = 0; iPoint < nPointTotal; iPoint++) {

    /*--- Index of the point ---*/
    restart_file << iPoint << "\t";

    /*--- Write the grid coordinates first ---*/
    for (iDim = 0; iDim < nDim; iDim++) {
      restart_file << scientific << Coords[iDim][iPoint] << "\t";
    }

    /*--- Loop over the variables and write the values to file ---*/
    for (iVar = 0; iVar < nVar_Total; iVar++) {
      restart_file << scientific << Data[iVar][iPoint] << "\t";
    }
    restart_file << "\n";
  }

  /*--- Write the general header and flow conditions ----*/

//  if (dual_time)
    restart_file <<"EXT_ITER= " << config->GetInnerIter() + 1 << endl;
//  else
//    restart_file <<"EXT_ITER= " << config->GetInnerIter() + config->GetInnerIter_OffSet() + 1 << endl;
  restart_file <<"AOA= " << config->GetAoA() - config->GetAoA_Offset() << endl;
  restart_file <<"SIDESLIP_ANGLE= " << config->GetAoS() - config->GetAoS_Offset() << endl;
  restart_file <<"INITIAL_BCTHRUST= " << config->GetInitial_BCThrust() << endl;
  restart_file <<"DCD_DCL_VALUE= " << config->GetdCD_dCL() << endl;
  restart_file <<"DCMX_DCL_VALUE= " << config->GetdCMx_dCL() << endl;
  restart_file <<"DCMY_DCL_VALUE= " << config->GetdCMy_dCL() << endl;
  restart_file <<"DCMZ_DCL_VALUE= " << config->GetdCMz_dCL() << endl;
  if (adjoint) restart_file << "SENS_AOA=" << solver[ADJFLOW_SOL]->GetTotal_Sens_AoA() * PI_NUMBER / 180.0 << endl;

  /*--- Close the data portion of the restart file. ---*/

  restart_file.close();

}

void COutputLegacy::DeallocateCoordinates(CConfig *config, CGeometry *geometry) {

  unsigned short iDim, nDim = geometry->GetnDim();

  /*--- The master node alone owns all data found in this routine. ---*/

  if (rank == MASTER_NODE) {

    /*--- Deallocate memory for coordinate data ---*/

    for (iDim = 0; iDim < nDim; iDim++) {
      delete [] Coords[iDim];
    }
    delete [] Coords;
  }
}

void COutputLegacy::DeallocateConnectivity(CConfig *config, CGeometry *geometry, bool surf_sol) {

  /*--- The master node alone owns all data found in this routine. ---*/
  if (rank == MASTER_NODE) {

    /*--- Deallocate memory for connectivity data ---*/
    if (surf_sol) {
      if (nGlobal_Line > 0      && Conn_Line      != NULL) delete [] Conn_Line;
      if (nGlobal_BoundTria > 0 && Conn_BoundTria != NULL) delete [] Conn_BoundTria;
      if (nGlobal_BoundQuad > 0 && Conn_BoundQuad != NULL) delete [] Conn_BoundQuad;
    }
    else {
      if (nGlobal_Tria > 0 && Conn_Tria != NULL) delete [] Conn_Tria;
      if (nGlobal_Quad > 0 && Conn_Quad != NULL) delete [] Conn_Quad;
      if (nGlobal_Tetr > 0 && Conn_Tetr != NULL) delete [] Conn_Tetr;
      if (nGlobal_Hexa > 0 && Conn_Hexa != NULL) delete [] Conn_Hexa;
      if (nGlobal_Pris > 0 && Conn_Pris != NULL) delete [] Conn_Pris;
      if (nGlobal_Pyra > 0 && Conn_Pyra != NULL) delete [] Conn_Pyra;

    }

  }
}

void COutputLegacy::DeallocateSolution(CConfig *config, CGeometry *geometry) {

  /*--- The master node alone owns all data found in this routine. ---*/
  if (rank == MASTER_NODE) {

    /*--- Deallocate memory for solution data ---*/
    for (unsigned short iVar = 0; iVar < nVar_Total; iVar++) {
      delete [] Data[iVar];
    }
    delete [] Data;

  }
}

void COutputLegacy::SetConvHistory_Header(ofstream *ConvHist_file, CConfig *config, unsigned short val_iZone, unsigned short val_iInst) {
  char cstr[200], turb_resid[1000], adj_turb_resid[1000];
  unsigned short iMarker_Monitoring;
  string Monitoring_Tag, monitoring_coeff, aeroelastic_coeff, turbo_coeff;

  bool rotating_frame = config->GetRotating_Frame();
  bool aeroelastic = config->GetAeroelastic_Simulation();
  bool equiv_area = config->GetEquivArea();
  bool buffet = (config->GetBuffet_Monitoring() || config->GetKind_ObjFunc() == BUFFET_SENSOR);
  bool engine        = ((config->GetnMarker_EngineInflow() != 0) || (config->GetnMarker_EngineExhaust() != 0));
  bool actuator_disk = ((config->GetnMarker_ActDiskInlet() != 0) || (config->GetnMarker_ActDiskOutlet() != 0));
  bool turbulent = ((config->GetKind_Solver() == RANS) || (config->GetKind_Solver() == ADJ_RANS) ||
                    (config->GetKind_Solver() == DISC_ADJ_RANS) || (config->GetKind_Solver() == DISC_ADJ_INC_RANS) ||
                    (config->GetKind_Solver() == INC_RANS));
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
  bool radiation            = config->AddRadiation();

  if (config->GetKind_Solver() == RANS || config->GetKind_Solver()  == NAVIER_STOKES ||
      config->GetKind_Solver() == INC_RANS || config->GetKind_Solver()  == INC_NAVIER_STOKES) {
    thermal = true;
  }

  /*--- Write file name with extension ---*/
  string filename = config->GetConv_FileName();
  string hist_ext = ".csv";
  if (config->GetTabular_FileFormat() == TAB_TECPLOT) hist_ext = ".dat";

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
    case SA:case SA_NEG:case SA_E: case SA_COMP: case SA_E_COMP:
      SPRINTF (turb_resid, ",\"Res_Turb[0]\"");
      break;
    case SST:case SST_SUST:
      SPRINTF (turb_resid, ",\"Res_Turb[0]\",\"Res_Turb[1]\"");
      break;
  }
  switch (config->GetKind_Turb_Model()) {
    case SA:case SA_NEG:case SA_E: case SA_COMP: case SA_E_COMP:
      SPRINTF (adj_turb_resid, ",\"Res_AdjTurb[0]\"");
      break;
    case SST:case SST_SUST:
      SPRINTF (adj_turb_resid, ",\"Res_AdjTurb[0]\",\"Res_AdjTurb[1]\"");
      break;
  }
  char fem_resid[]= ",\"Res_FEM[0]\",\"Res_FEM[1]\",\"Res_FEM[2]\"";
  char heat_resid[]= ",\"Res_Heat\"";
  char rad_resid[]= ",\"Res_P1-rad\"";

  /*--- End of the header ---*/

  char end[]= ",\"Linear_Solver_Iterations\",\"CFL_Number\",\"Time(min)\"\n";
  char endfea[]= ",\"Linear_Solver_Iterations\",\"Time(min)\"\n";

  if ((config->GetTabular_FileFormat() == TECPLOT) ||
      (config->GetTabular_FileFormat() == TECPLOT_BINARY)) {
    ConvHist_file[0] << "TITLE = \"SU2 Simulation\"" << endl;
    ConvHist_file[0] << "VARIABLES = ";
  }

  /*--- Write the header, case depending ---*/

  switch (config->GetKind_Solver()) {

    case EULER : case NAVIER_STOKES: case RANS :
    case INC_EULER : case INC_NAVIER_STOKES: case INC_RANS :
    case FEM_EULER : case FEM_NAVIER_STOKES: case FEM_RANS : case FEM_LES:
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

    case ADJ_EULER      : case ADJ_NAVIER_STOKES      : case ADJ_RANS:
    case DISC_ADJ_EULER: case DISC_ADJ_NAVIER_STOKES: case DISC_ADJ_RANS:
    case DISC_ADJ_INC_EULER: case DISC_ADJ_INC_NAVIER_STOKES: case DISC_ADJ_INC_RANS:
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

    case HEAT_EQUATION:
      ConvHist_file[0] << begin << heat_coeff;
      ConvHist_file[0] << heat_resid << end;
      break;

    case FEM_ELASTICITY:
      ConvHist_file[0] << begin << fem_coeff;
      if (incload) ConvHist_file[0] << fem_incload;
      ConvHist_file[0] << fem_resid << endfea;
      break;

    case DISC_ADJ_FEM:
      ConvHist_file[0] << begin << fem_coeff;
      ConvHist_file[0] << fem_resid << endfea;
      break;

  }

  if (config->GetTabular_FileFormat() == TECPLOT ||
      config->GetTabular_FileFormat() == TECPLOT_BINARY) {
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
  bool fea                  = ((config[val_iZone]->GetKind_Solver()== FEM_ELASTICITY)||(config[val_iZone]->GetKind_Solver()== DISC_ADJ_FEM));
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

  bool compressible = (config[val_iZone]->GetKind_Regime() == COMPRESSIBLE);
  bool incompressible = (config[val_iZone]->GetKind_Regime() == INCOMPRESSIBLE);

  if (!disc_adj && !cont_adj && !DualTime_Iteration) {

    if (fixed_cl &&
        (solver_container[val_iZone][val_iInst][MESH_0][FLOW_SOL]->GetStart_AoA_FD()) &&
        (iExtIter != solver_container[val_iZone][val_iInst][MESH_0][FLOW_SOL]->GetIter_Update_AoA())) {
      output_files = false;
    }

    if (fea || fluid_structure) output_files = false;

    /*--- We need to evaluate some of the objective functions to write the value on the history file ---*/

    if (((iExtIter % (config[val_iZone]->GetWrt_Sol_Freq())) == 0) ||
        (!fixed_cl && (iExtIter == (config[val_iZone]->GetnInner_Iter()-1))) ||
        /*--- If CL mode we need to compute the complete solution at two very particular iterations ---*/
        (fixed_cl && (iExtIter == (config[val_iZone]->GetnInner_Iter()-2) ||
          (solver_container[val_iZone][val_iInst][MESH_0][FLOW_SOL]->GetStart_AoA_FD() &&
          iExtIter == solver_container[val_iZone][val_iInst][MESH_0][FLOW_SOL]->GetIter_Update_AoA())))) {


      if ((rank == MASTER_NODE) && output_files) cout << endl << "------------------------ Evaluate Special Output ------------------------";

      switch (config[val_iZone]->GetKind_Solver()) {

        case EULER: case NAVIER_STOKES: case RANS:
        case INC_EULER: case INC_NAVIER_STOKES: case INC_RANS:
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

          /*--- For specific applications, evaluate and plot the equivalent area. ---*/

          if (config[val_iZone]->GetnMarker_NearFieldBound() != 0) {
            SpecialOutput_SonicBoom(solver_container[val_iZone][val_iInst][MESH_0][FLOW_SOL],
                                    geometry[val_iZone][val_iInst][MESH_0], config[val_iZone], output_files);
          }

          /*--- For specific applications, evaluate and plot the cp coefficent at different stations. ---*/

          if (config[val_iZone]->GetPlot_Section_Forces()) {
            SpecialOutput_SpanLoad(solver_container[val_iZone][val_iInst][MESH_0][FLOW_SOL],
                                   geometry[val_iZone][val_iInst][MESH_0], config[val_iZone], output_files);
          }

          break;
      }

      /*--- Output a file with the forces breakdown. ---*/

      if (config[val_iZone]->GetTime_Marching() == HARMONIC_BALANCE) {
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
      solver_container[val_iZone][val_iInst][FinestMesh][FLOW_SOL]->SetTotal_ComboObj(0.0);
      switch (config[val_iZone]->GetKind_Solver()) {
      case EULER:                   case NAVIER_STOKES:                   case RANS:
      case INC_EULER:               case INC_NAVIER_STOKES:               case INC_RANS:
        solver_container[val_iZone][val_iInst][FinestMesh][FLOW_SOL]->Evaluate_ObjFunc(config[val_iZone]);
        break;
      }
    }

    unsigned long ExtIter_OffSet = config[val_iZone]->GetExtIter_OffSet();
    if (config[val_iZone]->GetTime_Marching() == DT_STEPPING_1ST ||
        config[val_iZone]->GetTime_Marching() == DT_STEPPING_2ND)
      ExtIter_OffSet = 0;

    /*--- WARNING: These buffers have hard-coded lengths. Note that you
     may have to adjust them to be larger if adding more entries. ---*/

    char begin[1000], direct_coeff[1000], heat_coeff[1000], equivalent_area_coeff[1000], engine_coeff[1000], rotating_frame_coeff[1000], Cp_inverse_design[1000], Heat_inverse_design[1000], surface_coeff[1000], aeroelastic_coeff[1000], monitoring_coeff[10000], buffet_coeff[1000],
    adjoint_coeff[1000], flow_resid[1000], adj_flow_resid[1000], turb_resid[1000], trans_resid[1000],
    adj_turb_resid[1000], rad_resid[1000],
    begin_fem[1000], fem_coeff[1000], heat_resid[1000], combo_obj[1000],
    fem_resid[1000], end[1000], end_fem[1000], surface_outputs[1000], d_surface_outputs[1000], d_direct_coeff[1000], turbo_coeff[10000];


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
    bool thermal = (config[val_iZone]->GetKind_Solver() == RANS || config[val_iZone]->GetKind_Solver()  == NAVIER_STOKES ||
                    config[val_iZone]->GetKind_Solver() == INC_NAVIER_STOKES || config[val_iZone]->GetKind_Solver() == INC_RANS );
    bool turbulent = ((config[val_iZone]->GetKind_Solver() == RANS) || (config[val_iZone]->GetKind_Solver() == ADJ_RANS) ||
                      (config[val_iZone]->GetKind_Solver() == DISC_ADJ_RANS) || config[val_iZone]->GetKind_Solver()  == INC_RANS ||
                      config[val_iZone]->GetKind_Solver()  == DISC_ADJ_INC_RANS);
    bool adjoint =  cont_adj || disc_adj;
    bool frozen_visc = (cont_adj && config[val_iZone]->GetFrozen_Visc_Cont()) ||( disc_adj && config[val_iZone]->GetFrozen_Visc_Disc());
    bool heat =  ((config[val_iZone]->GetKind_Solver() == HEAT_EQUATION) || (config[val_iZone]->GetWeakly_Coupled_Heat()));
    bool weakly_coupled_heat = config[val_iZone]->GetWeakly_Coupled_Heat();
    bool flow = (config[val_iZone]->GetKind_Solver() == EULER) || (config[val_iZone]->GetKind_Solver() == NAVIER_STOKES) ||
                (config[val_iZone]->GetKind_Solver() == RANS) || (config[val_iZone]->GetKind_Solver() == FEM_EULER) ||
                (config[val_iZone]->GetKind_Solver() == FEM_NAVIER_STOKES) || (config[val_iZone]->GetKind_Solver() == FEM_RANS) ||
                (config[val_iZone]->GetKind_Solver() == FEM_LES) || (config[val_iZone]->GetKind_Solver() == ADJ_EULER) ||
                (config[val_iZone]->GetKind_Solver() == ADJ_NAVIER_STOKES) || (config[val_iZone]->GetKind_Solver() == ADJ_RANS) ||
                (config[val_iZone]->GetKind_Solver() == INC_EULER) || (config[val_iZone]->GetKind_Solver() == INC_NAVIER_STOKES) ||
                (config[val_iZone]->GetKind_Solver() == INC_RANS);
    bool buffet = (config[val_iZone]->GetBuffet_Monitoring() || config[val_iZone]->GetKind_ObjFunc() == BUFFET_SENSOR);

    bool fem = ((config[val_iZone]->GetKind_Solver() == FEM_ELASTICITY) ||          // FEM structural solver.
                (config[val_iZone]->GetKind_Solver() == DISC_ADJ_FEM));
    bool linear_analysis = (config[val_iZone]->GetGeometricConditions() == SMALL_DEFORMATIONS);  // Linear analysis.
    bool nonlinear_analysis = (config[val_iZone]->GetGeometricConditions() == LARGE_DEFORMATIONS);  // Nonlinear analysis.
    bool fsi = (config[val_iZone]->GetFSI_Simulation());          // FEM structural solver.
    bool discadj_fem = (config[val_iZone]->GetKind_Solver() == DISC_ADJ_FEM);

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
      else if ((config[ExtraHeatOutputZone]->GetKind_Solver() != HEAT_EQUATION)) {
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
    su2double *residual_flow         = NULL,
    *residual_turbulent    = NULL,
    *residual_transition   = NULL;
    su2double *residual_adjflow      = NULL,
    *residual_adjturbulent = NULL,
    *residual_adjheat = NULL;
    su2double *residual_fea          = NULL;
    su2double *residual_fem          = NULL;
    su2double *residual_heat         = NULL;
    su2double *residual_rad          = NULL;

    /*--- Coefficients Monitored arrays ---*/
    su2double *aeroelastic_plunge = NULL,
    *aeroelastic_pitch     = NULL,
    *Surface_CL            = NULL,
    *Surface_CD            = NULL,
    *Surface_CSF           = NULL,
    *Surface_CEff          = NULL,
    *Surface_CFx           = NULL,
    *Surface_CFy           = NULL,
    *Surface_CFz           = NULL,
    *Surface_CMx           = NULL,
    *Surface_CMy           = NULL,
    *Surface_CMz           = NULL,
    *Surface_Buffet_Metric = NULL;

    /*--- Initialize number of variables ---*/
    unsigned short nVar_Flow = 0, nVar_Turb = 0,
    nVar_Trans = 0, nVar_Heat = 0,
    nVar_AdjFlow = 0, nVar_AdjTurb = 0, nVar_AdjHeat = 0,
    nVar_FEM = 0, nVar_Rad = 0;

    /*--- Direct problem variables ---*/
    if (compressible) nVar_Flow = nDim+2; else nVar_Flow = nDim+2;
    if (turbulent) {
      switch (config[val_iZone]->GetKind_Turb_Model()) {
        case SA: case SA_NEG: case SA_E: case SA_E_COMP: case SA_COMP: nVar_Turb = 1; break;
        case SST: case SST_SUST: nVar_Turb = 2; break;
      }
    }
    if (transition) nVar_Trans = 2;
    if (heat) nVar_Heat = 1;

    if (fem) {
      if (linear_analysis) nVar_FEM = nDim;
      if (nonlinear_analysis) nVar_FEM = 3;

      if (config[val_iZone]->GetKind_Solver() == DISC_ADJ_FEM) nVar_FEM = nDim;

    }

    if (radiation) nVar_Rad = 1;

    /*--- Adjoint problem variables ---*/
    if (compressible) nVar_AdjFlow = nDim+2; else nVar_AdjFlow = nDim+2;
    if (turbulent) {
      switch (config[val_iZone]->GetKind_Turb_Model()) {
        case SA: case SA_NEG: case SA_E: case SA_E_COMP: case SA_COMP: nVar_AdjTurb = 1; break;
        case SST: case SST_SUST: nVar_AdjTurb = 2; break;
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

      case EULER:                   case NAVIER_STOKES:                   case RANS:
      case INC_EULER:               case INC_NAVIER_STOKES:               case INC_RANS:
      case FEM_EULER:               case FEM_NAVIER_STOKES:               case FEM_RANS:      case FEM_LES:
      case ADJ_EULER:               case ADJ_NAVIER_STOKES:               case ADJ_RANS:
      case DISC_ADJ_EULER:          case DISC_ADJ_NAVIER_STOKES:          case DISC_ADJ_RANS:
      case DISC_ADJ_INC_EULER:      case DISC_ADJ_INC_NAVIER_STOKES:      case DISC_ADJ_INC_RANS:

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


      case HEAT_EQUATION:

        /*--- Heat coefficients  ---*/

        Total_Heat         = solver_container[val_iZone][val_iInst][FinestMesh][HEAT_SOL]->GetTotal_HeatFlux();
        Total_MaxHeat      = solver_container[val_iZone][val_iInst][FinestMesh][HEAT_SOL]->GetTotal_MaxHeatFlux();
        Total_Temperature  = solver_container[val_iZone][val_iInst][FinestMesh][HEAT_SOL]->GetTotal_AvgTemperature();

        /*--- Heat Residuals ---*/

        for (iVar = 0; iVar < nVar_Heat; iVar++) {
          residual_heat[iVar] = solver_container[val_iZone][val_iInst][FinestMesh][HEAT_SOL]->GetRes_RMS(iVar);
        }

        break;

      case FEM_ELASTICITY:

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

      case DISC_ADJ_FEM:

        /*--- FEM coefficients -- As of now, this is the Von Mises Stress ---*/

        Total_VMStress = solver_container[val_iZone][val_iInst][FinestMesh][FEA_SOL]->GetTotal_CFEA();

        /*--- Residuals: ---*/
        /*--- Linear analysis: RMS of the displacements in the nDim coordinates ---*/
        /*--- Nonlinear analysis: UTOL, RTOL and DTOL (defined in the Postprocessing function) ---*/
         for (iVar = 0; iVar < nVar_FEM; iVar++) {
           residual_fem[iVar] = solver_container[val_iZone][val_iInst][FinestMesh][ADJFEA_SOL]->GetRes_RMS(iVar);
         }

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

    bool Unsteady = ((config[val_iZone]->GetTime_Marching() == DT_STEPPING_1ST) ||
                     (config[val_iZone]->GetTime_Marching() == DT_STEPPING_2ND));
    bool In_NoDualTime = (!DualTime_Iteration && (iExtIter % config[val_iZone]->GetWrt_Con_Freq() == 0));
    bool In_DualTime_0 = (DualTime_Iteration && (iIntIter % config[val_iZone]->GetWrt_Con_Freq_DualTime() == 0));
    bool In_DualTime_1 = (!DualTime_Iteration && Unsteady);
    bool In_DualTime_2 = (Unsteady && DualTime_Iteration && (iExtIter % config[val_iZone]->GetWrt_Con_Freq() == 0));
    bool In_DualTime_3 = (Unsteady && !DualTime_Iteration && (iExtIter % config[val_iZone]->GetWrt_Con_Freq() == 0));

    /*--- Header frequency: analogy for dynamic structural analysis ---*/
    /*--- DualTime_Iteration is a bool we receive, which is true if it comes from FEM_StructuralIteration and false from SU2_CFD ---*/
    /*--- We maintain the name, as it is an input of the function ---*/
    /*--- The function GetWrt_Con_Freq_DualTime should be modified to be able to define different frequencies ---*/
    /*--- dynamic determines if the problem is, or not, time dependent ---*/
    bool dynamic = (config[val_iZone]->GetTime_Domain());              // Dynamic simulations.
    bool In_NoDynamic = (!DualTime_Iteration && (iExtIter % config[val_iZone]->GetWrt_Con_Freq() == 0));
    bool In_Dynamic_0 = (DualTime_Iteration && (iIntIter % config[val_iZone]->GetWrt_Con_Freq_DualTime() == 0));
    bool In_Dynamic_1 = (!DualTime_Iteration && nonlinear_analysis);
    bool In_Dynamic_2 = (nonlinear_analysis && DualTime_Iteration && (iExtIter % config[val_iZone]->GetWrt_Con_Freq() == 0));
    bool In_Dynamic_3 = (nonlinear_analysis && !DualTime_Iteration && (iExtIter % config[val_iZone]->GetWrt_Con_Freq() == 0));

    bool write_heads;
    if (Unsteady) write_heads = (iIntIter == 0);
    else write_heads = (((iExtIter % (config[val_iZone]->GetWrt_Con_Freq()*40)) == 0));

    bool write_turbo = (((iExtIter % (config[val_iZone]->GetWrt_Con_Freq()*40)) == 0) || (iExtIter == (config[val_iZone]->GetnInner_Iter() -1)));

    /*--- Analogous for dynamic problems (as of now I separate the problems, it may be worthy to do all together later on ---*/
    bool write_heads_FEM;
    if (nonlinear_analysis) write_heads_FEM = (iIntIter == 0);
    else write_heads_FEM = (((iExtIter % (config[val_iZone]->GetWrt_Con_Freq()*40)) == 0));

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

          case EULER : case NAVIER_STOKES: case RANS:
          case INC_EULER : case INC_NAVIER_STOKES: case INC_RANS:
          case FEM_EULER : case FEM_NAVIER_STOKES: case FEM_RANS: case FEM_LES:
          case ADJ_EULER: case ADJ_NAVIER_STOKES: case ADJ_RANS: case DISC_ADJ_EULER:
          case DISC_ADJ_NAVIER_STOKES: case DISC_ADJ_RANS:
          case DISC_ADJ_INC_EULER: case DISC_ADJ_INC_NAVIER_STOKES: case DISC_ADJ_INC_RANS:

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

          case HEAT_EQUATION:

            SPRINTF (direct_coeff, ", %14.8e, %14.8e, %14.8e", Total_Heat, Total_MaxHeat, Total_Temperature);
            SPRINTF (heat_resid, ", %14.8e", log10 (residual_heat[0]));

            break;

          case FEM_ELASTICITY:

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

          case DISC_ADJ_FEM:

            SPRINTF (direct_coeff, ", %12.10f", Total_VMStress);
            if (nDim == 2) {
              SPRINTF (fem_resid, ", %14.8e, %14.8e, %14.8e, %14.8e, %14.8e", log10 (residual_fem[0]), log10 (residual_fem[1]), dummy, dummy, dummy);
            }
            else {
              SPRINTF (fem_resid, ", %14.8e, %14.8e, %14.8e, %14.8e, %14.8e", log10 (residual_fem[0]), log10 (residual_fem[1]), log10 (residual_fem[2]), dummy, dummy);
            }

            break;

        }
      }
      if ((val_iZone == 0 && val_iInst == 0)|| fluid_structure){
        /*--- Write the screen header---*/
        if (  (!fem && ((write_heads) && !(!DualTime_Iteration && Unsteady))) ||
            (fem && ((write_heads_FEM) && !(!DualTime_Iteration && nonlinear_analysis)))
        ) {

          if (!fem) {
            if (!Unsteady && (config[val_iZone]->GetTime_Marching() != TIME_STEPPING)) {
              switch (config[val_iZone]->GetKind_Solver()) {
              case EULER : case NAVIER_STOKES: case RANS:
              case INC_EULER : case INC_NAVIER_STOKES: case INC_RANS:
              case FEM_EULER : case FEM_NAVIER_STOKES: case FEM_RANS: case FEM_LES:
              case ADJ_EULER : case ADJ_NAVIER_STOKES: case ADJ_RANS:

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

              case DISC_ADJ_EULER: case DISC_ADJ_NAVIER_STOKES: case DISC_ADJ_RANS:
              case DISC_ADJ_INC_EULER: case DISC_ADJ_INC_NAVIER_STOKES: case DISC_ADJ_INC_RANS:

                cout << endl;
                cout << "------------------------ Discrete Adjoint Summary -----------------------" << endl;
                cout << "Total Geometry Sensitivity (updated every "  << config[val_iZone]->GetWrt_Sol_Freq() << " iterations): ";
                cout.precision(4);
                cout.setf(ios::scientific, ios::floatfield);
                cout << Total_Sens_Geo;
                cout << endl << "-------------------------------------------------------------------------" << endl;
                break;

              }
            }
            else {
              if (flow) {
                if ((config[val_iZone]->GetTime_Marching() == TIME_STEPPING) && (config[val_iZone]->GetUnst_CFL()== 0.0))
                {
                  cout << endl << "Min DT: " << solver_container[val_iZone][val_iInst][FinestMesh][FLOW_SOL]->GetMin_Delta_Time()<< ".Max DT: " << solver_container[val_iZone][val_iInst][FinestMesh][FLOW_SOL]->GetMax_Delta_Time() << ".Time step: " << config[val_iZone]->GetDelta_UnstTimeND() << ".";
                } else if ((config[val_iZone]->GetTime_Marching() == TIME_STEPPING) && (config[val_iZone]->GetUnst_CFL()!= 0.0)) {
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
          case EULER :                  case NAVIER_STOKES:
          case INC_EULER :              case INC_NAVIER_STOKES:
          case FEM_EULER : case FEM_NAVIER_STOKES: case FEM_LES:

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

          case RANS : case INC_RANS:

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
              case SST:	case SST_SUST: cout << "     Res[kine]" << "     Res[omega]"; break;
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

            case HEAT_EQUATION :
              if (!Unsteady) cout << endl << " Iter" << "    Time(s)";
              else cout << endl << " IntIter" << "  ExtIter";

              cout <<  "      Res[Heat]" << "   HFlux(Total)";
              break;

            case FEM_ELASTICITY :
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
           break;

            case ADJ_EULER :              case ADJ_NAVIER_STOKES :
            case DISC_ADJ_EULER:          case DISC_ADJ_NAVIER_STOKES:
            case DISC_ADJ_INC_EULER:      case DISC_ADJ_INC_NAVIER_STOKES:

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

            case ADJ_RANS : case DISC_ADJ_RANS: case DISC_ADJ_INC_RANS:

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
        case EULER : case NAVIER_STOKES:
        case INC_EULER : case INC_NAVIER_STOKES:
        case FEM_EULER : case FEM_NAVIER_STOKES: case FEM_LES:

          /*--- Write history file ---*/

          if ((!DualTime_Iteration) && (output_files)) {
            if (!turbo) {
              config[val_iZone]->GetHistFile()[0] << begin << direct_coeff;
              if (buffet) config[val_iZone]->GetHistFile()[0] << buffet_coeff;
              if (thermal) config[val_iZone]->GetHistFile()[0] << heat_coeff;
              if (equiv_area) config[val_iZone]->GetHistFile()[0] << equivalent_area_coeff;
              if (engine || actuator_disk) config[val_iZone]->GetHistFile()[0] << engine_coeff;
              if (inv_design) {
                config[val_iZone]->GetHistFile()[0] << Cp_inverse_design;
                if (thermal) config[val_iZone]->GetHistFile()[0] << Heat_inverse_design;
              }
              if (rotating_frame && !turbo) config[val_iZone]->GetHistFile()[0] << rotating_frame_coeff;
              config[val_iZone]->GetHistFile()[0] << flow_resid;
              if (weakly_coupled_heat) config[val_iZone]->GetHistFile()[0] << heat_resid;
              if (radiation) config[val_iZone]->GetHistFile()[0] << rad_resid;
            }
            else {
              config[val_iZone]->GetHistFile()[0] << begin << turbo_coeff << flow_resid;
            }

            if (aeroelastic) config[val_iZone]->GetHistFile()[0] << aeroelastic_coeff;
            if (output_per_surface) config[val_iZone]->GetHistFile()[0] << monitoring_coeff;
            if (output_surface) config[val_iZone]->GetHistFile()[0] << surface_outputs;
            if (direct_diff != NO_DERIVATIVE) {
              config[val_iZone]->GetHistFile()[0] << d_direct_coeff;
              if (output_surface) config[val_iZone]->GetHistFile()[0] << d_surface_outputs;
            }
            if (output_comboObj) config[val_iZone]->GetHistFile()[0] << combo_obj;
            config[val_iZone]->GetHistFile()[0] << end;
            config[val_iZone]->GetHistFile()[0].flush();
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

        case RANS : case INC_RANS:

          /*--- Write history file ---*/

          if ((!DualTime_Iteration) && (output_files)) {

            if (!turbo) {
              config[val_iZone]->GetHistFile()[0] << begin << direct_coeff;
              if (buffet) config[val_iZone]->GetHistFile()[0] << buffet_coeff;
              if (thermal) config[val_iZone]->GetHistFile()[0] << heat_coeff;
              if (equiv_area) config[val_iZone]->GetHistFile()[0] << equivalent_area_coeff;
              if (engine || actuator_disk) config[val_iZone]->GetHistFile()[0] << engine_coeff;
              if (inv_design) {
                config[val_iZone]->GetHistFile()[0] << Cp_inverse_design;
                if (thermal) config[val_iZone]->GetHistFile()[0] << Heat_inverse_design;
              }
              if (rotating_frame && !turbo) config[val_iZone]->GetHistFile()[0] << rotating_frame_coeff;
              config[val_iZone]->GetHistFile()[0] << flow_resid << turb_resid;
              if (weakly_coupled_heat) config[val_iZone]->GetHistFile()[0] << heat_resid;
            }
            else {
              config[val_iZone]->GetHistFile()[0] << begin << turbo_coeff << flow_resid << turb_resid;
            }

            if (aeroelastic) config[val_iZone]->GetHistFile()[0] << aeroelastic_coeff;
            if (output_per_surface) config[val_iZone]->GetHistFile()[0] << monitoring_coeff;
            if (output_surface) config[val_iZone]->GetHistFile()[0] << surface_outputs;
            if (direct_diff != NO_DERIVATIVE) {
              config[val_iZone]->GetHistFile()[0] << d_direct_coeff;
              if (output_surface) config[val_iZone]->GetHistFile()[0] << d_surface_outputs;
            }
            if (output_comboObj) config[val_iZone]->GetHistFile()[0] << combo_obj;
            config[val_iZone]->GetHistFile()[0] << end;
            config[val_iZone]->GetHistFile()[0].flush();
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
            config[val_iZone]->GetHistFile()[0] << begin << direct_coeff << heat_resid << end;
            config[val_iZone]->GetHistFile()[0].flush();
          }
          break;

        case FEM_ELASTICITY:

          if (!DualTime_Iteration) {
            config[val_iZone]->GetHistFile()[0] << begin << fem_coeff << fem_resid << end_fem;
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
          break;

        case DISC_ADJ_FEM:

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

        case ADJ_EULER :              case ADJ_NAVIER_STOKES :
        case DISC_ADJ_EULER:          case DISC_ADJ_NAVIER_STOKES:
        case DISC_ADJ_INC_EULER:      case DISC_ADJ_INC_NAVIER_STOKES:

          if (!DualTime_Iteration) {
            config[val_iZone]->GetHistFile()[0] << begin << adjoint_coeff << adj_flow_resid << end;
            config[val_iZone]->GetHistFile()[0].flush();
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

        case ADJ_RANS : case DISC_ADJ_RANS: case DISC_ADJ_INC_RANS:

          if (!DualTime_Iteration) {
            config[val_iZone]->GetHistFile()[0] << begin << adjoint_coeff << adj_flow_resid;
            if (!frozen_visc)
              config[val_iZone]->GetHistFile()[0] << adj_turb_resid;
            config[val_iZone]->GetHistFile()[0] << end;
            config[val_iZone]->GetHistFile()[0].flush();
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

void COutputLegacy::SpecialOutput_ForcesBreakdown(CSolver *****solver, CGeometry ****geometry, CConfig **config, unsigned short val_iZone, bool output) {

  char cstr[200];
  unsigned short iDim, iMarker_Monitoring;
  ofstream Breakdown_file;

  bool compressible       = (config[val_iZone]->GetKind_Regime() == COMPRESSIBLE);
  bool incompressible     = (config[val_iZone]->GetKind_Regime() == INCOMPRESSIBLE);
  bool unsteady           = (config[val_iZone]->GetTime_Marching() != NO);
  bool viscous            = config[val_iZone]->GetViscous();
  bool dynamic_grid       = config[val_iZone]->GetDynamic_Grid();
  bool gravity            = config[val_iZone]->GetGravityForce();
  bool turbulent          = config[val_iZone]->GetKind_Solver() == RANS;
  bool fixed_cl           = config[val_iZone]->GetFixed_CL_Mode();
  unsigned short Kind_Solver = config[val_iZone]->GetKind_Solver();
  unsigned short Kind_Turb_Model = config[val_iZone]->GetKind_Turb_Model();
  unsigned short Ref_NonDim = config[val_iZone]->GetRef_NonDim();

  unsigned short FinestMesh = config[val_iZone]->GetFinestMesh();
  unsigned short nDim = geometry[val_iZone][INST_0][FinestMesh]->GetnDim();
  bool flow = ((config[val_iZone]->GetKind_Solver() == EULER) ||
               (config[val_iZone]->GetKind_Solver() == NAVIER_STOKES) ||
               (config[val_iZone]->GetKind_Solver() == RANS) ||
               (config[val_iZone]->GetKind_Solver() == INC_EULER) ||
               (config[val_iZone]->GetKind_Solver() == INC_NAVIER_STOKES) ||
               (config[val_iZone]->GetKind_Solver() == INC_RANS));

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
    0.0, *Surface_CL = NULL, *Surface_CD = NULL,
    *Surface_CSF = NULL, *Surface_CEff = NULL, *Surface_CFx = NULL,
    *Surface_CFy = NULL, *Surface_CFz = NULL,
    *Surface_CMx = NULL, *Surface_CMy = NULL, *Surface_CMz = NULL,
    *Surface_CL_Inv = NULL,
    *Surface_CD_Inv = NULL, *Surface_CSF_Inv = NULL,
    *Surface_CEff_Inv = NULL, *Surface_CFx_Inv = NULL, *Surface_CFy_Inv =
    NULL, *Surface_CFz_Inv = NULL, *Surface_CMx_Inv = NULL,
    *Surface_CMy_Inv = NULL, *Surface_CMz_Inv = NULL,
    *Surface_CL_Visc = NULL,
    *Surface_CD_Visc = NULL, *Surface_CSF_Visc = NULL,
    *Surface_CEff_Visc = NULL, *Surface_CFx_Visc = NULL, *Surface_CFy_Visc =
    NULL, *Surface_CFz_Visc = NULL, *Surface_CMx_Visc = NULL,
    *Surface_CMy_Visc = NULL, *Surface_CMz_Visc = NULL,
    *Surface_CL_Mnt = NULL,
    *Surface_CD_Mnt = NULL, *Surface_CSF_Mnt = NULL,
    *Surface_CEff_Mnt = NULL, *Surface_CFx_Mnt = NULL, *Surface_CFy_Mnt =
    NULL, *Surface_CFz_Mnt = NULL, *Surface_CMx_Mnt = NULL,
    *Surface_CMy_Mnt = NULL, *Surface_CMz_Mnt = NULL;

    /*--- WARNING: when compiling on Windows, ctime() is not available. Comment out
     the two lines below that use the dt variable. ---*/
    //time_t now = time(0);
    //string dt = ctime(&now); dt[24] = '.';

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
    Breakdown_file <<"|   / __| | | |_  )   Release 6.1.0  \"Falcon\"                           |" << "\n";
    Breakdown_file <<"|   \\__ \\ |_| |/ /                                                      |" << "\n";
    Breakdown_file <<"|   |___/\\___//___|   Suite (Computational Fluid Dynamics Code)         |" << "\n";
    Breakdown_file << "|                                                                       |" << "\n";
    //Breakdown_file << "|   Local date and time: " << dt << "                      |" << "\n";
    Breakdown_file <<"-------------------------------------------------------------------------" << "\n";
    Breakdown_file << "| The current SU2 release has been coordinated by the                   |" << "\n";
    Breakdown_file << "| SU2 International Developers Society <www.su2devsociety.org>          |" << "\n";
    Breakdown_file << "| with selected contributions from the open-source community            |" << "\n";
    Breakdown_file <<"-------------------------------------------------------------------------" << "\n";
    Breakdown_file << "| The main research teams contributing to the current release are:      |" << "\n";
    Breakdown_file << "| - Prof. Juan J. Alonso's group at Stanford University.                |" << "\n";
    Breakdown_file << "| - Prof. Piero Colonna's group at Delft University of Technology.      |" << "\n";
    Breakdown_file << "| - Prof. Nicolas R. Gauger's group at Kaiserslautern U. of Technology. |" << "\n";
    Breakdown_file << "| - Prof. Alberto Guardone's group at Polytechnic University of Milan.  |" << "\n";
    Breakdown_file << "| - Prof. Rafael Palacios' group at Imperial College London.            |" << "\n";
    Breakdown_file << "| - Prof. Vincent Terrapon's group at the University of Liege.          |" << "\n";
    Breakdown_file << "| - Prof. Edwin van der Weide's group at the University of Twente.      |" << "\n";
    Breakdown_file << "| - Lab. of New Concepts in Aeronautics at Tech. Inst. of Aeronautics.  |" << "\n";
    Breakdown_file <<"-------------------------------------------------------------------------" << "\n";
    Breakdown_file << "| Copyright 2012-2018, Francisco D. Palacios, Thomas D. Economon,       |" << "\n";
    Breakdown_file << "|                      Tim Albring, and the SU2 contributors.           |" << "\n";
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
      case EULER: case INC_EULER:
        if (compressible) Breakdown_file << "Compressible Euler equations." << "\n";
        if (incompressible) Breakdown_file << "Incompressible Euler equations." << "\n";
        break;
      case NAVIER_STOKES: case INC_NAVIER_STOKES:
        if (compressible) Breakdown_file << "Compressible Laminar Navier-Stokes' equations." << "\n";
        if (incompressible) Breakdown_file << "Incompressible Laminar Navier-Stokes' equations." << "\n";
        break;
      case RANS: case INC_RANS:
        if (compressible) Breakdown_file << "Compressible RANS equations." << "\n";
        if (incompressible) Breakdown_file << "Incompressible RANS equations." << "\n";
        Breakdown_file << "Turbulence model: ";
        switch (Kind_Turb_Model) {
          case SA:        Breakdown_file << "Spalart Allmaras" << "\n"; break;
          case SA_NEG:    Breakdown_file << "Negative Spalart Allmaras" << "\n"; break;
          case SA_E:      Breakdown_file << "Edwards Spalart Allmaras" << "\n"; break;
          case SA_COMP:   Breakdown_file << "Compressibility Correction Spalart Allmaras" << "\n"; break;
          case SA_E_COMP: Breakdown_file << "Compressibility Correction Edwards Spalart Allmaras" << "\n"; break;
          case SST:       Breakdown_file << "Menter's SST"     << "\n"; break;
          case SST_SUST:  Breakdown_file << "Menter's SST with sustaining terms" << "\n"; break;
        }
        break;
    }


    /*--- Compressible version of console output ---*/

    if (compressible) {


    if ((compressible) && (Kind_Solver != FEM_ELASTICITY)) {
      Breakdown_file << "Mach number: " << config[val_iZone]->GetMach() <<"."<< "\n";
      Breakdown_file << "Angle of attack (AoA): " << config[val_iZone]->GetAoA() <<" deg, and angle of sideslip (AoS): " << config[val_iZone]->GetAoS() <<" deg."<< "\n";
      if ((Kind_Solver == NAVIER_STOKES) || (Kind_Solver == ADJ_NAVIER_STOKES) ||
          (Kind_Solver == RANS) || (Kind_Solver == ADJ_RANS))
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

          case CONSTANT_VISCOSITY:
            Breakdown_file << "Viscosity Model: CONSTANT_VISCOSITY  "<< "\n";
            Breakdown_file << "Laminar Viscosity: " << config[val_iZone]->GetMu_Constant();
            if (config[val_iZone]->GetSystemMeasurements() == SI) Breakdown_file << " N.s/m^2." << "\n";
            else if (config[val_iZone]->GetSystemMeasurements() == US) Breakdown_file << " lbf.s/ft^2." << "\n";
            Breakdown_file << "Laminar Viscosity (non-dim): " << config[val_iZone]->GetMu_ConstantND()<< "\n";
            break;

          case SUTHERLAND:
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

        }
        switch (config[val_iZone]->GetKind_ConductivityModel()) {

          case CONSTANT_PRANDTL:
            Breakdown_file << "Conductivity Model: CONSTANT_PRANDTL  "<< "\n";
            Breakdown_file << "Prandtl: " << config[val_iZone]->GetPrandtl_Lam()<< "\n";
            break;

          case CONSTANT_CONDUCTIVITY:
            Breakdown_file << "Conductivity Model: CONSTANT_CONDUCTIVITY "<< "\n";
            Breakdown_file << "Molecular Conductivity: " << config[val_iZone]->GetKt_Constant()<< " W/m^2.K." << "\n";
            Breakdown_file << "Molecular Conductivity (non-dim): " << config[val_iZone]->GetKt_ConstantND()<< "\n";
            break;

        }

        if ((Kind_Solver == RANS) || (Kind_Solver == ADJ_RANS) || (Kind_Solver == DISC_ADJ_RANS)) {
          switch (config[val_iZone]->GetKind_ConductivityModel_Turb()) {
            case CONSTANT_PRANDTL_TURB:
              Breakdown_file << "Turbulent Conductivity Model: CONSTANT_PRANDTL_TURB  "<< "\n";
              Breakdown_file << "Turbulent Prandtl: " << config[val_iZone]->GetPrandtl_Turb()<< "\n";
              break;
            case NO_CONDUCTIVITY_TURB:
              Breakdown_file << "Turbulent Conductivity Model: NO_CONDUCTIVITY_TURB "<< "\n";
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
        Breakdown_file << "Reference conductivity: " << config[val_iZone]->GetConductivity_Ref();
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
      bool boussinesq = (config[val_iZone]->GetKind_DensityModel() == BOUSSINESQ);

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

        case CONSTANT:
          if (energy) Breakdown_file << "Energy equation is active and decoupled." << "\n";
          else Breakdown_file << "No energy equation." << "\n";
          break;

        case BOUSSINESQ:
          if (energy) Breakdown_file << "Energy equation is active and coupled through Boussinesq approx." << "\n";
          break;

        case VARIABLE:
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

          case CONSTANT_VISCOSITY:
            Breakdown_file << "Viscosity Model: CONSTANT_VISCOSITY  "<< "\n";
            Breakdown_file << "Constant Laminar Viscosity: " << config[val_iZone]->GetMu_Constant();
            if (config[val_iZone]->GetSystemMeasurements() == SI) Breakdown_file << " N.s/m^2." << "\n";
            else if (config[val_iZone]->GetSystemMeasurements() == US) Breakdown_file << " lbf.s/ft^2." << "\n";
            Breakdown_file << "Laminar Viscosity (non-dim): " << config[val_iZone]->GetMu_ConstantND()<< "\n";
            break;

          case SUTHERLAND:
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

          case POLYNOMIAL_VISCOSITY:
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

            case CONSTANT_PRANDTL:
              Breakdown_file << "Conductivity Model: CONSTANT_PRANDTL  "<< "\n";
              Breakdown_file << "Prandtl (Laminar): " << config[val_iZone]->GetPrandtl_Lam()<< "\n";
              break;

            case CONSTANT_CONDUCTIVITY:
              Breakdown_file << "Conductivity Model: CONSTANT_CONDUCTIVITY "<< "\n";
              Breakdown_file << "Molecular Conductivity: " << config[val_iZone]->GetKt_Constant()<< " W/m^2.K." << "\n";
              Breakdown_file << "Molecular Conductivity (non-dim): " << config[val_iZone]->GetKt_ConstantND()<< "\n";
              break;

            case POLYNOMIAL_CONDUCTIVITY:
              Breakdown_file << "Viscosity Model: POLYNOMIAL_CONDUCTIVITY "<< endl;
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

          if ((Kind_Solver == RANS) || (Kind_Solver == ADJ_RANS) || (Kind_Solver == DISC_ADJ_RANS)) {
            switch (config[val_iZone]->GetKind_ConductivityModel_Turb()) {
              case CONSTANT_PRANDTL_TURB:
                Breakdown_file << "Turbulent Conductivity Model: CONSTANT_PRANDTL_TURB  "<< "\n";
                Breakdown_file << "Turbulent Prandtl: " << config[val_iZone]->GetPrandtl_Turb()<< "\n";
                break;
              case NO_CONDUCTIVITY_TURB:
                Breakdown_file << "Turbulent Conductivity Model: NO_CONDUCTIVITY_TURB "<< "\n";
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

void COutputLegacy::SetResult_Files(CSolver *****solver_container, CGeometry ****geometry, CConfig **config,
                              unsigned long iExtIter, unsigned short val_nZone) {

  unsigned short iZone;

  for (iZone = 0; iZone < val_nZone; iZone++) {

    /*--- Flags identifying the types of files to be written. ---*/

    bool Wrt_Vol = config[iZone]->GetWrt_Vol_Sol();
    bool Wrt_Srf = config[iZone]->GetWrt_Srf_Sol();
    bool Wrt_Csv = config[iZone]->GetWrt_Csv_Sol();

#ifdef HAVE_MPI
    /*--- Do not merge the volume solutions if we are running in parallel.
     Force the use of SU2_SOL to merge the volume sols in this case. ---*/

    SU2_MPI::Comm_size(MPI_COMM_WORLD, &size);
    if (size > SINGLE_NODE) {
      Wrt_Vol = false;
      Wrt_Srf = false;
    }
#endif

    if (rank == MASTER_NODE) cout << endl << "Writing comma-separated values (CSV) surface files." << endl;

    switch (config[iZone]->GetKind_Solver()) {

      case EULER : case NAVIER_STOKES : case RANS :
      case INC_EULER : case INC_NAVIER_STOKES : case INC_RANS :
      case FEM_EULER : case FEM_NAVIER_STOKES : case FEM_RANS : case FEM_LES :

        if (Wrt_Csv) SetSurfaceCSV_Flow(config[iZone], geometry[iZone][INST_0][MESH_0], solver_container[iZone][INST_0][MESH_0][FLOW_SOL], iExtIter, iZone, INST_0);
        break;


      case ADJ_EULER : case ADJ_NAVIER_STOKES : case ADJ_RANS : case DISC_ADJ_EULER: case DISC_ADJ_NAVIER_STOKES: case DISC_ADJ_RANS:
         case DISC_ADJ_INC_NAVIER_STOKES: case DISC_ADJ_INC_RANS: case DISC_ADJ_INC_EULER:
        if (Wrt_Csv) SetSurfaceCSV_Adjoint(config[iZone], geometry[iZone][INST_0][MESH_0], solver_container[iZone][INST_0][MESH_0][ADJFLOW_SOL], solver_container[iZone][INST_0][MESH_0][FLOW_SOL], iExtIter, iZone, INST_0);
        break;

    }

    /*--- Get the file output format ---*/

    unsigned short FileFormat = config[iZone]->GetTabular_FileFormat();

    /*--- Merge the node coordinates and connectivity, if necessary. This
     is only performed if a volume solution file is requested, and it
     is active by default. ---*/

    if (Wrt_Vol || Wrt_Srf) {
      if (rank == MASTER_NODE) cout << "Merging connectivities in the Master node." << endl;
      MergeConnectivity(config[iZone], geometry[iZone][INST_0][MESH_0], iZone);
    }

    /*--- Merge coordinates of all grid nodes (excluding ghost points).
     The grid coordinates are always merged and included first in the
     restart files. ---*/

    if (rank == MASTER_NODE) cout << "Merging coordinates in the Master node." << endl;
    MergeCoordinates(config[iZone], geometry[iZone][INST_0][MESH_0]);

    if ((rank == MASTER_NODE) && (Wrt_Vol || Wrt_Srf)) {
      if (FileFormat == TECPLOT_BINARY) {
        if (rank == MASTER_NODE) cout << "Writing Tecplot binary volume and surface mesh files." << endl;
//        SetTecplotBinary_DomainMesh(config[iZone], geometry[iZone][INST_0][MESH_0], iZone);
//        SetTecplotBinary_SurfaceMesh(config[iZone], geometry[iZone][INST_0][MESH_0], iZone);
        if (!wrote_base_file)
          DeallocateConnectivity(config[iZone], geometry[iZone][INST_0][MESH_0], false);
        if (!wrote_surf_file)
          DeallocateConnectivity(config[iZone], geometry[iZone][INST_0][MESH_0], wrote_surf_file);
      }
    }

    /*--- Merge the solution data needed for volume solutions and restarts ---*/

    if (rank == MASTER_NODE) cout << "Merging solution in the Master node." << endl;
    MergeSolution(config[iZone], geometry[iZone][INST_0][MESH_0], solver_container[iZone][INST_0][MESH_0], iZone);

    /*--- Write restart, or Tecplot files using the merged data.
     This data lives only on the master, and these routines are currently
     executed by the master proc alone (as if in serial). ---*/

    if (rank == MASTER_NODE) {

      /*--- Write a native restart file ---*/

      if (rank == MASTER_NODE) cout << "Writing SU2 native restart file." << endl;
      SetRestart(config[iZone], geometry[iZone][INST_0][MESH_0], solver_container[iZone][INST_0][MESH_0] , iZone);

      if (Wrt_Vol) {

        switch (FileFormat) {

          case TECPLOT:

            /*--- Write a Tecplot ASCII file ---*/

            if (rank == MASTER_NODE) cout << "Writing Tecplot ASCII file volume solution file." << endl;
//            SetTecplotASCII(config[iZone], geometry[iZone][INST_0][MESH_0], solver_container[iZone][INST_0][MESH_0], iZone, val_nZone, false);
            DeallocateConnectivity(config[iZone], geometry[iZone][INST_0][MESH_0], false);
            break;



          case TECPLOT_BINARY:

            /*--- Write a Tecplot binary solution file ---*/

            if (rank == MASTER_NODE) cout << "Writing Tecplot binary volume solution file." << endl;
//            SetTecplotBinary_DomainSolution(config[iZone], geometry[iZone][INST_0][MESH_0], iZone);
            break;


          case PARAVIEW:

            /*--- Write a Paraview ASCII file ---*/

            if (rank == MASTER_NODE) cout << "Writing Paraview ASCII volume solution file." << endl;
//            SetParaview_ASCII(config[iZone], geometry[iZone][INST_0][MESH_0], iZone, val_nZone, false);
            DeallocateConnectivity(config[iZone], geometry[iZone][INST_0][MESH_0], false);
            break;

          case PARAVIEW_BINARY:

            /*--- Write a ParaView ASCII file instead for now in serial. ---*/

            if (rank == MASTER_NODE) cout << "ParaView binary volume files not available in this mode." << endl;
            if (rank == MASTER_NODE) cout << " Writing ParaView ASCII volume solution file instead." << endl;
//            SetParaview_ASCII(config[iZone], geometry[iZone][INST_0][MESH_0], iZone, val_nZone, false);
            DeallocateConnectivity(config[iZone], geometry[iZone][INST_0][MESH_0], false);
            break;

          default:
            break;
        }

      }

      if (Wrt_Srf) {

        switch (FileFormat) {

          case TECPLOT:

            /*--- Write a Tecplot ASCII file ---*/

            if (rank == MASTER_NODE) cout << "Writing Tecplot ASCII surface solution file." << endl;
//            SetTecplotASCII(config[iZone], geometry[iZone][INST_0][MESH_0], solver_container[iZone][INST_0][MESH_0] , iZone, val_nZone, true);
            DeallocateConnectivity(config[iZone], geometry[iZone][INST_0][MESH_0], true);
            break;

          case TECPLOT_BINARY:

            /*--- Write a Tecplot binary solution file ---*/

            if (rank == MASTER_NODE) cout << "Writing Tecplot binary surface solution file." << endl;
//            SetTecplotBinary_SurfaceSolution(config[iZone], geometry[iZone][INST_0][MESH_0], iZone);
            break;

          case PARAVIEW:

            /*--- Write a Paraview ASCII file ---*/

            if (rank == MASTER_NODE) cout << "Writing Paraview ASCII surface solution file." << endl;
//            SetParaview_ASCII(config[iZone], geometry[iZone][INST_0][MESH_0], iZone, val_nZone, true);
            DeallocateConnectivity(config[iZone], geometry[iZone][INST_0][MESH_0], true);
            break;

          case PARAVIEW_BINARY:

            /*--- Write a ParaView ASCII file instead for now in serial. ---*/

            if (rank == MASTER_NODE) cout << "ParaView binary surface files not available in this mode." << endl;
            if (rank == MASTER_NODE) cout << " Writing ParaView ASCII surface solution file instead." << endl;
//            SetParaview_ASCII(config[iZone], geometry[iZone][INST_0][MESH_0], iZone, val_nZone, true);
            DeallocateConnectivity(config[iZone], geometry[iZone][INST_0][MESH_0], true);
            break;

          default:
            break;
        }

      }

      /*--- Release memory needed for merging the solution data. ---*/

      DeallocateCoordinates(config[iZone], geometry[iZone][INST_0][MESH_0]);
      DeallocateSolution(config[iZone], geometry[iZone][INST_0][MESH_0]);

    }

    /*--- Final broadcast (informing other procs that the base output
     file was written). ---*/

#ifdef HAVE_MPI
    SU2_MPI::Bcast(&wrote_base_file, 1, MPI_UNSIGNED_SHORT, MASTER_NODE, MPI_COMM_WORLD);
    SU2_MPI::Bcast(&wrote_surf_file, 1, MPI_UNSIGNED_SHORT, MASTER_NODE, MPI_COMM_WORLD);
#endif

  }
}

void COutputLegacy::SetBaselineResult_Files(CSolver ***solver, CGeometry ***geometry, CConfig **config,
                                      unsigned long iExtIter, unsigned short val_nZone) {

  unsigned short iZone, iInst, nInst;

  for (iZone = 0; iZone < val_nZone; iZone++) {

    nInst = config[iZone]->GetnTimeInstances();

    for (iInst = 0; iInst < nInst; iInst++) {

      config[iZone]->SetiInst(iInst);

      /*--- Flags identifying the types of files to be written. ---*/

      bool Wrt_Vol = config[iZone]->GetWrt_Vol_Sol();
      if (config[iZone]->GetKind_SU2() == SU2_DOT) { Wrt_Vol = false; }
      bool Wrt_Srf = config[iZone]->GetWrt_Srf_Sol();

      /*--- Get the file output format ---*/

      unsigned short FileFormat = config[iZone]->GetTabular_FileFormat();

      /*--- Merge the node coordinates and connectivity if necessary. This
     is only performed if a volume solution file is requested, and it
     is active by default. ---*/

      if ((Wrt_Vol || Wrt_Srf)) {
        if (rank == MASTER_NODE) cout << "Merging connectivities in the Master node." << endl;
        MergeConnectivity(config[iZone], geometry[iZone][iInst], iZone);
      }

      /*--- Merge the solution data needed for volume solutions and restarts ---*/

      if ((Wrt_Vol || Wrt_Srf)) {
        if (rank == MASTER_NODE) cout << "Merging solution in the Master node." << endl;
        MergeBaselineSolution(config[iZone], geometry[iZone][iInst], solver[iZone][iInst], iZone);
      }

      /*--- Write restart, Tecplot or Paraview files using the merged data.
     This data lives only on the master, and these routines are currently
     executed by the master proc alone (as if in serial). ---*/


      if (rank == MASTER_NODE) {

        if (Wrt_Vol) {

          switch (FileFormat) {

          case TECPLOT:

            /*--- Write a Tecplot ASCII file ---*/

            if (rank == MASTER_NODE) cout << "Writing Tecplot ASCII file (volume grid)." << endl;
//            SetTecplotASCII(config[iZone], geometry[iZone][iInst], &solver[iZone][iInst], iZone, val_nZone, false);
            DeallocateConnectivity(config[iZone], geometry[iZone][iInst], false);
            break;

          case TECPLOT_BINARY:

            /*--- Write a Tecplot binary solution file ---*/

            if (rank == MASTER_NODE) cout << "Writing Tecplot Binary file (volume grid)." << endl;
//            SetTecplotBinary_DomainMesh(config[iZone], geometry[iZone][iInst], iZone);
//            SetTecplotBinary_DomainSolution(config[iZone], geometry[iZone][iInst], iZone);
            break;


          case PARAVIEW:

            /*--- Write a Paraview ASCII file ---*/

            if (rank == MASTER_NODE) cout << "Writing Paraview ASCII file (volume grid)." << endl;
//            SetParaview_ASCII(config[iZone], geometry[iZone][iInst], iZone, val_nZone, false);
            DeallocateConnectivity(config[iZone], geometry[iZone][iInst], false);
            break;

          case PARAVIEW_BINARY:

            /*--- Write a ParaView ASCII file instead for now in serial. ---*/

            if (rank == MASTER_NODE) cout << "ParaView binary volume files not available in this mode." << endl;
            if (rank == MASTER_NODE) cout << " Writing ParaView ASCII volume solution file instead." << endl;
//            SetParaview_ASCII(config[iZone], geometry[iZone][iInst], iZone, val_nZone, false);
            DeallocateConnectivity(config[iZone], geometry[iZone][iInst], false);
            break;

          default:
            break;
          }

        }

        if (Wrt_Srf) {

          switch (FileFormat) {

          case TECPLOT:

            /*--- Write a Tecplot ASCII file ---*/

            if (rank == MASTER_NODE) cout << "Writing Tecplot ASCII file (surface grid)." << endl;
//            SetTecplotASCII(config[iZone], geometry[iZone][iInst], &solver[iZone][iInst], iZone, val_nZone, true);
            DeallocateConnectivity(config[iZone], geometry[iZone][iInst], true);
            break;

          case TECPLOT_BINARY:

            /*--- Write a Tecplot binary solution file ---*/

            if (rank == MASTER_NODE) cout << "Writing Tecplot Binary file (surface grid)." << endl;
//            SetTecplotBinary_SurfaceMesh(config[iZone], geometry[iZone][iInst], iZone);
//            SetTecplotBinary_SurfaceSolution(config[iZone], geometry[iZone][iInst], iZone);
            break;

          case PARAVIEW:

            /*--- Write a Paraview ASCII file ---*/

            if (rank == MASTER_NODE) cout << "Writing Paraview ASCII file (surface grid)." << endl;
//            SetParaview_ASCII(config[iZone], geometry[iZone][iInst], iZone, val_nZone, true);
            DeallocateConnectivity(config[iZone], geometry[iZone][iInst], true);
            break;

          case PARAVIEW_BINARY:

            /*--- Write a ParaView ASCII file instead for now in serial. ---*/

            if (rank == MASTER_NODE) cout << "ParaView binary surface files not available in this mode." << endl;
            if (rank == MASTER_NODE) cout << " Writing ParaView ASCII surface solution file instead." << endl;
//            SetParaview_ASCII(config[iZone], geometry[iZone][iInst], iZone, val_nZone, true);
            DeallocateConnectivity(config[iZone], geometry[iZone][iInst], true);
            break;

          default:
            break;
          }
        }

        if (FileFormat == TECPLOT_BINARY) {
          if (!wrote_base_file)
            DeallocateConnectivity(config[iZone], geometry[iZone][iInst], false);
          if (!wrote_surf_file)
            DeallocateConnectivity(config[iZone], geometry[iZone][iInst], wrote_surf_file);
        }

        if (Wrt_Vol || Wrt_Srf)
          DeallocateSolution(config[iZone], geometry[iZone][iInst]);
      }



      /*--- Final broadcast (informing other procs that the base output
     file was written). ---*/

#ifdef HAVE_MPI
      SU2_MPI::Bcast(&wrote_base_file, 1, MPI_UNSIGNED_SHORT, MASTER_NODE, MPI_COMM_WORLD);
#endif

    }

  }
}

void COutputLegacy::SetMesh_Files(CGeometry **geometry, CConfig **config, unsigned short val_nZone, bool new_file, bool su2_file) {

  char cstr[MAX_STRING_SIZE], out_file[MAX_STRING_SIZE];
  unsigned short iZone;
  ofstream output_file;
  string str;

  /*--- Read the name of the output and input file ---*/

  if (su2_file) {
    if (rank == MASTER_NODE) {
      str = config[ZONE_0]->GetMesh_Out_FileName();
      strcpy (out_file, str.c_str());
      strcpy (cstr, out_file);
      output_file.precision(15);
      output_file.open(cstr, ios::out);
      if (val_nZone > 1) {
        output_file << "NZONE= " << val_nZone << endl;
      }
    }
  }

  for (iZone = 0; iZone < val_nZone; iZone++) {

    /*--- Flags identifying the types of files to be written. ---*/

    bool Wrt_Vol = config[iZone]->GetVisualize_Volume_Def();
    bool Wrt_Srf = config[iZone]->GetVisualize_Surface_Def();
    //bool Wrt_Crd = config[iZone]->GetWrt_Crd_Sol();

    /*--- Merge the node coordinates and connectivity if necessary. This
     is only performed if a volume solution file is requested, and it
     is active by default. ---*/

    if (rank == MASTER_NODE) cout <<"Merging grid connectivity." << endl;
    MergeConnectivity(config[iZone], geometry[iZone], iZone);

    /*--- Merge coordinates of all grid nodes (excluding ghost points).
     The grid coordinates are always merged and included first in the
     restart files. ---*/

    if (rank == MASTER_NODE) cout <<"Merging grid coordinates." << endl;
    MergeCoordinates(config[iZone], geometry[iZone]);

    /*--- Write restart, Tecplot or Paraview files using the merged data.
     This data lives only on the master, and these routines are currently
     executed by the master proc alone (as if in serial). ---*/

    if (rank == MASTER_NODE) {

      if (Wrt_Vol) {

        if (rank == MASTER_NODE) cout <<"Writing volume mesh file." << endl;

//        /*--- Write a Tecplot ASCII file ---*/

////        if (config[iZone]->GetOutput_FileFormat() == PARAVIEW) SetParaview_MeshASCII(config[iZone], geometry[iZone], iZone,  val_nZone, false, new_file);
//        else if (config[iZone]->GetTabular_FileFormat() == PARAVIEW_BINARY) {
//          if (rank == MASTER_NODE) cout <<"Writing ASCII paraview volume mesh file by default." << endl;
////          SetParaview_MeshASCII(config[iZone], geometry[iZone], iZone,  val_nZone, false, new_file);
//        }
////        else if (config[iZone]->GetOutput_FileFormat() == TECPLOT) SetTecplotASCII_Mesh(config[iZone], geometry[iZone], iZone, false, new_file);
//        else if (config[iZone]->GetTabular_FileFormat() == TECPLOT_BINARY) {
//          if (rank == MASTER_NODE) cout <<"Writing ASCII tecplot volume mesh file by default." << endl;
////          SetTecplotASCII_Mesh(config[iZone], geometry[iZone], iZone, false, new_file);
//        }

      }

      if (Wrt_Srf) {

        if (rank == MASTER_NODE) cout <<"Writing surface mesh file." << endl;

        /*--- Write a Tecplot ASCII file ---*/

//        if (config[iZone]->GetOutput_FileFormat() == PARAVIEW) SetParaview_MeshASCII(config[iZone], geometry[iZone], iZone,  val_nZone, true, new_file);
//        else if (config[iZone]->GetOutput_FileFormat() == PARAVIEW_BINARY) {
//          if (rank == MASTER_NODE) cout <<"Writing ASCII paraview surface mesh file by default." << endl;
//          SetParaview_MeshASCII(config[iZone], geometry[iZone], iZone,  val_nZone, true, new_file);
//        }
//        else if (config[iZone]->GetOutput_FileFormat() == TECPLOT) SetTecplotASCII_Mesh(config[iZone], geometry[iZone], iZone, true, new_file);
//        else if (config[iZone]->GetOutput_FileFormat() == TECPLOT_BINARY) {
//          if (rank == MASTER_NODE) cout <<"Writing ASCII tecplot surface mesh file by default." << endl;
//          SetTecplotASCII_Mesh(config[iZone], geometry[iZone], iZone, true, new_file);
//        }

      }

      /*--- Write a .su2 ASCII file ---*/

      if (su2_file) {

        if (rank == MASTER_NODE) cout <<"Writing .su2 file." << endl;

//        SetSU2_MeshASCII(config[iZone], geometry[iZone], iZone, output_file);

        /*--- Write an stl surface file ---*/

        if (rank == MASTER_NODE) cout <<"Writing .stl surface file." << endl;

//        SetSTL_MeshASCII(config[iZone], geometry[iZone]);

      }

      /*--- Write a binary file with the grid coordinates alone. ---*/

//      if (Wrt_Crd) {
//        if (rank == MASTER_NODE) cout <<"Writing .dat binary coordinates file." << endl;
//        WriteCoordinates_Binary(config[iZone], geometry[iZone], iZone);
//      }


      /*--- Deallocate connectivity ---*/

      DeallocateConnectivity(config[iZone], geometry[iZone], true);
      DeallocateConnectivity(config[iZone], geometry[iZone], false);
      DeallocateCoordinates(config[iZone], geometry[iZone]);

    }

    /*--- Final broadcast (informing other procs that the base output
     file was written). ---*/

#ifdef HAVE_MPI
    SU2_MPI::Bcast(&wrote_base_file, 1, MPI_UNSIGNED_SHORT, MASTER_NODE, MPI_COMM_WORLD);
#endif

    /*--- Write an csv surface file, done in parallel ---*/

    if (rank == MASTER_NODE) cout <<"Writing .csv surface file." << endl;

//    if (su2_file) SetCSV_MeshASCII(config[iZone], geometry[iZone]);

  }

  if (rank == MASTER_NODE) {
    if (su2_file){
      output_file.close();
    }
  }
}

void COutputLegacy::SpecialOutput_SpanLoad(CSolver *solver, CGeometry *geometry, CConfig *config, bool output) {

  short iSection, nSection;
  unsigned long iVertex, iPoint, Trailing_Point;
  su2double *Plane_P0, *Plane_P0_, *Plane_Normal, *Plane_Normal_, *CPressure,
  Force[3], ForceInviscid[3], MomentInviscid[3] =
  { 0.0, 0.0, 0.0 }, MomentDist[3] = { 0.0, 0.0, 0.0 }, RefDensity,
  RefPressure, RefArea, *Velocity_Inf, Gas_Constant, Mach2Vel,
  Mach_Motion, Gamma, RefVel2 = 0.0, factor, NDPressure, *Origin,
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
  Origin = config->GetRefOriginMoment(0);
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

          MomentDist[0] = 0.5 	* (Xcoord_Airfoil[iVertex] + Xcoord_Airfoil[iVertex + 1]) - Origin[0];
          MomentDist[1] = 0.5 * (Ycoord_Airfoil[iVertex] + Ycoord_Airfoil[iVertex + 1]) - Origin[1];
          MomentDist[2] = 0.5 	* (Zcoord_Airfoil[iVertex] + Zcoord_Airfoil[iVertex + 1]) - Origin[3];

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

            if ((config->GetTabular_FileFormat() == TAB_CSV)) {
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
            if ((config->GetTabular_FileFormat() == TAB_CSV)) Load_File.open("load_distribution.csv", ios::app);
            else Load_File.open("load_distribution.dat", ios::app);
          }


          /*--- CL and spanload ---*/

          if ((config->GetTabular_FileFormat() == TAB_CSV))
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
  su2double XCoord, YCoord, ZCoord, Pressure, PressureCoeff = 0, Cp, CpTarget, *Normal = NULL, Area, PressDiff;
  bool *PointInDomain;
  string text_line, surfCp_filename;
  ifstream Surface_file;
  char buffer[50], cstr[200];


  nPointLocal = geometry->GetnPoint();
#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&nPointLocal, &nPointGlobal, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
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
        iPoint = geometry->node[geometry->vertex[iMarker][iVertex]->GetNode()]->GetGlobalIndex();
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

  if ((config->GetTime_Marching() && config->GetTime_Domain()) ||
      (config->GetTime_Marching() == HARMONIC_BALANCE)) {
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
  SU2_MPI::Allreduce(&MyPressDiff, &PressDiff, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif

  /*--- Update the total Cp difference coeffient ---*/

  solver_container->SetTotal_CpDiff(PressDiff);

  delete [] Point2Vertex;
  delete [] PointInDomain;

}

void COutputLegacy::SetHeatFlux_InverseDesign(CSolver *solver_container, CGeometry *geometry, CConfig *config, unsigned long iExtIter) {

  unsigned short iMarker, icommas, Boundary, iDim;
  unsigned long iVertex, iPoint, (*Point2Vertex)[2], nPointLocal = 0, nPointGlobal = 0;
  su2double XCoord, YCoord, ZCoord, PressureCoeff, HeatFlux = 0.0, HeatFluxDiff, HeatFluxTarget, *Normal = NULL, Area,
  Pressure, Cf;
  bool *PointInDomain;
  string text_line, surfHeatFlux_filename;
  ifstream Surface_file;
  char buffer[50], cstr[200];


  nPointLocal = geometry->GetnPoint();
#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&nPointLocal, &nPointGlobal, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
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
        iPoint = geometry->node[geometry->vertex[iMarker][iVertex]->GetNode()]->GetGlobalIndex();
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

  if ((config->GetTime_Marching() && config->GetTime_Domain()) ||
      (config->GetTime_Marching() == HARMONIC_BALANCE)) {
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
  SU2_MPI::Allreduce(&MyHeatFluxDiff, &HeatFluxDiff, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif

  /*--- Update the total HeatFlux difference coeffient ---*/

  solver_container->SetTotal_HeatFluxDiff(HeatFluxDiff);

  delete [] Point2Vertex;
  delete [] PointInDomain;

}

void COutputLegacy::SpecialOutput_SonicBoom(CSolver *solver, CGeometry *geometry, CConfig *config, bool output) {

  ofstream EquivArea_file, FuncGrad_file;
  unsigned short iMarker = 0, iDim;
  short *AzimuthalAngle = NULL;
  su2double Gamma, auxXCoord, auxYCoord, auxZCoord, InverseDesign = 0.0, DeltaX, Coord_i, Coord_j, jp1Coord, *Coord = NULL, MeanFuntion,
  *Face_Normal = NULL, auxArea, auxPress, Mach, Beta, R_Plane, Pressure_Inf,
  ModVelocity_Inf, Velocity_Inf[3], factor, *Xcoord = NULL, *Ycoord = NULL, *Zcoord = NULL,
  *Pressure = NULL, *FaceArea = NULL, *EquivArea = NULL, *TargetArea = NULL, *NearFieldWeight = NULL,
  *Weight = NULL, jFunction, jp1Function;
  unsigned long jVertex, iVertex, iPoint, nVertex_NearField = 0, auxPoint,
  *IdPoint = NULL, *IdDomain = NULL, auxDomain;
  unsigned short iPhiAngle;
  ofstream NearFieldEA_file; ifstream TargetEA_file;

  su2double XCoordBegin_OF = config->GetEA_IntLimit(0);
  su2double XCoordEnd_OF = config->GetEA_IntLimit(1);

  unsigned short nDim = geometry->GetnDim();
  su2double AoA = -(config->GetAoA()*PI_NUMBER/180.0);
  su2double EAScaleFactor = config->GetEA_ScaleFactor(); // The EA Obj. Func. should be ~ force based Obj. Func.

  Mach  = config->GetMach();
  Gamma = config->GetGamma();
  Beta = sqrt(Mach*Mach-1.0);
  R_Plane = fabs(config->GetEA_IntLimit(2));
  Pressure_Inf = config->GetPressure_FreeStreamND();
  Velocity_Inf[0] = config->GetVelocity_FreeStreamND()[0];
  Velocity_Inf[1] = config->GetVelocity_FreeStreamND()[1];
  Velocity_Inf[2] = config->GetVelocity_FreeStreamND()[2];
  ModVelocity_Inf = 0;
  for (iDim = 0; iDim < 3; iDim++)
    ModVelocity_Inf += Velocity_Inf[iDim] * Velocity_Inf[iDim];

  factor = 4.0*sqrt(2.0*Beta*R_Plane) / (Gamma*Pressure_Inf*Mach*Mach);

  if (rank == MASTER_NODE) cout << endl << "Writing Equivalent Area files.";

#ifndef HAVE_MPI

  /*--- Compute the total number of points on the near-field ---*/

  nVertex_NearField = 0;
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
    if (config->GetMarker_All_KindBC(iMarker) == NEARFIELD_BOUNDARY)
      for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        Face_Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
        Coord = geometry->node[iPoint]->GetCoord();

        /*--- Using Face_Normal(z), and Coord(z) we identify only a surface,
         note that there are 2 NEARFIELD_BOUNDARY surfaces ---*/

        if ((Face_Normal[nDim-1] > 0.0) && (Coord[nDim-1] < 0.0)) nVertex_NearField ++;
      }

  /*--- Create an array with all the coordinates, points, pressures, face area,
   equivalent area, and nearfield weight ---*/

  Xcoord = new su2double[nVertex_NearField];
  Ycoord = new su2double[nVertex_NearField];
  Zcoord = new su2double[nVertex_NearField];
  AzimuthalAngle = new short[nVertex_NearField];
  IdPoint = new unsigned long[nVertex_NearField];
  IdDomain = new unsigned long[nVertex_NearField];
  Pressure = new su2double[nVertex_NearField];
  FaceArea = new su2double[nVertex_NearField];
  EquivArea = new su2double[nVertex_NearField];
  TargetArea = new su2double[nVertex_NearField];
  NearFieldWeight = new su2double[nVertex_NearField];
  Weight = new su2double[nVertex_NearField];

  /*--- Copy the boundary information to an array ---*/

  nVertex_NearField = 0;
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
    if (config->GetMarker_All_KindBC(iMarker) == NEARFIELD_BOUNDARY)
      for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        Face_Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
        Coord = geometry->node[iPoint]->GetCoord();

        if ((Face_Normal[nDim-1] > 0.0) && (Coord[nDim-1] < 0.0)) {

          IdPoint[nVertex_NearField] = iPoint;
          Xcoord[nVertex_NearField] = geometry->node[iPoint]->GetCoord(0);
          Ycoord[nVertex_NearField] = geometry->node[iPoint]->GetCoord(1);

          if (nDim ==2) {
            AzimuthalAngle[nVertex_NearField] = 0;
          }

          if (nDim == 3) {
            Zcoord[nVertex_NearField] = geometry->node[iPoint]->GetCoord(2);

            /*--- Rotate the nearfield cylinder (AoA) only 3D ---*/

            su2double YcoordRot = Ycoord[nVertex_NearField];
            su2double ZcoordRot = Xcoord[nVertex_NearField]*sin(AoA) + Zcoord[nVertex_NearField]*cos(AoA);

            /*--- Compute the Azimuthal angle (resolution of degress in the Azimuthal angle)---*/

            su2double AngleDouble; short AngleInt;
            AngleDouble = fabs(atan(-YcoordRot/ZcoordRot)*180.0/PI_NUMBER);

            /*--- Fix an azimuthal line due to misalignments of the near-field ---*/

            su2double FixAzimuthalLine = config->GetFixAzimuthalLine();

            if ((AngleDouble >= FixAzimuthalLine - 0.1) && (AngleDouble <= FixAzimuthalLine + 0.1)) AngleDouble = FixAzimuthalLine - 0.1;

            AngleInt = SU2_TYPE::Short(floor(AngleDouble + 0.5));
            if (AngleInt >= 0) AzimuthalAngle[nVertex_NearField] = AngleInt;
            else AzimuthalAngle[nVertex_NearField] = 180 + AngleInt;
          }

          if (AzimuthalAngle[nVertex_NearField] <= 60) {
            Pressure[nVertex_NearField] = solver->GetNodes()->GetPressure(iPoint);
            FaceArea[nVertex_NearField] = fabs(Face_Normal[nDim-1]);
            nVertex_NearField ++;
          }

        }
      }

#else

  int nProcessor;
  SU2_MPI::Comm_size(MPI_COMM_WORLD, &nProcessor);

  unsigned long nLocalVertex_NearField = 0, MaxLocalVertex_NearField = 0;
  int iProcessor;

  unsigned long *Buffer_Receive_nVertex = NULL;
  if (rank == MASTER_NODE) {
    Buffer_Receive_nVertex = new unsigned long [nProcessor];
  }

  /*--- Compute the total number of points of the near-field ghost nodes ---*/

  nLocalVertex_NearField = 0;
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
    if (config->GetMarker_All_KindBC(iMarker) == NEARFIELD_BOUNDARY)
      for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        Face_Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
        Coord = geometry->node[iPoint]->GetCoord();

        if (geometry->node[iPoint]->GetDomain())
          if ((Face_Normal[nDim-1] > 0.0) && (Coord[nDim-1] < 0.0))
            nLocalVertex_NearField ++;
      }

  unsigned long *Buffer_Send_nVertex = new unsigned long [1];
  Buffer_Send_nVertex[0] = nLocalVertex_NearField;

  /*--- Send Near-Field vertex information --*/

  SU2_MPI::Allreduce(&nLocalVertex_NearField, &nVertex_NearField, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&nLocalVertex_NearField, &MaxLocalVertex_NearField, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
  SU2_MPI::Gather(Buffer_Send_nVertex, 1, MPI_UNSIGNED_LONG, Buffer_Receive_nVertex, 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
  delete [] Buffer_Send_nVertex;

  su2double *Buffer_Send_Xcoord = new su2double[MaxLocalVertex_NearField];
  su2double *Buffer_Send_Ycoord = new su2double[MaxLocalVertex_NearField];
  su2double *Buffer_Send_Zcoord = new su2double[MaxLocalVertex_NearField];
  unsigned long *Buffer_Send_IdPoint = new unsigned long [MaxLocalVertex_NearField];
  su2double *Buffer_Send_Pressure = new su2double [MaxLocalVertex_NearField];
  su2double *Buffer_Send_FaceArea = new su2double[MaxLocalVertex_NearField];

  su2double *Buffer_Receive_Xcoord = NULL;
  su2double *Buffer_Receive_Ycoord = NULL;
  su2double *Buffer_Receive_Zcoord = NULL;
  unsigned long *Buffer_Receive_IdPoint = NULL;
  su2double *Buffer_Receive_Pressure = NULL;
  su2double *Buffer_Receive_FaceArea = NULL;

  if (rank == MASTER_NODE) {
    Buffer_Receive_Xcoord = new su2double[nProcessor*MaxLocalVertex_NearField];
    Buffer_Receive_Ycoord = new su2double[nProcessor*MaxLocalVertex_NearField];
    Buffer_Receive_Zcoord = new su2double[nProcessor*MaxLocalVertex_NearField];
    Buffer_Receive_IdPoint = new unsigned long[nProcessor*MaxLocalVertex_NearField];
    Buffer_Receive_Pressure = new su2double[nProcessor*MaxLocalVertex_NearField];
    Buffer_Receive_FaceArea = new su2double[nProcessor*MaxLocalVertex_NearField];
  }

  unsigned long nBuffer_Xcoord = MaxLocalVertex_NearField;
  unsigned long nBuffer_Ycoord = MaxLocalVertex_NearField;
  unsigned long nBuffer_Zcoord = MaxLocalVertex_NearField;
  unsigned long nBuffer_IdPoint = MaxLocalVertex_NearField;
  unsigned long nBuffer_Pressure = MaxLocalVertex_NearField;
  unsigned long nBuffer_FaceArea = MaxLocalVertex_NearField;

  for (iVertex = 0; iVertex < MaxLocalVertex_NearField; iVertex++) {
    Buffer_Send_IdPoint[iVertex] = 0; Buffer_Send_Pressure[iVertex] = 0.0;
    Buffer_Send_FaceArea[iVertex] = 0.0; Buffer_Send_Xcoord[iVertex] = 0.0;
    Buffer_Send_Ycoord[iVertex] = 0.0; Buffer_Send_Zcoord[iVertex] = 0.0;
  }

  /*--- Copy coordinates, index points, and pressures to the auxiliar vector --*/

  nLocalVertex_NearField = 0;
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
    if (config->GetMarker_All_KindBC(iMarker) == NEARFIELD_BOUNDARY)
      for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        Face_Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
        Coord = geometry->node[iPoint]->GetCoord();

        if (geometry->node[iPoint]->GetDomain())
          if ((Face_Normal[nDim-1] > 0.0) && (Coord[nDim-1] < 0.0)) {
            Buffer_Send_IdPoint[nLocalVertex_NearField] = iPoint;
            Buffer_Send_Xcoord[nLocalVertex_NearField] = geometry->node[iPoint]->GetCoord(0);
            Buffer_Send_Ycoord[nLocalVertex_NearField] = geometry->node[iPoint]->GetCoord(1);
            Buffer_Send_Zcoord[nLocalVertex_NearField] = geometry->node[iPoint]->GetCoord(2);
            Buffer_Send_Pressure[nLocalVertex_NearField] = solver->GetNodes()->GetPressure(iPoint);
            Buffer_Send_FaceArea[nLocalVertex_NearField] = fabs(Face_Normal[nDim-1]);
            nLocalVertex_NearField++;
          }
      }

  /*--- Send all the information --*/

  SU2_MPI::Gather(Buffer_Send_Xcoord, nBuffer_Xcoord, MPI_DOUBLE, Buffer_Receive_Xcoord, nBuffer_Xcoord, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  SU2_MPI::Gather(Buffer_Send_Ycoord, nBuffer_Ycoord, MPI_DOUBLE, Buffer_Receive_Ycoord, nBuffer_Ycoord, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  SU2_MPI::Gather(Buffer_Send_Zcoord, nBuffer_Zcoord, MPI_DOUBLE, Buffer_Receive_Zcoord, nBuffer_Zcoord, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  SU2_MPI::Gather(Buffer_Send_IdPoint, nBuffer_IdPoint, MPI_UNSIGNED_LONG, Buffer_Receive_IdPoint, nBuffer_IdPoint, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
  SU2_MPI::Gather(Buffer_Send_Pressure, nBuffer_Pressure, MPI_DOUBLE, Buffer_Receive_Pressure, nBuffer_Pressure, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  SU2_MPI::Gather(Buffer_Send_FaceArea, nBuffer_FaceArea, MPI_DOUBLE, Buffer_Receive_FaceArea, nBuffer_FaceArea, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  delete [] Buffer_Send_Xcoord;
  delete [] Buffer_Send_Ycoord;
  delete [] Buffer_Send_Zcoord;
  delete [] Buffer_Send_IdPoint;
  delete [] Buffer_Send_Pressure;
  delete [] Buffer_Send_FaceArea;

  if (rank == MASTER_NODE) {

    Xcoord = new su2double[nVertex_NearField];
    Ycoord = new su2double[nVertex_NearField];
    Zcoord = new su2double[nVertex_NearField];
    AzimuthalAngle = new short[nVertex_NearField];
    IdPoint = new unsigned long[nVertex_NearField];
    IdDomain = new unsigned long[nVertex_NearField];
    Pressure = new su2double[nVertex_NearField];
    FaceArea = new su2double[nVertex_NearField];
    EquivArea = new su2double[nVertex_NearField];
    TargetArea = new su2double[nVertex_NearField];
    NearFieldWeight = new su2double[nVertex_NearField];
    Weight = new su2double[nVertex_NearField];

    nVertex_NearField = 0;
    for (iProcessor = 0; iProcessor < nProcessor; iProcessor++)
      for (iVertex = 0; iVertex < Buffer_Receive_nVertex[iProcessor]; iVertex++) {
        Xcoord[nVertex_NearField] = Buffer_Receive_Xcoord[iProcessor*MaxLocalVertex_NearField+iVertex];
        Ycoord[nVertex_NearField] = Buffer_Receive_Ycoord[iProcessor*MaxLocalVertex_NearField+iVertex];

        if (nDim == 2) {
          AzimuthalAngle[nVertex_NearField] = 0;
        }

        if (nDim == 3) {
          Zcoord[nVertex_NearField] = Buffer_Receive_Zcoord[iProcessor*MaxLocalVertex_NearField+iVertex];

          /*--- Rotate the nearfield cylinder  ---*/

          su2double YcoordRot = Ycoord[nVertex_NearField];
          su2double ZcoordRot = Xcoord[nVertex_NearField]*sin(AoA) + Zcoord[nVertex_NearField]*cos(AoA);

          /*--- Compute the Azimuthal angle ---*/

          su2double AngleDouble; short AngleInt;
          AngleDouble = fabs(atan(-YcoordRot/ZcoordRot)*180.0/PI_NUMBER);

          /*--- Fix an azimuthal line due to misalignments of the near-field ---*/

          su2double FixAzimuthalLine = config->GetFixAzimuthalLine();

          if ((AngleDouble >= FixAzimuthalLine - 0.1) && (AngleDouble <= FixAzimuthalLine + 0.1))
            AngleDouble = FixAzimuthalLine - 0.1;

          AngleInt = SU2_TYPE::Short(floor(AngleDouble + 0.5));

          if (AngleInt >= 0) AzimuthalAngle[nVertex_NearField] = AngleInt;
          else AzimuthalAngle[nVertex_NearField] = 180 + AngleInt;
        }

        if (AzimuthalAngle[nVertex_NearField] <= 60) {
          IdPoint[nVertex_NearField] = Buffer_Receive_IdPoint[iProcessor*MaxLocalVertex_NearField+iVertex];
          Pressure[nVertex_NearField] = Buffer_Receive_Pressure[iProcessor*MaxLocalVertex_NearField+iVertex];
          FaceArea[nVertex_NearField] = Buffer_Receive_FaceArea[iProcessor*MaxLocalVertex_NearField+iVertex];
          IdDomain[nVertex_NearField] = iProcessor;
          nVertex_NearField++;
        }

      }

    delete [] Buffer_Receive_nVertex;

    delete [] Buffer_Receive_Xcoord;
    delete [] Buffer_Receive_Ycoord;
    delete [] Buffer_Receive_Zcoord;
    delete [] Buffer_Receive_IdPoint;
    delete [] Buffer_Receive_Pressure;
    delete [] Buffer_Receive_FaceArea;

  }

#endif

  if (rank == MASTER_NODE) {

    vector<short> PhiAngleList;
    vector<short>::iterator IterPhiAngleList;

    for (iVertex = 0; iVertex < nVertex_NearField; iVertex++)
      PhiAngleList.push_back(AzimuthalAngle[iVertex]);

    sort( PhiAngleList.begin(), PhiAngleList.end());
    IterPhiAngleList = unique( PhiAngleList.begin(), PhiAngleList.end());
    PhiAngleList.resize( IterPhiAngleList - PhiAngleList.begin() );

    /*--- Create vectors and distribute the values among the different PhiAngle queues ---*/

    vector<vector<su2double> > Xcoord_PhiAngle; Xcoord_PhiAngle.resize(PhiAngleList.size());
    vector<vector<su2double> > Ycoord_PhiAngle; Ycoord_PhiAngle.resize(PhiAngleList.size());
    vector<vector<su2double> > Zcoord_PhiAngle; Zcoord_PhiAngle.resize(PhiAngleList.size());
    vector<vector<unsigned long> > IdPoint_PhiAngle; IdPoint_PhiAngle.resize(PhiAngleList.size());
    vector<vector<unsigned long> > IdDomain_PhiAngle; IdDomain_PhiAngle.resize(PhiAngleList.size());
    vector<vector<su2double> > Pressure_PhiAngle; Pressure_PhiAngle.resize(PhiAngleList.size());
    vector<vector<su2double> > FaceArea_PhiAngle; FaceArea_PhiAngle.resize(PhiAngleList.size());
    vector<vector<su2double> > EquivArea_PhiAngle; EquivArea_PhiAngle.resize(PhiAngleList.size());
    vector<vector<su2double> > TargetArea_PhiAngle; TargetArea_PhiAngle.resize(PhiAngleList.size());
    vector<vector<su2double> > NearFieldWeight_PhiAngle; NearFieldWeight_PhiAngle.resize(PhiAngleList.size());
    vector<vector<su2double> > Weight_PhiAngle; Weight_PhiAngle.resize(PhiAngleList.size());

    /*--- Distribute the values among the different PhiAngles ---*/

    for (iVertex = 0; iVertex < nVertex_NearField; iVertex++)
      for (iPhiAngle = 0; iPhiAngle < PhiAngleList.size(); iPhiAngle++)
        if (AzimuthalAngle[iVertex] == PhiAngleList[iPhiAngle]) {
          Xcoord_PhiAngle[iPhiAngle].push_back(Xcoord[iVertex]);
          Ycoord_PhiAngle[iPhiAngle].push_back(Ycoord[iVertex]);
          Zcoord_PhiAngle[iPhiAngle].push_back(Zcoord[iVertex]);
          IdPoint_PhiAngle[iPhiAngle].push_back(IdPoint[iVertex]);
          IdDomain_PhiAngle[iPhiAngle].push_back(IdDomain[iVertex]);
          Pressure_PhiAngle[iPhiAngle].push_back(Pressure[iVertex]);
          FaceArea_PhiAngle[iPhiAngle].push_back(FaceArea[iVertex]);
          EquivArea_PhiAngle[iPhiAngle].push_back(EquivArea[iVertex]);
          TargetArea_PhiAngle[iPhiAngle].push_back(TargetArea[iVertex]);
          NearFieldWeight_PhiAngle[iPhiAngle].push_back(NearFieldWeight[iVertex]);
          Weight_PhiAngle[iPhiAngle].push_back(Weight[iVertex]);
        }

    /*--- Order the arrays (x Coordinate, Pressure, Point, and Domain) ---*/

    for (iPhiAngle = 0; iPhiAngle < PhiAngleList.size(); iPhiAngle++)
      for (iVertex = 0; iVertex < Xcoord_PhiAngle[iPhiAngle].size(); iVertex++)
        for (jVertex = 0; jVertex < Xcoord_PhiAngle[iPhiAngle].size() - 1 - iVertex; jVertex++)
          if (Xcoord_PhiAngle[iPhiAngle][jVertex] > Xcoord_PhiAngle[iPhiAngle][jVertex+1]) {
            auxXCoord = Xcoord_PhiAngle[iPhiAngle][jVertex]; Xcoord_PhiAngle[iPhiAngle][jVertex] = Xcoord_PhiAngle[iPhiAngle][jVertex+1]; Xcoord_PhiAngle[iPhiAngle][jVertex+1] = auxXCoord;
            auxYCoord = Ycoord_PhiAngle[iPhiAngle][jVertex]; Ycoord_PhiAngle[iPhiAngle][jVertex] = Ycoord_PhiAngle[iPhiAngle][jVertex+1]; Ycoord_PhiAngle[iPhiAngle][jVertex+1] = auxYCoord;
            auxZCoord = Zcoord_PhiAngle[iPhiAngle][jVertex]; Zcoord_PhiAngle[iPhiAngle][jVertex] = Zcoord_PhiAngle[iPhiAngle][jVertex+1]; Zcoord_PhiAngle[iPhiAngle][jVertex+1] = auxZCoord;
            auxPress = Pressure_PhiAngle[iPhiAngle][jVertex]; Pressure_PhiAngle[iPhiAngle][jVertex] = Pressure_PhiAngle[iPhiAngle][jVertex+1]; Pressure_PhiAngle[iPhiAngle][jVertex+1] = auxPress;
            auxArea = FaceArea_PhiAngle[iPhiAngle][jVertex]; FaceArea_PhiAngle[iPhiAngle][jVertex] = FaceArea_PhiAngle[iPhiAngle][jVertex+1]; FaceArea_PhiAngle[iPhiAngle][jVertex+1] = auxArea;
            auxPoint = IdPoint_PhiAngle[iPhiAngle][jVertex]; IdPoint_PhiAngle[iPhiAngle][jVertex] = IdPoint_PhiAngle[iPhiAngle][jVertex+1]; IdPoint_PhiAngle[iPhiAngle][jVertex+1] = auxPoint;
            auxDomain = IdDomain_PhiAngle[iPhiAngle][jVertex]; IdDomain_PhiAngle[iPhiAngle][jVertex] = IdDomain_PhiAngle[iPhiAngle][jVertex+1]; IdDomain_PhiAngle[iPhiAngle][jVertex+1] = auxDomain;
          }


    /*--- Check that all the azimuth lists have the same size ---*/

    unsigned short nVertex = Xcoord_PhiAngle[0].size();
    for (iPhiAngle = 0; iPhiAngle < PhiAngleList.size(); iPhiAngle++) {
      unsigned short nVertex_aux = Xcoord_PhiAngle[iPhiAngle].size();
      if (nVertex_aux != nVertex) cout <<"Be careful!!! one azimuth list is shorter than the other"<< endl;
      nVertex = min(nVertex, nVertex_aux);
    }

    /*--- Compute equivalent area distribution at each azimuth angle ---*/

    for (iPhiAngle = 0; iPhiAngle < PhiAngleList.size(); iPhiAngle++) {
      EquivArea_PhiAngle[iPhiAngle][0] = 0.0;
      for (iVertex = 1; iVertex < EquivArea_PhiAngle[iPhiAngle].size(); iVertex++) {
        EquivArea_PhiAngle[iPhiAngle][iVertex] = 0.0;

        Coord_i = Xcoord_PhiAngle[iPhiAngle][iVertex]*cos(AoA) - Zcoord_PhiAngle[iPhiAngle][iVertex]*sin(AoA);

        for (jVertex = 0; jVertex < iVertex-1; jVertex++) {

          Coord_j = Xcoord_PhiAngle[iPhiAngle][jVertex]*cos(AoA) - Zcoord_PhiAngle[iPhiAngle][jVertex]*sin(AoA);
          jp1Coord = Xcoord_PhiAngle[iPhiAngle][jVertex+1]*cos(AoA) - Zcoord_PhiAngle[iPhiAngle][jVertex+1]*sin(AoA);

          jFunction = factor*(Pressure_PhiAngle[iPhiAngle][jVertex] - Pressure_Inf)*sqrt(Coord_i-Coord_j);
          jp1Function = factor*(Pressure_PhiAngle[iPhiAngle][jVertex+1] - Pressure_Inf)*sqrt(Coord_i-jp1Coord);

          DeltaX = (jp1Coord-Coord_j);
          MeanFuntion = 0.5*(jp1Function + jFunction);
          EquivArea_PhiAngle[iPhiAngle][iVertex] += DeltaX * MeanFuntion;
        }
      }
    }

    /*--- Create a file with the equivalent area distribution at each azimuthal angle ---*/

    NearFieldEA_file.precision(15);

    if (output) {

      NearFieldEA_file.open("Equivalent_Area.dat", ios::out);
      NearFieldEA_file << "TITLE = \"Equivalent Area evaluation at each azimuthal angle\"" << "\n";

      if (config->GetSystemMeasurements() == US)
        NearFieldEA_file << "VARIABLES = \"Height (in) at r="<< R_Plane*12.0 << " in. (cyl. coord. system)\"";
      else
        NearFieldEA_file << "VARIABLES = \"Height (m) at r="<< R_Plane << " m. (cylindrical coordinate system)\"";

      for (iPhiAngle = 0; iPhiAngle < PhiAngleList.size(); iPhiAngle++) {
        if (config->GetSystemMeasurements() == US)
          NearFieldEA_file << ", \"Equivalent Area (ft<sup>2</sup>), <greek>F</greek>= " << PhiAngleList[iPhiAngle] << " deg.\"";
        else
          NearFieldEA_file << ", \"Equivalent Area (m<sup>2</sup>), <greek>F</greek>= " << PhiAngleList[iPhiAngle] << " deg.\"";
      }

      NearFieldEA_file << "\n";
      for (iVertex = 0; iVertex < EquivArea_PhiAngle[0].size(); iVertex++) {

        su2double XcoordRot = Xcoord_PhiAngle[0][iVertex]*cos(AoA) - Zcoord_PhiAngle[0][iVertex]*sin(AoA);
        su2double XcoordRot_init = Xcoord_PhiAngle[0][0]*cos(AoA) - Zcoord_PhiAngle[0][0]*sin(AoA);

        if (config->GetSystemMeasurements() == US)
          NearFieldEA_file << scientific << (XcoordRot - XcoordRot_init) * 12.0;
        else
          NearFieldEA_file << scientific << (XcoordRot - XcoordRot_init);

        for (iPhiAngle = 0; iPhiAngle < PhiAngleList.size(); iPhiAngle++) {
          NearFieldEA_file << scientific << ", " << EquivArea_PhiAngle[iPhiAngle][iVertex];
        }

        NearFieldEA_file << "\n";

      }
      NearFieldEA_file.close();

    }

    /*--- Read target equivalent area from the configuration file,
     this first implementation requires a complete table (same as the original
     EA table). so... no interpolation. ---*/

    vector<vector<su2double> > TargetArea_PhiAngle_Trans;
    TargetEA_file.open("TargetEA.dat", ios::in);

    if (TargetEA_file.fail()) {
      /*--- Set the table to 0 ---*/
      for (iPhiAngle = 0; iPhiAngle < PhiAngleList.size(); iPhiAngle++)
        for (iVertex = 0; iVertex < TargetArea_PhiAngle[iPhiAngle].size(); iVertex++)
          TargetArea_PhiAngle[iPhiAngle][iVertex] = 0.0;
    }
    else {

      /*--- skip header lines ---*/

      string line;
      getline(TargetEA_file, line);
      getline(TargetEA_file, line);

      while (TargetEA_file) {

        string line;
        getline(TargetEA_file, line);
        istringstream is(line);
        vector<su2double> row;
        unsigned short iter = 0;

        while (is.good()) {
          string token;
          getline(is, token,',');

          istringstream js(token);

          su2double data;
          js >> data;

          /*--- The first element in the table is the coordinate (in or m)---*/

          if (iter != 0) row.push_back(data);
          iter++;

        }
        TargetArea_PhiAngle_Trans.push_back(row);
      }

      for (iPhiAngle = 0; iPhiAngle < PhiAngleList.size(); iPhiAngle++)
        for (iVertex = 0; iVertex < EquivArea_PhiAngle[iPhiAngle].size(); iVertex++)
          TargetArea_PhiAngle[iPhiAngle][iVertex] = TargetArea_PhiAngle_Trans[iVertex][iPhiAngle];

    }

    /*--- Divide by the number of Phi angles in the nearfield ---*/

    su2double PhiFactor = 1.0/su2double(PhiAngleList.size());

    /*--- Evaluate the objective function ---*/

    InverseDesign = 0;
    for (iPhiAngle = 0; iPhiAngle < PhiAngleList.size(); iPhiAngle++)
      for (iVertex = 0; iVertex < EquivArea_PhiAngle[iPhiAngle].size(); iVertex++) {
        Weight_PhiAngle[iPhiAngle][iVertex] = 1.0;
        Coord_i = Xcoord_PhiAngle[iPhiAngle][iVertex];

        su2double Difference = EquivArea_PhiAngle[iPhiAngle][iVertex]-TargetArea_PhiAngle[iPhiAngle][iVertex];
        su2double percentage = fabs(Difference)*100/fabs(TargetArea_PhiAngle[iPhiAngle][iVertex]);

        if ((percentage < 0.1) || (Coord_i < XCoordBegin_OF) || (Coord_i > XCoordEnd_OF)) Difference = 0.0;

        InverseDesign += EAScaleFactor*PhiFactor*Weight_PhiAngle[iPhiAngle][iVertex]*Difference*Difference;

      }

    /*--- Evaluate the weight of the nearfield pressure (adjoint input) ---*/

    for (iPhiAngle = 0; iPhiAngle < PhiAngleList.size(); iPhiAngle++)
      for (iVertex = 0; iVertex < EquivArea_PhiAngle[iPhiAngle].size(); iVertex++) {
        Coord_i = Xcoord_PhiAngle[iPhiAngle][iVertex];
        NearFieldWeight_PhiAngle[iPhiAngle][iVertex] = 0.0;
        for (jVertex = iVertex; jVertex < EquivArea_PhiAngle[iPhiAngle].size(); jVertex++) {
          Coord_j = Xcoord_PhiAngle[iPhiAngle][jVertex];
          Weight_PhiAngle[iPhiAngle][iVertex] = 1.0;

          su2double Difference = EquivArea_PhiAngle[iPhiAngle][jVertex]-TargetArea_PhiAngle[iPhiAngle][jVertex];
          su2double percentage = fabs(Difference)*100/fabs(TargetArea_PhiAngle[iPhiAngle][jVertex]);

          if ((percentage < 0.1) || (Coord_j < XCoordBegin_OF) || (Coord_j > XCoordEnd_OF)) Difference = 0.0;

          NearFieldWeight_PhiAngle[iPhiAngle][iVertex] += EAScaleFactor*PhiFactor*Weight_PhiAngle[iPhiAngle][iVertex]*2.0*Difference*factor*sqrt(Coord_j-Coord_i);
        }
      }

    /*--- Write the Nearfield pressure at each Azimuthal PhiAngle ---*/

    EquivArea_file.precision(15);

    if (output) {

      EquivArea_file.open("nearfield_flow.dat", ios::out);
      EquivArea_file << "TITLE = \"Equivalent Area evaluation at each azimuthal angle\"" << "\n";

      if (config->GetSystemMeasurements() == US)
        EquivArea_file << "VARIABLES = \"Height (in) at r="<< R_Plane*12.0 << " in. (cyl. coord. system)\",\"Equivalent Area (ft<sup>2</sup>)\",\"Target Equivalent Area (ft<sup>2</sup>)\",\"Cp\"" << "\n";
      else
        EquivArea_file << "VARIABLES = \"Height (m) at r="<< R_Plane << " m. (cylindrical coordinate system)\",\"Equivalent Area (m<sup>2</sup>)\",\"Target Equivalent Area (m<sup>2</sup>)\",\"Cp\"" << "\n";

      for (iPhiAngle = 0; iPhiAngle < PhiAngleList.size(); iPhiAngle++) {
        EquivArea_file << fixed << "ZONE T= \"<greek>F</greek>=" << PhiAngleList[iPhiAngle] << " deg.\"" << "\n";
        for (iVertex = 0; iVertex < Xcoord_PhiAngle[iPhiAngle].size(); iVertex++) {

          su2double XcoordRot = Xcoord_PhiAngle[0][iVertex]*cos(AoA) - Zcoord_PhiAngle[0][iVertex]*sin(AoA);
          su2double XcoordRot_init = Xcoord_PhiAngle[0][0]*cos(AoA) - Zcoord_PhiAngle[0][0]*sin(AoA);

          if (config->GetSystemMeasurements() == US)
            EquivArea_file << scientific << (XcoordRot - XcoordRot_init) * 12.0;
          else
            EquivArea_file << scientific << (XcoordRot - XcoordRot_init);

          EquivArea_file << scientific << ", " << EquivArea_PhiAngle[iPhiAngle][iVertex]
          << ", " << TargetArea_PhiAngle[iPhiAngle][iVertex] << ", " << (Pressure_PhiAngle[iPhiAngle][iVertex]-Pressure_Inf)/Pressure_Inf << "\n";
        }
      }

      EquivArea_file.close();

    }

    /*--- Write Weight file for adjoint computation ---*/

    FuncGrad_file.precision(15);

    if (output) {

      FuncGrad_file.open("WeightNF.dat", ios::out);

      FuncGrad_file << scientific << "-1.0";
      for (iPhiAngle = 0; iPhiAngle < PhiAngleList.size(); iPhiAngle++)
        FuncGrad_file << scientific << "\t" << PhiAngleList[iPhiAngle];
      FuncGrad_file << "\n";

      for (iVertex = 0; iVertex < NearFieldWeight_PhiAngle[0].size(); iVertex++) {
        su2double XcoordRot = Xcoord_PhiAngle[0][iVertex]*cos(AoA) - Zcoord_PhiAngle[0][iVertex]*sin(AoA);
        FuncGrad_file << scientific << XcoordRot;
        for (iPhiAngle = 0; iPhiAngle < PhiAngleList.size(); iPhiAngle++)
          FuncGrad_file << scientific << "\t" << NearFieldWeight_PhiAngle[iPhiAngle][iVertex];
        FuncGrad_file << "\n";
      }
      FuncGrad_file.close();

    }

    /*--- Delete structures ---*/

    delete [] Xcoord; delete [] Ycoord; delete [] Zcoord;
    delete [] AzimuthalAngle; delete [] IdPoint; delete [] IdDomain;
    delete [] Pressure; delete [] FaceArea;
    delete [] EquivArea; delete [] TargetArea;
    delete [] NearFieldWeight; delete [] Weight;

  }

#ifndef HAVE_MPI

  /*--- Store the value of the NearField coefficient ---*/

  solver->SetTotal_CEquivArea(InverseDesign);

#else

  /*--- Send the value of the NearField coefficient to all the processors ---*/

  SU2_MPI::Bcast(&InverseDesign, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);

  /*--- Store the value of the NearField coefficient ---*/

  solver->SetTotal_CEquivArea(InverseDesign);

#endif

}

void COutputLegacy::SpecialOutput_Distortion(CSolver *solver, CGeometry *geometry, CConfig *config, bool output) {

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
  unsigned long Buffer_Send_nVertex[1], *Buffer_Recv_nVertex = NULL;
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

    if ((config->GetTabular_FileFormat() == TAB_CSV)) strcpy (cstr, "surface_analysis.vtk");
    else strcpy (cstr, "surface_analysis.dat");

    SurfFlow_file.precision(15);

    SurfFlow_file.open(cstr, ios::out);

    if ((config->GetTabular_FileFormat() == TAB_CSV)) {
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
          if (geometry->node[iPoint]->GetDomain()) nLocalVertex_Surface++;
        }
      }
    }

    /*--- Communicate the number of local vertices on each partition
     to the master node ---*/

    Buffer_Send_nVertex[0] = nLocalVertex_Surface;
    if (rank == MASTER_NODE) Buffer_Recv_nVertex = new unsigned long [nProcessor];

#ifdef HAVE_MPI
    SU2_MPI::Allreduce(&nLocalVertex_Surface, &MaxLocalVertex_Surface, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
    SU2_MPI::Gather(&Buffer_Send_nVertex, 1, MPI_UNSIGNED_LONG, Buffer_Recv_nVertex, 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
#else
    MaxLocalVertex_Surface = nLocalVertex_Surface;
    Buffer_Recv_nVertex[MASTER_NODE] = Buffer_Send_nVertex[MASTER_NODE];
#endif

    /*--- Send and Recv buffers ---*/

    su2double *Buffer_Send_Coord_x = NULL, *Buffer_Recv_Coord_x = NULL;
    Buffer_Send_Coord_x = new su2double [MaxLocalVertex_Surface];

    su2double *Buffer_Send_Coord_y = NULL, *Buffer_Recv_Coord_y = NULL;
    Buffer_Send_Coord_y = new su2double [MaxLocalVertex_Surface];

    su2double *Buffer_Send_Coord_z = NULL, *Buffer_Recv_Coord_z = NULL;
    if (nDim == 3)  Buffer_Send_Coord_z = new su2double [MaxLocalVertex_Surface];

    su2double *Buffer_Send_PT = NULL, *Buffer_Recv_PT = NULL;
    Buffer_Send_PT = new su2double [MaxLocalVertex_Surface];

    su2double *Buffer_Send_TT = NULL, *Buffer_Recv_TT = NULL;
    Buffer_Send_TT = new su2double [MaxLocalVertex_Surface];

    su2double *Buffer_Send_P = NULL, *Buffer_Recv_P = NULL;
    Buffer_Send_P = new su2double [MaxLocalVertex_Surface];

    su2double *Buffer_Send_T = NULL, *Buffer_Recv_T = NULL;
    Buffer_Send_T = new su2double [MaxLocalVertex_Surface];

    su2double *Buffer_Send_Mach = NULL, *Buffer_Recv_Mach = NULL;
    Buffer_Send_Mach = new su2double [MaxLocalVertex_Surface];

    su2double *Buffer_Send_Vel_x = NULL, *Buffer_Recv_Vel_x = NULL;
    Buffer_Send_Vel_x = new su2double [MaxLocalVertex_Surface];

    su2double *Buffer_Send_Vel_y = NULL, *Buffer_Recv_Vel_y = NULL;
    Buffer_Send_Vel_y = new su2double [MaxLocalVertex_Surface];

    su2double *Buffer_Send_Vel_z = NULL, *Buffer_Recv_Vel_z = NULL;
    if (nDim == 3) Buffer_Send_Vel_z = new su2double [MaxLocalVertex_Surface];

    su2double *Buffer_Send_q = NULL, *Buffer_Recv_q = NULL;
    Buffer_Send_q = new su2double [MaxLocalVertex_Surface];

    su2double *Buffer_Send_Area = NULL, *Buffer_Recv_Area = NULL;
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

          if (geometry->node[iPoint]->GetDomain()) {

            Buffer_Send_Coord_x[nVertex_Surface] = geometry->node[iPoint]->GetCoord(0);
            Buffer_Send_Coord_y[nVertex_Surface] = geometry->node[iPoint]->GetCoord(1);
            if (nDim == 3) { Buffer_Send_Coord_z[nVertex_Surface] = geometry->node[iPoint]->GetCoord(2); }

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

    SU2_MPI::Gather(Buffer_Send_Coord_x, MaxLocalVertex_Surface, MPI_DOUBLE, Buffer_Recv_Coord_x, MaxLocalVertex_Surface, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    SU2_MPI::Gather(Buffer_Send_Coord_y, MaxLocalVertex_Surface, MPI_DOUBLE, Buffer_Recv_Coord_y, MaxLocalVertex_Surface, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    if (nDim == 3) SU2_MPI::Gather(Buffer_Send_Coord_z, MaxLocalVertex_Surface, MPI_DOUBLE, Buffer_Recv_Coord_z, MaxLocalVertex_Surface, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    SU2_MPI::Gather(Buffer_Send_PT, MaxLocalVertex_Surface, MPI_DOUBLE, Buffer_Recv_PT, MaxLocalVertex_Surface, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    SU2_MPI::Gather(Buffer_Send_TT, MaxLocalVertex_Surface, MPI_DOUBLE, Buffer_Recv_TT, MaxLocalVertex_Surface, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    SU2_MPI::Gather(Buffer_Send_P, MaxLocalVertex_Surface, MPI_DOUBLE, Buffer_Recv_P, MaxLocalVertex_Surface, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    SU2_MPI::Gather(Buffer_Send_T, MaxLocalVertex_Surface, MPI_DOUBLE, Buffer_Recv_T, MaxLocalVertex_Surface, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    SU2_MPI::Gather(Buffer_Send_Mach, MaxLocalVertex_Surface, MPI_DOUBLE, Buffer_Recv_Mach, MaxLocalVertex_Surface, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    SU2_MPI::Gather(Buffer_Send_Vel_x, MaxLocalVertex_Surface, MPI_DOUBLE, Buffer_Recv_Vel_x, MaxLocalVertex_Surface, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    SU2_MPI::Gather(Buffer_Send_Vel_y, MaxLocalVertex_Surface, MPI_DOUBLE, Buffer_Recv_Vel_y, MaxLocalVertex_Surface, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    if (nDim == 3) SU2_MPI::Gather(Buffer_Send_Vel_z, MaxLocalVertex_Surface, MPI_DOUBLE, Buffer_Recv_Vel_z, MaxLocalVertex_Surface, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    SU2_MPI::Gather(Buffer_Send_q, MaxLocalVertex_Surface, MPI_DOUBLE, Buffer_Recv_q, MaxLocalVertex_Surface, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    SU2_MPI::Gather(Buffer_Send_Area, MaxLocalVertex_Surface, MPI_DOUBLE, Buffer_Recv_Area, MaxLocalVertex_Surface, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);

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

        if (config->GetTabular_FileFormat() == TAB_CSV) {

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

void COutputLegacy::SetSensitivity_Files(CGeometry ***geometry, CConfig **config, unsigned short val_nZone) {

  unsigned short iMarker,iDim, nDim, iVar, nMarker, nVar;
  unsigned long iVertex, iPoint, nPoint, nVertex;
  su2double *Normal, Prod, Sens = 0.0, SensDim, Area;

  unsigned short iZone;

  CSolver ***solver = new CSolver**[val_nZone];
  for (iZone = 0; iZone < val_nZone; iZone++) {
    solver[iZone] = new CSolver*[1];
  }

  for (iZone = 0; iZone < val_nZone; iZone++) {

    nPoint = geometry[iZone][INST_0]->GetnPoint();
    nDim   = geometry[iZone][INST_0]->GetnDim();
    nMarker = config[iZone]->GetnMarker_All();
    nVar = nDim + 1;

    /*--- We create a baseline solver to easily merge the sensitivity information ---*/

    vector<string> fieldnames;
    fieldnames.push_back("\"Point\"");
    fieldnames.push_back("\"x\"");
    fieldnames.push_back("\"y\"");
    if (nDim == 3) {
      fieldnames.push_back("\"z\"");
    }
    fieldnames.push_back("\"Sensitivity_x\"");
    fieldnames.push_back("\"Sensitivity_y\"");
    if (nDim == 3) {
      fieldnames.push_back("\"Sensitivity_z\"");
    }
    fieldnames.push_back("\"Surface_Sensitivity\"");

    solver[iZone][INST_0] = new CBaselineSolver(geometry[iZone][INST_0], config[iZone], nVar+nDim, fieldnames);

    for (iPoint = 0; iPoint < nPoint; iPoint++) {
      for (iDim = 0; iDim < nDim; iDim++) {
        solver[iZone][INST_0]->GetNodes()->SetSolution(iPoint, iDim, geometry[iZone][INST_0]->node[iPoint]->GetCoord(iDim));
      }
      for (iVar = 0; iVar < nDim; iVar++) {
        solver[iZone][INST_0]->GetNodes()->SetSolution(iPoint, iVar+nDim, geometry[iZone][INST_0]->GetSensitivity(iPoint, iVar));
      }
    }

    /*--- Compute the sensitivity in normal direction ---*/

    for (iMarker = 0; iMarker < nMarker; iMarker++) {

      if((config[iZone]->GetMarker_All_KindBC(iMarker) == HEAT_FLUX ) ||
         (config[iZone]->GetMarker_All_KindBC(iMarker) == EULER_WALL ) ||
         (config[iZone]->GetMarker_All_KindBC(iMarker) == ISOTHERMAL ) ||
         (config[iZone]->GetMarker_All_KindBC(iMarker) == CHT_WALL_INTERFACE )) {


        nVertex = geometry[iZone][INST_0]->GetnVertex(iMarker);

        for (iVertex = 0; iVertex < nVertex; iVertex++) {
          iPoint = geometry[iZone][INST_0]->vertex[iMarker][iVertex]->GetNode();
          Normal = geometry[iZone][INST_0]->vertex[iMarker][iVertex]->GetNormal();
          Prod = 0.0;
          Area = 0.0;
          for (iDim = 0; iDim < nDim; iDim++) {

            /*--- Retrieve the gradient calculated with discrete adjoint method ---*/

            SensDim = geometry[iZone][INST_0]->GetSensitivity(iPoint, iDim);

            /*--- Calculate scalar product for projection onto the normal vector ---*/

            Prod += Normal[iDim]*SensDim;

            Area += Normal[iDim]*Normal[iDim];
          }

          Area = sqrt(Area);

          /*--- Projection of the gradient onto the normal vector of the surface ---*/

          Sens = Prod/Area;

          solver[iZone][INST_0]->GetNodes()->SetSolution(iPoint, 2*nDim, Sens);

        }
      }
    }
  }

  /*--- Merge the information and write the output files ---*/

  SetBaselineResult_Files(solver, geometry, config, 0, val_nZone);

  for (iZone = 0; iZone < val_nZone; iZone++) {
    delete solver[iZone][0];
    delete solver[iZone];
  }
  delete [] solver;
}

void COutputLegacy::WriteTurboPerfConvHistory(CConfig *config){

  unsigned short iMarker_Monitoring;
  string inMarker_Tag, outMarker_Tag, inMarkerTag_Mix;
  unsigned short nZone       = config->GetnZone();
  bool turbulent = ((config->GetKind_Solver() == RANS) || (config->GetKind_Solver() == DISC_ADJ_RANS));
  bool menter_sst = (config->GetKind_Turb_Model() == SST) || (config->GetKind_Turb_Model() == SST_SUST);

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

void COutputLegacy::SpecialOutput_HarmonicBalance(CSolver *****solver, CGeometry ****geometry, CConfig **config, unsigned short iInst, unsigned short val_nInst, bool output) {

  /*--- Write file with flow quantities for harmonic balance HB ---*/
  ofstream HB_output_file;
  ofstream mean_HB_file;

  /*--- MPI Send/Recv buffers ---*/
  su2double *sbuf_var = NULL,  *rbuf_var = NULL;

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

void COutputLegacy::SetResult_Files_Parallel(CSolver *****solver_container,
                                       CGeometry ****geometry,
                                       CConfig **config,
                                       unsigned long iExtIter,
                                       unsigned short val_nZone) {

  unsigned short iZone, iVar, iInst;
  unsigned long iPoint;
  unsigned short nInst = 1;
  bool compressible = true;

  for (iZone = 0; iZone < val_nZone; iZone++) {

    /*--- Bool to distinguish between the FVM and FEM solvers. ---*/
    const bool fem_solver = config[iZone]->GetFEMSolver();

    /*--- Get the file output format ---*/

    unsigned short FileFormat = config[iZone]->GetTabular_FileFormat();
    nInst = config[iZone]->GetnTimeInstances();

    for (iInst = 0; iInst < nInst; iInst++){

      bool cont_adj = config[iZone]->GetContinuous_Adjoint();
      bool disc_adj = config[iZone]->GetDiscrete_Adjoint();

      /*--- Flags identifying the types of files to be written. ---*/
      /*--- For now, we are disabling the parallel writers for Tecplot
          ASCII until we have parallel versions of all file formats
          available. SU2_SOL will remain intact for writing files
          until this capability is completed. ---*/

      bool Wrt_Vol = config[iZone]->GetWrt_Vol_Sol();
      bool Wrt_Srf = config[iZone]->GetWrt_Srf_Sol();
      bool Wrt_Csv = config[iZone]->GetWrt_Csv_Sol();

#ifdef HAVE_MPI
      /*--- Do not merge the connectivity or write the visualization files
     if we are running in parallel, unless we are using ParaView binary.
       Force the use of SU2_SOL to merge and write the viz. files in this
       case to save overhead. ---*/

      if ((size > SINGLE_NODE) && (FileFormat != PARAVIEW_BINARY)) {
        Wrt_Vol = false;
        Wrt_Srf = false;
      }
#endif

    /*--- Check for compressible/incompressible flow problems. ---*/

    compressible = (config[iZone]->GetKind_Regime() == COMPRESSIBLE);

    /*--- First, prepare the offsets needed throughout below. ---*/

    PrepareOffsets(config[iZone], geometry[iZone][iInst][MESH_0]);

    /*--- Write out CSV files in parallel for flow and adjoint. ---*/

    if (rank == MASTER_NODE) cout << endl << "Writing comma-separated values (CSV) surface files." << endl;

    switch (config[iZone]->GetKind_Solver()) {
      case EULER : case NAVIER_STOKES : case RANS :
      case INC_EULER : case INC_NAVIER_STOKES : case INC_RANS :
        if (Wrt_Csv) SetSurfaceCSV_Flow(config[iZone], geometry[iZone][iInst][MESH_0],
            solver_container[iZone][iInst][MESH_0][FLOW_SOL], iExtIter, iZone, iInst);
        break;
      case ADJ_EULER : case ADJ_NAVIER_STOKES : case ADJ_RANS :
      case DISC_ADJ_EULER: case DISC_ADJ_NAVIER_STOKES: case DISC_ADJ_RANS:
      case DISC_ADJ_INC_EULER: case DISC_ADJ_INC_NAVIER_STOKES: case DISC_ADJ_INC_RANS:
        if (Wrt_Csv) SetSurfaceCSV_Adjoint(config[iZone], geometry[iZone][iInst][MESH_0],
            solver_container[iZone][iInst][MESH_0][ADJFLOW_SOL],
            solver_container[iZone][iInst][MESH_0][FLOW_SOL], iExtIter, iZone, iInst);
        break;
      default: break;
    }

    /*--- This switch statement will become a call to a virtual function
     defined within each of the "physics" output child classes that loads
     the local data for that particular problem alone. ---*/

      if (rank == MASTER_NODE)
        cout << "Loading solution output data locally on each rank." << endl;

      switch (config[iZone]->GetKind_Solver()) {
        case EULER : case NAVIER_STOKES : case RANS :
        case INC_EULER : case INC_NAVIER_STOKES : case INC_RANS :
        if (compressible)
          LoadLocalData_Flow(config[iZone], geometry[iZone][iInst][MESH_0], solver_container[iZone][iInst][MESH_0], iZone);
        else
          LoadLocalData_IncFlow(config[iZone], geometry[iZone][iInst][MESH_0], solver_container[iZone][iInst][MESH_0], iZone);
        break;
        case ADJ_EULER : case ADJ_NAVIER_STOKES : case ADJ_RANS :
        case DISC_ADJ_EULER: case DISC_ADJ_NAVIER_STOKES: case DISC_ADJ_RANS:
        case DISC_ADJ_INC_EULER: case DISC_ADJ_INC_NAVIER_STOKES: case DISC_ADJ_INC_RANS:
        LoadLocalData_AdjFlow(config[iZone], geometry[iZone][iInst][MESH_0], solver_container[iZone][iInst][MESH_0], iZone);
        break;
      case FEM_ELASTICITY: case DISC_ADJ_FEM:
        LoadLocalData_Elasticity(config[iZone], geometry[iZone][iInst][MESH_0], solver_container[iZone][iInst][MESH_0], iZone);
        break;
      case HEAT_EQUATION:
        LoadLocalData_Base(config[iZone], geometry[iZone][iInst][MESH_0], solver_container[iZone][iInst][MESH_0], iZone);
        break;
      case FEM_EULER: case FEM_NAVIER_STOKES: case FEM_RANS: case FEM_LES:
        LoadLocalData_FEM(config[iZone], geometry[iZone][iInst][MESH_0], solver_container[iZone][iInst][MESH_0], iZone);
      default: break;
    }

    /*--- Store the solution to be used on the final iteration with cte. lift mode. ---*/

    if ((!cont_adj) && (!disc_adj) && (config[iZone]->GetFixed_CL_Mode()) &&
        (solver_container[iZone][iInst][MESH_0][FLOW_SOL]->GetStart_AoA_FD()) &&
        iExtIter == solver_container[iZone][iInst][MESH_0][FLOW_SOL]->GetIter_Update_AoA()) {
        //(config[iZone]->GetnExtIter()-config[iZone]->GetIter_dCL_dAlpha() -1 == iExtIter)) {

      if (rank == MASTER_NODE)
        cout << "Storing solution output data locally on each rank (cte. CL mode)." << endl;

      Local_Data_Copy = new su2double*[geometry[iZone][iInst][MESH_0]->GetnPoint()];
      for (iPoint = 0; iPoint < geometry[iZone][iInst][MESH_0]->GetnPoint(); iPoint++) {
        Local_Data_Copy[iPoint] = new su2double[nVar_Par];
      }

      for (iPoint = 0; iPoint < geometry[iZone][iInst][MESH_0]->GetnPoint(); iPoint++) {
        for (iVar = 0; iVar < nVar_Par; iVar++) {
          Local_Data_Copy[iPoint][iVar] = Local_Data[iPoint][iVar];
        }
      }

    }

    /*--- Recover the solution to be used on the final iteration with cte. lift mode. ---*/

    if ((!cont_adj) && (!disc_adj) && (config[iZone]->GetFixed_CL_Mode()) &&
        (solver_container[iZone][iInst][MESH_0][FLOW_SOL]->GetEnd_AoA_FD()) &&
        //(iExtIter - solver_container[iZone][iInst][MESH_0][FLOW_SOL]->GetIter_Update_AoA() == config[iZone]->GetIter_dCL_dAlpha()) &&
        (Local_Data_Copy != NULL)) {

      if (rank == MASTER_NODE)
        cout << "Recovering solution output data locally on each rank (cte. CL mode)." << endl;

      for (iPoint = 0; iPoint < geometry[iZone][iInst][MESH_0]->GetnPoint(); iPoint++) {
        for (iVar = 0; iVar < nVar_Par; iVar++) {
          Local_Data[iPoint][iVar] = Local_Data_Copy[iPoint][iVar];
        }
      }

      for (iPoint = 0; iPoint < geometry[iZone][iInst][MESH_0]->GetnPoint(); iPoint++)
        delete [] Local_Data_Copy[iPoint];
      delete [] Local_Data_Copy;

    }

    /*--- After loading the data local to a processor, we perform a sorting,
     i.e., a linear partitioning of the data across all ranks in the communicator. ---*/

    if (rank == MASTER_NODE) cout << "Sorting output data across all ranks." << endl;

    if (fem_solver)
      SortOutputData_FEM(config[iZone], geometry[iZone][iInst][MESH_0]);
    else
      SortOutputData(config[iZone], geometry[iZone][iInst][MESH_0]);

    /*--- Write either a binary or ASCII restart file in parallel. ---*/

    if (config[iZone]->GetWrt_Binary_Restart()) {
      if (rank == MASTER_NODE) cout << "Writing binary SU2 native restart file." << endl;
      WriteRestart_Parallel_Binary(config[iZone], geometry[iZone][iInst][MESH_0], solver_container[iZone][iInst][MESH_0], iZone, iInst);
    } else {
      if (rank == MASTER_NODE) cout << "Writing ASCII SU2 native restart file." << endl;
      WriteRestart_Parallel_ASCII(config[iZone], geometry[iZone][iInst][MESH_0], solver_container[iZone][iInst][MESH_0], iZone, iInst);
    }

    /*--- Write a slice on a structured mesh if requested. ---*/

    if (config[iZone]->GetWrt_Slice()) {
      WriteCSV_Slice(config[iZone], geometry[iZone][iInst][MESH_0],
                     solver_container[iZone][iInst][MESH_0][FLOW_SOL], iExtIter, iZone, 0);
      WriteCSV_Slice(config[iZone], geometry[iZone][iInst][MESH_0],
                     solver_container[iZone][iInst][MESH_0][FLOW_SOL], iExtIter, iZone, 1);
    }

    /*--- Write the solution files if they are requested and we are executing
     with a single rank (all data on one proc and no comm. overhead). Once we
     have parallel binary versions of Tecplot / ParaView / CGNS / etc., we
     can allow the write of the viz. files as well. ---*/

    if (Wrt_Vol || Wrt_Srf) {

      /*--- First, sort all connectivity into linearly partitioned chunks of elements. ---*/

      if (rank == MASTER_NODE)
        cout << "Preparing element connectivity across all ranks." << endl;

      if (fem_solver)
        SortConnectivity_FEM(config[iZone], geometry[iZone][iInst][MESH_0], iZone);
      else
        SortConnectivity(config[iZone], geometry[iZone][iInst][MESH_0], iZone);

      /*--- Sort the surface data and renumber if for writing. ---*/
      if (Wrt_Srf){
        if (fem_solver)
          SortOutputData_Surface_FEM(config[iZone], geometry[iZone][iInst][MESH_0]);
        else
          SortOutputData_Surface(config[iZone], geometry[iZone][iInst][MESH_0]);
      }
      /*--- Write Tecplot/ParaView ASCII files for the volume and/or surface solutions. ---*/

      if (Wrt_Vol) {

        switch (FileFormat) {

          case TECPLOT:

            /*--- Write a Tecplot ASCII file ---*/

            if (rank == MASTER_NODE) cout << "Writing Tecplot ASCII file volume solution file." << endl;
//            WriteTecplotASCII_Parallel(config[iZone], geometry[iZone][iInst][MESH_0],
//                solver_container[iZone][iInst][MESH_0], iZone, val_nZone, iInst, nInst, false);
            break;



          case TECPLOT_BINARY:

            /*--- Write a Tecplot ASCII file instead for now in serial. ---*/

            if (rank == MASTER_NODE) cout << "Tecplot binary volume files not available in serial with SU2_CFD." << endl;
            if (rank == MASTER_NODE) cout << "  Run SU2_SOL to generate Tecplot binary." << endl;
            if (rank == MASTER_NODE) cout << "Writing Tecplot ASCII file volume solution file instead." << endl;
//            WriteTecplotASCII_Parallel(config[iZone], geometry[iZone][iInst][MESH_0],
//                solver_container[iZone][iInst][MESH_0], iZone, val_nZone, iInst, nInst, false);
            break;



          case PARAVIEW:

            /*--- Write a Paraview ASCII file ---*/

            if (rank == MASTER_NODE) cout << "Writing Paraview ASCII volume solution file." << endl;
//            WriteParaViewASCII_Parallel(config[iZone], geometry[iZone][iInst][MESH_0],
//                solver_container[iZone][iInst][MESH_0], iZone, val_nZone, iInst, nInst, false);
            break;

          case PARAVIEW_BINARY:

            /*--- Write a Paraview binary file ---*/

            if (rank == MASTER_NODE) cout << "Writing Paraview binary volume solution file." << endl;
//            WriteParaViewBinary_Parallel(config[iZone], geometry[iZone][iInst][MESH_0],
//                                        solver_container[iZone][iInst][MESH_0], iZone, val_nZone, false);
            break;

          default:
            break;
          }

        }

        if (Wrt_Srf) {

          switch (FileFormat) {

          case TECPLOT:

            /*--- Write a Tecplot ASCII file ---*/

            if (rank == MASTER_NODE) cout << "Writing Tecplot ASCII file surface solution file." << endl;
//            WriteTecplotASCII_Parallel(config[iZone], geometry[iZone][iInst][MESH_0],
//                solver_container[iZone][iInst][MESH_0], iZone, val_nZone, iInst, nInst, true);
            break;

          case TECPLOT_BINARY:

            /*--- Write a Tecplot ASCII file instead for now in serial. ---*/

            if (rank == MASTER_NODE) cout << "Tecplot binary surface files not available in serial with SU2_CFD." << endl;
            if (rank == MASTER_NODE) cout << "  Run SU2_SOL to generate Tecplot binary." << endl;
            if (rank == MASTER_NODE) cout << "Writing Tecplot ASCII file surface solution file instead." << endl;
//            WriteTecplotASCII_Parallel(config[iZone], geometry[iZone][iInst][MESH_0],
//                solver_container[iZone][iInst][MESH_0], iZone, val_nZone, iInst, nInst, true);
            break;

          case PARAVIEW:

            /*--- Write a Paraview ASCII file ---*/

            if (rank == MASTER_NODE) cout << "Writing Paraview ASCII surface solution file." << endl;
//            WriteParaViewASCII_Parallel(config[iZone], geometry[iZone][iInst][MESH_0],
//                solver_container[iZone][iInst][MESH_0], iZone, val_nZone, iInst, nInst, true);
            break;

          case PARAVIEW_BINARY:

            /*--- Write a Paraview binary file ---*/

            if (rank == MASTER_NODE) cout << "Writing Paraview binary surface solution file." << endl;
//            WriteParaViewBinary_Parallel(config[iZone], geometry[iZone][iInst][MESH_0],
//                                         solver_container[iZone][iInst][MESH_0], iZone, val_nZone, true);
            break;


          default:
            break;
          }

        }

        /*--- Clean up the connectivity data that was allocated for output. ---*/

        DeallocateConnectivity_Parallel(config[iZone], geometry[iZone][iInst][MESH_0], false);
        if (Wrt_Srf) DeallocateConnectivity_Parallel(config[iZone], geometry[iZone][iInst][MESH_0], true);

        /*--- Clean up the surface data that was only needed for output. ---*/

        if (Wrt_Srf) DeallocateSurfaceData_Parallel(config[iZone], geometry[iZone][iInst][MESH_0]);

      }

      /*--- Deallocate the nodal data needed for writing restarts. ---*/

      DeallocateData_Parallel(config[iZone], geometry[iZone][iInst][MESH_0]);

      /*--- Clear the variable names list. ---*/

      Variable_Names.clear();

    }
  }
}

void COutputLegacy::LoadLocalData_Flow(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned short val_iZone) {

  unsigned short iDim;
  unsigned short Kind_Solver = config->GetKind_Solver();
  unsigned short nDim = geometry->GetnDim();

  unsigned long iVar, jVar;
  unsigned long iPoint, jPoint, FirstIndex = NONE, SecondIndex = NONE, iMarker, iVertex;
  unsigned long nVar_First = 0, nVar_Second = 0, nVar_Consv_Par = 0;

  su2double RefArea = config->GetRefArea();
  su2double Gamma = config->GetGamma();
  su2double RefVel2;
  su2double Gas_Constant, Mach2Vel, Mach_Motion, RefDensity, RefPressure = 0.0, factor = 0.0;
  su2double *Aux_Frict_x = NULL, *Aux_Frict_y = NULL, *Aux_Frict_z = NULL, *Aux_Heat = NULL, *Aux_yPlus = NULL, *Aux_Buffet = NULL;
  su2double *Grid_Vel = NULL;
  su2double Q, Grad_Vel[3][3] = {{0.0, 0.0, 0.0},{0.0, 0.0, 0.0},{0.0, 0.0, 0.0}};

  bool transition           = (config->GetKind_Trans_Model() == BC);
  bool dynamic_grid         = config->GetDynamic_Grid();
  bool rotating_frame       = config->GetRotating_Frame();
  bool Wrt_Halo             = config->GetWrt_Halo(), isPeriodic;

  int *Local_Halo = NULL;

  stringstream varname;

  /*--- Set the non-dimensionalization for coefficients. ---*/

  if (dynamic_grid) {
    Gas_Constant = config->GetGas_ConstantND();
    Mach2Vel = sqrt(Gamma*Gas_Constant*config->GetTemperature_FreeStreamND());
    Mach_Motion = config->GetMach_Motion();
    RefVel2 = (Mach_Motion*Mach2Vel)*(Mach_Motion*Mach2Vel);
  }
  else {
    RefVel2 = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      RefVel2  += solver[FLOW_SOL]->GetVelocity_Inf(iDim)*solver[FLOW_SOL]->GetVelocity_Inf(iDim);
  }
  RefDensity  = solver[FLOW_SOL]->GetDensity_Inf();
  RefPressure = solver[FLOW_SOL]->GetPressure_Inf();
  factor = 1.0 / (0.5*RefDensity*RefArea*RefVel2);

  /*--- Use a switch statement to decide how many solver containers we have
   in this zone for output. ---*/

  switch (config->GetKind_Solver()) {
    case EULER : case NAVIER_STOKES: FirstIndex = FLOW_SOL; SecondIndex = NONE; break;
    case RANS : FirstIndex = FLOW_SOL; SecondIndex = TURB_SOL; break;
    default: SecondIndex = NONE; break;
  }

  nVar_First = solver[FirstIndex]->GetnVar();
  if (SecondIndex != NONE) nVar_Second = solver[SecondIndex]->GetnVar();
  nVar_Consv_Par = nVar_First + nVar_Second;

  /*--------------------------------------------------------------------------*/
  /*--- Step 1: Register the variables that will be output. To register a  ---*/
  /*---         variable, two things are required. First, increment the    ---*/
  /*---         counter for the number of variables (nVar_Par), which      ---*/
  /*---         controls the size of the data structure allocation, i.e.,  ---*/
  /*---         the number of columns in an nPoint x nVar structure.       ---*/
  /*---         Second, add a name for the variable to the vector that     ---*/
  /*---         holds the string names.                                    ---*/
  /*--------------------------------------------------------------------------*/

  /*--- All output files first need the grid coordinates. ---*/

  nVar_Par  = 1; Variable_Names.push_back("x");
  nVar_Par += 1; Variable_Names.push_back("y");
  if (geometry->GetnDim() == 3) {
    nVar_Par += 1; Variable_Names.push_back("z");
  }

  /*--- At a mininum, the restarts and visualization files need the
   conservative variables, so these follow next. ---*/

  nVar_Par += nVar_Consv_Par;

  Variable_Names.push_back("Density");
  Variable_Names.push_back("Momentum_x");
  Variable_Names.push_back("Momentum_y");
  if (geometry->GetnDim() == 3) Variable_Names.push_back("Momentum_z");
  Variable_Names.push_back("Energy");

  if (SecondIndex != NONE) {
    if ((config->GetKind_Turb_Model() == SST) || (config->GetKind_Turb_Model() == SST_SUST)) {
      Variable_Names.push_back("TKE");
      Variable_Names.push_back("Omega");
    } else {
      /*--- S-A variants ---*/
      Variable_Names.push_back("Nu_Tilde");
    }
  }

  /*--- If requested, register the limiter and residuals for all of the
   equations in the current flow problem. ---*/

  /*--- Add the limiters ---*/

  if (config->GetWrt_Limiters()) {
    nVar_Par += nVar_Consv_Par;

    Variable_Names.push_back("Limiter_Density");
    Variable_Names.push_back("Limiter_Momentum_x");
    Variable_Names.push_back("Limiter_Momentum_y");
    if (geometry->GetnDim() == 3) Variable_Names.push_back("Limiter_Momentum_z");
    Variable_Names.push_back("Limiter_Energy");

    if (SecondIndex != NONE) {
      if ((config->GetKind_Turb_Model() == SST) || (config->GetKind_Turb_Model() == SST_SUST)) {
        Variable_Names.push_back("Limiter_TKE");
        Variable_Names.push_back("Limiter_Omega");
      } else {
        /*--- S-A variants ---*/
        Variable_Names.push_back("Limiter_Nu_Tilde");
      }
    }
  }

  /*--- Add the residuals ---*/

  if (config->GetWrt_Residuals()) {
    nVar_Par += nVar_Consv_Par;

    Variable_Names.push_back("Residual_Density");
    Variable_Names.push_back("Residual_Momentum_x");
    Variable_Names.push_back("Residual_Momentum_y");
    if (geometry->GetnDim() == 3) Variable_Names.push_back("Residual_Momentum_z");
    Variable_Names.push_back("Residual_Energy");

    if (SecondIndex != NONE) {
      if ((config->GetKind_Turb_Model() == SST) || (config->GetKind_Turb_Model() == SST_SUST)) {
        Variable_Names.push_back("Residual_TKE");
        Variable_Names.push_back("Residual_Omega");
      } else {
        /*--- S-A variants ---*/
        Variable_Names.push_back("Residual_Nu_Tilde");
      }
    }
  }

  /*--- Add the grid velocity. ---*/

  if (dynamic_grid) {
    if (geometry->GetnDim() == 2) nVar_Par += 2;
    else if (geometry->GetnDim() == 3) nVar_Par += 3;

    Variable_Names.push_back("Grid_Velocity_x");
    Variable_Names.push_back("Grid_Velocity_y");
    if (geometry->GetnDim() == 3) Variable_Names.push_back("Grid_Velocity_z");
  }


  /*--- Add Pressure, Temperature, Cp, Mach. ---*/

  nVar_Par += 1;
  Variable_Names.push_back("Pressure");

  nVar_Par += 2;
  Variable_Names.push_back("Temperature");
  Variable_Names.push_back("Mach");

  nVar_Par += 1;
  if ((config->GetTabular_FileFormat() == TAB_CSV)){
    Variable_Names.push_back("Pressure_Coefficient");
  } else {
    Variable_Names.push_back("C<sub>p</sub>");
  }

  /*--- Add Laminar Viscosity, Skin Friction, Heat Flux, & yPlus to the restart file ---*/

  if ((Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) {
    if ((config->GetTabular_FileFormat() == TAB_CSV)){
      nVar_Par += 1; Variable_Names.push_back("Laminar_Viscosity");
      nVar_Par += 2;
      Variable_Names.push_back("Skin_Friction_Coefficient_x");
      Variable_Names.push_back("Skin_Friction_Coefficient_y");
      if (geometry->GetnDim() == 3) {
        nVar_Par += 1; Variable_Names.push_back("Skin_Friction_Coefficient_z");
      }
      nVar_Par += 1;
      Variable_Names.push_back("Heat_Flux");
    } else {
      nVar_Par += 1; Variable_Names.push_back("<greek>m</greek>");
      nVar_Par += 2;
      Variable_Names.push_back("C<sub>f</sub>_x");
      Variable_Names.push_back("C<sub>f</sub>_y");
      if (geometry->GetnDim() == 3) {
        nVar_Par += 1; Variable_Names.push_back("C<sub>f</sub>_z");
      }
      if (config->GetBuffet_Monitoring() || config->GetKind_ObjFunc() == BUFFET_SENSOR){
        Variable_Names.push_back("Buffet_Sensor");
        nVar_Par += 1;
      }
      nVar_Par += 1;
      Variable_Names.push_back("h");
    }
  }

  /*--- Add Eddy Viscosity. ---*/

  if (Kind_Solver == RANS) {
    nVar_Par += 2;
    if ((config->GetTabular_FileFormat() == TAB_CSV)){
      Variable_Names.push_back("Y_Plus");
      Variable_Names.push_back("Eddy_Viscosity");
    } else {
      Variable_Names.push_back("y<sup>+</sup>");
      Variable_Names.push_back("<greek>m</greek><sub>t</sub>");
    }
  }

  /*--- Add the distance to the nearest sharp edge if requested. ---*/

  if (config->GetWrt_SharpEdges()) {
    nVar_Par += 1;
    Variable_Names.push_back("Sharp_Edge_Dist");
  }

  /*--- Add the intermittency for the BC trans. model. ---*/

  if (transition) {
    nVar_Par += 1;
    if ((config->GetTabular_FileFormat() == TAB_CSV)){
      Variable_Names.push_back("gamma_BC");
    } else {
      Variable_Names.push_back("<greek>g</greek><sub>BC</sub>");
    }
  }

  if (config->GetKind_HybridRANSLES()!=NO_HYBRIDRANSLES){
    nVar_Par +=1;
    Variable_Names.push_back("DES_LengthScale");
    nVar_Par +=1;
    Variable_Names.push_back("Wall_Distance");
  }

  if (config->GetKind_RoeLowDiss() != NO_ROELOWDISS){
    nVar_Par +=1;
    Variable_Names.push_back("Roe_Dissipation");
  }

  if (solver[FLOW_SOL]->VerificationSolution) {
    if (solver[FLOW_SOL]->VerificationSolution->ExactSolutionKnown()) {
      nVar_Par += 2*nVar_Consv_Par;
      Variable_Names.push_back("Verification_Density");
      Variable_Names.push_back("Verification_Momentum_x");
      Variable_Names.push_back("Verification_Momentum_y");
      if (geometry->GetnDim() == 3) Variable_Names.push_back("Verification_Momentum_z");
      Variable_Names.push_back("Verification_Energy");
      Variable_Names.push_back("Error_Density");
      Variable_Names.push_back("Error_Momentum_x");
      Variable_Names.push_back("Error_Momentum_y");
      if (geometry->GetnDim() == 3) Variable_Names.push_back("Error_Momentum_z");
      Variable_Names.push_back("Error_Energy");
    }
  }

  /*--- New variables get registered here before the end of the loop. ---*/

  if ((Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) {
    nVar_Par += 2;
    Variable_Names.push_back("Vorticity_x");
    Variable_Names.push_back("Vorticity_y");
    if (geometry->GetnDim() == 3) {
      nVar_Par += 1; Variable_Names.push_back("Vorticity_z");
    }

    if (geometry->GetnDim() == 3) {
      nVar_Par +=1;
      Variable_Names.push_back("Q_Criterion");
    }

  }

  if (rotating_frame) {
    nVar_Par += 2;
    Variable_Names.push_back("Relative_Velocity_x");
    Variable_Names.push_back("Relative_Velocity_y");
    if (geometry->GetnDim() == 3) {
      nVar_Par += 1; Variable_Names.push_back("Relative_Velocity_z");
    }
  }

  if (config->GetWrt_MeshQuality()) {
    nVar_Par +=1;
    Variable_Names.push_back("Orthogonality");
    nVar_Par +=1;
    Variable_Names.push_back("Aspect_Ratio");
    nVar_Par +=1;
    Variable_Names.push_back("Volume_Ratio");
  }



  /*--- Auxiliary vectors for variables defined on surfaces only. ---*/

  if ((Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) {
    Aux_Frict_x = new su2double[geometry->GetnPoint()];
    Aux_Frict_y = new su2double[geometry->GetnPoint()];
    Aux_Frict_z = new su2double[geometry->GetnPoint()];
    Aux_Heat    = new su2double[geometry->GetnPoint()];
    Aux_yPlus   = new su2double[geometry->GetnPoint()];
    Aux_Buffet  = new su2double[geometry->GetnPoint()];

    /*--- First, loop through the mesh in order to find and store the
     value of the viscous coefficients at any surface nodes. They
     will be placed in an auxiliary vector and then communicated like
     all other volumetric variables. ---*/

    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
      Aux_Frict_x[iPoint] = 0.0;
      Aux_Frict_y[iPoint] = 0.0;
      Aux_Frict_z[iPoint] = 0.0;
      Aux_Heat[iPoint]    = 0.0;
      Aux_yPlus[iPoint]   = 0.0;
      Aux_Buffet[iPoint]  = 0.0;
    }
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if (config->GetMarker_All_Plotting(iMarker) == YES) {
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          Aux_Frict_x[iPoint] = solver[FLOW_SOL]->GetCSkinFriction(iMarker, iVertex, 0);
          Aux_Frict_y[iPoint] = solver[FLOW_SOL]->GetCSkinFriction(iMarker, iVertex, 1);
          if (geometry->GetnDim() == 3) Aux_Frict_z[iPoint] = solver[FLOW_SOL]->GetCSkinFriction(iMarker, iVertex, 2);
          Aux_Heat[iPoint] = solver[FLOW_SOL]->GetHeatFlux(iMarker, iVertex);
          Aux_yPlus[iPoint] = solver[FLOW_SOL]->GetYPlus(iMarker, iVertex);
          if (config->GetBuffet_Monitoring() || config->GetKind_ObjFunc() == BUFFET_SENSOR) Aux_Buffet[iPoint] = solver[FLOW_SOL]->GetBuffetSensor(iMarker, iVertex);
        }
      }
    }
  }

  /*--- Allocate the local data structure now that we know how many
   variables are in the output. ---*/

  Local_Data = new su2double*[geometry->GetnPoint()];
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
    Local_Data[iPoint] = new su2double[nVar_Par];
  }

  Local_Halo = new int[geometry->GetnPoint()];
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
    Local_Halo[iPoint] = !geometry->node[iPoint]->GetDomain();

  /*--- Search all send/recv boundaries on this partition for any periodic
   nodes that were part of the original domain. We want to recover these
   for visualization purposes. ---*/

  if (!Wrt_Halo) {
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) {

        /*--- Checking for less than or equal to the rank, because there may
         be some periodic halo nodes that send info to the same rank. ---*/

        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          isPeriodic = ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() > 0) &&
                        (geometry->vertex[iMarker][iVertex]->GetRotation_Type() % 2 == 1));
          if (isPeriodic) Local_Halo[iPoint] = false;
        }
      }
    }
  }

  /*--------------------------------------------------------------------------*/
  /*--- Step 2: Loop over all grid nodes and load up the desired data for  ---*/
  /*---         the restart and vizualization files. Note that we need to  ---*/
  /*---         increment the iVar variable after each variable load.      ---*/
  /*---         The idea is that we're filling up the columns of field     ---*/
  /*---         data for each iPoint (row) of the data structure.          ---*/
  /*---         This data will then be sorted, communicated, and written   ---*/
  /*---         to files automatically after this routine. Note that the   ---*/
  /*---         ordering of the data loading MUST match the order of the   ---*/
  /*---         variable registration above for the files to be correct.   ---*/
  /*--------------------------------------------------------------------------*/

  jPoint = 0;

  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {

    /*--- Check for halos & write only if requested ---*/

    if (!Local_Halo[iPoint] || Wrt_Halo) {

      /*--- Restart the column index with each new point. ---*/

      iVar = 0;

      /*--- Load the grid node coordinate values. ---*/

      for (iDim = 0; iDim < geometry->GetnDim(); iDim++) {
        Local_Data[jPoint][iVar] = geometry->node[iPoint]->GetCoord(iDim);
        if (config->GetSystemMeasurements() == US)
          Local_Data[jPoint][iVar] *= 12.0;
        iVar++;
      }

      /*--- Load the conservative variable states for the mean flow variables. ---*/

      for (jVar = 0; jVar < nVar_First; jVar++) {
        Local_Data[jPoint][iVar] = solver[FirstIndex]->GetNodes()->GetSolution(iPoint, jVar);
        iVar++;
      }

      /*--- If this is RANS, i.e., the second solver container is not empty,
       then load data for the conservative turbulence variables. ---*/

      if (SecondIndex != NONE) {
        for (jVar = 0; jVar < nVar_Second; jVar++) {
          Local_Data[jPoint][iVar] = solver[SecondIndex]->GetNodes()->GetSolution(iPoint, jVar);
          iVar++;
        }
      }

      /*--- If limiters and/or residuals are requested. ---*/

        /*--- Limiters ---*/
        if (config->GetWrt_Limiters()) {
          /*--- Mean Flow Limiters ---*/
          for (jVar = 0; jVar < nVar_First; jVar++) {
            Local_Data[jPoint][iVar] = solver[FirstIndex]->GetNodes()->GetLimiter_Primitive(iPoint, jVar);
            iVar++;
          }
          /*--- RANS Limiters ---*/
          if (SecondIndex != NONE) {
            for (jVar = 0; jVar < nVar_Second; jVar++) {
              Local_Data[jPoint][iVar] = solver[SecondIndex]->GetNodes()->GetLimiter_Primitive(iPoint, jVar);
              iVar++;
            }
          }
        }

        /*--- Residuals ---*/
        if (config->GetWrt_Residuals()) {
          /*--- Mean Flow Residuals ---*/
          for (jVar = 0; jVar < nVar_First; jVar++) {
            Local_Data[jPoint][iVar] = solver[FirstIndex]->LinSysRes(iPoint, jVar);
            iVar++;
          }
          /*--- RANS Residuals ---*/
          if (SecondIndex != NONE) {
            for (jVar = 0; jVar < nVar_Second; jVar++) {
              Local_Data[jPoint][iVar] = solver[SecondIndex]->LinSysRes(iPoint, jVar);
              iVar++;
            }
          }
        }



        /*--- Load buffers with the three grid velocity components. ---*/

        if (dynamic_grid) {
          Grid_Vel = geometry->node[iPoint]->GetGridVel();
          Local_Data[jPoint][iVar] = Grid_Vel[0]; iVar++;
          Local_Data[jPoint][iVar] = Grid_Vel[1]; iVar++;
          if (geometry->GetnDim() == 3) {
            Local_Data[jPoint][iVar] = Grid_Vel[2];
            iVar++;
          }
        }

        /*--- Load data for the pressure, temperature, Cp, and Mach variables. ---*/

        Local_Data[jPoint][iVar] = solver[FLOW_SOL]->GetNodes()->GetPressure(iPoint); iVar++;
        Local_Data[jPoint][iVar] = solver[FLOW_SOL]->GetNodes()->GetTemperature(iPoint); iVar++;
        Local_Data[jPoint][iVar] = sqrt(solver[FLOW_SOL]->GetNodes()->GetVelocity2(iPoint))/solver[FLOW_SOL]->GetNodes()->GetSoundSpeed(iPoint); iVar++;
        Local_Data[jPoint][iVar] = (solver[FLOW_SOL]->GetNodes()->GetPressure(iPoint) - RefPressure)*factor*RefArea; iVar++;

        if ((Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) {

          /*--- Load data for the laminar viscosity. ---*/

          Local_Data[jPoint][iVar] = solver[FLOW_SOL]->GetNodes()->GetLaminarViscosity(iPoint); iVar++;

          /*--- Load data for the skin friction, heat flux, buffet, and y-plus. ---*/

          Local_Data[jPoint][iVar] = Aux_Frict_x[iPoint]; iVar++;
          Local_Data[jPoint][iVar] = Aux_Frict_y[iPoint]; iVar++;
          if (geometry->GetnDim() == 3) {
            Local_Data[jPoint][iVar] = Aux_Frict_z[iPoint];
            iVar++;
          }
          if (config->GetBuffet_Monitoring() || config->GetKind_ObjFunc() == BUFFET_SENSOR) {
            Local_Data[jPoint][iVar] = Aux_Buffet[iPoint];
            iVar++;
          }
          Local_Data[jPoint][iVar] = Aux_Heat[iPoint]; iVar++;

        }

        /*--- Load data for the Eddy viscosity for RANS. ---*/

        if (Kind_Solver == RANS) {
          Local_Data[jPoint][iVar] = Aux_yPlus[iPoint]; iVar++;
          Local_Data[jPoint][iVar] = solver[FLOW_SOL]->GetNodes()->GetEddyViscosity(iPoint); iVar++;
        }

        /*--- Load data for the distance to the nearest sharp edge. ---*/

        if (config->GetWrt_SharpEdges()) {
          Local_Data[jPoint][iVar] = geometry->node[iPoint]->GetSharpEdge_Distance(); iVar++;
        }

        /*--- Load data for the intermittency of the BC trans. model. ---*/

        if (transition) {
          Local_Data[jPoint][iVar] = solver[TURB_SOL]->GetNodes()->GetGammaBC(iPoint); iVar++;
        }

        if (config->GetKind_HybridRANSLES()!=NO_HYBRIDRANSLES){
          Local_Data[jPoint][iVar] = solver[FLOW_SOL]->GetNodes()->GetDES_LengthScale(iPoint); iVar++;
          Local_Data[jPoint][iVar] = geometry->node[iPoint]->GetWall_Distance(); iVar++;
        }

        if (config->GetKind_RoeLowDiss() != NO_ROELOWDISS){
          Local_Data[jPoint][iVar] = solver[FLOW_SOL]->GetNodes()->GetRoe_Dissipation(iPoint); iVar++;
        }

        if (solver[FLOW_SOL]->VerificationSolution) {
          if (solver[FLOW_SOL]->VerificationSolution->ExactSolutionKnown()) {

            /*--- Get the physical time if necessary. ---*/
            su2double time = 0.0;
            if (config->GetTime_Marching()) time = config->GetPhysicalTime();

            /* Set the pointers to the coordinates and solution of this DOF. */
            const su2double *coor = geometry->node[iPoint]->GetCoord();
            su2double *solDOF     = solver[FLOW_SOL]->GetNodes()->GetSolution(iPoint);
            su2double mmsSol[5]   = {0.0,0.0,0.0,0.0,0.0};
            su2double error[5]    = {0.0,0.0,0.0,0.0,0.0};

            /* Get the verification solution. */
            solver[FLOW_SOL]->VerificationSolution->GetSolution(coor, time, mmsSol);
            for (jVar = 0; jVar < nVar_First; jVar++) {
              Local_Data[jPoint][iVar] = mmsSol[jVar];
              iVar++;
            }

            /* Get local error from the verification solution class. */
            solver[FLOW_SOL]->VerificationSolution->GetLocalError(coor, time, solDOF, error);
            for (jVar = 0; jVar < nVar_First; jVar++) {
              Local_Data[jPoint][iVar] = error[jVar];
              iVar++;
            }

          }
        }

        /*--- New variables can be loaded to the Local_Data structure here,
         assuming they were registered above correctly. ---*/

        if ((Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) {
          Local_Data[jPoint][iVar] = solver[FLOW_SOL]->GetNodes()->GetVorticity(iPoint)[0]; iVar++;
          Local_Data[jPoint][iVar] = solver[FLOW_SOL]->GetNodes()->GetVorticity(iPoint)[1]; iVar++;
          if (geometry->GetnDim() == 3) {
            Local_Data[jPoint][iVar] = solver[FLOW_SOL]->GetNodes()->GetVorticity(iPoint)[2];
            iVar++;
          }

          if (nDim == 3){
            for (iDim = 0; iDim < nDim; iDim++) {
              for (unsigned short jDim = 0; jDim < nDim; jDim++) {
                Grad_Vel[iDim][jDim] = solver[FLOW_SOL]->GetNodes()->GetGradient_Primitive(iPoint, iDim+1, jDim);
              }
            }

            /*--- Q Criterion Eq 1.2 of---*/
            /*--- HALLER, G. (2005). An objective definition of a vortex. Journal of Fluid Mechanics, 525, 1-26. doi:10.1017/S0022112004002526 ---*/

            su2double s11 = Grad_Vel[0][0];
            su2double s12 = 0.5 * (Grad_Vel[0][1] + Grad_Vel[1][0]);
            su2double s13 = 0.5 * (Grad_Vel[0][2] + Grad_Vel[2][0]);
            su2double s22 = Grad_Vel[1][1];
            su2double s23 = 0.5 * (Grad_Vel[1][2] + Grad_Vel[2][1]);
            su2double s33 = Grad_Vel[2][2];
            su2double omega12 = 0.5 * (Grad_Vel[0][1] - Grad_Vel[1][0]);
            su2double omega13 = 0.5 * (Grad_Vel[0][2] - Grad_Vel[2][0]);
            su2double omega23 = 0.5 * (Grad_Vel[1][2] - Grad_Vel[2][1]);

            Q = 2. * pow( omega12, 2.) + 2. * pow( omega13, 2.) + 2. * pow( omega23, 2.) - \
                pow( s11, 2.) - pow( s22, 2.) - pow( s33, 2.0) - 2. * pow( s12, 2.) - 2. * pow( s13, 2.) - 2. * pow( s23, 2.0);
            Local_Data[jPoint][iVar] = Q; iVar++;
          }
        }

        /*--- For rotating frame problems, compute the relative velocity. ---*/

        if (rotating_frame) {
          Grid_Vel = geometry->node[iPoint]->GetGridVel();
          su2double *Solution = solver[FLOW_SOL]->GetNodes()->GetSolution(iPoint);
          Local_Data[jPoint][iVar] = Solution[1]/Solution[0] - Grid_Vel[0]; iVar++;
          Local_Data[jPoint][iVar] = Solution[2]/Solution[0] - Grid_Vel[1]; iVar++;
          if (geometry->GetnDim() == 3) {
            Local_Data[jPoint][iVar] = Solution[3]/Solution[0] - Grid_Vel[2];
            iVar++;
          }
        }

        /*--- Output the mesh quality metrics. ---*/

        if (config->GetWrt_MeshQuality()) {
          Local_Data[jPoint][iVar] = geometry->Orthogonality[iPoint]; iVar++;
          Local_Data[jPoint][iVar] = geometry->Aspect_Ratio[iPoint];  iVar++;
          Local_Data[jPoint][iVar] = geometry->Volume_Ratio[iPoint];  iVar++;
        }

    }

    /*--- Increment the point counter, as there may have been halos we
       skipped over during the data loading. ---*/

    jPoint++;

  }


  /*--- Free memory for auxiliary vectors. ---*/

  if ((Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) {
    delete [] Aux_Frict_x;
    delete [] Aux_Frict_y;
    delete [] Aux_Frict_z;
    delete [] Aux_Heat;
    delete [] Aux_yPlus;
    delete [] Aux_Buffet;
  }

  delete [] Local_Halo;

}

void COutputLegacy::LoadLocalData_IncFlow(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned short val_iZone) {

  unsigned short iDim;
  unsigned short Kind_Solver = config->GetKind_Solver();
  unsigned short nDim = geometry->GetnDim();

  unsigned long iVar, jVar;
  unsigned long iPoint, jPoint, FirstIndex = NONE, SecondIndex = NONE, iMarker, iVertex;
  unsigned long nVar_First = 0, nVar_Second = 0, nVar_Consv_Par = 0;

  su2double RefArea = config->GetRefArea();
  su2double RefVel2 = 0.0;
  su2double factor = 0.0;
  su2double RefDensity = 0.0, RefPressure = 0.0;
  su2double *Aux_Frict_x = NULL, *Aux_Frict_y = NULL, *Aux_Frict_z = NULL, *Aux_Heat = NULL, *Aux_yPlus = NULL;
  su2double *Grid_Vel = NULL;
  su2double Q, Grad_Vel[3][3] = {{0.0, 0.0, 0.0},{0.0, 0.0, 0.0},{0.0, 0.0, 0.0}};

  bool transition           = (config->GetKind_Trans_Model() == BC);
  bool dynamic_grid         = config->GetDynamic_Grid();
  bool Wrt_Halo             = config->GetWrt_Halo(), isPeriodic;
  bool variable_density     = (config->GetKind_DensityModel() == VARIABLE);
  bool energy               = config->GetEnergy_Equation();
  bool weakly_coupled_heat  = config->GetWeakly_Coupled_Heat();
  bool wrt_cp           = (variable_density &&
                           (config->GetKind_FluidModel() == INC_IDEAL_GAS_POLY));
  bool wrt_kt           = ((config->GetKind_ConductivityModel() != CONSTANT_CONDUCTIVITY) &&
                           (config->GetViscous()));
  bool p1_radiation = (config->GetKind_RadiationModel() == P1_MODEL);
  int *Local_Halo = NULL;

  stringstream varname;

  /*--- Evaluate reference values for non-dimensionalization.
   For dimensional or non-dim based on initial values, use
   the far-field state (inf). For a custom non-dim based
   on user-provided reference values, use the ref values
   to compute the forces. ---*/

  if ((config->GetRef_Inc_NonDim() == DIMENSIONAL) ||
      (config->GetRef_Inc_NonDim() == INITIAL_VALUES)) {
    RefDensity  = solver[FLOW_SOL]->GetDensity_Inf();
    RefVel2 = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      RefVel2  += solver[FLOW_SOL]->GetVelocity_Inf(iDim)*solver[FLOW_SOL]->GetVelocity_Inf(iDim);
  }
  else if (config->GetRef_Inc_NonDim() == REFERENCE_VALUES) {
    RefDensity = config->GetInc_Density_Ref();
    RefVel2    = config->GetInc_Velocity_Ref()*config->GetInc_Velocity_Ref();
  }

  /*--- Reference pressure is always the far-field value (0.0). ---*/

  RefPressure = solver[FLOW_SOL]->GetPressure_Inf();

  /*--- Compute factor for force coefficients. ---*/

  factor = 1.0 / (0.5*RefDensity*RefArea*RefVel2);

  /*--- Use a switch statement to decide how many solver containers we have
   in this zone for output. ---*/

  switch (config->GetKind_Solver()) {
    case INC_EULER : case INC_NAVIER_STOKES: FirstIndex = FLOW_SOL; SecondIndex = NONE; break;
    case INC_RANS : FirstIndex = FLOW_SOL; SecondIndex = TURB_SOL; break;
    default: SecondIndex = NONE; break;
  }

  nVar_First = solver[FirstIndex]->GetnVar();
  if ((!energy) && (!weakly_coupled_heat)) nVar_First--;
  if (SecondIndex != NONE) nVar_Second = solver[SecondIndex]->GetnVar();
  nVar_Consv_Par = nVar_First + nVar_Second;

  /*--------------------------------------------------------------------------*/
  /*--- Step 1: Register the variables that will be output. To register a  ---*/
  /*---         variable, two things are required. First, increment the    ---*/
  /*---         counter for the number of variables (nVar_Par), which      ---*/
  /*---         controls the size of the data structure allocation, i.e.,  ---*/
  /*---         the number of columns in an nPoint x nVar structure.       ---*/
  /*---         Second, add a name for the variable to the vector that     ---*/
  /*---         holds the string names.                                    ---*/
  /*--------------------------------------------------------------------------*/

  /*--- All output files first need the grid coordinates. ---*/

  nVar_Par  = 1; Variable_Names.push_back("x");
  nVar_Par += 1; Variable_Names.push_back("y");
  if (geometry->GetnDim() == 3) {
    nVar_Par += 1; Variable_Names.push_back("z");
  }

  /*--- At a mininum, the restarts and visualization files need the
   conservative variables, so these follow next. ---*/

  nVar_Par += nVar_Consv_Par;

  /*--- The incompressible solver uses primitives as the working variables. ---*/

  Variable_Names.push_back("Pressure");
  Variable_Names.push_back("Velocity_x");
  Variable_Names.push_back("Velocity_y");
  if (geometry->GetnDim() == 3) Variable_Names.push_back("Velocity_z");
  if (energy || weakly_coupled_heat) Variable_Names.push_back("Temperature");

  if (SecondIndex != NONE) {
    if ((config->GetKind_Turb_Model() == SST) || (config->GetKind_Turb_Model() == SST_SUST)) {
      Variable_Names.push_back("TKE");
      Variable_Names.push_back("Omega");
    } else {
      /*--- S-A variants ---*/
      Variable_Names.push_back("Nu_Tilde");
    }
  }

  /*--- Add restart data for P1 radiative model. ---*/
  if (p1_radiation){
    nVar_Par += 1;
    Variable_Names.push_back("Radiative_Energy(P1)");
  }

  /*--- If requested, register the limiter and residuals for all of the
   equations in the current flow problem. ---*/


  /*--- Add the limiters ---*/

  if (config->GetWrt_Limiters()) {
    nVar_Par += nVar_Consv_Par;

    Variable_Names.push_back("Limiter_Pressure");
    Variable_Names.push_back("Limiter_Velocity_x");
    Variable_Names.push_back("Limiter_Velocity_y");
    if (geometry->GetnDim() == 3) Variable_Names.push_back("Limiter_Velocity_z");
    if (energy || weakly_coupled_heat)
      Variable_Names.push_back("Limiter_Temperature");

    if (SecondIndex != NONE) {
      if ((config->GetKind_Turb_Model() == SST) || (config->GetKind_Turb_Model() == SST_SUST)) {
        Variable_Names.push_back("Limiter_TKE");
        Variable_Names.push_back("Limiter_Omega");
      } else {
        /*--- S-A variants ---*/
        Variable_Names.push_back("Limiter_Nu_Tilde");
      }
    }
  }

  /*--- Add the residuals ---*/

  if (config->GetWrt_Residuals()) {
    nVar_Par += nVar_Consv_Par;

    Variable_Names.push_back("Residual_Pressure");
    Variable_Names.push_back("Residual_Velocity_x");
    Variable_Names.push_back("Residual_Velocity_y");
    if (geometry->GetnDim() == 3) Variable_Names.push_back("Residual_Velocity_z");
    if (energy || weakly_coupled_heat)
      Variable_Names.push_back("Residual_Temperature");

    if (SecondIndex != NONE) {
      if ((config->GetKind_Turb_Model() == SST) || (config->GetKind_Turb_Model() == SST_SUST)) {
        Variable_Names.push_back("Residual_TKE");
        Variable_Names.push_back("Residual_Omega");
      } else {
        /*--- S-A variants ---*/
        Variable_Names.push_back("Residual_Nu_Tilde");
      }
    }
  }

  /*--- Add the grid velocity. ---*/

  if (dynamic_grid) {
    if (geometry->GetnDim() == 2) nVar_Par += 2;
    else if (geometry->GetnDim() == 3) nVar_Par += 3;

    Variable_Names.push_back("Grid_Velocity_x");
    Variable_Names.push_back("Grid_Velocity_y");
    if (geometry->GetnDim() == 3) Variable_Names.push_back("Grid_Velocity_z");
  }

  /*--- Add Cp. ---*/

  nVar_Par += 1;
  if ((config->GetTabular_FileFormat() == TAB_CSV)){
    Variable_Names.push_back("Pressure_Coefficient");
  } else {
    Variable_Names.push_back("C<sub>p</sub>");
  }

  /*--- Add Laminar Viscosity, Skin Friction, and Heat Flux to the restart file ---*/

  if ((Kind_Solver == INC_NAVIER_STOKES) || (Kind_Solver == INC_RANS)) {
    if ((config->GetTabular_FileFormat() == TAB_CSV)){
      nVar_Par += 1; Variable_Names.push_back("Laminar_Viscosity");
      nVar_Par += 2;
      Variable_Names.push_back("Skin_Friction_Coefficient_x");
      Variable_Names.push_back("Skin_Friction_Coefficient_y");
      if (geometry->GetnDim() == 3) {
        nVar_Par += 1; Variable_Names.push_back("Skin_Friction_Coefficient_z");
      }
      if (energy || weakly_coupled_heat) {
        nVar_Par += 1;
        Variable_Names.push_back("Heat_Flux");
      }
    } else {
      nVar_Par += 1; Variable_Names.push_back("<greek>m</greek>");
      nVar_Par += 2;
      Variable_Names.push_back("C<sub>f</sub>_x");
      Variable_Names.push_back("C<sub>f</sub>_y");
      if (geometry->GetnDim() == 3) {
        nVar_Par += 1; Variable_Names.push_back("C<sub>f</sub>_z");
      }
      if (energy || weakly_coupled_heat) {
        nVar_Par += 1;
        Variable_Names.push_back("h");
      }
    }
  }

  /*--- Add Y+ and Eddy Viscosity. ---*/

  if (Kind_Solver == INC_RANS) {
    nVar_Par += 2;
    if ((config->GetTabular_FileFormat() == TAB_CSV)){
      Variable_Names.push_back("Y_Plus");
      Variable_Names.push_back("Eddy_Viscosity");
    } else {
      Variable_Names.push_back("y<sup>+</sup>");
      Variable_Names.push_back("<greek>m</greek><sub>t</sub>");
    }
  }

  /*--- Add the distance to the nearest sharp edge if requested. ---*/

  if (config->GetWrt_SharpEdges()) {
    nVar_Par += 1;
    Variable_Names.push_back("Sharp_Edge_Dist");
  }

  /*--- Add the intermittency for the BC trans. model. ---*/

  if (transition) {
    nVar_Par += 1;
    if ((config->GetTabular_FileFormat() == TAB_CSV)){
      Variable_Names.push_back("gamma_BC");
    } else {
      Variable_Names.push_back("<greek>g</greek><sub>BC</sub>");
    }
  }

  if (variable_density) {
    nVar_Par += 1;
    Variable_Names.push_back("Density");
  }

  if (wrt_cp) {
    nVar_Par += 1;
    Variable_Names.push_back("Specific_Heat");
  }

  if (wrt_kt) {
    nVar_Par += 1;
    Variable_Names.push_back("Thermal_Conductivity");
  }

  if (solver[FLOW_SOL]->VerificationSolution) {
    if (solver[FLOW_SOL]->VerificationSolution->ExactSolutionKnown()) {
      nVar_Par += 2*nVar_Consv_Par;
      Variable_Names.push_back("Verification_Pressure");
      Variable_Names.push_back("Verification_Velocity_x");
      Variable_Names.push_back("Verification_Velocity_y");
      if (geometry->GetnDim() == 3) Variable_Names.push_back("Verification_Velocity_z");
      if (energy || weakly_coupled_heat) Variable_Names.push_back("Verification_Temperature");
      Variable_Names.push_back("Error_Pressure");
      Variable_Names.push_back("Error_Velocity_x");
      Variable_Names.push_back("Error_Velocity_y");
      if (geometry->GetnDim() == 3) Variable_Names.push_back("Error_Velocity_z");
      if (energy || weakly_coupled_heat) Variable_Names.push_back("Error_Temperature");
    }
  }

  if ((Kind_Solver == INC_NAVIER_STOKES) || (Kind_Solver == INC_RANS)) {
    nVar_Par += 2;
    Variable_Names.push_back("Vorticity_x");
    Variable_Names.push_back("Vorticity_y");
    nVar_Par += 1; Variable_Names.push_back("Vorticity_z");

    if (geometry->GetnDim() == 3) {
      nVar_Par +=1;
      Variable_Names.push_back("Q_Criterion");
    }
  }

  if (config->GetWrt_MeshQuality()) {
    nVar_Par +=1;
    Variable_Names.push_back("Orthogonality");
    nVar_Par +=1;
    Variable_Names.push_back("Aspect_Ratio");
    nVar_Par +=1;
    Variable_Names.push_back("Volume_Ratio");
  }

  /*--- New variables get registered here before the end of the loop. ---*/



  /*--- Auxiliary vectors for variables defined on surfaces only. ---*/

  if ((Kind_Solver == INC_NAVIER_STOKES) || (Kind_Solver == INC_RANS)) {
    Aux_Frict_x = new su2double[geometry->GetnPoint()];
    Aux_Frict_y = new su2double[geometry->GetnPoint()];
    Aux_Frict_z = new su2double[geometry->GetnPoint()];
    Aux_Heat    = new su2double[geometry->GetnPoint()];
    Aux_yPlus   = new su2double[geometry->GetnPoint()];

    /*--- First, loop through the mesh in order to find and store the
     value of the viscous coefficients at any surface nodes. They
     will be placed in an auxiliary vector and then communicated like
     all other volumetric variables. ---*/

    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
      Aux_Frict_x[iPoint] = 0.0;
      Aux_Frict_y[iPoint] = 0.0;
      Aux_Frict_z[iPoint] = 0.0;
      Aux_Heat[iPoint]    = 0.0;
      Aux_yPlus[iPoint]   = 0.0;
    }
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if (config->GetMarker_All_Plotting(iMarker) == YES) {
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          Aux_Frict_x[iPoint] = solver[FLOW_SOL]->GetCSkinFriction(iMarker, iVertex, 0);
          Aux_Frict_y[iPoint] = solver[FLOW_SOL]->GetCSkinFriction(iMarker, iVertex, 1);
          if (geometry->GetnDim() == 3) Aux_Frict_z[iPoint] = solver[FLOW_SOL]->GetCSkinFriction(iMarker, iVertex, 2);
          Aux_Heat[iPoint] = solver[FLOW_SOL]->GetHeatFlux(iMarker, iVertex);
          Aux_yPlus[iPoint] = solver[FLOW_SOL]->GetYPlus(iMarker, iVertex);
        }
      }
    }
  }

  /*--- Allocate the local data structure now that we know how many
   variables are in the output. ---*/

  Local_Data = new su2double*[geometry->GetnPoint()];
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
    Local_Data[iPoint] = new su2double[nVar_Par];
  }

  Local_Halo = new int[geometry->GetnPoint()];
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
    Local_Halo[iPoint] = !geometry->node[iPoint]->GetDomain();

  /*--- Search all send/recv boundaries on this partition for any periodic
   nodes that were part of the original domain. We want to recover these
   for visualization purposes. ---*/

  if (!Wrt_Halo) {
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) {

        /*--- Checking for less than or equal to the rank, because there may
         be some periodic halo nodes that send info to the same rank. ---*/

        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          isPeriodic = ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() > 0) &&
                        (geometry->vertex[iMarker][iVertex]->GetRotation_Type() % 2 == 1));
          if (isPeriodic) Local_Halo[iPoint] = false;
        }
      }
    }
  }

  /*--------------------------------------------------------------------------*/
  /*--- Step 2: Loop over all grid nodes and load up the desired data for  ---*/
  /*---         the restart and vizualization files. Note that we need to  ---*/
  /*---         increment the iVar variable after each variable load.      ---*/
  /*---         The idea is that we're filling up the columns of field     ---*/
  /*---         data for each iPoint (row) of the data structure.          ---*/
  /*---         This data will then be sorted, communicated, and written   ---*/
  /*---         to files automatically after this routine. Note that the   ---*/
  /*---         ordering of the data loading MUST match the order of the   ---*/
  /*---         variable registration above for the files to be correct.   ---*/
  /*--------------------------------------------------------------------------*/

  jPoint = 0;

  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {

    /*--- Check for halos & write only if requested ---*/

    if (!Local_Halo[iPoint] || Wrt_Halo) {

      /*--- Restart the column index with each new point. ---*/

      iVar = 0;

      /*--- Load the grid node coordinate values. ---*/

      for (iDim = 0; iDim < geometry->GetnDim(); iDim++) {
        Local_Data[jPoint][iVar] = geometry->node[iPoint]->GetCoord(iDim); iVar++;
      }

      /*--- Load the conservative variable states for the mean flow variables. ---*/

      Local_Data[jPoint][iVar] = solver[FirstIndex]->GetNodes()->GetSolution(iPoint, 0); iVar++;
      Local_Data[jPoint][iVar] = solver[FirstIndex]->GetNodes()->GetSolution(iPoint, 1); iVar++;
      Local_Data[jPoint][iVar] = solver[FirstIndex]->GetNodes()->GetSolution(iPoint, 2); iVar++;
      if (nDim == 3) {
        Local_Data[jPoint][iVar] = solver[FirstIndex]->GetNodes()->GetSolution(iPoint, 3); iVar++;
      }
      if (weakly_coupled_heat) {
        Local_Data[jPoint][iVar] = solver[HEAT_SOL]->GetNodes()->GetSolution(iPoint, 0);
        iVar++;
      } else if (energy) {
        Local_Data[jPoint][iVar] = solver[FirstIndex]->GetNodes()->GetSolution(iPoint, nDim+1);
        iVar++;
      }

      /*--- If this is RANS, i.e., the second solver container is not empty,
       then load data for the conservative turbulence variables. ---*/

      if (SecondIndex != NONE) {
        for (jVar = 0; jVar < nVar_Second; jVar++) {
          Local_Data[jPoint][iVar] = solver[SecondIndex]->GetNodes()->GetSolution(iPoint, jVar); iVar++;
        }
      }

      /*--- Add restart data for P1 radiative model. ---*/
      if (p1_radiation){
        Local_Data[jPoint][iVar] = solver[RAD_SOL]->GetNodes()->GetSolution(iPoint,0); iVar++;
      }

      /*--- If limiters and/or residuals are requested. ---*/

      /*--- Limiters ---*/
      if (config->GetWrt_Limiters()) {
        /*--- Mean Flow Limiters ---*/
        for (jVar = 0; jVar < nVar_First; jVar++) {
          Local_Data[jPoint][iVar] = solver[FirstIndex]->GetNodes()->GetLimiter_Primitive(iPoint, jVar);
          iVar++;
        }
        /*--- RANS Limiters ---*/
        if (SecondIndex != NONE) {
          for (jVar = 0; jVar < nVar_Second; jVar++) {
            Local_Data[jPoint][iVar] = solver[SecondIndex]->GetNodes()->GetLimiter_Primitive(iPoint, jVar);
            iVar++;
          }
        }
      }

      /*--- Residuals ---*/
      if (config->GetWrt_Residuals()) {
        /*--- Mean Flow Residuals ---*/
        for (jVar = 0; jVar < nVar_First; jVar++) {
          Local_Data[jPoint][iVar] = solver[FirstIndex]->LinSysRes(iPoint, jVar);
          iVar++;
        }

        /*--- RANS Residuals ---*/
        if (SecondIndex != NONE) {
          for (jVar = 0; jVar < nVar_Second; jVar++) {
            Local_Data[jPoint][iVar] = solver[SecondIndex]->LinSysRes(iPoint, jVar);
            iVar++;
          }
        }
      }


      /*--- Load buffers with the three grid velocity components. ---*/

      if (dynamic_grid) {
        Grid_Vel = geometry->node[iPoint]->GetGridVel();
        Local_Data[jPoint][iVar] = Grid_Vel[0]; iVar++;
        Local_Data[jPoint][iVar] = Grid_Vel[1]; iVar++;
        if (geometry->GetnDim() == 3) {
          Local_Data[jPoint][iVar] = Grid_Vel[2];
          iVar++;
        }
      }

      /*--- Load data for Cp and Mach variables. ---*/

      Local_Data[jPoint][iVar] = (solver[FLOW_SOL]->GetNodes()->GetPressure(iPoint) - RefPressure)*factor*RefArea; iVar++;

      if ((Kind_Solver == INC_NAVIER_STOKES) || (Kind_Solver == INC_RANS)) {

        /*--- Load data for the laminar viscosity. ---*/

        Local_Data[jPoint][iVar] = solver[FLOW_SOL]->GetNodes()->GetLaminarViscosity(iPoint); iVar++;

        /*--- Load data for the skin friction, heat flux, and y-plus. ---*/

        Local_Data[jPoint][iVar] = Aux_Frict_x[iPoint]; iVar++;
        Local_Data[jPoint][iVar] = Aux_Frict_y[iPoint]; iVar++;
        if (geometry->GetnDim() == 3) {
          Local_Data[jPoint][iVar] = Aux_Frict_z[iPoint];
          iVar++;
        }

        if (energy || weakly_coupled_heat) {
          Local_Data[jPoint][iVar] = Aux_Heat[iPoint]; iVar++;
        }

      }

      /*--- Load data for the Eddy viscosity for RANS. ---*/

      if (Kind_Solver == INC_RANS) {
        Local_Data[jPoint][iVar] = Aux_yPlus[iPoint]; iVar++;
        Local_Data[jPoint][iVar] = solver[FLOW_SOL]->GetNodes()->GetEddyViscosity(iPoint); iVar++;
      }

      /*--- Load data for the distance to the nearest sharp edge. ---*/

      if (config->GetWrt_SharpEdges()) {
        Local_Data[jPoint][iVar] = geometry->node[iPoint]->GetSharpEdge_Distance(); iVar++;
      }

      /*--- Load data for the intermittency of the BC trans. model. ---*/

      if (transition) {
        Local_Data[jPoint][iVar] = solver[TURB_SOL]->GetNodes()->GetGammaBC(iPoint); iVar++;
      }

      /*--- Load density if we are solving a variable density problem. ---*/

      if (variable_density) {
        Local_Data[jPoint][iVar] = solver[FLOW_SOL]->GetNodes()->GetDensity(iPoint); iVar++;
      }

      /*--- Load Cp and conductivity if they are temperature-dependent. ---*/
      if (wrt_cp) {
        Local_Data[jPoint][iVar] = solver[FLOW_SOL]->GetNodes()->GetSpecificHeatCp(iPoint); iVar++;
      }

      if (wrt_kt) {
        Local_Data[jPoint][iVar] = solver[FLOW_SOL]->GetNodes()->GetThermalConductivity(iPoint); iVar++;
      }

      if (solver[FLOW_SOL]->VerificationSolution) {
        if (solver[FLOW_SOL]->VerificationSolution->ExactSolutionKnown()) {

          /*--- Get the physical time if necessary. ---*/
          su2double time = 0.0;
          if (config->GetTime_Marching()) time = config->GetPhysicalTime();

          /* Set the pointers to the coordinates and solution of this DOF. */
          const su2double *coor = geometry->node[iPoint]->GetCoord();
          su2double *solDOF     = solver[FLOW_SOL]->GetNodes()->GetSolution(iPoint);
          su2double mmsSol[5]   = {0.0,0.0,0.0,0.0,0.0};
          su2double error[5]    = {0.0,0.0,0.0,0.0,0.0};

          /* Get the verification solution. */
          solver[FLOW_SOL]->VerificationSolution->GetSolution(coor, time, mmsSol);
          for (jVar = 0; jVar < nVar_First; jVar++) {
            Local_Data[jPoint][iVar] = mmsSol[jVar];
            iVar++;
          }

          /* Get local error from the verification solution class. */
          solver[FLOW_SOL]->VerificationSolution->GetLocalError(coor, time, solDOF, error);
          for (jVar = 0; jVar < nVar_First; jVar++) {
            Local_Data[jPoint][iVar] = error[jVar];
            iVar++;
          }

        }
      }

      if ((Kind_Solver == INC_NAVIER_STOKES) || (Kind_Solver == INC_RANS)) {

        Local_Data[jPoint][iVar] = solver[FLOW_SOL]->GetNodes()->GetVorticity(iPoint)[0]; iVar++;
        Local_Data[jPoint][iVar] = solver[FLOW_SOL]->GetNodes()->GetVorticity(iPoint)[1]; iVar++;
        Local_Data[jPoint][iVar] = solver[FLOW_SOL]->GetNodes()->GetVorticity(iPoint)[2]; iVar++;

        if (nDim == 3){
          for (iDim = 0; iDim < nDim; iDim++) {
            for (unsigned short jDim = 0; jDim < nDim; jDim++) {
              Grad_Vel[iDim][jDim] = solver[FLOW_SOL]->GetNodes()->GetGradient_Primitive(iPoint, iDim+1, jDim);
            }
          }

          /*--- Q Criterion Eq 1.2 of---*/
          /*--- HALLER, G. (2005). An objective definition of a vortex. Journal of Fluid Mechanics, 525, 1-26. doi:10.1017/S0022112004002526 ---*/

          su2double s11 = Grad_Vel[0][0];
          su2double s12 = 0.5 * (Grad_Vel[0][1] + Grad_Vel[1][0]);
          su2double s13 = 0.5 * (Grad_Vel[0][2] + Grad_Vel[2][0]);
          su2double s22 = Grad_Vel[1][1];
          su2double s23 = 0.5 * (Grad_Vel[1][2] + Grad_Vel[2][1]);
          su2double s33 = Grad_Vel[2][2];
          su2double omega12 = 0.5 * (Grad_Vel[0][1] - Grad_Vel[1][0]);
          su2double omega13 = 0.5 * (Grad_Vel[0][2] - Grad_Vel[2][0]);
          su2double omega23 = 0.5 * (Grad_Vel[1][2] - Grad_Vel[2][1]);

          Q = 2. * pow( omega12, 2.) + 2. * pow( omega13, 2.) + 2. * pow( omega23, 2.) - \
              pow( s11, 2.) - pow( s22, 2.) - pow( s33, 2.0) - 2. * pow( s12, 2.) - 2. * pow( s13, 2.) - 2. * pow( s23, 2.0);

          Local_Data[jPoint][iVar] = Q; iVar++;
        }
      }

      /*--- Output the mesh quality metrics. ---*/

      if (config->GetWrt_MeshQuality()) {
        Local_Data[jPoint][iVar] = geometry->Orthogonality[iPoint]; iVar++;
        Local_Data[jPoint][iVar] = geometry->Aspect_Ratio[iPoint];  iVar++;
        Local_Data[jPoint][iVar] = geometry->Volume_Ratio[iPoint];  iVar++;
      }

      /*--- New variables can be loaded to the Local_Data structure here,
         assuming they were registered above correctly. ---*/

    }

    /*--- Increment the point counter, as there may have been halos we
       skipped over during the data loading. ---*/

    jPoint++;

  }


  /*--- Free memory for auxiliary vectors. ---*/

  if ((Kind_Solver == INC_NAVIER_STOKES) || (Kind_Solver == INC_RANS)) {
    delete [] Aux_Frict_x;
    delete [] Aux_Frict_y;
    delete [] Aux_Frict_z;
    delete [] Aux_Heat;
    delete [] Aux_yPlus;
  }

  delete [] Local_Halo;

}

void COutputLegacy::LoadLocalData_AdjFlow(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned short val_iZone) {

  unsigned short iDim;
  unsigned short Kind_Solver = config->GetKind_Solver();
  unsigned short nDim = geometry->GetnDim();

  unsigned long iVar, jVar;
  unsigned long iPoint, jPoint, FirstIndex = NONE, SecondIndex = NONE, iMarker, iVertex;
  unsigned long nVar_First = 0, nVar_Second = 0, nVar_Consv_Par = 0;

  su2double *Aux_Sens = NULL;
  su2double *Grid_Vel = NULL;
  su2double *Normal, Area;

  bool incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  bool dynamic_grid   = config->GetDynamic_Grid();

  bool Wrt_Halo       = config->GetWrt_Halo(), isPeriodic;

  int *Local_Halo;

  stringstream varname;

  /*--- Use a switch statement to decide how many solver containers we have
   in this zone for output. ---*/

  switch (config->GetKind_Solver()) {
    case ADJ_EULER : case ADJ_NAVIER_STOKES : FirstIndex = ADJFLOW_SOL; SecondIndex = NONE; break;
    case ADJ_RANS : FirstIndex = ADJFLOW_SOL; if (config->GetFrozen_Visc_Cont()) SecondIndex = NONE; else SecondIndex = ADJTURB_SOL; break;
    case DISC_ADJ_EULER: case DISC_ADJ_NAVIER_STOKES:
    case DISC_ADJ_INC_EULER: case DISC_ADJ_INC_NAVIER_STOKES:
      FirstIndex = ADJFLOW_SOL; SecondIndex = NONE;  break;
    case DISC_ADJ_RANS: case DISC_ADJ_INC_RANS: FirstIndex = ADJFLOW_SOL; if (config->GetFrozen_Visc_Disc()) SecondIndex = NONE; else SecondIndex = ADJTURB_SOL;  break;
  }

  nVar_First = solver[FirstIndex]->GetnVar();
  if (SecondIndex != NONE) nVar_Second = solver[SecondIndex]->GetnVar();
  nVar_Consv_Par = nVar_First + nVar_Second;

  /*--------------------------------------------------------------------------*/
  /*--- Step 1: Register the variables that will be output. To register a  ---*/
  /*---         variable, two things are required. First, increment the    ---*/
  /*---         counter for the number of variables (nVar_Par), which      ---*/
  /*---         controls the size of the data structure allocation, i.e.,  ---*/
  /*---         the number of columns in an nPoint x nVar structure.       ---*/
  /*---         Second, add a name for the variable to the vector that     ---*/
  /*---         holds the string names.                                    ---*/
  /*--------------------------------------------------------------------------*/

  /*--- All output files first need the grid coordinates. ---*/

  nVar_Par  = 1; Variable_Names.push_back("x");
  nVar_Par += 1; Variable_Names.push_back("y");
  if (geometry->GetnDim() == 3) {
    nVar_Par += 1; Variable_Names.push_back("z");
  }

  /*--- At a mininum, the restarts and visualization files need the
   conservative variables, so these follow next. ---*/

  nVar_Par += nVar_Consv_Par;

  /*--- For now, leave the names as "Conservative_", etc., in order
   to avoid confusion with the serial version, which still prints these
   names. Names can be set alternatively by using the commented code
   below. ---*/

  if (incompressible) {
    Variable_Names.push_back("Adjoint_Pressure");
    Variable_Names.push_back("Adjoint_Velocity_x");
    Variable_Names.push_back("Adjoint_Velocity_y");
    if (geometry->GetnDim() == 3) Variable_Names.push_back("Adjoint_Velocity_z");
    Variable_Names.push_back("Adjoint_Temperature");
  } else {
    Variable_Names.push_back("Adjoint_Density");
    Variable_Names.push_back("Adjoint_Momentum_x");
    Variable_Names.push_back("Adjoint_Momentum_y");
    if (geometry->GetnDim() == 3)
      Variable_Names.push_back("Adjoint_Momentum_z");
    Variable_Names.push_back("Adjoint_Energy");
  }
  if (SecondIndex != NONE) {
    if ((config->GetKind_Turb_Model() == SST) || (config->GetKind_Turb_Model() == SST_SUST)) {
      Variable_Names.push_back("Adjoint_TKE");
      Variable_Names.push_back("Adjoint_Omega");
    } else {
      /*--- S-A variants ---*/
      Variable_Names.push_back("Adjoint_Nu_Tilde");
    }
  }

  /*--- For the discrete adjoint, we have the full field of sensitivity
   in each coordinate direction. ---*/

  if ((Kind_Solver == DISC_ADJ_EULER)         ||
      (Kind_Solver == DISC_ADJ_NAVIER_STOKES) ||
      (Kind_Solver == DISC_ADJ_RANS) ||
      (Kind_Solver == DISC_ADJ_INC_EULER)         ||
      (Kind_Solver == DISC_ADJ_INC_NAVIER_STOKES) ||
      (Kind_Solver == DISC_ADJ_INC_RANS)) {
    nVar_Par += nDim;
    Variable_Names.push_back("Sensitivity_x");
    Variable_Names.push_back("Sensitivity_y");
    if (geometry->GetnDim()== 3)
      Variable_Names.push_back("Sensitivity_z");
  }

  /*--- If requested, register the limiter and residuals for all of the
   equations in the current flow problem. ---*/

  /*--- Add the limiters ---*/

  if (config->GetWrt_Limiters()) {
    nVar_Par += nVar_Consv_Par;
    if (incompressible) {
      Variable_Names.push_back("Limiter_Adjoint_Pressure");
      Variable_Names.push_back("Limiter_Adjoint_Velocity_x");
      Variable_Names.push_back("Limiter_Adjoint_Velocity_y");
      if (geometry->GetnDim() == 3) Variable_Names.push_back("Limiter_Adjoint_Velocity_z");
      Variable_Names.push_back("Limiter_Adjoint_Temperature");
    } else {
      Variable_Names.push_back("Limiter_Adjoint_Density");
      Variable_Names.push_back("Limiter_Adjoint_Momentum_x");
      Variable_Names.push_back("Limiter_Adjoint_Momentum_y");
      if (geometry->GetnDim() == 3)
        Variable_Names.push_back("Limiter_Adjoint_Momentum_z");
      Variable_Names.push_back("Limiter_Adjoint_Energy");
    }
    if (SecondIndex != NONE) {
      if ((config->GetKind_Turb_Model() == SST) || (config->GetKind_Turb_Model() == SST_SUST)) {
        Variable_Names.push_back("Limiter_Adjoint_TKE");
        Variable_Names.push_back("Limiter_Adjoint_Omega");
      } else {
        /*--- S-A variants ---*/
        Variable_Names.push_back("Limiter_Adjoint_Nu_Tilde");
      }
    }
  }

  /*--- Add the residuals ---*/

  if (config->GetWrt_Residuals()) {
    nVar_Par += nVar_Consv_Par;
    if (incompressible) {
      Variable_Names.push_back("Residual_Adjoint_Pressure");
      Variable_Names.push_back("Residual_Adjoint_Velocity_x");
      Variable_Names.push_back("Residual_Adjoint_Velocity_y");
      if (geometry->GetnDim() == 3) Variable_Names.push_back("Residual_Adjoint_Velocity_z");
      Variable_Names.push_back("Residual_Adjoint_Temperature");
    } else {
      Variable_Names.push_back("Residual_Adjoint_Density");
      Variable_Names.push_back("Residual_Adjoint_Momentum_x");
      Variable_Names.push_back("Residual_Adjoint_Momentum_y");
      if (geometry->GetnDim() == 3)
        Variable_Names.push_back("Residual_Adjoint_Momentum_z");
      Variable_Names.push_back("Residual_Adjoint_Energy");
    }
    if (SecondIndex != NONE) {
      if ((config->GetKind_Turb_Model() == SST) || (config->GetKind_Turb_Model() == SST_SUST)) {
        Variable_Names.push_back("Residual_Adjoint_TKE");
        Variable_Names.push_back("Residual_Adjoint_Omega");
      } else {
        /*--- S-A variants ---*/
        Variable_Names.push_back("Residual_Adjoint_Nu_Tilde");
      }
    }
  }

  /*--- Add the grid velocity. ---*/

  if (dynamic_grid) {
    if (geometry->GetnDim() == 2) nVar_Par += 2;
    else if (geometry->GetnDim() == 3) nVar_Par += 3;
    Variable_Names.push_back("Grid_Velocity_x");
    Variable_Names.push_back("Grid_Velocity_y");
    if (geometry->GetnDim() == 3) Variable_Names.push_back("Grid_Velocity_z");
  }

  /*--- All adjoint solvers write the surface sensitivity. ---*/

  nVar_Par += 1; Variable_Names.push_back("Surface_Sensitivity");

  /*--- For the continouus adjoint, we write either convective scheme's
     dissipation sensor (centered) or limiter (uwpind) for adj. density. ---*/

  if (( Kind_Solver == ADJ_EULER              ) ||
      ( Kind_Solver == ADJ_NAVIER_STOKES      ) ||
      ( Kind_Solver == ADJ_RANS               )) {
    nVar_Par += 1;
    if (config->GetKind_ConvNumScheme() == SPACE_CENTERED) {
      Variable_Names.push_back("Dissipation_Sensor");
    } else {
      Variable_Names.push_back("Limiter_Adjoint_Density");
    }
  }

  /*--- New variables get registered here before the end of the loop. ---*/



  /*--- Auxiliary vectors for variables defined on surfaces only. ---*/

  Aux_Sens = new su2double[geometry->GetnPoint()];

  /*--- First, loop through the mesh in order to find and store the
   value of the viscous coefficients at any surface nodes. They
   will be placed in an auxiliary vector and then communicated like
   all other volumetric variables. ---*/

  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
    Aux_Sens[iPoint] = 0.0;
  }
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
    if (config->GetMarker_All_Plotting(iMarker) == YES) {
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
        Area = 0.0;
        for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim];
        Area = sqrt (Area);
        Aux_Sens[iPoint] = solver[ADJFLOW_SOL]->GetCSensitivity(iMarker, iVertex)/Area;
      }
    }

  /*--- Allocate the local data structure now that we know how many
   variables are in the output. ---*/

  Local_Data = new su2double*[geometry->GetnPoint()];
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
    Local_Data[iPoint] = new su2double[nVar_Par];
  }

  Local_Halo = new int[geometry->GetnPoint()];
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
    Local_Halo[iPoint] = !geometry->node[iPoint]->GetDomain();

  /*--- Search all send/recv boundaries on this partition for any periodic
   nodes that were part of the original domain. We want to recover these
   for visualization purposes. ---*/

  if (!Wrt_Halo) {
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) {

        /*--- Checking for less than or equal to the rank, because there may
         be some periodic halo nodes that send info to the same rank. ---*/

        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          isPeriodic = ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() > 0) &&
                        (geometry->vertex[iMarker][iVertex]->GetRotation_Type() % 2 == 1));
          if (isPeriodic) Local_Halo[iPoint] = false;
        }
      }
    }
  }

  /*--------------------------------------------------------------------------*/
  /*--- Step 2: Loop over all grid nodes and load up the desired data for  ---*/
  /*---         the restart and vizualization files. Note that we need to  ---*/
  /*---         increment the iVar variable after each variable load.      ---*/
  /*---         The idea is that we're filling up the columns of field     ---*/
  /*---         data for each iPoint (row) of the data structure. This     ---*/
  /*---         This data will then be sorted, communicated, and written   ---*/
  /*---         to files automatically after this routine. Note that the   ---*/
  /*---         ordering of the data loading MUST match the order of the   ---*/
  /*---         variable registration above for the files to be correct.   ---*/
  /*--------------------------------------------------------------------------*/

  jPoint = 0;

  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {

    /*--- Check for halos & write only if requested ---*/

    if (!Local_Halo[iPoint] || Wrt_Halo) {

      /*--- Restart the column index with each new point. ---*/

      iVar = 0;

      /*--- Load the grid node coordinate values. ---*/

      for (iDim = 0; iDim < geometry->GetnDim(); iDim++) {
        Local_Data[jPoint][iVar] = geometry->node[iPoint]->GetCoord(iDim);
        if (config->GetSystemMeasurements() == US)
        	Local_Data[jPoint][iVar] *= 12.0;
        iVar++;
      }

      /*--- Load the conservative variable states for the mean flow variables.
       If requested, load the limiters and residuals as well. ---*/

      for (jVar = 0; jVar < nVar_First; jVar++) {
        Local_Data[jPoint][iVar] = solver[FirstIndex]->GetNodes()->GetSolution(iPoint, jVar);
        iVar++;
      }


      /*--- If this is Adj. RANS, i.e., the second solver container is not empty,
       then load data for the conservative turbulence variables and the
       limiters / residuals (if requested). ----*/

      if (SecondIndex != NONE) {
        for (jVar = 0; jVar < nVar_Second; jVar++) {
          Local_Data[jPoint][iVar] = solver[SecondIndex]->GetNodes()->GetSolution(iPoint, jVar);
          iVar++;
        }
      }

      /*--- Load data for the discrete sensitivities. ---*/

      if ((Kind_Solver == DISC_ADJ_EULER)         ||
          (Kind_Solver == DISC_ADJ_NAVIER_STOKES) ||
          (Kind_Solver == DISC_ADJ_RANS) ||
          (Kind_Solver == DISC_ADJ_INC_EULER)         ||
          (Kind_Solver == DISC_ADJ_INC_NAVIER_STOKES) ||
          (Kind_Solver == DISC_ADJ_INC_RANS)) {
        Local_Data[jPoint][iVar] = solver[ADJFLOW_SOL]->GetNodes()->GetSensitivity(iPoint, 0); iVar++;
        Local_Data[jPoint][iVar] = solver[ADJFLOW_SOL]->GetNodes()->GetSensitivity(iPoint, 1); iVar++;
        if (geometry->GetnDim()== 3) {
          Local_Data[jPoint][iVar] = solver[ADJFLOW_SOL]->GetNodes()->GetSensitivity(iPoint, 2);
          iVar++;
        }
      }

      if (config->GetWrt_Limiters()) {
        for (jVar = 0; jVar < nVar_First; jVar++) {
          Local_Data[jPoint][iVar] = solver[FirstIndex]->GetNodes()->GetLimiter(iPoint, jVar);
          iVar++;
        }
        if (SecondIndex != NONE) {
          for (jVar = 0; jVar < nVar_Second; jVar++) {
            Local_Data[jPoint][iVar] = solver[SecondIndex]->GetNodes()->GetLimiter(iPoint, jVar);
            iVar++;
          }
        }
      }

      if (config->GetWrt_Residuals()) {
        for (jVar = 0; jVar < nVar_First; jVar++) {
          if (!config->GetDiscrete_Adjoint()) {
            Local_Data[jPoint][iVar] = solver[FirstIndex]->LinSysRes(iPoint, jVar);
          } else {
            Local_Data[jPoint][iVar] = solver[FirstIndex]->GetNodes()->GetSolution(iPoint, jVar) -
                                       solver[FirstIndex]->GetNodes()->GetSolution_Old(iPoint, jVar);
          }
          iVar++;
        }
        if (SecondIndex != NONE) {
          for (jVar = 0; jVar < nVar_Second; jVar++) {
            if (!config->GetDiscrete_Adjoint()) {
              Local_Data[jPoint][iVar] = solver[SecondIndex]->LinSysRes(iPoint, jVar);
            } else {
              Local_Data[jPoint][iVar] = solver[SecondIndex]->GetNodes()->GetSolution(iPoint, jVar) -
                                         solver[SecondIndex]->GetNodes()->GetSolution_Old(iPoint, jVar);
            }
            iVar++;
          }
        }
      }

      /*--- Load buffers with the three grid velocity components. ---*/

      if (dynamic_grid) {
        Grid_Vel = geometry->node[iPoint]->GetGridVel();
        Local_Data[jPoint][iVar] = Grid_Vel[0]; iVar++;
        Local_Data[jPoint][iVar] = Grid_Vel[1]; iVar++;
        if (geometry->GetnDim() == 3) {
          Local_Data[jPoint][iVar] = Grid_Vel[2];
          iVar++;
        }
      }

      /*--- Load data for the surface sensitivity. ---*/

      Local_Data[iPoint][iVar] = Aux_Sens[iPoint]; iVar++;

      /*--- Load data for the convective scheme sensor. ---*/

      if (( Kind_Solver == ADJ_EULER              ) ||
          ( Kind_Solver == ADJ_NAVIER_STOKES      ) ||
          ( Kind_Solver == ADJ_RANS               )) {
        if (config->GetKind_ConvNumScheme() == SPACE_CENTERED) {
          Local_Data[jPoint][iVar] = solver[ADJFLOW_SOL]->GetNodes()->GetSensor(iPoint); iVar++;
        } else {
          Local_Data[jPoint][iVar] = solver[ADJFLOW_SOL]->GetNodes()->GetLimiter(iPoint, 0); iVar++;
        }
      }

      /*--- New variables can be loaded to the Local_Data structure here,
         assuming they were registered above correctly. ---*/

    }

    /*--- Increment the point counter, as there may have been halos we
       skipped over during the data loading. ---*/

    jPoint++;
  }


  /*--- Free memory for auxiliary vectors. ---*/

  delete [] Aux_Sens;
  delete [] Local_Halo;

}

void COutputLegacy::LoadLocalData_Elasticity(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned short val_iZone) {

  unsigned short iDim;

  unsigned long iVar, jVar;
  unsigned long iPoint, jPoint, FirstIndex = NONE, iMarker, iVertex;
  unsigned long nVar_First = 0, nVar_Consv_Par = 0;

  su2double *Node_Vel = NULL, *Node_Accel = NULL;
  const su2double *Stress = NULL;

  bool Wrt_Halo   = config->GetWrt_Halo(), isPeriodic;
  int *Local_Halo = NULL;

  stringstream varname;

  /*--- Use a switch statement to decide how many solver containers we have
   in this zone for output. ---*/

  switch (config->GetKind_Solver()) {
    case FEM_ELASTICITY: FirstIndex = FEA_SOL; break;
    case DISC_ADJ_FEM: FirstIndex = ADJFEA_SOL; break;
  }

  nVar_First = solver[FirstIndex]->GetnVar();
  nVar_Consv_Par = nVar_First;

  /*--------------------------------------------------------------------------*/
  /*--- Step 1: Register the variables that will be output. To register a  ---*/
  /*---         variable, two things are required. First, increment the    ---*/
  /*---         counter for the number of variables (nVar_Par), which      ---*/
  /*---         controls the size of the data structure allocation, i.e.,  ---*/
  /*---         the number of columns in an nPoint x nVar structure.       ---*/
  /*---         Second, add a name for the variable to the vector that     ---*/
  /*---         holds the string names.                                    ---*/
  /*--------------------------------------------------------------------------*/

  /*--- All output files first need the grid coordinates. ---*/

  nVar_Par  = 1; Variable_Names.push_back("x");
  nVar_Par += 1; Variable_Names.push_back("y");
  if (geometry->GetnDim() == 3) {
    nVar_Par += 1; Variable_Names.push_back("z");
  }

  /*--- At a mininum, the restarts and visualization files need the
   conservative variables, so these follow next. ---*/

  nVar_Par += nVar_Consv_Par;

  /*--- For now, leave the names as "Conservative_", etc., in order
   to avoid confusion with the serial version, which still prints these
   names. Names can be set alternatively by using the commented code
   below. ---*/

  Variable_Names.push_back("Displacement_x");
  Variable_Names.push_back("Displacement_y");
  if (geometry->GetnDim() == 3)
    Variable_Names.push_back("Displacement_z");

  /*--- If requested, register the limiter and residuals for all of the
   equations in the current flow problem. ---*/


  /*--- Add the residuals ---*/

  if (config->GetWrt_Residuals()) {
    nVar_Par += nVar_Consv_Par;
    Variable_Names.push_back("Residual_Displacement_x");
    Variable_Names.push_back("Residual_Displacement_y");
    if (geometry->GetnDim() == 3)
      Variable_Names.push_back("Residual_Displacement_z");
  }

  /*--- If the analysis is dynamic... ---*/
  if (config->GetTime_Domain()) {

    /*--- Velocities ---*/
    nVar_Par += 2;
    Variable_Names.push_back("Velocity_x");
    Variable_Names.push_back("Velocity_y");
    if (geometry->GetnDim() == 3) {
      nVar_Par += 1;
      Variable_Names.push_back("Velocity_z");
    }

    /*--- Accelerations ---*/
    nVar_Par += 2;
    Variable_Names.push_back("Acceleration_x");
    Variable_Names.push_back("Acceleration_y");
    if (geometry->GetnDim() == 3) {
      nVar_Par += 1;
      Variable_Names.push_back("Acceleration_z");
    }
  }

  if (!(config->GetDiscrete_Adjoint())) {

    /*--- Add the stresses. ---*/

    nVar_Par += 3;
    Variable_Names.push_back("Sxx");
    Variable_Names.push_back("Syy");
    Variable_Names.push_back("Sxy");
    if (geometry->GetnDim() == 3) {
      nVar_Par += 3;
      Variable_Names.push_back("Szz");
      Variable_Names.push_back("Sxz");
      Variable_Names.push_back("Syz");
    }

    /*--- Add the Von Mises Stress. ---*/

    nVar_Par += 1;
    Variable_Names.push_back("Von_Mises_Stress");

  }

  /*--- New variables get registered here before the end of the loop. ---*/



  /*--- Allocate the local data structure now that we know how many
   variables are in the output. ---*/

  Local_Data = new su2double*[geometry->GetnPoint()];
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
    Local_Data[iPoint] = new su2double[nVar_Par];
  }

  Local_Halo = new int[geometry->GetnPoint()];
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
    Local_Halo[iPoint] = !geometry->node[iPoint]->GetDomain();

  /*--- Search all send/recv boundaries on this partition for any periodic
   nodes that were part of the original domain. We want to recover these
   for visualization purposes. ---*/

  if (!Wrt_Halo) {
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) {

        /*--- Checking for less than or equal to the rank, because there may
         be some periodic halo nodes that send info to the same rank. ---*/

        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          isPeriodic = ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() > 0) &&
                        (geometry->vertex[iMarker][iVertex]->GetRotation_Type() % 2 == 1));
          if (isPeriodic) Local_Halo[iPoint] = false;
        }
      }
    }
  }

  /*--------------------------------------------------------------------------*/
  /*--- Step 2: Loop over all grid nodes and load up the desired data for  ---*/
  /*---         the restart and vizualization files. Note that we need to  ---*/
  /*---         increment the iVar variable after each variable load.      ---*/
  /*---         The idea is that we're filling up the columns of field     ---*/
  /*---         data for each iPoint (row) of the data structure. This     ---*/
  /*---         This data will then be sorted, communicated, and written   ---*/
  /*---         to files automatically after this routine. Note that the   ---*/
  /*---         ordering of the data loading MUST match the order of the   ---*/
  /*---         variable registration above for the files to be correct.   ---*/
  /*--------------------------------------------------------------------------*/

  jPoint = 0;

  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {

    /*--- Check for halos & write only if requested ---*/

    if (!Local_Halo[iPoint] || Wrt_Halo) {

      /*--- Restart the column index with each new point. ---*/

      iVar = 0;

      /*--- Load the grid node coordinate values. ---*/

      for (iDim = 0; iDim < geometry->GetnDim(); iDim++) {
        Local_Data[jPoint][iVar] = geometry->node[iPoint]->GetCoord(iDim);
        if (config->GetSystemMeasurements() == US)
        	Local_Data[jPoint][iVar] *= 12.0;
        iVar++;
      }

      /*--- Load the conservative variable states for the mean flow variables.
       If requested, load the limiters and residuals as well. ---*/

      for (jVar = 0; jVar < nVar_First; jVar++) {
        Local_Data[jPoint][iVar] = solver[FirstIndex]->GetNodes()->GetSolution(iPoint, jVar);
        iVar++;
      }

      if (config->GetWrt_Residuals()) {
        for (jVar = 0; jVar < nVar_First; jVar++) {
          Local_Data[jPoint][iVar] = solver[FirstIndex]->LinSysRes(iPoint, jVar);
          iVar++;
        }
      }



      /*--- Load the velocities and accelerations (dynamic calculations). ---*/

      if (config->GetTime_Domain()) {

        /*--- Velocities ---*/

        Node_Vel = solver[FEA_SOL]->GetNodes()->GetSolution_Vel(iPoint);
        Local_Data[jPoint][iVar] = Node_Vel[0]; iVar++;
        Local_Data[jPoint][iVar] = Node_Vel[1]; iVar++;
        if (geometry->GetnDim() == 3) {
          Local_Data[jPoint][iVar] = Node_Vel[2];
          iVar++;
        }

        /*--- Accelerations ---*/

        Node_Accel = solver[FEA_SOL]->GetNodes()->GetSolution_Accel(iPoint);
        Local_Data[jPoint][iVar] = Node_Accel[0]; iVar++;
        Local_Data[jPoint][iVar] = Node_Accel[1]; iVar++;
        if (geometry->GetnDim() == 3) {
          Local_Data[jPoint][iVar] = Node_Accel[2];
          iVar++;
        }
      }

      if (!(config->GetDiscrete_Adjoint())) {

        /*--- Add the stresses. ---*/

        Stress = solver[FEA_SOL]->GetNodes()->GetStress_FEM(iPoint);

        /*--- Sigma xx ---*/
        Local_Data[jPoint][iVar] = Stress[0]; iVar++;
        /*--- Sigma yy ---*/
        Local_Data[jPoint][iVar] = Stress[1]; iVar++;
        /*--- Sigma xy ---*/
        Local_Data[jPoint][iVar] = Stress[2]; iVar++;

        if (geometry->GetnDim() == 3) {
          /*--- Sigma zz ---*/
          Local_Data[jPoint][iVar] = Stress[3]; iVar++;
          /*--- Sigma xz ---*/
          Local_Data[jPoint][iVar] = Stress[4]; iVar++;
          /*--- Sigma yz ---*/
          Local_Data[jPoint][iVar] = Stress[5]; iVar++;
        }

        /*--- Add the Von Mises Stress. ---*/

        Local_Data[iPoint][iVar] = solver[FEA_SOL]->GetNodes()->GetVonMises_Stress(iPoint); iVar++;


        /*--- New variables can be loaded to the Local_Data structure here,
           assuming they were registered above correctly. ---*/

      }



      /*--- Increment the point counter, as there may have been halos we
       skipped over during the data loading. ---*/

      jPoint++;
    }
  }

  /*--- Free memory for auxiliary vectors. ---*/

  delete [] Local_Halo;

}

void COutputLegacy::LoadLocalData_Base(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned short val_iZone) {

  unsigned short iDim;

  unsigned long iVar, jVar;
  unsigned long iPoint, jPoint, FirstIndex = NONE, iMarker, iVertex;
  unsigned long nVar_First = 0, nVar_Consv_Par = 0;

  bool Wrt_Halo = config->GetWrt_Halo(), isPeriodic;

  int *Local_Halo;

  stringstream varname;

  /*--- Use a switch statement to decide how many solver containers we have
   in this zone for output. ---*/

  switch (config->GetKind_Solver()) {
    case HEAT_EQUATION: FirstIndex = HEAT_SOL;     break;
  }

  nVar_First = solver[FirstIndex]->GetnVar();
  nVar_Consv_Par = nVar_First;

  /*--------------------------------------------------------------------------*/
  /*--- Step 1: Register the variables that will be output. To register a  ---*/
  /*---         variable, two things are required. First, increment the    ---*/
  /*---         counter for the number of variables (nVar_Par), which      ---*/
  /*---         controls the size of the data structure allocation, i.e.,  ---*/
  /*---         the number of columns in an nPoint x nVar structure.       ---*/
  /*---         Second, add a name for the variable to the vector that     ---*/
  /*---         holds the string names.                                    ---*/
  /*--------------------------------------------------------------------------*/

  /*--- All output files first need the grid coordinates. ---*/

  nVar_Par  = 1; Variable_Names.push_back("x");
  nVar_Par += 1; Variable_Names.push_back("y");
  if (geometry->GetnDim() == 3) {
    nVar_Par += 1; Variable_Names.push_back("z");
  }

  /*--- At a mininum, the restarts and visualization files need the
   conservative variables, so these follow next. ---*/

  nVar_Par += nVar_Consv_Par;
  for (iVar = 0; iVar < nVar_Consv_Par; iVar++) {
    varname << "Conservative_" << iVar+1;
    Variable_Names.push_back(varname.str());
    varname.str("");
  }

  /*--- If requested, register the residuals for all of the
   equations in the current problem. ---*/


  /*--- Add the residuals ---*/

  if (config->GetWrt_Residuals()) {
    nVar_Par += nVar_Consv_Par;
    for (iVar = 0; iVar < nVar_Consv_Par; iVar++) {
      varname << "Residual_" << iVar+1;
      Variable_Names.push_back(varname.str());
      varname.str("");
    }
  }

  /*--- New variables get registered here before the end of the loop. ---*/



  /*--- Allocate the local data structure now that we know how many
   variables are in the output. ---*/

  Local_Data = new su2double*[geometry->GetnPoint()];
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
    Local_Data[iPoint] = new su2double[nVar_Par];
  }

  Local_Halo = new int[geometry->GetnPoint()];
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
    Local_Halo[iPoint] = !geometry->node[iPoint]->GetDomain();

  /*--- Search all send/recv boundaries on this partition for any periodic
   nodes that were part of the original domain. We want to recover these
   for visualization purposes. ---*/

  if (!Wrt_Halo) {
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) {

        /*--- Checking for less than or equal to the rank, because there may
         be some periodic halo nodes that send info to the same rank. ---*/

        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          isPeriodic = ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() > 0) &&
                        (geometry->vertex[iMarker][iVertex]->GetRotation_Type() % 2 == 1));
          if (isPeriodic) Local_Halo[iPoint] = false;
        }
      }
    }
  }

  /*--------------------------------------------------------------------------*/
  /*--- Step 2: Loop over all grid nodes and load up the desired data for  ---*/
  /*---         the restart and vizualization files. Note that we need to  ---*/
  /*---         increment the iVar variable after each variable load.      ---*/
  /*---         The idea is that we're filling up the columns of field     ---*/
  /*---         data for each iPoint (row) of the data structure. This     ---*/
  /*---         This data will then be sorted, communicated, and written   ---*/
  /*---         to files automatically after this routine. Note that the   ---*/
  /*---         ordering of the data loading MUST match the order of the   ---*/
  /*---         variable registration above for the files to be correct.   ---*/
  /*--------------------------------------------------------------------------*/

  jPoint = 0;

  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {

    /*--- Check for halos & write only if requested ---*/

    if (!Local_Halo[iPoint] || Wrt_Halo) {

      /*--- Restart the column index with each new point. ---*/

      iVar = 0;

      /*--- Load the grid node coordinate values. ---*/

      for (iDim = 0; iDim < geometry->GetnDim(); iDim++) {
        Local_Data[jPoint][iVar] = geometry->node[iPoint]->GetCoord(iDim);
        if (config->GetSystemMeasurements() == US)
        	Local_Data[jPoint][iVar] *= 12.0;
        iVar++;
      }

      /*--- Load the conservative variable states for the mean flow variables.
       If requested, load the limiters and residuals as well. ---*/

      for (jVar = 0; jVar < nVar_First; jVar++) {
        Local_Data[jPoint][iVar] = solver[FirstIndex]->GetNodes()->GetSolution(iPoint, jVar);
        iVar++;
      }

      if (config->GetWrt_Residuals()) {
        for (jVar = 0; jVar < nVar_First; jVar++) {
          Local_Data[jPoint][iVar] = solver[FirstIndex]->LinSysRes(iPoint, jVar);
          iVar++;
        }
      }


      /*--- New variables can be loaded to the Local_Data structure here,
       assuming they were registered above correctly. ---*/

      /*--- Increment the point counter, as there may have been halos we
       skipped over during the data loading. ---*/

      jPoint++;

    }
  }

  /*--- Free memory for auxiliary vectors. ---*/

  delete [] Local_Halo;

}

void COutputLegacy::SortConnectivity(CConfig *config, CGeometry *geometry, unsigned short val_iZone) {

  /*--- Flags identifying the types of files to be written. ---*/

  bool Wrt_Vol = config->GetWrt_Vol_Sol();
  bool Wrt_Srf = config->GetWrt_Srf_Sol();

  /*--- Sort connectivity for each type of element (excluding halos). Note
   In these routines, we sort the connectivity into a linear partitioning
   across all processors based on the global index of the grid nodes. ---*/

  /*--- Sort volumetric grid connectivity. ---*/

  if (Wrt_Vol) {

    if ((rank == MASTER_NODE) && (size != SINGLE_NODE))
      cout <<"Sorting volumetric grid connectivity." << endl;

    SortVolumetricConnectivity(config, geometry, TRIANGLE     );
    SortVolumetricConnectivity(config, geometry, QUADRILATERAL);
    SortVolumetricConnectivity(config, geometry, TETRAHEDRON  );
    SortVolumetricConnectivity(config, geometry, HEXAHEDRON   );
    SortVolumetricConnectivity(config, geometry, PRISM        );
    SortVolumetricConnectivity(config, geometry, PYRAMID      );

  }

  /*--- Sort surface grid connectivity. ---*/

  if (Wrt_Srf) {

    if ((rank == MASTER_NODE) && (size != SINGLE_NODE))
      cout <<"Sorting surface grid connectivity." << endl;

    SortSurfaceConnectivity(config, geometry, LINE         );
    SortSurfaceConnectivity(config, geometry, TRIANGLE     );
    SortSurfaceConnectivity(config, geometry, QUADRILATERAL);

  }

  /*--- Reduce the total number of cells we will be writing in the output files. ---*/

  unsigned long nTotal_Elem = nParallel_Tria + nParallel_Quad + nParallel_Tetr + nParallel_Hexa + nParallel_Pris + nParallel_Pyra;
  unsigned long nTotal_Surf_Elem = nParallel_Line + nParallel_BoundTria + nParallel_BoundQuad;
#ifndef HAVE_MPI
  nGlobal_Elem_Par = nTotal_Elem;
  nSurf_Elem_Par   = nTotal_Surf_Elem;
#else
  SU2_MPI::Allreduce(&nTotal_Elem, &nGlobal_Elem_Par, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&nTotal_Surf_Elem, &nSurf_Elem_Par, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#endif

}

void COutputLegacy::SortVolumetricConnectivity(CConfig *config, CGeometry *geometry, unsigned short Elem_Type) {

  unsigned long iProcessor;
  unsigned short NODES_PER_ELEMENT = 0;
  unsigned long iPoint, jPoint, kPoint, nLocalPoint, nTotalPoint;
  unsigned long nElem_Total = 0, Global_Index;

  unsigned long iVertex, iMarker;
  int SendRecv, RecvFrom;

  bool notPeriodic, notHalo, addedPeriodic, isPeriodic;

  int *Local_Halo = NULL;
  int *Conn_Elem  = NULL;

#ifdef HAVE_MPI
  SU2_MPI::Request *send_req, *recv_req;
  SU2_MPI::Status status;
  int ind;
#endif

  /*--- Store the local number of this element type and the number of nodes
   per this element type. In serial, this will be the total number of this
   element type in the entire mesh. In parallel, it is the number on only
   the current partition. ---*/

  switch (Elem_Type) {
    case TRIANGLE:
      NODES_PER_ELEMENT = N_POINTS_TRIANGLE;
      break;
    case QUADRILATERAL:
      NODES_PER_ELEMENT = N_POINTS_QUADRILATERAL;
      break;
    case TETRAHEDRON:
      NODES_PER_ELEMENT = N_POINTS_TETRAHEDRON;
      break;
    case HEXAHEDRON:
      NODES_PER_ELEMENT = N_POINTS_HEXAHEDRON;
      break;
    case PRISM:
      NODES_PER_ELEMENT = N_POINTS_PRISM;
      break;
    case PYRAMID:
      NODES_PER_ELEMENT = N_POINTS_PYRAMID;
      break;
    default:
      SU2_MPI::Error("Unrecognized element type", CURRENT_FUNCTION);
  }

  /*--- Force the removal of all added periodic elements (use global index).
   First, we isolate and create a list of all added periodic points, excluding
   those that were part of the original domain (we want these to be in the
   output files). ---*/

  vector<unsigned long> Added_Periodic;
  Added_Periodic.clear();

  if (config->GetKind_SU2() != SU2_DEF) {
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) {
        SendRecv = config->GetMarker_All_SendRecv(iMarker);
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

          if ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() > 0) &&
              (geometry->vertex[iMarker][iVertex]->GetRotation_Type() % 2 == 0) &&
              (SendRecv < 0)) {
            Added_Periodic.push_back(geometry->node[iPoint]->GetGlobalIndex());
          }
        }
      }
    }
  }

  /*--- Now we communicate this information to all processors, so that they
   can force the removal of these particular nodes by flagging them as halo
   points. In general, this should be a small percentage of the total mesh,
   so the communication/storage costs here shouldn't be prohibitive. ---*/

  /*--- First communicate the number of points that each rank has found. ---*/

  unsigned long nAddedPeriodic = 0, maxAddedPeriodic = 0;
  unsigned long Buffer_Send_nAddedPeriodic[1], *Buffer_Recv_nAddedPeriodic = NULL;
  Buffer_Recv_nAddedPeriodic = new unsigned long[size];

  nAddedPeriodic = Added_Periodic.size();
  Buffer_Send_nAddedPeriodic[0] = nAddedPeriodic;

#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&nAddedPeriodic, &maxAddedPeriodic, 1, MPI_UNSIGNED_LONG,
                     MPI_MAX, MPI_COMM_WORLD);
  SU2_MPI::Allgather(&Buffer_Send_nAddedPeriodic, 1, MPI_UNSIGNED_LONG,
                     Buffer_Recv_nAddedPeriodic,  1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
#else
  maxAddedPeriodic = nAddedPeriodic;
  Buffer_Recv_nAddedPeriodic[0] = Buffer_Send_nAddedPeriodic[0];
#endif

  /*--- Communicate the global index values of all added periodic nodes. ---*/
  unsigned long *Buffer_Send_AddedPeriodic = new unsigned long[maxAddedPeriodic];
  unsigned long *Buffer_Recv_AddedPeriodic = new unsigned long[size*maxAddedPeriodic];

  for (iPoint = 0; iPoint < Added_Periodic.size(); iPoint++) {
    Buffer_Send_AddedPeriodic[iPoint] = Added_Periodic[iPoint];
  }

  /*--- Gather the element connectivity information. All processors will now
   have a copy of the global index values for all added periodic points. ---*/

#ifdef HAVE_MPI
  SU2_MPI::Allgather(Buffer_Send_AddedPeriodic, maxAddedPeriodic, MPI_UNSIGNED_LONG,
                     Buffer_Recv_AddedPeriodic, maxAddedPeriodic, MPI_UNSIGNED_LONG,
                     MPI_COMM_WORLD);
#else
  for (iPoint = 0; iPoint < maxAddedPeriodic; iPoint++)
    Buffer_Recv_AddedPeriodic[iPoint] = Buffer_Send_AddedPeriodic[iPoint];
#endif

  /*--- Search all send/recv boundaries on this partition for halo cells. In
   particular, consider only the recv conditions (these are the true halo
   nodes). Check the ranks of the processors that are communicating and
   choose to keep only the halo cells from the higher rank processor. Here,
   we are also choosing to keep periodic nodes that were part of the original
   domain. We will check the communicated list of added periodic points. ---*/

  Local_Halo = new int[geometry->GetnPoint()];
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
    Local_Halo[iPoint] = !geometry->node[iPoint]->GetDomain();

  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) {
      SendRecv = config->GetMarker_All_SendRecv(iMarker);
      RecvFrom = abs(SendRecv)-1;

      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        Global_Index = geometry->node[iPoint]->GetGlobalIndex();

        /*--- We need to keep one copy of overlapping halo cells. ---*/

        notHalo = ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() == 0) &&
                   (SendRecv < 0) && (rank > RecvFrom));

        /*--- We want to keep the periodic nodes that were part of the original domain.
         For SU2_DEF we want to keep all periodic nodes. ---*/

        if (config->GetKind_SU2() == SU2_DEF) {
          isPeriodic = ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() > 0));
        }else {
          isPeriodic = ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() > 0) &&
                        (geometry->vertex[iMarker][iVertex]->GetRotation_Type() % 2 == 1));
        }

        notPeriodic = (isPeriodic && (SendRecv < 0));

        /*--- Lastly, check that this isn't an added periodic point that
         we will forcibly remove. Use the communicated list of these points. ---*/

        addedPeriodic = false; kPoint = 0;
        for (iProcessor = 0; iProcessor < (unsigned long)size; iProcessor++) {
          for (jPoint = 0; jPoint < Buffer_Recv_nAddedPeriodic[iProcessor]; jPoint++) {
            if (Global_Index == Buffer_Recv_AddedPeriodic[kPoint+jPoint])
              addedPeriodic = true;
          }

          /*--- Adjust jNode to index of next proc's data in the buffers. ---*/

          kPoint = (iProcessor+1)*maxAddedPeriodic;

        }

        /*--- If we found either of these types of nodes, flag them to be kept. ---*/

        if ((notHalo || notPeriodic) && !addedPeriodic) {
          Local_Halo[iPoint] = false;
        }

      }
    }
  }

  /*--- Now that we've done the gymnastics to find any periodic points,
   compute the total number of local and global points for the output. ---*/

  nLocalPoint = 0;
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
      if (Local_Halo[iPoint] == false)
        nLocalPoint++;

#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&nLocalPoint, &nTotalPoint, 1,
                     MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#else
  nTotalPoint = nLocalPoint;
#endif

  /*--- Compute the number of points that will be on each processor.
   This is a linear partitioning with the addition of a simple load
   balancing for any remainder points. ---*/

  unsigned long *npoint_procs  = new unsigned long[size];
  unsigned long *starting_node = new unsigned long[size];
  unsigned long *ending_node   = new unsigned long[size];
  unsigned long *nPoint_Linear = new unsigned long[size+1];

  unsigned long total_pt_accounted = 0;
  for (int ii = 0; ii < size; ii++) {
    npoint_procs[ii] = nTotalPoint/size;
    total_pt_accounted = total_pt_accounted + npoint_procs[ii];
  }

  /*--- Get the number of remainder points after the even division. ---*/

  unsigned long rem_points = nTotalPoint-total_pt_accounted;
  for (unsigned long ii = 0; ii < rem_points; ii++) {
    npoint_procs[ii]++;
  }

  /*--- Store the local number of nodes and the beginning/end index ---*/

  starting_node[0] = 0;
  ending_node[0]   = starting_node[0] + npoint_procs[0];
  nPoint_Linear[0] = 0;
  for (int ii = 1; ii < size; ii++) {
    starting_node[ii] = ending_node[ii-1];
    ending_node[ii]   = starting_node[ii] + npoint_procs[ii];
    nPoint_Linear[ii] = nPoint_Linear[ii-1] + npoint_procs[ii-1];
  }
  nPoint_Linear[size] = nTotalPoint;

  /*--- We start with the connectivity distributed across all procs with
   no particular ordering assumed. We need to loop through our local partition
   and decide how many elements we must send to each other rank in order to
   have all elements sorted according to a linear partitioning of the grid
   nodes, i.e., rank 0 holds the first nPoint()/nProcessors nodes.
   First, initialize a counter and flag. ---*/

  int *nElem_Send = new int[size+1]; nElem_Send[0] = 0;
  int *nElem_Recv = new int[size+1]; nElem_Recv[0] = 0;
  int *nElem_Flag = new int[size];

  for (int ii=0; ii < size; ii++) {
    nElem_Send[ii] = 0;
    nElem_Recv[ii] = 0;
    nElem_Flag[ii]= -1;
  }
  nElem_Send[size] = 0; nElem_Recv[size] = 0;

  for (int ii = 0; ii < (int)geometry->GetnElem(); ii++ ) {
    if (geometry->elem[ii]->GetVTK_Type() == Elem_Type) {
      for ( int jj = 0; jj < NODES_PER_ELEMENT; jj++ ) {

        /*--- Get the index of the current point. ---*/

        iPoint = geometry->elem[ii]->GetNode(jj);
        Global_Index = geometry->node[iPoint]->GetGlobalIndex();

        /*--- Search for the lowest global index in this element. We
         send the element to the processor owning the range that includes
         the lowest global index value. ---*/

        for (int kk = 0; kk < NODES_PER_ELEMENT; kk++) {
          jPoint = geometry->elem[ii]->GetNode(kk);
          unsigned long newID = geometry->node[jPoint]->GetGlobalIndex();
          if (newID < Global_Index) Global_Index = newID;
        }

        /*--- Search for the processor that owns this point ---*/

        iProcessor = Global_Index/npoint_procs[0];
        if (iProcessor >= (unsigned long)size)
          iProcessor = (unsigned long)size-1;
        if (Global_Index >= nPoint_Linear[iProcessor])
          while(Global_Index >= nPoint_Linear[iProcessor+1]) iProcessor++;
        else
          while(Global_Index <  nPoint_Linear[iProcessor])   iProcessor--;

        /*--- If we have not visited this element yet, increment our
         number of elements that must be sent to a particular proc. ---*/

        if ((nElem_Flag[iProcessor] != ii)) {
          nElem_Flag[iProcessor] = ii;
          nElem_Send[iProcessor+1]++;
        }

      }
    }
  }

  /*--- Communicate the number of cells to be sent/recv'd amongst
   all processors. After this communication, each proc knows how
   many cells it will receive from each other processor. ---*/

#ifdef HAVE_MPI
  SU2_MPI::Alltoall(&(nElem_Send[1]), 1, MPI_INT,
                    &(nElem_Recv[1]), 1, MPI_INT, MPI_COMM_WORLD);
#else
  nElem_Recv[1] = nElem_Send[1];
#endif

  /*--- Prepare to send connectivities. First check how many
   messages we will be sending and receiving. Here we also put
   the counters into cumulative storage format to make the
   communications simpler. ---*/

  int nSends = 0, nRecvs = 0;
  for (int ii=0; ii < size; ii++) nElem_Flag[ii] = -1;

  for (int ii = 0; ii < size; ii++) {
    if ((ii != rank) && (nElem_Send[ii+1] > 0)) nSends++;
    if ((ii != rank) && (nElem_Recv[ii+1] > 0)) nRecvs++;

    nElem_Send[ii+1] += nElem_Send[ii];
    nElem_Recv[ii+1] += nElem_Recv[ii];
  }

  /*--- Allocate memory to hold the connectivity that we are
   sending. ---*/

  unsigned long *connSend = NULL;
  connSend = new unsigned long[NODES_PER_ELEMENT*nElem_Send[size]];
  for (int ii = 0; ii < NODES_PER_ELEMENT*nElem_Send[size]; ii++)
    connSend[ii] = 0;

  /*--- Allocate arrays for storing halo flags. ---*/

  unsigned short *haloSend = new unsigned short[nElem_Send[size]];
  for (int ii = 0; ii < nElem_Send[size]; ii++)
    haloSend[ii] = false;

  /*--- Create an index variable to keep track of our index
   position as we load up the send buffer. ---*/

  unsigned long *index = new unsigned long[size];
  for (int ii=0; ii < size; ii++) index[ii] = NODES_PER_ELEMENT*nElem_Send[ii];

  unsigned long *haloIndex = new unsigned long[size];
  for (int ii=0; ii < size; ii++) haloIndex[ii] = nElem_Send[ii];

  /*--- Loop through our elements and load the elems and their
   additional data that we will send to the other procs. ---*/

  for (int ii = 0; ii < (int)geometry->GetnElem(); ii++) {
    if (geometry->elem[ii]->GetVTK_Type() == Elem_Type) {
      for ( int jj = 0; jj < NODES_PER_ELEMENT; jj++ ) {

        /*--- Get the index of the current point. ---*/

        iPoint = geometry->elem[ii]->GetNode(jj);
        Global_Index = geometry->node[iPoint]->GetGlobalIndex();

        /*--- Search for the lowest global index in this element. We
         send the element to the processor owning the range that includes
         the lowest global index value. ---*/

        for (int kk = 0; kk < NODES_PER_ELEMENT; kk++) {
          jPoint = geometry->elem[ii]->GetNode(kk);
          unsigned long newID = geometry->node[jPoint]->GetGlobalIndex();
          if (newID < Global_Index) Global_Index = newID;
        }

        /*--- Search for the processor that owns this point ---*/

        iProcessor = Global_Index/npoint_procs[0];
        if (iProcessor >= (unsigned long)size)
          iProcessor = (unsigned long)size-1;
        if (Global_Index >= nPoint_Linear[iProcessor])
          while(Global_Index >= nPoint_Linear[iProcessor+1]) iProcessor++;
        else
          while(Global_Index <  nPoint_Linear[iProcessor])   iProcessor--;

        /*--- Load connectivity into the buffer for sending ---*/

        if (nElem_Flag[iProcessor] != ii) {

          nElem_Flag[iProcessor] = ii;
          unsigned long nn = index[iProcessor];
          unsigned long mm = haloIndex[iProcessor];

          /*--- Load the connectivity values. ---*/

          for (int kk = 0; kk < NODES_PER_ELEMENT; kk++) {
            iPoint = geometry->elem[ii]->GetNode(kk);
            connSend[nn] = geometry->node[iPoint]->GetGlobalIndex(); nn++;

            /*--- Check if this is a halo node. If so, flag this element
             as a halo cell. We will use this later to sort and remove
             any duplicates from the connectivity list. ---*/

            if (Local_Halo[iPoint]) haloSend[mm] = true;

          }

          /*--- Increment the index by the message length ---*/

          index[iProcessor]    += NODES_PER_ELEMENT;
          haloIndex[iProcessor]++;

        }
      }
    }
  }

  /*--- Free memory after loading up the send buffer. ---*/

  delete [] index;
  delete [] haloIndex;

  /*--- Allocate the memory that we need for receiving the conn
   values and then cue up the non-blocking receives. Note that
   we do not include our own rank in the communications. We will
   directly copy our own data later. ---*/

  unsigned long *connRecv = NULL;
  connRecv = new unsigned long[NODES_PER_ELEMENT*nElem_Recv[size]];
  for (int ii = 0; ii < NODES_PER_ELEMENT*nElem_Recv[size]; ii++)
    connRecv[ii] = 0;

  unsigned short *haloRecv = new unsigned short[nElem_Recv[size]];
  for (int ii = 0; ii < nElem_Recv[size]; ii++)
    haloRecv[ii] = false;

#ifdef HAVE_MPI
  /*--- We need double the number of messages to send both the conn.
   and the flags for the halo cells. ---*/

  send_req = new SU2_MPI::Request[2*nSends];
  recv_req = new SU2_MPI::Request[2*nRecvs];

  /*--- Launch the non-blocking recv's for the connectivity. ---*/

  unsigned long iMessage = 0;
  for (int ii=0; ii<size; ii++) {
    if ((ii != rank) && (nElem_Recv[ii+1] > nElem_Recv[ii])) {
      int ll     = NODES_PER_ELEMENT*nElem_Recv[ii];
      int kk     = nElem_Recv[ii+1] - nElem_Recv[ii];
      int count  = NODES_PER_ELEMENT*kk;
      int source = ii;
      int tag    = ii + 1;
      SU2_MPI::Irecv(&(connRecv[ll]), count, MPI_UNSIGNED_LONG, source, tag,
                     MPI_COMM_WORLD, &(recv_req[iMessage]));
      iMessage++;
    }
  }

  /*--- Launch the non-blocking sends of the connectivity. ---*/

  iMessage = 0;
  for (int ii=0; ii<size; ii++) {
    if ((ii != rank) && (nElem_Send[ii+1] > nElem_Send[ii])) {
      int ll = NODES_PER_ELEMENT*nElem_Send[ii];
      int kk = nElem_Send[ii+1] - nElem_Send[ii];
      int count  = NODES_PER_ELEMENT*kk;
      int dest = ii;
      int tag    = rank + 1;
      SU2_MPI::Isend(&(connSend[ll]), count, MPI_UNSIGNED_LONG, dest, tag,
                     MPI_COMM_WORLD, &(send_req[iMessage]));
      iMessage++;
    }
  }

  /*--- Repeat the process to communicate the halo flags. ---*/

  iMessage = 0;
  for (int ii=0; ii<size; ii++) {
    if ((ii != rank) && (nElem_Recv[ii+1] > nElem_Recv[ii])) {
      int ll     = nElem_Recv[ii];
      int kk     = nElem_Recv[ii+1] - nElem_Recv[ii];
      int count  = kk;
      int source = ii;
      int tag    = ii + 1;
      SU2_MPI::Irecv(&(haloRecv[ll]), count, MPI_UNSIGNED_SHORT, source, tag,
                     MPI_COMM_WORLD, &(recv_req[iMessage+nRecvs]));
      iMessage++;
    }
  }

  /*--- Launch the non-blocking sends of the halo flags. ---*/

  iMessage = 0;
  for (int ii=0; ii<size; ii++) {
    if ((ii != rank) && (nElem_Send[ii+1] > nElem_Send[ii])) {
      int ll = nElem_Send[ii];
      int kk = nElem_Send[ii+1] - nElem_Send[ii];
      int count  = kk;
      int dest   = ii;
      int tag    = rank + 1;
      SU2_MPI::Isend(&(haloSend[ll]), count, MPI_UNSIGNED_SHORT, dest, tag,
                     MPI_COMM_WORLD, &(send_req[iMessage+nSends]));
      iMessage++;
    }
  }
#endif

  /*--- Copy my own rank's data into the recv buffer directly. ---*/

  int mm = NODES_PER_ELEMENT*nElem_Recv[rank];
  int ll = NODES_PER_ELEMENT*nElem_Send[rank];
  int kk = NODES_PER_ELEMENT*nElem_Send[rank+1];

  for (int nn=ll; nn<kk; nn++, mm++) connRecv[mm] = connSend[nn];

  mm = nElem_Recv[rank];
  ll = nElem_Send[rank];
  kk = nElem_Send[rank+1];

  for (int nn=ll; nn<kk; nn++, mm++) haloRecv[mm] = haloSend[nn];

  /*--- Wait for the non-blocking sends and recvs to complete. ---*/

#ifdef HAVE_MPI
  int number = 2*nSends;
  for (int ii = 0; ii < number; ii++)
    SU2_MPI::Waitany(number, send_req, &ind, &status);

  number = 2*nRecvs;
  for (int ii = 0; ii < number; ii++)
    SU2_MPI::Waitany(number, recv_req, &ind, &status);

  delete [] send_req;
  delete [] recv_req;
#endif

  /*--- Store the connectivity for this rank in the proper data
   structure before post-processing below. Note that we add 1 here
   to the connectivity for vizualization packages. First, allocate
   appropriate amount of memory for this section. ---*/

  if (nElem_Recv[size] > 0) Conn_Elem = new int[NODES_PER_ELEMENT*nElem_Recv[size]];
  int count = 0; nElem_Total = 0;
  for (int ii = 0; ii < nElem_Recv[size]; ii++) {
    if (!haloRecv[ii]) {
      nElem_Total++;
      for (int jj = 0; jj < NODES_PER_ELEMENT; jj++) {
        Conn_Elem[count] = (int)connRecv[ii*NODES_PER_ELEMENT+jj] + 1;
        count++;
      }
    }
  }

  /*--- Store the particular global element count in the class data,
   and set the class data pointer to the connectivity array. ---*/

  switch (Elem_Type) {
    case TRIANGLE:
      nParallel_Tria = nElem_Total;
      if (nParallel_Tria > 0) Conn_Tria_Par = Conn_Elem;
      break;
    case QUADRILATERAL:
      nParallel_Quad = nElem_Total;
      if (nParallel_Quad > 0) Conn_Quad_Par = Conn_Elem;
      break;
    case TETRAHEDRON:
      nParallel_Tetr = nElem_Total;
      if (nParallel_Tetr > 0) Conn_Tetr_Par = Conn_Elem;
      break;
    case HEXAHEDRON:
      nParallel_Hexa = nElem_Total;
      if (nParallel_Hexa > 0) Conn_Hexa_Par = Conn_Elem;
      break;
    case PRISM:
      nParallel_Pris = nElem_Total;
      if (nParallel_Pris > 0) Conn_Pris_Par = Conn_Elem;
      break;
    case PYRAMID:
      nParallel_Pyra = nElem_Total;
      if (nParallel_Pyra > 0) Conn_Pyra_Par = Conn_Elem;
      break;
    default:
      SU2_MPI::Error("Unrecognized element type", CURRENT_FUNCTION);
      break;
  }

  /*--- Free temporary memory from communications ---*/

  delete [] connSend;
  delete [] connRecv;
  delete [] haloSend;
  delete [] haloRecv;
  delete [] Local_Halo;
  delete [] nElem_Recv;
  delete [] nElem_Send;
  delete [] nElem_Flag;
  delete [] Buffer_Recv_nAddedPeriodic;
  delete [] Buffer_Send_AddedPeriodic;
  delete [] Buffer_Recv_AddedPeriodic;
  delete [] npoint_procs;
  delete [] starting_node;
  delete [] ending_node;
  delete [] nPoint_Linear;

}

void COutputLegacy::SortSurfaceConnectivity(CConfig *config, CGeometry *geometry, unsigned short Elem_Type) {

  unsigned long iProcessor;
  unsigned short NODES_PER_ELEMENT;
  unsigned long iPoint, jPoint, kPoint, nLocalPoint, nTotalPoint;
  unsigned long nElem_Total = 0, Global_Index;

  unsigned long iVertex, iMarker;
  int SendRecv, RecvFrom;

  bool notPeriodic, notHalo, addedPeriodic, isPeriodic;

  int *Local_Halo = NULL;
  int *Conn_Elem  = NULL;

#ifdef HAVE_MPI
  SU2_MPI::Request *send_req, *recv_req;
  SU2_MPI::Status status;
  int ind;
#endif

  /*--- Store the local number of this element type and the number of nodes
   per this element type. In serial, this will be the total number of this
   element type in the entire mesh. In parallel, it is the number on only
   the current partition. ---*/

  switch (Elem_Type) {
    case LINE:
      NODES_PER_ELEMENT = N_POINTS_LINE;
      break;
    case TRIANGLE:
      NODES_PER_ELEMENT = N_POINTS_TRIANGLE;
      break;
    case QUADRILATERAL:
      NODES_PER_ELEMENT = N_POINTS_QUADRILATERAL;
      break;
    default:
      SU2_MPI::Error("Unrecognized element type", CURRENT_FUNCTION);
      NODES_PER_ELEMENT = 0;
      break;
  }

  /*--- Force the removal of all added periodic elements (use global index).
   First, we isolate and create a list of all added periodic points, excluding
   those that were part of the original domain (we want these to be in the
   output files). ---*/

  vector<unsigned long> Added_Periodic;
  Added_Periodic.clear();

  if (config->GetKind_SU2() != SU2_DEF) {
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) {
        SendRecv = config->GetMarker_All_SendRecv(iMarker);
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

          if ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() > 0) &&
              (geometry->vertex[iMarker][iVertex]->GetRotation_Type() % 2 == 0) &&
              (SendRecv < 0)) {
            Added_Periodic.push_back(geometry->node[iPoint]->GetGlobalIndex());
          }
        }
      }
    }
  }

  /*--- Now we communicate this information to all processors, so that they
   can force the removal of these particular nodes by flagging them as halo
   points. In general, this should be a small percentage of the total mesh,
   so the communication/storage costs here shouldn't be prohibitive. ---*/

  /*--- First communicate the number of points that each rank has found. ---*/

  unsigned long nAddedPeriodic = 0, maxAddedPeriodic = 0;
  unsigned long Buffer_Send_nAddedPeriodic[1], *Buffer_Recv_nAddedPeriodic = NULL;
  Buffer_Recv_nAddedPeriodic = new unsigned long[size];

  nAddedPeriodic = Added_Periodic.size();
  Buffer_Send_nAddedPeriodic[0] = nAddedPeriodic;

#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&nAddedPeriodic, &maxAddedPeriodic, 1, MPI_UNSIGNED_LONG,
                     MPI_MAX, MPI_COMM_WORLD);
  SU2_MPI::Allgather(&Buffer_Send_nAddedPeriodic, 1, MPI_UNSIGNED_LONG,
                     Buffer_Recv_nAddedPeriodic,  1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
#else
  maxAddedPeriodic = nAddedPeriodic;
  Buffer_Recv_nAddedPeriodic[0] = Buffer_Send_nAddedPeriodic[0];
#endif

  /*--- Communicate the global index values of all added periodic nodes. ---*/
  unsigned long *Buffer_Send_AddedPeriodic = new unsigned long[maxAddedPeriodic];
  unsigned long *Buffer_Recv_AddedPeriodic = new unsigned long[size*maxAddedPeriodic];

  for (iPoint = 0; iPoint < Added_Periodic.size(); iPoint++) {
    Buffer_Send_AddedPeriodic[iPoint] = Added_Periodic[iPoint];
  }

  /*--- Gather the element connectivity information. All processors will now
   have a copy of the global index values for all added periodic points. ---*/

#ifdef HAVE_MPI
  SU2_MPI::Allgather(Buffer_Send_AddedPeriodic, maxAddedPeriodic, MPI_UNSIGNED_LONG,
                     Buffer_Recv_AddedPeriodic, maxAddedPeriodic, MPI_UNSIGNED_LONG,
                     MPI_COMM_WORLD);
#else
  for (iPoint = 0; iPoint < maxAddedPeriodic; iPoint++)
    Buffer_Recv_AddedPeriodic[iPoint] = Buffer_Send_AddedPeriodic[iPoint];
#endif

  /*--- Search all send/recv boundaries on this partition for halo cells. In
   particular, consider only the recv conditions (these are the true halo
   nodes). Check the ranks of the processors that are communicating and
   choose to keep only the halo cells from the higher rank processor. Here,
   we are also choosing to keep periodic nodes that were part of the original
   domain. We will check the communicated list of added periodic points. ---*/

  Local_Halo = new int[geometry->GetnPoint()];
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
    Local_Halo[iPoint] = !geometry->node[iPoint]->GetDomain();

  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) {
      SendRecv = config->GetMarker_All_SendRecv(iMarker);
      RecvFrom = abs(SendRecv)-1;

      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        Global_Index = geometry->node[iPoint]->GetGlobalIndex();

        /*--- We need to keep one copy of overlapping halo cells. ---*/

        notHalo = ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() == 0) &&
                   (SendRecv < 0) && (rank > RecvFrom));

        /*--- We want to keep the periodic nodes that were part of the original domain.
         For SU2_DEF we want to keep all periodic nodes. ---*/

        if (config->GetKind_SU2() == SU2_DEF) {
          isPeriodic = ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() > 0));
        }else {
          isPeriodic = ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() > 0) &&
                        (geometry->vertex[iMarker][iVertex]->GetRotation_Type() % 2 == 1));
        }

        notPeriodic = (isPeriodic && (SendRecv < 0));

        /*--- Lastly, check that this isn't an added periodic point that
         we will forcibly remove. Use the communicated list of these points. ---*/

        addedPeriodic = false; kPoint = 0;
        for (iProcessor = 0; iProcessor < (unsigned long)size; iProcessor++) {
          for (jPoint = 0; jPoint < Buffer_Recv_nAddedPeriodic[iProcessor]; jPoint++) {
            if (Global_Index == Buffer_Recv_AddedPeriodic[kPoint+jPoint])
              addedPeriodic = true;
          }

          /*--- Adjust jNode to index of next proc's data in the buffers. ---*/

          kPoint = (iProcessor+1)*maxAddedPeriodic;

        }

        /*--- If we found either of these types of nodes, flag them to be kept. ---*/

        if ((notHalo || notPeriodic) && !addedPeriodic) {
          Local_Halo[iPoint] = false;
        }

      }
    }
  }

  /*--- Now that we've done the gymnastics to find any periodic points,
   compute the total number of local and global points for the output. ---*/

  nLocalPoint = 0;
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
    if (Local_Halo[iPoint] == false)
      nLocalPoint++;

#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&nLocalPoint, &nTotalPoint, 1,
                     MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#else
  nTotalPoint = nLocalPoint;
#endif

  /*--- Compute the number of points that will be on each processor.
   This is a linear partitioning with the addition of a simple load
   balancing for any remainder points. ---*/

  unsigned long *npoint_procs  = new unsigned long[size];
  unsigned long *starting_node = new unsigned long[size];
  unsigned long *ending_node   = new unsigned long[size];
  unsigned long *nPoint_Linear = new unsigned long[size+1];

  unsigned long total_pt_accounted = 0;
  for (int ii = 0; ii < size; ii++) {
    npoint_procs[ii] = nTotalPoint/size;
    total_pt_accounted = total_pt_accounted + npoint_procs[ii];
  }

  /*--- Get the number of remainder points after the even division. ---*/

  unsigned long rem_points = nTotalPoint-total_pt_accounted;
  for (unsigned long ii = 0; ii < rem_points; ii++) {
    npoint_procs[ii]++;
  }

  /*--- Store the local number of nodes and the beginning/end index ---*/

  starting_node[0] = 0;
  ending_node[0]   = starting_node[0] + npoint_procs[0];
  nPoint_Linear[0] = 0;
  for (int ii = 1; ii < size; ii++) {
    starting_node[ii] = ending_node[ii-1];
    ending_node[ii]   = starting_node[ii] + npoint_procs[ii];
    nPoint_Linear[ii] = nPoint_Linear[ii-1] + npoint_procs[ii-1];
  }
  nPoint_Linear[size] = nTotalPoint;

  /*--- We start with the connectivity distributed across all procs with
   no particular ordering assumed. We need to loop through our local partition
   and decide how many elements we must send to each other rank in order to
   have all elements sorted according to a linear partitioning of the grid
   nodes, i.e., rank 0 holds the first nPoint()/nProcessors nodes.
   First, initialize a counter and flag. ---*/

  int *nElem_Send = new int[size+1]; nElem_Send[0] = 0;
  int *nElem_Recv = new int[size+1]; nElem_Recv[0] = 0;
  int *nElem_Flag = new int[size];

  for (int ii=0; ii < size; ii++) {
    nElem_Send[ii] = 0;
    nElem_Recv[ii] = 0;
    nElem_Flag[ii]= -1;
  }
  nElem_Send[size] = 0; nElem_Recv[size] = 0;

  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_Plotting(iMarker) == YES) {

      for (int ii = 0; ii < (int)geometry->GetnElem_Bound(iMarker); ii++) {

        if (geometry->bound[iMarker][ii]->GetVTK_Type() == Elem_Type) {
          for ( int jj = 0; jj < NODES_PER_ELEMENT; jj++ ) {

            /*--- Get the index of the current point. ---*/

            iPoint = geometry->bound[iMarker][ii]->GetNode(jj);
            Global_Index = geometry->node[iPoint]->GetGlobalIndex();

            /*--- Search for the lowest global index in this element. We
             send the element to the processor owning the range that includes
             the lowest global index value. ---*/

            for (int kk = 0; kk < NODES_PER_ELEMENT; kk++) {
              jPoint = geometry->bound[iMarker][ii]->GetNode(kk);
              unsigned long newID = geometry->node[jPoint]->GetGlobalIndex();
              if (newID < Global_Index) Global_Index = newID;
            }

            /*--- Search for the processor that owns this point ---*/

            iProcessor = Global_Index/npoint_procs[0];
            if (iProcessor >= (unsigned long)size)
              iProcessor = (unsigned long)size-1;
            if (Global_Index >= nPoint_Linear[iProcessor])
              while(Global_Index >= nPoint_Linear[iProcessor+1]) iProcessor++;
            else
              while(Global_Index <  nPoint_Linear[iProcessor])   iProcessor--;

            /*--- If we have not visited this element yet, increment our
             number of elements that must be sent to a particular proc. ---*/

            if ((nElem_Flag[iProcessor] != ii)) {
              nElem_Flag[iProcessor] = ii;
              nElem_Send[iProcessor+1]++;
            }

          }
        }
      }
    }
  }

  /*--- Communicate the number of cells to be sent/recv'd amongst
   all processors. After this communication, each proc knows how
   many cells it will receive from each other processor. ---*/

#ifdef HAVE_MPI
  SU2_MPI::Alltoall(&(nElem_Send[1]), 1, MPI_INT,
                    &(nElem_Recv[1]), 1, MPI_INT, MPI_COMM_WORLD);
#else
  nElem_Recv[1] = nElem_Send[1];
#endif

  /*--- Prepare to send connectivities. First check how many
   messages we will be sending and receiving. Here we also put
   the counters into cumulative storage format to make the
   communications simpler. ---*/

  int nSends = 0, nRecvs = 0;
  for (int ii=0; ii < size; ii++) nElem_Flag[ii] = -1;

  for (int ii = 0; ii < size; ii++) {
    if ((ii != rank) && (nElem_Send[ii+1] > 0)) nSends++;
    if ((ii != rank) && (nElem_Recv[ii+1] > 0)) nRecvs++;

    nElem_Send[ii+1] += nElem_Send[ii];
    nElem_Recv[ii+1] += nElem_Recv[ii];
  }

  /*--- Allocate memory to hold the connectivity that we are
   sending. ---*/

  unsigned long *connSend = NULL;
  connSend = new unsigned long[NODES_PER_ELEMENT*nElem_Send[size]];
  for (int ii = 0; ii < NODES_PER_ELEMENT*nElem_Send[size]; ii++)
    connSend[ii] = 0;

  /*--- Allocate arrays for storing halo flags. ---*/

  unsigned short *haloSend = new unsigned short[nElem_Send[size]];
  for (int ii = 0; ii < nElem_Send[size]; ii++)
    haloSend[ii] = false;

  /*--- Create an index variable to keep track of our index
   position as we load up the send buffer. ---*/

  unsigned long *index = new unsigned long[size];
  for (int ii=0; ii < size; ii++) index[ii] = NODES_PER_ELEMENT*nElem_Send[ii];

  unsigned long *haloIndex = new unsigned long[size];
  for (int ii=0; ii < size; ii++) haloIndex[ii] = nElem_Send[ii];

  /*--- Loop through our elements and load the elems and their
   additional data that we will send to the other procs. ---*/

  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_Plotting(iMarker) == YES) {

      for (int ii = 0; ii < (int)geometry->GetnElem_Bound(iMarker); ii++) {

        if (geometry->bound[iMarker][ii]->GetVTK_Type() == Elem_Type) {
          for ( int jj = 0; jj < NODES_PER_ELEMENT; jj++ ) {

            /*--- Get the index of the current point. ---*/

            iPoint = geometry->bound[iMarker][ii]->GetNode(jj);
            Global_Index = geometry->node[iPoint]->GetGlobalIndex();

            /*--- Search for the lowest global index in this element. We
             send the element to the processor owning the range that includes
             the lowest global index value. ---*/

            for (int kk = 0; kk < NODES_PER_ELEMENT; kk++) {
              jPoint = geometry->bound[iMarker][ii]->GetNode(kk);
              unsigned long newID = geometry->node[jPoint]->GetGlobalIndex();
              if (newID < Global_Index) Global_Index = newID;
            }

            /*--- Search for the processor that owns this point ---*/

            iProcessor = Global_Index/npoint_procs[0];
            if (iProcessor >= (unsigned long)size)
              iProcessor = (unsigned long)size-1;
            if (Global_Index >= nPoint_Linear[iProcessor])
              while(Global_Index >= nPoint_Linear[iProcessor+1]) iProcessor++;
            else
              while(Global_Index <  nPoint_Linear[iProcessor])   iProcessor--;

            /*--- Load connectivity into the buffer for sending ---*/

            if (nElem_Flag[iProcessor] != ii) {

              nElem_Flag[iProcessor] = ii;
              unsigned long nn = index[iProcessor];
              unsigned long mm = haloIndex[iProcessor];

              /*--- Load the connectivity values. ---*/

              for (int kk = 0; kk < NODES_PER_ELEMENT; kk++) {
                iPoint = geometry->bound[iMarker][ii]->GetNode(kk);
                connSend[nn] = geometry->node[iPoint]->GetGlobalIndex(); nn++;

                /*--- Check if this is a halo node. If so, flag this element
                 as a halo cell. We will use this later to sort and remove
                 any duplicates from the connectivity list. ---*/

                if (Local_Halo[iPoint]) haloSend[mm] = true;

              }

              /*--- Increment the index by the message length ---*/

              index[iProcessor]    += NODES_PER_ELEMENT;
              haloIndex[iProcessor]++;

            }
          }
        }
      }
    }
  }

  /*--- Free memory after loading up the send buffer. ---*/

  delete [] index;
  delete [] haloIndex;

  /*--- Allocate the memory that we need for receiving the conn
   values and then cue up the non-blocking receives. Note that
   we do not include our own rank in the communications. We will
   directly copy our own data later. ---*/

  unsigned long *connRecv = NULL;
  connRecv = new unsigned long[NODES_PER_ELEMENT*nElem_Recv[size]];
  for (int ii = 0; ii < NODES_PER_ELEMENT*nElem_Recv[size]; ii++)
    connRecv[ii] = 0;

  unsigned short *haloRecv = new unsigned short[nElem_Recv[size]];
  for (int ii = 0; ii < nElem_Recv[size]; ii++)
    haloRecv[ii] = false;

#ifdef HAVE_MPI
  /*--- We need double the number of messages to send both the conn.
   and the flags for the halo cells. ---*/

  send_req = new SU2_MPI::Request[2*nSends];
  recv_req = new SU2_MPI::Request[2*nRecvs];

  /*--- Launch the non-blocking recv's for the connectivity. ---*/

  unsigned long iMessage = 0;
  for (int ii=0; ii<size; ii++) {
    if ((ii != rank) && (nElem_Recv[ii+1] > nElem_Recv[ii])) {
      int ll     = NODES_PER_ELEMENT*nElem_Recv[ii];
      int kk     = nElem_Recv[ii+1] - nElem_Recv[ii];
      int count  = NODES_PER_ELEMENT*kk;
      int source = ii;
      int tag    = ii + 1;
      SU2_MPI::Irecv(&(connRecv[ll]), count, MPI_UNSIGNED_LONG, source, tag,
                     MPI_COMM_WORLD, &(recv_req[iMessage]));
      iMessage++;
    }
  }

  /*--- Launch the non-blocking sends of the connectivity. ---*/

  iMessage = 0;
  for (int ii=0; ii<size; ii++) {
    if ((ii != rank) && (nElem_Send[ii+1] > nElem_Send[ii])) {
      int ll = NODES_PER_ELEMENT*nElem_Send[ii];
      int kk = nElem_Send[ii+1] - nElem_Send[ii];
      int count  = NODES_PER_ELEMENT*kk;
      int dest = ii;
      int tag    = rank + 1;
      SU2_MPI::Isend(&(connSend[ll]), count, MPI_UNSIGNED_LONG, dest, tag,
                     MPI_COMM_WORLD, &(send_req[iMessage]));
      iMessage++;
    }
  }

  /*--- Repeat the process to communicate the halo flags. ---*/

  iMessage = 0;
  for (int ii=0; ii<size; ii++) {
    if ((ii != rank) && (nElem_Recv[ii+1] > nElem_Recv[ii])) {
      int ll     = nElem_Recv[ii];
      int kk     = nElem_Recv[ii+1] - nElem_Recv[ii];
      int count  = kk;
      int source = ii;
      int tag    = ii + 1;
      SU2_MPI::Irecv(&(haloRecv[ll]), count, MPI_UNSIGNED_SHORT, source, tag,
                     MPI_COMM_WORLD, &(recv_req[iMessage+nRecvs]));
      iMessage++;
    }
  }

  /*--- Launch the non-blocking sends of the halo flags. ---*/

  iMessage = 0;
  for (int ii=0; ii<size; ii++) {
    if ((ii != rank) && (nElem_Send[ii+1] > nElem_Send[ii])) {
      int ll = nElem_Send[ii];
      int kk = nElem_Send[ii+1] - nElem_Send[ii];
      int count  = kk;
      int dest   = ii;
      int tag    = rank + 1;
      SU2_MPI::Isend(&(haloSend[ll]), count, MPI_UNSIGNED_SHORT, dest, tag,
                     MPI_COMM_WORLD, &(send_req[iMessage+nSends]));
      iMessage++;
    }
  }
#endif

  /*--- Copy my own rank's data into the recv buffer directly. ---*/

  int mm = NODES_PER_ELEMENT*nElem_Recv[rank];
  int ll = NODES_PER_ELEMENT*nElem_Send[rank];
  int kk = NODES_PER_ELEMENT*nElem_Send[rank+1];

  for (int nn=ll; nn<kk; nn++, mm++) connRecv[mm] = connSend[nn];

  mm = nElem_Recv[rank];
  ll = nElem_Send[rank];
  kk = nElem_Send[rank+1];

  for (int nn=ll; nn<kk; nn++, mm++) haloRecv[mm] = haloSend[nn];

  /*--- Wait for the non-blocking sends and recvs to complete. ---*/

#ifdef HAVE_MPI
  int number = 2*nSends;
  for (int ii = 0; ii < number; ii++)
    SU2_MPI::Waitany(number, send_req, &ind, &status);

  number = 2*nRecvs;
  for (int ii = 0; ii < number; ii++)
    SU2_MPI::Waitany(number, recv_req, &ind, &status);

  delete [] send_req;
  delete [] recv_req;
#endif

  /*--- Store the connectivity for this rank in the proper data
   structure before post-processing below. Note that we add 1 here
   to the connectivity for vizualization packages. First, allocate
   appropriate amount of memory for this section. ---*/

  if (nElem_Recv[size] > 0) Conn_Elem = new int[NODES_PER_ELEMENT*nElem_Recv[size]];
  int count = 0; nElem_Total = 0;
  for (int ii = 0; ii < nElem_Recv[size]; ii++) {
    if (!haloRecv[ii]) {
      nElem_Total++;
      for (int jj = 0; jj < NODES_PER_ELEMENT; jj++) {
        Conn_Elem[count] = (int)connRecv[ii*NODES_PER_ELEMENT+jj] + 1;
        count++;
      }
    }
  }

  /*--- Store the particular global element count in the class data,
   and set the class data pointer to the connectivity array. ---*/

  switch (Elem_Type) {
    case LINE:
      nParallel_Line = nElem_Total;
      Conn_BoundLine_Par = Conn_Elem;
      break;
    case TRIANGLE:
      nParallel_BoundTria = nElem_Total;
      Conn_BoundTria_Par = Conn_Elem;
      break;
    case QUADRILATERAL:
      nParallel_BoundQuad = nElem_Total;
      Conn_BoundQuad_Par = Conn_Elem;
      break;
    default:
      SU2_MPI::Error("Unrecognized element type", CURRENT_FUNCTION);
      break;
  }

  /*--- Free temporary memory from communications ---*/

  delete [] connSend;
  delete [] connRecv;
  delete [] haloSend;
  delete [] haloRecv;
  delete [] Local_Halo;
  delete [] nElem_Recv;
  delete [] nElem_Send;
  delete [] nElem_Flag;
  delete [] Buffer_Recv_nAddedPeriodic;
  delete [] Buffer_Send_AddedPeriodic;
  delete [] Buffer_Recv_AddedPeriodic;
  delete [] npoint_procs;
  delete [] starting_node;
  delete [] ending_node;
  delete [] nPoint_Linear;

}

void COutputLegacy::SortOutputData(CConfig *config, CGeometry *geometry) {

  unsigned short iMarker;
  unsigned long iProcessor;
  unsigned long iPoint, Global_Index, nLocalPoint, nTotalPoint, iVertex;

  int VARS_PER_POINT = nVar_Par;
  int *Local_Halo = NULL;

  bool isPeriodic;

#ifdef HAVE_MPI
  SU2_MPI::Request *send_req, *recv_req;
  SU2_MPI::Status status;
  int ind;
#endif

  /*--- Search all send/recv boundaries on this partition for any periodic
   nodes that were part of the original domain. We want to recover these
   for visualization purposes. ---*/

  Local_Halo = new int[geometry->GetnPoint()];
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
    Local_Halo[iPoint] = !geometry->node[iPoint]->GetDomain();

  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) {

      /*--- Checking for less than or equal to the rank, because there may
       be some periodic halo nodes that send info to the same rank. ---*/

      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        isPeriodic = ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() > 0) &&
                      (geometry->vertex[iMarker][iVertex]->GetRotation_Type() % 2 == 1));
        if (isPeriodic) Local_Halo[iPoint] = false;
      }
    }
  }

  /*--- Sum total number of nodes that belong to the domain ---*/

  nLocalPoint = 0;
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
    if (Local_Halo[iPoint] == false)
      nLocalPoint++;

#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&nLocalPoint, &nTotalPoint, 1,
                     MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#else
  nTotalPoint = nLocalPoint;
#endif

  /*--- Now that we know the actual number of points we need to output,
   compute the number of points that will be on each processor.
   This is a linear partitioning with the addition of a simple load
   balancing for any remainder points. ---*/

  unsigned long *npoint_procs  = new unsigned long[size];
  unsigned long *starting_node = new unsigned long[size];
  unsigned long *ending_node   = new unsigned long[size];
  unsigned long *nPoint_Linear = new unsigned long[size+1];

  unsigned long total_pt_accounted = 0;
  for (int ii = 0; ii < size; ii++) {
    npoint_procs[ii] = nTotalPoint/size;
    total_pt_accounted = total_pt_accounted + npoint_procs[ii];
  }

  /*--- Get the number of remainder points after the even division. ---*/

  unsigned long rem_points = nTotalPoint-total_pt_accounted;
  for (unsigned long ii = 0; ii < rem_points; ii++) {
    npoint_procs[ii]++;
  }

  /*--- Store the local number of nodes and the beginning/end index ---*/

  starting_node[0] = 0;
  ending_node[0]   = starting_node[0] + npoint_procs[0];
  nPoint_Linear[0] = 0;
  for (int ii = 1; ii < size; ii++) {
    starting_node[ii] = ending_node[ii-1];
    ending_node[ii]   = starting_node[ii] + npoint_procs[ii];
    nPoint_Linear[ii] = nPoint_Linear[ii-1] + npoint_procs[ii-1];
  }
  nPoint_Linear[size] = nTotalPoint;

  /*--- We start with the grid nodes distributed across all procs with
   no particular ordering assumed. We need to loop through our local partition
   and decide how many nodes we must send to each other rank in order to
   have all nodes sorted according to a linear partitioning of the grid
   nodes, i.e., rank 0 holds the first ~ nGlobalPoint()/nProcessors nodes.
   First, initialize a counter and flag. ---*/

  int *nPoint_Send = new int[size+1]; nPoint_Send[0] = 0;
  int *nPoint_Recv = new int[size+1]; nPoint_Recv[0] = 0;
  int *nPoint_Flag = new int[size];

  for (int ii=0; ii < size; ii++) {
    nPoint_Send[ii] = 0;
    nPoint_Recv[ii] = 0;
    nPoint_Flag[ii]= -1;
  }
  nPoint_Send[size] = 0; nPoint_Recv[size] = 0;

  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++ ) {

    /*--- We only write interior points and recovered periodic points. ---*/

    if (!Local_Halo[iPoint]) {

      /*--- Get the global index of the current point. ---*/

      Global_Index = geometry->node[iPoint]->GetGlobalIndex();

      /*--- Search for the processor that owns this point ---*/

      iProcessor = Global_Index/npoint_procs[0];
      if (iProcessor >= (unsigned long)size)
        iProcessor = (unsigned long)size-1;
      if (Global_Index >= nPoint_Linear[iProcessor])
        while(Global_Index >= nPoint_Linear[iProcessor+1]) iProcessor++;
      else
        while(Global_Index <  nPoint_Linear[iProcessor])   iProcessor--;

      /*--- If we have not visited this node yet, increment our
       number of elements that must be sent to a particular proc. ---*/

      if (nPoint_Flag[iProcessor] != (int)iPoint) {
        nPoint_Flag[iProcessor] = (int)iPoint;
        nPoint_Send[iProcessor+1]++;
      }

    }
  }

  /*--- Communicate the number of nodes to be sent/recv'd amongst
   all processors. After this communication, each proc knows how
   many cells it will receive from each other processor. ---*/

#ifdef HAVE_MPI
  SU2_MPI::Alltoall(&(nPoint_Send[1]), 1, MPI_INT,
                    &(nPoint_Recv[1]), 1, MPI_INT, MPI_COMM_WORLD);
#else
  nPoint_Recv[1] = nPoint_Send[1];
#endif

  /*--- Prepare to send coordinates. First check how many
   messages we will be sending and receiving. Here we also put
   the counters into cumulative storage format to make the
   communications simpler. ---*/

  int nSends = 0, nRecvs = 0;
  for (int ii=0; ii < size; ii++) nPoint_Flag[ii] = -1;

  for (int ii = 0; ii < size; ii++) {
    if ((ii != rank) && (nPoint_Send[ii+1] > 0)) nSends++;
    if ((ii != rank) && (nPoint_Recv[ii+1] > 0)) nRecvs++;

    nPoint_Send[ii+1] += nPoint_Send[ii];
    nPoint_Recv[ii+1] += nPoint_Recv[ii];
  }

  /*--- Allocate memory to hold the connectivity that we are
   sending. ---*/

  su2double *connSend = NULL;
  connSend = new su2double[VARS_PER_POINT*nPoint_Send[size]];
  for (int ii = 0; ii < VARS_PER_POINT*nPoint_Send[size]; ii++)
    connSend[ii] = 0;

  /*--- Allocate arrays for sending the global ID. ---*/

  unsigned long *idSend = new unsigned long[nPoint_Send[size]];
  for (int ii = 0; ii < nPoint_Send[size]; ii++)
    idSend[ii] = 0;

  /*--- Create an index variable to keep track of our index
   positions as we load up the send buffer. ---*/

  unsigned long *index = new unsigned long[size];
  for (int ii=0; ii < size; ii++) index[ii] = VARS_PER_POINT*nPoint_Send[ii];

  unsigned long *idIndex = new unsigned long[size];
  for (int ii=0; ii < size; ii++) idIndex[ii] = nPoint_Send[ii];

  /*--- Loop through our elements and load the elems and their
   additional data that we will send to the other procs. ---*/

  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {

    /*--- We only write interior points and recovered periodic points. ---*/

    if (!Local_Halo[iPoint]) {

      /*--- Get the index of the current point. ---*/

      Global_Index = geometry->node[iPoint]->GetGlobalIndex();

      /*--- Search for the processor that owns this point. ---*/

      iProcessor = Global_Index/npoint_procs[0];
      if (iProcessor >= (unsigned long)size)
        iProcessor = (unsigned long)size-1;
      if (Global_Index >= nPoint_Linear[iProcessor])
        while(Global_Index >= nPoint_Linear[iProcessor+1]) iProcessor++;
      else
        while(Global_Index <  nPoint_Linear[iProcessor])   iProcessor--;

      /*--- Load node coordinates into the buffer for sending. ---*/

      if (nPoint_Flag[iProcessor] != (int)iPoint) {

        nPoint_Flag[iProcessor] = (int)iPoint;
        unsigned long nn = index[iProcessor];

        /*--- Load the data values. ---*/

        for (unsigned short kk = 0; kk < VARS_PER_POINT; kk++) {
          connSend[nn] = Local_Data[iPoint][kk]; nn++;
        }

        /*--- Load the global ID (minus offset) for sorting the
         points once they all reach the correct processor. ---*/

        nn = idIndex[iProcessor];
        idSend[nn] = Global_Index - starting_node[iProcessor];

        /*--- Increment the index by the message length ---*/

        index[iProcessor]  += VARS_PER_POINT;
        idIndex[iProcessor]++;

      }
    }
  }

  /*--- Free memory after loading up the send buffer. ---*/

  delete [] index;
  delete [] idIndex;

  /*--- Allocate the memory that we need for receiving the conn
   values and then cue up the non-blocking receives. Note that
   we do not include our own rank in the communications. We will
   directly copy our own data later. ---*/

  su2double *connRecv = NULL;
  connRecv = new su2double[VARS_PER_POINT*nPoint_Recv[size]];
  for (int ii = 0; ii < VARS_PER_POINT*nPoint_Recv[size]; ii++)
    connRecv[ii] = 0;

  unsigned long *idRecv = new unsigned long[nPoint_Recv[size]];
  for (int ii = 0; ii < nPoint_Recv[size]; ii++)
    idRecv[ii] = 0;

#ifdef HAVE_MPI
  /*--- We need double the number of messages to send both the conn.
   and the global IDs. ---*/

  send_req = new SU2_MPI::Request[2*nSends];
  recv_req = new SU2_MPI::Request[2*nRecvs];

  unsigned long iMessage = 0;
  for (int ii=0; ii<size; ii++) {
    if ((ii != rank) && (nPoint_Recv[ii+1] > nPoint_Recv[ii])) {
      int ll     = VARS_PER_POINT*nPoint_Recv[ii];
      int kk     = nPoint_Recv[ii+1] - nPoint_Recv[ii];
      int count  = VARS_PER_POINT*kk;
      int source = ii;
      int tag    = ii + 1;
      SU2_MPI::Irecv(&(connRecv[ll]), count, MPI_DOUBLE, source, tag,
                     MPI_COMM_WORLD, &(recv_req[iMessage]));
      iMessage++;
    }
  }

  /*--- Launch the non-blocking sends of the connectivity. ---*/

  iMessage = 0;
  for (int ii=0; ii<size; ii++) {
    if ((ii != rank) && (nPoint_Send[ii+1] > nPoint_Send[ii])) {
      int ll = VARS_PER_POINT*nPoint_Send[ii];
      int kk = nPoint_Send[ii+1] - nPoint_Send[ii];
      int count  = VARS_PER_POINT*kk;
      int dest = ii;
      int tag    = rank + 1;
      SU2_MPI::Isend(&(connSend[ll]), count, MPI_DOUBLE, dest, tag,
                     MPI_COMM_WORLD, &(send_req[iMessage]));
      iMessage++;
    }
  }

  /*--- Repeat the process to communicate the global IDs. ---*/

  iMessage = 0;
  for (int ii=0; ii<size; ii++) {
    if ((ii != rank) && (nPoint_Recv[ii+1] > nPoint_Recv[ii])) {
      int ll     = nPoint_Recv[ii];
      int kk     = nPoint_Recv[ii+1] - nPoint_Recv[ii];
      int count  = kk;
      int source = ii;
      int tag    = ii + 1;
      SU2_MPI::Irecv(&(idRecv[ll]), count, MPI_UNSIGNED_LONG, source, tag,
                     MPI_COMM_WORLD, &(recv_req[iMessage+nRecvs]));
      iMessage++;
    }
  }

  /*--- Launch the non-blocking sends of the global IDs. ---*/

  iMessage = 0;
  for (int ii=0; ii<size; ii++) {
    if ((ii != rank) && (nPoint_Send[ii+1] > nPoint_Send[ii])) {
      int ll = nPoint_Send[ii];
      int kk = nPoint_Send[ii+1] - nPoint_Send[ii];
      int count  = kk;
      int dest   = ii;
      int tag    = rank + 1;
      SU2_MPI::Isend(&(idSend[ll]), count, MPI_UNSIGNED_LONG, dest, tag,
                     MPI_COMM_WORLD, &(send_req[iMessage+nSends]));
      iMessage++;
    }
  }
#endif

  /*--- Copy my own rank's data into the recv buffer directly. ---*/

  int mm = VARS_PER_POINT*nPoint_Recv[rank];
  int ll = VARS_PER_POINT*nPoint_Send[rank];
  int kk = VARS_PER_POINT*nPoint_Send[rank+1];

  for (int nn=ll; nn<kk; nn++, mm++) connRecv[mm] = connSend[nn];

  mm = nPoint_Recv[rank];
  ll = nPoint_Send[rank];
  kk = nPoint_Send[rank+1];

  for (int nn=ll; nn<kk; nn++, mm++) idRecv[mm] = idSend[nn];

  /*--- Wait for the non-blocking sends and recvs to complete. ---*/

#ifdef HAVE_MPI
  int number = 2*nSends;
  for (int ii = 0; ii < number; ii++)
    SU2_MPI::Waitany(number, send_req, &ind, &status);

  number = 2*nRecvs;
  for (int ii = 0; ii < number; ii++)
    SU2_MPI::Waitany(number, recv_req, &ind, &status);

  delete [] send_req;
  delete [] recv_req;
#endif

  /*--- Store the connectivity for this rank in the proper data
   structure before post-processing below. First, allocate the
   appropriate amount of memory for this section. ---*/

  Parallel_Data = new su2double*[VARS_PER_POINT];
  for (int jj = 0; jj < VARS_PER_POINT; jj++) {
    Parallel_Data[jj] = new su2double[nPoint_Recv[size]];
    for (int ii = 0; ii < nPoint_Recv[size]; ii++) {
      Parallel_Data[jj][idRecv[ii]] = connRecv[ii*VARS_PER_POINT+jj];
    }
  }

  /*--- Store the total number of local points my rank has for
   the current section after completing the communications. ---*/

  nParallel_Poin = nPoint_Recv[size];

  /*--- Reduce the total number of points we will write in the output files. ---*/

#ifndef HAVE_MPI
  nGlobal_Poin_Par = nParallel_Poin;
#else
  SU2_MPI::Allreduce(&nParallel_Poin, &nGlobal_Poin_Par, 1,
                     MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#endif

  /*--- Free temporary memory from communications ---*/

  delete [] connSend;
  delete [] connRecv;
  delete [] idSend;
  delete [] idRecv;
  delete [] nPoint_Recv;
  delete [] nPoint_Send;
  delete [] nPoint_Flag;

  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
    delete [] Local_Data[iPoint];
  delete [] Local_Data;

  delete [] Local_Halo;
  delete [] npoint_procs;
  delete [] starting_node;
  delete [] ending_node;
  delete [] nPoint_Linear;

}

void COutputLegacy::SortOutputData_Surface(CConfig *config, CGeometry *geometry) {

  unsigned short iMarker;
  unsigned long iProcessor;
  unsigned long iPoint, jPoint, kPoint, iElem;
  unsigned long Global_Index, nLocalPoint, nTotalPoint, iVertex;

  int VARS_PER_POINT = nVar_Par;
  int *Local_Halo = NULL;
  int iNode, count;
  int SendRecv, RecvFrom;

  bool notPeriodic, notHalo, addedPeriodic, isPeriodic;

#ifdef HAVE_MPI
  SU2_MPI::Request *send_req, *recv_req;
  SU2_MPI::Status status;
  int ind;
#endif

  /*--------------------------------------------------------------------------*/
  /*--- Step 1: We already have the surface connectivity spread out in     ---*/
  /*---         linear partitions across all procs and the output data     ---*/
  /*---         for the entire field is similarly linearly partitioned.    ---*/
  /*---         We need to identify which nodes in the volume data are     ---*/
  /*---         also surface points. Our first step is to loop over all    ---*/
  /*---         of the sorted surface connectivity and create a data       ---*/
  /*---         structure on each proc that can identify the local surf    ---*/
  /*---         points. Note that the linear partitioning is slightly      ---*/
  /*---         different between the nodes and elements, so we will       ---*/
  /*---         have to move between the two systems in this routine.      ---*/
  /*--------------------------------------------------------------------------*/

  /*--- Search all send/recv boundaries on this partition for any periodic
   nodes that were part of the original domain. We want to recover these
   for visualization purposes. This is the linear partitioning for nodes. ---*/

  Local_Halo = new int[geometry->GetnPoint()];
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
    Local_Halo[iPoint] = !geometry->node[iPoint]->GetDomain();

  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) {

      /*--- Checking for less than or equal to the rank, because there may
       be some periodic halo nodes that send info to the same rank. ---*/

      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        isPeriodic = ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() > 0) &&
                      (geometry->vertex[iMarker][iVertex]->GetRotation_Type() % 2 == 1));
        if (isPeriodic) Local_Halo[iPoint] = false;
      }
    }
  }

  /*--- Sum total number of nodes that belong to the domain ---*/

  nLocalPoint = 0;
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
    if (Local_Halo[iPoint] == false)
      nLocalPoint++;

#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&nLocalPoint, &nTotalPoint, 1,
                     MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#else
  nTotalPoint = nLocalPoint;
#endif

  /*--- Now that we know the actual number of points we need to output,
   compute the number of points that will be on each processor.
   This is a linear partitioning with the addition of a simple load
   balancing for any remainder points. ---*/

  unsigned long *npoint_procs  = new unsigned long[size];
  unsigned long *starting_node = new unsigned long[size];
  unsigned long *ending_node   = new unsigned long[size];

  unsigned long *nPoint_Linear_Nodes = new unsigned long[size+1];
  unsigned long *nPoint_Linear_Elems = new unsigned long[size+1];

  unsigned long total_pt_accounted = 0;
  for (int ii = 0; ii < size; ii++) {
    npoint_procs[ii] = nTotalPoint/size;
    total_pt_accounted = total_pt_accounted + npoint_procs[ii];
  }

  /*--- Get the number of remainder points after the even division. ---*/

  unsigned long rem_points = nTotalPoint-total_pt_accounted;
  for (unsigned long ii = 0; ii < rem_points; ii++) {
    npoint_procs[ii]++;
  }

  /*--- Store the local number of nodes and the beginning/end index ---*/

  starting_node[0] = 0;
  ending_node[0]   = starting_node[0] + npoint_procs[0];
  nPoint_Linear_Nodes[0] = 0;
  for (int ii = 1; ii < size; ii++) {
    starting_node[ii] = ending_node[ii-1];
    ending_node[ii]   = starting_node[ii] + npoint_procs[ii];
    nPoint_Linear_Nodes[ii] = nPoint_Linear_Nodes[ii-1] + npoint_procs[ii-1];
  }
  nPoint_Linear_Nodes[size] = nTotalPoint;

  /*--- Prepare to check and communicate the nodes that each proc has
   locally from the surface connectivity. ---*/

  int *nElem_Send = new int[size+1]; nElem_Send[0] = 0;
  int *nElem_Recv = new int[size+1]; nElem_Recv[0] = 0;
  int *nElem_Flag = new int[size];

  for (int ii=0; ii < size; ii++) {
    nElem_Send[ii] = 0;
    nElem_Recv[ii] = 0;
    nElem_Flag[ii]= -1;
  }
  nElem_Send[size] = 0; nElem_Recv[size] = 0;

  /*--- Loop through our local line elements and check where each
   of the grid nodes resides based on global index. ---*/

  for (int ii = 0; ii < (int)nParallel_Line; ii++) {
    for ( int jj = 0; jj < N_POINTS_LINE; jj++ ) {

      /*--- Get global index. Note the -1 since it was 1-based for viz. ---*/

      iNode = ii*N_POINTS_LINE+jj;
      Global_Index = Conn_BoundLine_Par[iNode]-1;

      /*--- Search for the processor that owns this point ---*/

      iProcessor = Global_Index/npoint_procs[0];
      if (iProcessor >= (unsigned long)size)
        iProcessor = (unsigned long)size-1;
      if (Global_Index >= nPoint_Linear_Nodes[iProcessor])
        while(Global_Index >= nPoint_Linear_Nodes[iProcessor+1]) iProcessor++;
      else
        while(Global_Index <  nPoint_Linear_Nodes[iProcessor])   iProcessor--;

      /*--- If we have not visited this element yet, increment our
       number of elements that must be sent to a particular proc. ---*/

      if ((nElem_Flag[iProcessor] != iNode)) {
        nElem_Flag[iProcessor] = iNode;
        nElem_Send[iProcessor+1]++;
      }

    }
  }

  /*--- Reset out flags and then loop through our local triangle surface
   elements performing the same check for where each grid node resides. ---*/

  for (int ii=0; ii < size; ii++) nElem_Flag[ii]= -1;

  for (int ii = 0; ii < (int)nParallel_BoundTria; ii++) {
    for ( int jj = 0; jj < N_POINTS_TRIANGLE; jj++ ) {

      /*--- Get global index. Note the -1 since it was 1-based for viz. ---*/

      iNode = ii*N_POINTS_TRIANGLE + jj;
      Global_Index = Conn_BoundTria_Par[iNode]-1;

      /*--- Search for the processor that owns this point ---*/

      iProcessor = Global_Index/npoint_procs[0];
      if (iProcessor >= (unsigned long)size)
        iProcessor = (unsigned long)size-1;
      if (Global_Index >= nPoint_Linear_Nodes[iProcessor])
        while(Global_Index >= nPoint_Linear_Nodes[iProcessor+1]) iProcessor++;
      else
        while(Global_Index <  nPoint_Linear_Nodes[iProcessor])   iProcessor--;

      /*--- If we have not visited this element yet, increment our
       number of elements that must be sent to a particular proc. ---*/

      if ((nElem_Flag[iProcessor] != iNode)) {
        nElem_Flag[iProcessor] = iNode;
        nElem_Send[iProcessor+1]++;
      }

    }
  }

  /*--- Reset out flags and then loop through our local quad surface
   elements performing the same check for where each grid node resides. ---*/

  for (int ii=0; ii < size; ii++) nElem_Flag[ii]= -1;

  for (int ii = 0; ii < (int)nParallel_BoundQuad; ii++) {
    for ( int jj = 0; jj < N_POINTS_QUADRILATERAL; jj++ ) {

      /*--- Get global index. Note the -1 since it was 1-based for viz. ---*/

      iNode = ii*N_POINTS_QUADRILATERAL+jj;
      Global_Index = Conn_BoundQuad_Par[iNode]-1;

      /*--- Search for the processor that owns this point ---*/

      iProcessor = Global_Index/npoint_procs[0];
      if (iProcessor >= (unsigned long)size)
        iProcessor = (unsigned long)size-1;
      if (Global_Index >= nPoint_Linear_Nodes[iProcessor])
        while(Global_Index >= nPoint_Linear_Nodes[iProcessor+1]) iProcessor++;
      else
        while(Global_Index <  nPoint_Linear_Nodes[iProcessor])   iProcessor--;

      /*--- If we have not visited this element yet, increment our
       number of elements that must be sent to a particular proc. ---*/

      if ((nElem_Flag[iProcessor] != iNode)) {
        nElem_Flag[iProcessor] = iNode;
        nElem_Send[iProcessor+1]++;
      }

    }
  }

  /*--- Communicate the number of nodes to be sent/recv'd amongst
   all processors. After this communication, each proc knows how
   many nodes it will receive from each other processor. ---*/

#ifdef HAVE_MPI
  SU2_MPI::Alltoall(&(nElem_Send[1]), 1, MPI_INT,
                    &(nElem_Recv[1]), 1, MPI_INT, MPI_COMM_WORLD);
#else
  nElem_Recv[1] = nElem_Send[1];
#endif

  /*--- Prepare to send. First check how many
   messages we will be sending and receiving. Here we also put
   the counters into cumulative storage format to make the
   communications simpler. ---*/

  int nSends = 0, nRecvs = 0;
  for (int ii=0; ii < size; ii++) nElem_Flag[ii] = -1;

  for (int ii = 0; ii < size; ii++) {
    if ((ii != rank) && (nElem_Send[ii+1] > 0)) nSends++;
    if ((ii != rank) && (nElem_Recv[ii+1] > 0)) nRecvs++;

    nElem_Send[ii+1] += nElem_Send[ii];
    nElem_Recv[ii+1] += nElem_Recv[ii];
  }

  /*--- Allocate arrays for sending the global ID. ---*/

  unsigned long *idSend = new unsigned long[nElem_Send[size]];
  for (int ii = 0; ii < nElem_Send[size]; ii++) idSend[ii] = 0;

  /*--- Create an index variable to keep track of our index
   positions as we load up the send buffer. ---*/

  unsigned long *idIndex = new unsigned long[size];
  for (int ii=0; ii < size; ii++) idIndex[ii] = nElem_Send[ii];

  /*--- Now loop back through the local connectivities for the surface
   elements and load up the global IDs for sending to their home proc. ---*/

  for (int ii = 0; ii < (int)nParallel_Line; ii++) {
    for ( int jj = 0; jj < N_POINTS_LINE; jj++ ) {

      /*--- Get global index. Note the -1 since it was 1-based for viz. ---*/

      iNode = ii*N_POINTS_LINE+jj;
      Global_Index = Conn_BoundLine_Par[iNode]-1;

      /*--- Search for the processor that owns this point ---*/

      iProcessor = Global_Index/npoint_procs[0];
      if (iProcessor >= (unsigned long)size)
        iProcessor = (unsigned long)size-1;
      if (Global_Index >= nPoint_Linear_Nodes[iProcessor])
        while(Global_Index >= nPoint_Linear_Nodes[iProcessor+1]) iProcessor++;
      else
        while(Global_Index <  nPoint_Linear_Nodes[iProcessor])   iProcessor--;

      /*--- Load global ID into the buffer for sending ---*/

      if (nElem_Flag[iProcessor] != iNode) {

        nElem_Flag[iProcessor] = iNode;
        unsigned long nn = idIndex[iProcessor];

        /*--- Load the connectivity values. ---*/

        idSend[nn] = Global_Index; nn++;

        /*--- Increment the index by the message length ---*/

        idIndex[iProcessor]++;

      }

    }
  }

  for (int ii=0; ii < size; ii++) nElem_Flag[ii]= -1;

  for (int ii = 0; ii < (int)nParallel_BoundTria; ii++) {
    for ( int jj = 0; jj < N_POINTS_TRIANGLE; jj++ ) {

      /*--- Get global index. Note the -1 since it was 1-based for viz. ---*/

      iNode = ii*N_POINTS_TRIANGLE + jj;
      Global_Index = Conn_BoundTria_Par[iNode]-1;

      /*--- Search for the processor that owns this point ---*/

      iProcessor = Global_Index/npoint_procs[0];
      if (iProcessor >= (unsigned long)size)
        iProcessor = (unsigned long)size-1;
      if (Global_Index >= nPoint_Linear_Nodes[iProcessor])
        while(Global_Index >= nPoint_Linear_Nodes[iProcessor+1]) iProcessor++;
      else
        while(Global_Index <  nPoint_Linear_Nodes[iProcessor])   iProcessor--;

      /*--- Load global ID into the buffer for sending ---*/

      if (nElem_Flag[iProcessor] != iNode) {

        nElem_Flag[iProcessor] = iNode;
        unsigned long nn = idIndex[iProcessor];

        /*--- Load the connectivity values. ---*/

        idSend[nn] = Global_Index; nn++;

        /*--- Increment the index by the message length ---*/

        idIndex[iProcessor]++;

      }

    }
  }

  for (int ii=0; ii < size; ii++) nElem_Flag[ii]= -1;

  for (int ii = 0; ii < (int)nParallel_BoundQuad; ii++) {
    for ( int jj = 0; jj < N_POINTS_QUADRILATERAL; jj++ ) {

      /*--- Get global index. Note the -1 since it was 1-based for viz. ---*/

      iNode = ii*N_POINTS_QUADRILATERAL+jj;
      Global_Index = Conn_BoundQuad_Par[iNode]-1;

      /*--- Search for the processor that owns this point ---*/

      iProcessor = Global_Index/npoint_procs[0];
      if (iProcessor >= (unsigned long)size)
        iProcessor = (unsigned long)size-1;
      if (Global_Index >= nPoint_Linear_Nodes[iProcessor])
        while(Global_Index >= nPoint_Linear_Nodes[iProcessor+1]) iProcessor++;
      else
        while(Global_Index <  nPoint_Linear_Nodes[iProcessor])   iProcessor--;

      /*--- Load global ID into the buffer for sending ---*/

      if (nElem_Flag[iProcessor] != iNode) {

        nElem_Flag[iProcessor] = iNode;
        unsigned long nn = idIndex[iProcessor];

        /*--- Load the connectivity values. ---*/

        idSend[nn] = Global_Index; nn++;

        /*--- Increment the index by the message length ---*/

        idIndex[iProcessor]++;

      }

    }
  }

  /*--- Allocate the memory that we need for receiving the global IDs
   values and then cue up the non-blocking receives. Note that
   we do not include our own rank in the communications. We will
   directly copy our own data later. ---*/

  unsigned long *idRecv = NULL;
  idRecv = new unsigned long[nElem_Recv[size]];
  for (int ii = 0; ii < nElem_Recv[size]; ii++)
    idRecv[ii] = 0;

#ifdef HAVE_MPI
  /*--- We need double the number of messages to send both the conn.
   and the flags for the halo cells. ---*/

  send_req = new SU2_MPI::Request[nSends];
  recv_req = new SU2_MPI::Request[nRecvs];

  /*--- Launch the non-blocking recv's for the global IDs. ---*/

  unsigned long iMessage = 0;
  for (int ii=0; ii<size; ii++) {
    if ((ii != rank) && (nElem_Recv[ii+1] > nElem_Recv[ii])) {
      int ll     = nElem_Recv[ii];
      int kk     = nElem_Recv[ii+1] - nElem_Recv[ii];
      int count  = kk;
      int source = ii;
      int tag    = ii + 1;
      SU2_MPI::Irecv(&(idRecv[ll]), count, MPI_UNSIGNED_LONG, source, tag,
                     MPI_COMM_WORLD, &(recv_req[iMessage]));
      iMessage++;
    }
  }

  /*--- Launch the non-blocking sends of the global IDs. ---*/

  iMessage = 0;
  for (int ii=0; ii<size; ii++) {
    if ((ii != rank) && (nElem_Send[ii+1] > nElem_Send[ii])) {
      int ll = nElem_Send[ii];
      int kk = nElem_Send[ii+1] - nElem_Send[ii];
      int count  = kk;
      int dest = ii;
      int tag    = rank + 1;
      SU2_MPI::Isend(&(idSend[ll]), count, MPI_UNSIGNED_LONG, dest, tag,
                     MPI_COMM_WORLD, &(send_req[iMessage]));
      iMessage++;
    }
  }
#endif

  /*--- Copy my own rank's data into the recv buffer directly. ---*/

  int mm = nElem_Recv[rank];
  int ll = nElem_Send[rank];
  int kk = nElem_Send[rank+1];

  for (int nn=ll; nn<kk; nn++, mm++) idRecv[mm] = idSend[nn];

  /*--- Wait for the non-blocking sends and recvs to complete. ---*/

#ifdef HAVE_MPI
  int number = nSends;
  for (int ii = 0; ii < number; ii++)
    SU2_MPI::Waitany(number, send_req, &ind, &status);

  number = nRecvs;
  for (int ii = 0; ii < number; ii++)
    SU2_MPI::Waitany(number, recv_req, &ind, &status);

  delete [] send_req;
  delete [] recv_req;
#endif

  /*--------------------------------------------------------------------------*/
  /*--- Step 2: Each proc now knows which is its local grid nodes from     ---*/
  /*---         the entire volume solution are part of the surface. We     ---*/
  /*---         now apply a mask to extract just those points on the       ---*/
  /*---         surface. We also need to perform a renumbering so that     ---*/
  /*---         the surface data (nodes and connectivity) have their       ---*/
  /*---         own global numbering. This is important for writing        ---*/
  /*---         output files in a later routine.                           ---*/
  /*--------------------------------------------------------------------------*/

  /*--- Create a local data structure that acts as a mask to extract the
   set of points within the local set that are on the surface. ---*/

  int *surfPoint = new int[nParallel_Poin];
  for (iPoint = 0; iPoint < nParallel_Poin; iPoint++) surfPoint[iPoint] = -1;

  for (int ii = 0; ii < nElem_Recv[size]; ii++) {
    surfPoint[(int)idRecv[ii]- starting_node[rank]] = (int)idRecv[ii];
  }

  /*--- First, add up the number of surface points I have on my rank. ---*/

  nSurf_Poin_Par = 0;
  for (iPoint = 0; iPoint < nParallel_Poin; iPoint++) {
    if (surfPoint[iPoint] != -1) {
      nSurf_Poin_Par++;
    }
  }

  /*--- Communicate this number of local surface points to all other
   processors so that it can be used to create offsets for the new
   global numbering for the surface points. ---*/

  int *nPoint_Send = new int[size+1]; nPoint_Send[0] = 0;
  int *nPoint_Recv = new int[size+1]; nPoint_Recv[0] = 0;

  for (int ii=1; ii < size+1; ii++) nPoint_Send[ii]= (int)nSurf_Poin_Par;

#ifdef HAVE_MPI
  SU2_MPI::Alltoall(&(nPoint_Send[1]), 1, MPI_INT,
                    &(nPoint_Recv[1]), 1, MPI_INT, MPI_COMM_WORLD);
#else
  nPoint_Recv[1] = nPoint_Send[1];
#endif

  /*--- Go to cumulative storage format to compute the offsets. ---*/

  for (int ii = 0; ii < size; ii++) {
    nPoint_Send[ii+1] += nPoint_Send[ii];
    nPoint_Recv[ii+1] += nPoint_Recv[ii];
  }

  /*--- Now that we know the number of local surface points that we have,
   we can allocate the new data structure to hold these points alone. Here,
   we also copy the data for those points from our volume data structure. ---*/

  Parallel_Surf_Data = new su2double*[VARS_PER_POINT];
  for (int jj = 0; jj < VARS_PER_POINT; jj++) {
    Parallel_Surf_Data[jj] = new su2double[nSurf_Poin_Par];
    count = 0;
    for (int ii = 0; ii < (int)nParallel_Poin; ii++) {
      if (surfPoint[ii] !=-1) {
        Parallel_Surf_Data[jj][count] = Parallel_Data[jj][ii];
        count++;
      }
    }
  }

  /*--- Reduce the total number of surf points we have. This will be
   needed for writing the surface solution files later. ---*/

#ifndef HAVE_MPI
  nGlobal_Surf_Poin = nSurf_Poin_Par;
#else
  SU2_MPI::Allreduce(&nSurf_Poin_Par, &nGlobal_Surf_Poin, 1,
                     MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#endif

  /*--- Now that we know every proc's global offset for the number of
   surface points, we can create the new global numbering. Here, we
   create a new mapping using two arrays, which will need to be
   communicated. We use our mask again here.  ---*/

  unsigned long *globalP = new unsigned long[nSurf_Poin_Par];
  unsigned long *renumbP = new unsigned long[nSurf_Poin_Par];

  count = 0;
  for (iPoint = 0; iPoint < nParallel_Poin; iPoint++) {
    if (surfPoint[iPoint] != -1) {
      globalP[count] = surfPoint[iPoint];
      renumbP[count] = count + nPoint_Recv[rank];
      count++;
    }
  }

  /*--------------------------------------------------------------------------*/
  /*--- Step 3: Communicate the arrays with the new global surface point   ---*/
  /*---         numbering to the procs that hold the connectivity for      ---*/
  /*---         each element. This will be done in two phases. First,      ---*/
  /*---         we send the arrays around to the other procs based on      ---*/
  /*---         the linear partitioning for the elems. This gets us        ---*/
  /*---         most of the way there, however, due to the type of         ---*/
  /*---         linear partitioning for the elements, there may exist      ---*/
  /*---         elements that have nodes outside of the linear part.       ---*/
  /*---         bounds. This is because the elems are distributed based    ---*/
  /*---         on the node with the smallest global ID.                   ---*/
  /*--------------------------------------------------------------------------*/

  /*--- First, we perform the linear partitioning again as it is done
   for elements, which is slightly different than for nodes (above). ---*/

  /*--- Force the removal of all added periodic elements (use global index).
   First, we isolate and create a list of all added periodic points, excluding
   those that were part of the original domain (we want these to be in the
   output files). ---*/

  vector<unsigned long> Added_Periodic;
  Added_Periodic.clear();

  if (config->GetKind_SU2() != SU2_DEF) {
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) {
        SendRecv = config->GetMarker_All_SendRecv(iMarker);
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

          if ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() > 0) &&
              (geometry->vertex[iMarker][iVertex]->GetRotation_Type() % 2 == 0) &&
              (SendRecv < 0)) {
            Added_Periodic.push_back(geometry->node[iPoint]->GetGlobalIndex());
          }
        }
      }
    }
  }

  /*--- Now we communicate this information to all processors, so that they
   can force the removal of these particular nodes by flagging them as halo
   points. In general, this should be a small percentage of the total mesh,
   so the communication/storage costs here shouldn't be prohibitive. ---*/

  /*--- First communicate the number of points that each rank has found. ---*/

  unsigned long nAddedPeriodic = 0, maxAddedPeriodic = 0;
  unsigned long Buffer_Send_nAddedPeriodic[1], *Buffer_Recv_nAddedPeriodic = NULL;
  Buffer_Recv_nAddedPeriodic = new unsigned long[size];

  nAddedPeriodic = Added_Periodic.size();
  Buffer_Send_nAddedPeriodic[0] = nAddedPeriodic;

#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&nAddedPeriodic, &maxAddedPeriodic, 1, MPI_UNSIGNED_LONG,
                     MPI_MAX, MPI_COMM_WORLD);
  SU2_MPI::Allgather(&Buffer_Send_nAddedPeriodic, 1, MPI_UNSIGNED_LONG,
                     Buffer_Recv_nAddedPeriodic,  1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
#else
  maxAddedPeriodic = nAddedPeriodic;
  Buffer_Recv_nAddedPeriodic[0] = Buffer_Send_nAddedPeriodic[0];
#endif

  /*--- Communicate the global index values of all added periodic nodes. ---*/
  unsigned long *Buffer_Send_AddedPeriodic = new unsigned long[maxAddedPeriodic];
  unsigned long *Buffer_Recv_AddedPeriodic = new unsigned long[size*maxAddedPeriodic];

  for (iPoint = 0; iPoint < Added_Periodic.size(); iPoint++) {
    Buffer_Send_AddedPeriodic[iPoint] = Added_Periodic[iPoint];
  }

  /*--- Gather the element connectivity information. All processors will now
   have a copy of the global index values for all added periodic points. ---*/

#ifdef HAVE_MPI
  SU2_MPI::Allgather(Buffer_Send_AddedPeriodic, maxAddedPeriodic, MPI_UNSIGNED_LONG,
                     Buffer_Recv_AddedPeriodic, maxAddedPeriodic, MPI_UNSIGNED_LONG,
                     MPI_COMM_WORLD);
#else
  for (iPoint = 0; iPoint < maxAddedPeriodic; iPoint++)
    Buffer_Recv_AddedPeriodic[iPoint] = Buffer_Send_AddedPeriodic[iPoint];
#endif

  /*--- Search all send/recv boundaries on this partition for halo cells. In
   particular, consider only the recv conditions (these are the true halo
   nodes). Check the ranks of the processors that are communicating and
   choose to keep only the halo cells from the higher rank processor. Here,
   we are also choosing to keep periodic nodes that were part of the original
   domain. We will check the communicated list of added periodic points. ---*/

  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
    Local_Halo[iPoint] = !geometry->node[iPoint]->GetDomain();

  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) {
      SendRecv = config->GetMarker_All_SendRecv(iMarker);
      RecvFrom = abs(SendRecv)-1;

      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        Global_Index = geometry->node[iPoint]->GetGlobalIndex();

        /*--- We need to keep one copy of overlapping halo cells. ---*/

        notHalo = ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() == 0) &&
                   (SendRecv < 0) && (rank > RecvFrom));

        /*--- We want to keep the periodic nodes that were part of the original domain.
         For SU2_DEF we want to keep all periodic nodes. ---*/

        if (config->GetKind_SU2() == SU2_DEF) {
          isPeriodic = ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() > 0));
        }else {
          isPeriodic = ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() > 0) &&
                        (geometry->vertex[iMarker][iVertex]->GetRotation_Type() % 2 == 1));
        }

        notPeriodic = (isPeriodic && (SendRecv < 0));

        /*--- Lastly, check that this isn't an added periodic point that
         we will forcibly remove. Use the communicated list of these points. ---*/

        addedPeriodic = false; kPoint = 0;
        for (iProcessor = 0; iProcessor < (unsigned long)size; iProcessor++) {
          for (jPoint = 0; jPoint < Buffer_Recv_nAddedPeriodic[iProcessor]; jPoint++) {
            if (Global_Index == Buffer_Recv_AddedPeriodic[kPoint+jPoint])
              addedPeriodic = true;
          }

          /*--- Adjust jNode to index of next proc's data in the buffers. ---*/

          kPoint = (iProcessor+1)*maxAddedPeriodic;

        }

        /*--- If we found either of these types of nodes, flag them to be kept. ---*/

        if ((notHalo || notPeriodic) && !addedPeriodic) {
          Local_Halo[iPoint] = false;
        }

      }
    }
  }

  /*--- Now that we've done the gymnastics to find any periodic points,
   compute the total number of local and global points for the output. ---*/

  nLocalPoint = 0;
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
    if (Local_Halo[iPoint] == false)
      nLocalPoint++;

#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&nLocalPoint, &nTotalPoint, 1,
                     MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#else
  nTotalPoint = nLocalPoint;
#endif

  /*--- Compute the number of points that will be on each processor.
   This is a linear partitioning with the addition of a simple load
   balancing for any remainder points. ---*/

  total_pt_accounted = 0;
  for (int ii = 0; ii < size; ii++) {
    npoint_procs[ii] = nTotalPoint/size;
    total_pt_accounted = total_pt_accounted + npoint_procs[ii];
  }

  /*--- Get the number of remainder points after the even division. ---*/

  rem_points = nTotalPoint-total_pt_accounted;
  for (unsigned long ii = 0; ii < rem_points; ii++) {
    npoint_procs[ii]++;
  }

  /*--- Store the local number of nodes and the beginning/end index ---*/

  starting_node[0] = 0;
  ending_node[0]   = starting_node[0] + npoint_procs[0];
  nPoint_Linear_Elems[0] = 0;
  for (int ii = 1; ii < size; ii++) {
    starting_node[ii] = ending_node[ii-1];
    ending_node[ii]   = starting_node[ii] + npoint_procs[ii];
    nPoint_Linear_Elems[ii] = nPoint_Linear_Elems[ii-1] + npoint_procs[ii-1];
  }
  nPoint_Linear_Elems[size] = nTotalPoint;

  /*--- Reset our flags and counters ---*/

  for (int ii=0; ii < size; ii++) {
    nElem_Send[ii] = 0;
    nElem_Recv[ii] = 0;
    nElem_Flag[ii]= -1;
  }
  nElem_Send[size] = 0; nElem_Recv[size] = 0;

  /*--- Loop through my local surface nodes, find which proc the global
   value lives on, then communicate the global ID and remumbered value. ---*/

  for (int ii = 0; ii < (int)nSurf_Poin_Par; ii++) {

    Global_Index = globalP[ii];

    /*--- Search for the processor that owns this point ---*/

    iProcessor = Global_Index/npoint_procs[0];
    if (iProcessor >= (unsigned long)size)
      iProcessor = (unsigned long)size-1;
    if (Global_Index >= nPoint_Linear_Elems[iProcessor])
      while(Global_Index >= nPoint_Linear_Elems[iProcessor+1]) iProcessor++;
    else
      while(Global_Index <  nPoint_Linear_Elems[iProcessor])   iProcessor--;

    /*--- If we have not visited this element yet, increment our
     number of elements that must be sent to a particular proc. ---*/

    if ((nElem_Flag[iProcessor] != ii)) {
      nElem_Flag[iProcessor] = ii;
      nElem_Send[iProcessor+1]++;
    }

  }

  /*--- Communicate the number of cells to be sent/recv'd amongst
   all processors. After this communication, each proc knows how
   many cells it will receive from each other processor. ---*/

#ifdef HAVE_MPI
  SU2_MPI::Alltoall(&(nElem_Send[1]), 1, MPI_INT,
                    &(nElem_Recv[1]), 1, MPI_INT, MPI_COMM_WORLD);
#else
  nElem_Recv[1] = nElem_Send[1];
#endif

  /*--- Prepare to send. First check how many
   messages we will be sending and receiving. Here we also put
   the counters into cumulative storage format to make the
   communications simpler. ---*/

  nSends = 0; nRecvs = 0;
  for (int ii=0; ii < size; ii++) nElem_Flag[ii] = -1;

  for (int ii = 0; ii < size; ii++) {
    if ((ii != rank) && (nElem_Send[ii+1] > 0)) nSends++;
    if ((ii != rank) && (nElem_Recv[ii+1] > 0)) nRecvs++;

    nElem_Send[ii+1] += nElem_Send[ii];
    nElem_Recv[ii+1] += nElem_Recv[ii];
  }

  /*--- Allocate memory to hold the globals that we are
   sending. ---*/

  unsigned long *globalSend = NULL;
  globalSend = new unsigned long[nElem_Send[size]];
  for (int ii = 0; ii < nElem_Send[size]; ii++)
    globalSend[ii] = 0;

  /*--- Allocate memory to hold the renumbering that we are
   sending. ---*/

  unsigned long *renumbSend = NULL;
  renumbSend = new unsigned long[nElem_Send[size]];
  for (int ii = 0; ii < nElem_Send[size]; ii++)
    renumbSend[ii] = 0;

  /*--- Create an index variable to keep track of our index
   position as we load up the send buffer. ---*/

  unsigned long *index = new unsigned long[size];
  for (int ii=0; ii < size; ii++) index[ii] = nElem_Send[ii];

  /*--- Loop back through and load up the buffers for the global IDs
   and their new renumbering values. ---*/

  for (int ii = 0; ii < (int)nSurf_Poin_Par; ii++) {

    Global_Index = globalP[ii];

    /*--- Search for the processor that owns this point ---*/

    iProcessor = Global_Index/npoint_procs[0];
    if (iProcessor >= (unsigned long)size)
      iProcessor = (unsigned long)size-1;
    if (Global_Index >= nPoint_Linear_Elems[iProcessor])
      while(Global_Index >= nPoint_Linear_Elems[iProcessor+1]) iProcessor++;
    else
      while(Global_Index <  nPoint_Linear_Elems[iProcessor])   iProcessor--;


    if (nElem_Flag[iProcessor] != ii) {

      nElem_Flag[iProcessor] = ii;
      unsigned long nn = index[iProcessor];

      globalSend[nn] = Global_Index;
      renumbSend[nn] = renumbP[ii];

      /*--- Increment the index by the message length ---*/

      index[iProcessor]++;

    }
  }

  /*--- Free memory after loading up the send buffer. ---*/

  delete [] index;

  /*--- Allocate the memory that we need for receiving the
   values and then cue up the non-blocking receives. Note that
   we do not include our own rank in the communications. We will
   directly copy our own data later. ---*/

  unsigned long *globalRecv = NULL;
  globalRecv = new unsigned long[nElem_Recv[size]];
  for (int ii = 0; ii < nElem_Recv[size]; ii++)
    globalRecv[ii] = 0;

  unsigned long *renumbRecv = NULL;
  renumbRecv = new unsigned long[nElem_Recv[size]];
  for (int ii = 0; ii < nElem_Recv[size]; ii++)
    renumbRecv[ii] = 0;

#ifdef HAVE_MPI
  /*--- We need double the number of messages to send both the conn.
   and the flags for the halo cells. ---*/

  send_req = new SU2_MPI::Request[2*nSends];
  recv_req = new SU2_MPI::Request[2*nRecvs];

  /*--- Launch the non-blocking recv's for the global ID. ---*/

  iMessage = 0;
  for (int ii=0; ii<size; ii++) {
    if ((ii != rank) && (nElem_Recv[ii+1] > nElem_Recv[ii])) {
      int ll     = nElem_Recv[ii];
      int kk     = nElem_Recv[ii+1] - nElem_Recv[ii];
      int count  = kk;
      int source = ii;
      int tag    = ii + 1;
      SU2_MPI::Irecv(&(globalRecv[ll]), count, MPI_UNSIGNED_LONG, source, tag,
                     MPI_COMM_WORLD, &(recv_req[iMessage]));
      iMessage++;
    }
  }

  /*--- Launch the non-blocking sends of the global ID. ---*/

  iMessage = 0;
  for (int ii=0; ii<size; ii++) {
    if ((ii != rank) && (nElem_Send[ii+1] > nElem_Send[ii])) {
      int ll = nElem_Send[ii];
      int kk = nElem_Send[ii+1] - nElem_Send[ii];
      int count  = kk;
      int dest = ii;
      int tag    = rank + 1;
      SU2_MPI::Isend(&(globalSend[ll]), count, MPI_UNSIGNED_LONG, dest, tag,
                     MPI_COMM_WORLD, &(send_req[iMessage]));
      iMessage++;
    }
  }

  /*--- Launch the non-blocking recv's for the renumbered ID. ---*/

  iMessage = 0;
  for (int ii=0; ii<size; ii++) {
    if ((ii != rank) && (nElem_Recv[ii+1] > nElem_Recv[ii])) {
      int ll     = nElem_Recv[ii];
      int kk     = nElem_Recv[ii+1] - nElem_Recv[ii];
      int count  = kk;
      int source = ii;
      int tag    = ii + 1;
      SU2_MPI::Irecv(&(renumbRecv[ll]), count, MPI_UNSIGNED_LONG, source, tag,
                     MPI_COMM_WORLD, &(recv_req[iMessage+nRecvs]));
      iMessage++;
    }
  }

  /*--- Launch the non-blocking sends of the renumbered ID. ---*/

  iMessage = 0;
  for (int ii=0; ii<size; ii++) {
    if ((ii != rank) && (nElem_Send[ii+1] > nElem_Send[ii])) {
      int ll = nElem_Send[ii];
      int kk = nElem_Send[ii+1] - nElem_Send[ii];
      int count  = kk;
      int dest = ii;
      int tag    = rank + 1;
      SU2_MPI::Isend(&(renumbSend[ll]), count, MPI_UNSIGNED_LONG, dest, tag,
                     MPI_COMM_WORLD, &(send_req[iMessage+nSends]));
      iMessage++;
    }
  }

#endif

  /*--- Load our own procs data into the buffers directly. ---*/

  mm = nElem_Recv[rank];
  ll = nElem_Send[rank];
  kk = nElem_Send[rank+1];

  for (int nn=ll; nn<kk; nn++, mm++) globalRecv[mm] = globalSend[nn];

  mm = nElem_Recv[rank];
  ll = nElem_Send[rank];
  kk = nElem_Send[rank+1];

  for (int nn=ll; nn<kk; nn++, mm++) renumbRecv[mm] = renumbSend[nn];

  /*--- Wait for the non-blocking sends and recvs to complete. ---*/

#ifdef HAVE_MPI
  number = 2*nSends;
  for (int ii = 0; ii < number; ii++)
    SU2_MPI::Waitany(number, send_req, &ind, &status);

  number = 2*nRecvs;
  for (int ii = 0; ii < number; ii++)
    SU2_MPI::Waitany(number, recv_req, &ind, &status);

  delete [] send_req;
  delete [] recv_req;
#endif

  /*-- Now update my local connectivitiy for the surface with the new
   numbering. Create a new mapping for global -> renumber for nodes. Note
   the adding of 1 back in here for the eventual viz. purposes. ---*/

  map<unsigned long,unsigned long> Global2Renumber;
  for (int ii = 0; ii < nElem_Recv[size]; ii++) {
    Global2Renumber[globalRecv[ii]] = renumbRecv[ii] + 1;
  }


  /*--- The final step is one last pass over all elements to check
   for points outside of the linear partitions of the elements. Again,
   note that elems were distributed based on their smallest global ID,
   so some nodes of the elem may have global IDs lying outside of the
   linear partitioning. We need to recover the mapping for these
   outliers. We loop over all local surface elements to find these. ---*/

  vector<unsigned long>::iterator it;
  vector<unsigned long> outliers;

  for (int ii = 0; ii < (int)nParallel_Line; ii++) {
    for ( int jj = 0; jj < N_POINTS_LINE; jj++ ) {

      iNode = ii*N_POINTS_LINE+jj;
      Global_Index = Conn_BoundLine_Par[iNode]-1;

      /*--- Search for the processor that owns this point ---*/

      iProcessor = Global_Index/npoint_procs[0];
      if (iProcessor >= (unsigned long)size)
        iProcessor = (unsigned long)size-1;
      if (Global_Index >= nPoint_Linear_Elems[iProcessor])
        while(Global_Index >= nPoint_Linear_Elems[iProcessor+1]) iProcessor++;
      else
        while(Global_Index <  nPoint_Linear_Elems[iProcessor])   iProcessor--;

      /*--- Store the global ID if it is outside our own linear partition. ---*/

      if ((iProcessor != (unsigned long)rank)) {
        outliers.push_back(Global_Index);
      }

    }
  }

  for (int ii=0; ii < size; ii++) nElem_Flag[ii]= -1;

  for (int ii = 0; ii < (int)nParallel_BoundTria; ii++) {
    for ( int jj = 0; jj < N_POINTS_TRIANGLE; jj++ ) {

      iNode = ii*N_POINTS_TRIANGLE + jj;
      Global_Index = Conn_BoundTria_Par[iNode]-1;

      /*--- Search for the processor that owns this point ---*/

      iProcessor = Global_Index/npoint_procs[0];
      if (iProcessor >= (unsigned long)size)
        iProcessor = (unsigned long)size-1;
      if (Global_Index >= nPoint_Linear_Elems[iProcessor])
        while(Global_Index >= nPoint_Linear_Elems[iProcessor+1]) iProcessor++;
      else
        while(Global_Index <  nPoint_Linear_Elems[iProcessor])   iProcessor--;

      /*--- Store the global ID if it is outside our own linear partition. ---*/

      if ((iProcessor != (unsigned long)rank)) {
        outliers.push_back(Global_Index);
      }

    }
  }

  for (int ii=0; ii < size; ii++) nElem_Flag[ii]= -1;

  for (int ii = 0; ii < (int)nParallel_BoundQuad; ii++) {
    for ( int jj = 0; jj < N_POINTS_QUADRILATERAL; jj++ ) {

      iNode = ii*N_POINTS_QUADRILATERAL+jj;
      Global_Index = Conn_BoundQuad_Par[iNode]-1;

      /*--- Search for the processor that owns this point ---*/

      iProcessor = Global_Index/npoint_procs[0];
      if (iProcessor >= (unsigned long)size)
        iProcessor = (unsigned long)size-1;
      if (Global_Index >= nPoint_Linear_Elems[iProcessor])
        while(Global_Index >= nPoint_Linear_Elems[iProcessor+1]) iProcessor++;
      else
        while(Global_Index <  nPoint_Linear_Elems[iProcessor])   iProcessor--;

      /*--- Store the global ID if it is outside our own linear partition. ---*/

      if ((iProcessor != (unsigned long)rank)) {
        outliers.push_back(Global_Index);
      }

    }
  }

  /*--- Create a unique list of global IDs that fall outside of our procs
   linear partition. ---*/

  sort(outliers.begin(), outliers.end());
  it = unique(outliers.begin(), outliers.end());
  outliers.resize(it - outliers.begin());

  /*--- Now loop over the outliers and communicate to those procs that
   hold the new numbering for our outlier points. We need to ask for the
   new numbering from these procs. ---*/

  for (int ii=0; ii < size; ii++) {
    nElem_Send[ii] = 0;
    nElem_Recv[ii] = 0;
    nElem_Flag[ii]= -1;
  }
  nElem_Send[size] = 0; nElem_Recv[size] = 0;

  for (int ii = 0; ii < (int)outliers.size(); ii++) {

    Global_Index = outliers[ii];

    /*--- Search for the processor that owns this point ---*/

    iProcessor = Global_Index/npoint_procs[0];
    if (iProcessor >= (unsigned long)size)
      iProcessor = (unsigned long)size-1;
    if (Global_Index >= nPoint_Linear_Nodes[iProcessor])
      while(Global_Index >= nPoint_Linear_Nodes[iProcessor+1]) iProcessor++;
    else
      while(Global_Index <  nPoint_Linear_Nodes[iProcessor])   iProcessor--;

    /*--- If we have not visited this element yet, increment our
     number of elements that must be sent to a particular proc. ---*/

    if ((nElem_Flag[iProcessor] != ii)) {
      nElem_Flag[iProcessor] = ii;
      nElem_Send[iProcessor+1]++;
    }

  }

  /*--- Communicate the number of cells to be sent/recv'd amongst
   all processors. After this communication, each proc knows how
   many cells it will receive from each other processor. ---*/

#ifdef HAVE_MPI
  SU2_MPI::Alltoall(&(nElem_Send[1]), 1, MPI_INT,
                    &(nElem_Recv[1]), 1, MPI_INT, MPI_COMM_WORLD);
#else
  nElem_Recv[1] = nElem_Send[1];
#endif

  /*--- Prepare to send connectivities. First check how many
   messages we will be sending and receiving. Here we also put
   the counters into cumulative storage format to make the
   communications simpler. ---*/

  nSends = 0; nRecvs = 0;
  for (int ii=0; ii < size; ii++) nElem_Flag[ii] = -1;

  for (int ii = 0; ii < size; ii++) {
    if ((ii != rank) && (nElem_Send[ii+1] > 0)) nSends++;
    if ((ii != rank) && (nElem_Recv[ii+1] > 0)) nRecvs++;

    nElem_Send[ii+1] += nElem_Send[ii];
    nElem_Recv[ii+1] += nElem_Recv[ii];
  }

  delete [] idSend;
  idSend = new unsigned long[nElem_Send[size]];
  for (int ii = 0; ii < nElem_Send[size]; ii++)
    idSend[ii] = 0;

  /*--- Reset our index variable for reuse. ---*/

  for (int ii=0; ii < size; ii++) idIndex[ii] = nElem_Send[ii];

  /*--- Loop over the outliers again and load up the global IDs. ---*/

  for (int ii = 0; ii < (int)outliers.size(); ii++) {

    Global_Index = outliers[ii];

    /*--- Search for the processor that owns this point ---*/

    iProcessor = Global_Index/npoint_procs[0];
    if (iProcessor >= (unsigned long)size)
      iProcessor = (unsigned long)size-1;
    if (Global_Index >= nPoint_Linear_Nodes[iProcessor])
      while(Global_Index >= nPoint_Linear_Nodes[iProcessor+1]) iProcessor++;
    else
      while(Global_Index <  nPoint_Linear_Nodes[iProcessor])   iProcessor--;

    /*--- If we have not visited this element yet, increment our
     number of elements that must be sent to a particular proc. ---*/

    if ((nElem_Flag[iProcessor] != ii)) {

      nElem_Flag[iProcessor] = ii;
      unsigned long nn = idIndex[iProcessor];

      /*--- Load the global ID values. ---*/

      idSend[nn] = Global_Index; nn++;

      /*--- Increment the index by the message length ---*/

      idIndex[iProcessor]++;

    }
  }

  /*--- Allocate the memory that we need for receiving the
   values and then cue up the non-blocking receives. Note that
   we do not include our own rank in the communications. We will
   directly copy our own data later. ---*/

  delete [] idRecv;
  idRecv = new unsigned long[nElem_Recv[size]];
  for (int ii = 0; ii < nElem_Recv[size]; ii++)
    idRecv[ii] = 0;

#ifdef HAVE_MPI
  /*--- We need double the number of messages to send both the conn.
   and the flags for the halo cells. ---*/

  send_req = new SU2_MPI::Request[nSends];
  recv_req = new SU2_MPI::Request[nRecvs];

  /*--- Launch the non-blocking recv's for the connectivity. ---*/

  iMessage = 0;
  for (int ii=0; ii<size; ii++) {
    if ((ii != rank) && (nElem_Recv[ii+1] > nElem_Recv[ii])) {
      int ll     = nElem_Recv[ii];
      int kk     = nElem_Recv[ii+1] - nElem_Recv[ii];
      int count  = kk;
      int source = ii;
      int tag    = ii + 1;
      SU2_MPI::Irecv(&(idRecv[ll]), count, MPI_UNSIGNED_LONG, source, tag,
                     MPI_COMM_WORLD, &(recv_req[iMessage]));
      iMessage++;
    }
  }

  /*--- Launch the non-blocking sends of the connectivity. ---*/

  iMessage = 0;
  for (int ii=0; ii<size; ii++) {
    if ((ii != rank) && (nElem_Send[ii+1] > nElem_Send[ii])) {
      int ll = nElem_Send[ii];
      int kk = nElem_Send[ii+1] - nElem_Send[ii];
      int count  = kk;
      int dest = ii;
      int tag    = rank + 1;
      SU2_MPI::Isend(&(idSend[ll]), count, MPI_UNSIGNED_LONG, dest, tag,
                     MPI_COMM_WORLD, &(send_req[iMessage]));
      iMessage++;
    }
  }
#endif

  /*--- Copy my own rank's data into the recv buffer directly. ---*/

  mm = nElem_Recv[rank];
  ll = nElem_Send[rank];
  kk = nElem_Send[rank+1];

  for (int nn=ll; nn<kk; nn++, mm++) idRecv[mm] = idSend[nn];

  /*--- Wait for the non-blocking sends and recvs to complete. ---*/

#ifdef HAVE_MPI
  number = nSends;
  for (int ii = 0; ii < number; ii++)
    SU2_MPI::Waitany(number, send_req, &ind, &status);

  number = nRecvs;
  for (int ii = 0; ii < number; ii++)
    SU2_MPI::Waitany(number, recv_req, &ind, &status);

  delete [] send_req;
  delete [] recv_req;
#endif

  /*--- The procs holding the outlier grid nodes now have the global IDs
   that they need to have their renumbering shared. ---*/

  for (int ii = 0; ii < nElem_Recv[size]; ii++) {
    for (iPoint = 0; iPoint < nSurf_Poin_Par; iPoint++) {
      if (idRecv[ii] == globalP[iPoint]) {
        idRecv[ii] = renumbP[iPoint];
      }
    }
  }

  /*--- Now simply reverse the last communication to give the renumbered IDs
   back to the owner of the outlier points. Note everything is flipped. ---*/

#ifdef HAVE_MPI
  /*--- We need double the number of messages to send both the conn.
   and the flags for the halo cells. ---*/

  send_req = new SU2_MPI::Request[nRecvs];
  recv_req = new SU2_MPI::Request[nSends];

  /*--- Launch the non-blocking sends of the connectivity. ---*/

  iMessage = 0;
  for (int ii=0; ii<size; ii++) {
    if ((ii != rank) && (nElem_Send[ii+1] > nElem_Send[ii])) {
      int ll = nElem_Send[ii];
      int kk = nElem_Send[ii+1] - nElem_Send[ii];
      int count  = kk;
      int dest = ii;
      int tag    = ii + 1;
      SU2_MPI::Irecv(&(idSend[ll]), count, MPI_UNSIGNED_LONG, dest, tag,
                     MPI_COMM_WORLD, &(recv_req[iMessage]));
      iMessage++;
    }
  }

  /*--- Launch the non-blocking recv's for the connectivity. ---*/

  iMessage = 0;
  for (int ii=0; ii<size; ii++) {
    if ((ii != rank) && (nElem_Recv[ii+1] > nElem_Recv[ii])) {
      int ll     = nElem_Recv[ii];
      int kk     = nElem_Recv[ii+1] - nElem_Recv[ii];
      int count  = kk;
      int source = ii;
      int tag    = rank + 1;
      SU2_MPI::Isend(&(idRecv[ll]), count, MPI_UNSIGNED_LONG, source, tag,
                     MPI_COMM_WORLD, &(send_req[iMessage]));
      iMessage++;
    }
  }
#endif

  /*--- Copy my own rank's data into the recv buffer directly. ---*/

  mm = nElem_Send[rank];
  ll = nElem_Recv[rank];
  kk = nElem_Recv[rank+1];

  for (int nn=ll; nn<kk; nn++, mm++) idSend[mm] = idRecv[nn];

  /*--- Wait for the non-blocking sends and recvs to complete. ---*/

#ifdef HAVE_MPI
  number = nRecvs;
  for (int ii = 0; ii < number; ii++)
    SU2_MPI::Waitany(number, send_req, &ind, &status);

  number = nSends;
  for (int ii = 0; ii < number; ii++)
    SU2_MPI::Waitany(number, recv_req, &ind, &status);

  delete [] send_req;
  delete [] recv_req;
#endif

  /*--- Add the renumbering for the outliers to the map from before carrying
   the global -> renumber transformation. Note that by construction,
   nElem_Send[ii] == outliers.size(). We also add in the 1 for viz. here. ---*/

  for (int ii = 0; ii < nElem_Send[size]; ii++) {
    Global2Renumber[outliers[ii]] = idSend[ii] + 1;
  }

  /*--- We can now overwrite the local connectivity for our surface elems
   using our completed map with the new global renumbering. Whew!! Note
   the -1 when accessing the conn from the map. ---*/

  for (iElem = 0; iElem < nParallel_Line; iElem++) {
    iNode = (int)iElem*N_POINTS_LINE;
    Conn_BoundLine_Par[iNode+0] = (int)Global2Renumber[Conn_BoundLine_Par[iNode+0]-1];
    Conn_BoundLine_Par[iNode+1] = (int)Global2Renumber[Conn_BoundLine_Par[iNode+1]-1];
  }

  for (iElem = 0; iElem < nParallel_BoundTria; iElem++) {
    iNode = (int)iElem*N_POINTS_TRIANGLE;
    Conn_BoundTria_Par[iNode+0] = (int)Global2Renumber[Conn_BoundTria_Par[iNode+0]-1];
    Conn_BoundTria_Par[iNode+1] = (int)Global2Renumber[Conn_BoundTria_Par[iNode+1]-1];
    Conn_BoundTria_Par[iNode+2] = (int)Global2Renumber[Conn_BoundTria_Par[iNode+2]-1];
  }

  for (iElem = 0; iElem < nParallel_BoundQuad; iElem++) {
    iNode = (int)iElem*N_POINTS_QUADRILATERAL;
    Conn_BoundQuad_Par[iNode+0] = (int)Global2Renumber[Conn_BoundQuad_Par[iNode+0]-1];
    Conn_BoundQuad_Par[iNode+1] = (int)Global2Renumber[Conn_BoundQuad_Par[iNode+1]-1];
    Conn_BoundQuad_Par[iNode+2] = (int)Global2Renumber[Conn_BoundQuad_Par[iNode+2]-1];
    Conn_BoundQuad_Par[iNode+3] = (int)Global2Renumber[Conn_BoundQuad_Par[iNode+3]-1];
  }

  /*--- Free temporary memory ---*/

  delete [] idIndex;
  delete [] surfPoint;
  delete [] globalP;
  delete [] renumbP;

  delete [] idSend;
  delete [] idRecv;
  delete [] globalSend;
  delete [] globalRecv;
  delete [] renumbSend;
  delete [] renumbRecv;
  delete [] nElem_Recv;
  delete [] nElem_Send;
  delete [] nElem_Flag;
  delete [] Local_Halo;
  delete [] Buffer_Recv_nAddedPeriodic;
  delete [] Buffer_Send_AddedPeriodic;
  delete [] Buffer_Recv_AddedPeriodic;
  delete [] npoint_procs;
  delete [] starting_node;
  delete [] ending_node;
  delete [] nPoint_Linear_Elems;
  delete [] nPoint_Linear_Nodes;
  delete [] nPoint_Send;
  delete [] nPoint_Recv;

}

void COutputLegacy::WriteRestart_Parallel_ASCII(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned short val_iZone, unsigned short val_iInst) {

  /*--- Local variables ---*/

  unsigned short nZone = geometry->GetnZone(), nInst = config->GetnTimeInstances();
  unsigned short iVar;
  unsigned long iPoint, iExtIter = config->GetInnerIter();
  bool fem       = (config->GetKind_Solver() == FEM_ELASTICITY);
  bool disc_adj_fem = (config->GetKind_Solver() == DISC_ADJ_FEM);
  bool adjoint   = (config->GetContinuous_Adjoint() ||
                    config->GetDiscrete_Adjoint());
  ofstream restart_file;
  string filename;

  int iProcessor;

  /*--- Retrieve filename from config ---*/

//  if ((config->GetContinuous_Adjoint()) || (config->GetDiscrete_Adjoint())) {
//    filename = config->GetRestart_AdjFileName();
//    filename = config->GetObjFunc_Extension(filename);
//  } else if (fem) {
//    filename = config->GetRestart_FEMFileName();
//  } else if (disc_adj_fem){
//    filename = config->GetRestart_AdjFEMFileName();
//  } else {
//    filename = config->GetRestart_FileName();
//  }

  /*--- Append the zone number if multizone problems ---*/
  if (nZone > 1)
    filename= config->GetMultizone_FileName(filename, val_iZone, ".dat");

  /*--- Append the zone number if multiple instance problems ---*/
  if (nInst > 1)
    filename= config->GetMultiInstance_FileName(filename, val_iInst, ".dat");

  /*--- Unsteady problems require an iteration number to be appended. ---*/
  if (config->GetTime_Marching() == HARMONIC_BALANCE) {
    filename = config->GetUnsteady_FileName(filename, SU2_TYPE::Int(val_iInst), ".dat");
  } else if (config->GetTime_Domain()) {
    filename = config->GetUnsteady_FileName(filename, SU2_TYPE::Int(iExtIter), ".dat");
  } else if ((fem || disc_adj_fem) && (config->GetWrt_Dynamic())) {
    filename = config->GetUnsteady_FileName(filename, SU2_TYPE::Int(iExtIter), ".dat");
  }

  /*--- Only the master node writes the header. ---*/

  if (rank == MASTER_NODE) {
    restart_file.open(filename.c_str(), ios::out);
    restart_file.precision(15);
    restart_file << "\"PointID\"";
    for (iVar = 0; iVar < Variable_Names.size()-1; iVar++)
      restart_file << "\t\"" << Variable_Names[iVar] << "\"";
    restart_file << "\t\"" << Variable_Names[Variable_Names.size()-1] << "\"" << endl;
    restart_file.close();
  }

#ifdef HAVE_MPI
  SU2_MPI::Barrier(MPI_COMM_WORLD);
#endif

  /*--- All processors open the file. ---*/

  restart_file.open(filename.c_str(), ios::out | ios::app);
  restart_file.precision(15);

  /*--- Write the restart file in parallel, processor by processor. ---*/

  unsigned long myPoint = 0, offset = 0, Global_Index;
  for (iProcessor = 0; iProcessor < size; iProcessor++) {
    if (rank == iProcessor) {
      for (iPoint = 0; iPoint < nParallel_Poin; iPoint++) {

        /*--- Global Index of the current point. (note outer loop over procs) ---*/

        Global_Index = iPoint + offset;

        /*--- Only write original domain points, i.e., exclude any periodic
         or halo nodes, even if they are output in the viz. files. ---*/

        if (Global_Index < nPoint_Restart) {

          /*--- Write global index. (note outer loop over procs) ---*/

          restart_file << Global_Index << "\t";
          myPoint++;

          /*--- Loop over the variables and write the values to file ---*/

          for (iVar = 0; iVar < nVar_Par; iVar++) {
            restart_file << scientific << Parallel_Data[iVar][iPoint] << "\t";
          }
          restart_file << "\n";
        }
      }
    }
    /*--- Flush the file and wait for all processors to arrive. ---*/
    restart_file.flush();
#ifdef HAVE_MPI
    SU2_MPI::Allreduce(&myPoint, &offset, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Barrier(MPI_COMM_WORLD);
#endif

  }

  /*--- Write the metadata (master rank alone) ----*/

  if (rank == MASTER_NODE) {
//    if (dual_time)
      restart_file <<"EXT_ITER= " << config->GetInnerIter() + 1 << endl;
//    else
//      restart_file <<"EXT_ITER= " << config->GetInnerIter() + config->GetInnerIter_OffSet() + 1 << endl;
    restart_file <<"AOA= " << config->GetAoA() - config->GetAoA_Offset() << endl;
    restart_file <<"SIDESLIP_ANGLE= " << config->GetAoS() - config->GetAoS_Offset() << endl;
    restart_file <<"INITIAL_BCTHRUST= " << config->GetInitial_BCThrust() << endl;
    restart_file <<"DCD_DCL_VALUE= " << config->GetdCD_dCL() << endl;
    restart_file <<"DCMX_DCL_VALUE= " << config->GetdCMx_dCL() << endl;
    restart_file <<"DCMY_DCL_VALUE= " << config->GetdCMy_dCL() << endl;
    restart_file <<"DCMZ_DCL_VALUE= " << config->GetdCMz_dCL() << endl;

    if (((config->GetKind_Solver() == DISC_ADJ_EULER)         ||
        (config->GetKind_Solver() == DISC_ADJ_NAVIER_STOKES) ||
        (config->GetKind_Solver() == DISC_ADJ_RANS) ||
        (config->GetKind_Solver() == DISC_ADJ_INC_EULER)         ||
        (config->GetKind_Solver() == DISC_ADJ_INC_NAVIER_STOKES) ||
        (config->GetKind_Solver() == DISC_ADJ_INC_RANS)) && adjoint) {
      restart_file << "SENS_AOA=" << solver[ADJFLOW_SOL]->GetTotal_Sens_AoA() * PI_NUMBER / 180.0 << endl;
    }
  }

  /*--- All processors close the file. ---*/

  restart_file.close();

}

void COutputLegacy::WriteRestart_Parallel_Binary(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned short val_iZone, unsigned short val_iInst) {

  /*--- Local variables ---*/

  unsigned short iVar, nZone = geometry->GetnZone(), nInst = config->GetnTimeInstances();
  unsigned long iPoint, iExtIter = config->GetInnerIter();
  bool fem       = (config->GetKind_Solver() == FEM_ELASTICITY);
  bool adjoint   = (config->GetContinuous_Adjoint() ||
                    config->GetDiscrete_Adjoint());
  bool dual_time = ((config->GetTime_Marching() == DT_STEPPING_1ST) ||
                    (config->GetTime_Marching() == DT_STEPPING_2ND));
  bool wrt_perf  = config->GetWrt_Performance();
  ofstream restart_file;
  string filename;
  char str_buf[CGNS_STRING_SIZE], fname[100];
  su2double file_size = 0.0, StartTime, StopTime, UsedTime, Bandwidth;

  /*--- Retrieve filename from config ---*/

//  if ((config->GetContinuous_Adjoint()) || (config->GetDiscrete_Adjoint())) {
//    filename = config->GetRestart_AdjFileName();
//    filename = config->GetObjFunc_Extension(filename);
//  } else if (fem) {
//    filename = config->GetRestart_FEMFileName();
//  } else {
//    filename = config->GetRestart_FileName();
//  }

  /*--- Append the zone number if multizone problems ---*/
  if (nZone > 1)
    filename= config->GetMultizone_FileName(filename, val_iZone, ".dat");

  /*--- Append the zone number if multiple instance problems ---*/
  if (nInst > 1)
    filename= config->GetMultiInstance_FileName(filename, val_iInst, ".dat");

  /*--- Unsteady problems require an iteration number to be appended. ---*/
  if (config->GetTime_Marching() == HARMONIC_BALANCE) {
    filename = config->GetUnsteady_FileName(filename, SU2_TYPE::Int(val_iInst), ".dat");
  } else if (config->GetTime_Domain()) {
    filename = config->GetUnsteady_FileName(filename, SU2_TYPE::Int(iExtIter), ".dat");
  } else if ((fem) && (config->GetWrt_Dynamic())) {
    filename = config->GetUnsteady_FileName(filename, SU2_TYPE::Int(iExtIter), ".dat");
  }

  strcpy(fname, filename.c_str());

  /*--- Prepare the first ints containing the counts. The first is a
   magic number that we can use to check for binary files (it is the hex
   representation for "SU2"). The second two values are number of variables
   and number of points (DoFs). The last two values are for metadata:
   one int for ExtIter and 8 su2doubles. ---*/

  int var_buf_size = 5;
  int var_buf[5] = {535532, nVar_Par, (int)nGlobalPoint_Sort, 1, 8};

  /*--- Prepare the 1D data buffer on this rank. ---*/

  passivedouble *buf = new passivedouble[nParallel_Poin*nVar_Par];

  /*--- For now, create a temp 1D buffer to load up the data for writing.
   This will be replaced with a derived data type most likely. ---*/

  for (iPoint = 0; iPoint < nParallel_Poin; iPoint++)
    for (iVar = 0; iVar < nVar_Par; iVar++)
      buf[iPoint*nVar_Par+iVar] = SU2_TYPE::GetValue(Parallel_Data[iVar][iPoint]);

  /*--- Prepare metadata. ---*/

  int Restart_ExtIter;
  if (dual_time)
    Restart_ExtIter= (int)config->GetInnerIter() + 1;
  else
    Restart_ExtIter = (int)config->GetInnerIter() + (int)config->GetExtIter_OffSet() + 1;

  passivedouble Restart_Metadata[8] = {
    SU2_TYPE::GetValue(config->GetAoA() - config->GetAoA_Offset()),
    SU2_TYPE::GetValue(config->GetAoS() - config->GetAoS_Offset()),
    SU2_TYPE::GetValue(config->GetInitial_BCThrust()),
    SU2_TYPE::GetValue(config->GetdCD_dCL()),
    SU2_TYPE::GetValue(config->GetdCMx_dCL()),
    SU2_TYPE::GetValue(config->GetdCMy_dCL()),
    SU2_TYPE::GetValue(config->GetdCMz_dCL()),
    0.0
  };

  if (((config->GetKind_Solver() == DISC_ADJ_EULER)         ||
      (config->GetKind_Solver() == DISC_ADJ_NAVIER_STOKES) ||
      (config->GetKind_Solver() == DISC_ADJ_RANS) ||
      (config->GetKind_Solver() == DISC_ADJ_INC_EULER)         ||
      (config->GetKind_Solver() == DISC_ADJ_INC_NAVIER_STOKES) ||
      (config->GetKind_Solver() == DISC_ADJ_INC_RANS)) && adjoint) {
    Restart_Metadata[4] = SU2_TYPE::GetValue(solver[ADJFLOW_SOL]->GetTotal_Sens_AoA() * PI_NUMBER / 180.0);
  }

  /*--- Set a timer for the binary file writing. ---*/

  StartTime = SU2_MPI::Wtime();

#ifndef HAVE_MPI

  FILE* fhw;
  fhw = fopen(fname, "wb");

  /*--- Error check for opening the file. ---*/

  if (!fhw) {
    SU2_MPI::Error(string("Unable to open SU2 restart file ") + string(fname), CURRENT_FUNCTION);
  }

  /*--- First, write the number of variables and points. ---*/

  fwrite(var_buf, var_buf_size, sizeof(int), fhw);
  file_size += (su2double)var_buf_size*sizeof(int);

  /*--- Write the variable names to the file. Note that we are adopting a
   fixed length of 33 for the string length to match with CGNS. This is
   needed for when we read the strings later. ---*/

  for (iVar = 0; iVar < nVar_Par; iVar++) {
    strncpy(str_buf, Variable_Names[iVar].c_str(), CGNS_STRING_SIZE);
    fwrite(str_buf, CGNS_STRING_SIZE, sizeof(char), fhw);
    file_size += (su2double)CGNS_STRING_SIZE*sizeof(char);
  }

  /*--- Call to write the entire restart file data in binary in one shot. ---*/

  fwrite(buf, nVar_Par*nParallel_Poin, sizeof(passivedouble), fhw);
  file_size += (su2double)nVar_Par*nParallel_Poin*sizeof(passivedouble);

  /*--- Write the external iteration. ---*/

  fwrite(&Restart_ExtIter, 1, sizeof(int), fhw);
  file_size += (su2double)sizeof(int);

  /*--- Write the metadata. ---*/

  fwrite(Restart_Metadata, 8, sizeof(passivedouble), fhw);
  file_size += (su2double)8*sizeof(passivedouble);

  /*--- Close the file. ---*/

  fclose(fhw);

#else

  /*--- Parallel binary output using MPI I/O. ---*/

  MPI_File fhw;
  SU2_MPI::Status status;
  MPI_Datatype etype, filetype;
  MPI_Offset disp;
  int ierr;

  /*--- We're writing only su2doubles in the data portion of the file. ---*/

  etype = MPI_DOUBLE;

  /*--- Define a derived datatype for this ranks contiguous chunk of data
   that will be placed in the restart (1D array size = num points * num vars). ---*/

  MPI_Type_contiguous(nVar_Par*nParallel_Poin, MPI_DOUBLE, &filetype);
  MPI_Type_commit(&filetype);

  /*--- All ranks open the file using MPI. Here, we try to open the file with
   exclusive so that an error is generated if the file exists. We always want
   to write a fresh restart file, so we delete any existing files and create
   a new one. ---*/

  ierr = MPI_File_open(MPI_COMM_WORLD, fname,
                       MPI_MODE_CREATE|MPI_MODE_EXCL|MPI_MODE_WRONLY,
                       MPI_INFO_NULL, &fhw);
  if (ierr != MPI_SUCCESS)  {
    MPI_File_close(&fhw);
    if (rank == 0)
      MPI_File_delete(fname, MPI_INFO_NULL);
    ierr = MPI_File_open(MPI_COMM_WORLD, fname,
                         MPI_MODE_CREATE|MPI_MODE_EXCL|MPI_MODE_WRONLY,
                         MPI_INFO_NULL, &fhw);
  }

  /*--- Error check opening the file. ---*/

  if (ierr) {
    SU2_MPI::Error(string("Unable to open SU2 restart file ") + string(fname), CURRENT_FUNCTION);
  }

  /*--- First, write the number of variables and points (i.e., cols and rows),
   which we will need in order to read the file later. Also, write the
   variable string names here. Only the master rank writes the header. ---*/

  if (rank == MASTER_NODE) {
    MPI_File_write(fhw, var_buf, var_buf_size, MPI_INT, MPI_STATUS_IGNORE);
    file_size += (su2double)var_buf_size*sizeof(int);

    /*--- Write the variable names to the file. Note that we are adopting a
     fixed length of 33 for the string length to match with CGNS. This is
     needed for when we read the strings later. ---*/

    for (iVar = 0; iVar < nVar_Par; iVar++) {
      disp = var_buf_size*sizeof(int) + iVar*CGNS_STRING_SIZE*sizeof(char);
      strncpy(str_buf, Variable_Names[iVar].c_str(), CGNS_STRING_SIZE);
      MPI_File_write_at(fhw, disp, str_buf, CGNS_STRING_SIZE, MPI_CHAR, MPI_STATUS_IGNORE);
      file_size += (su2double)CGNS_STRING_SIZE*sizeof(char);
    }
  }

  /*--- Compute the offset for this rank's linear partition of the data in bytes.
   After the calculations above, we have the partition sizes store in nPoint_Linear
   in cumulative storage format. ---*/

  disp = (var_buf_size*sizeof(int) + nVar_Par*CGNS_STRING_SIZE*sizeof(char) +
          nVar_Par*nPointCumulative[rank]*sizeof(passivedouble));

  /*--- Set the view for the MPI file write, i.e., describe the location in
   the file that this rank "sees" for writing its piece of the restart file. ---*/

  MPI_File_set_view(fhw, disp, etype, filetype, (char*)"native", MPI_INFO_NULL);

  /*--- Collective call for all ranks to write to their view simultaneously. ---*/

  MPI_File_write_all(fhw, buf, nVar_Par*nParallel_Poin, MPI_DOUBLE, &status);
  file_size += (su2double)nVar_Par*nParallel_Poin*sizeof(passivedouble);

  /*--- Free the derived datatype. ---*/

  MPI_Type_free(&filetype);

  /*--- Reset the file view before writing the metadata. ---*/

  MPI_File_set_view(fhw, 0, MPI_BYTE, MPI_BYTE, (char*)"native", MPI_INFO_NULL);

  /*--- Finally, the master rank writes the metadata. ---*/

  if (rank == MASTER_NODE) {

    /*--- External iteration. ---*/

    disp = (var_buf_size*sizeof(int) + nVar_Par*CGNS_STRING_SIZE*sizeof(char) +
            nVar_Par*nGlobalPoint_Sort*sizeof(passivedouble));
    MPI_File_write_at(fhw, disp, &Restart_ExtIter, 1, MPI_INT, MPI_STATUS_IGNORE);
    file_size += (su2double)sizeof(int);

    /*--- Additional doubles for AoA, AoS, etc. ---*/

    disp = (var_buf_size*sizeof(int) + nVar_Par*CGNS_STRING_SIZE*sizeof(char) +
            nVar_Par*nGlobalPoint_Sort*sizeof(passivedouble) + 1*sizeof(int));
    MPI_File_write_at(fhw, disp, Restart_Metadata, 8, MPI_DOUBLE, MPI_STATUS_IGNORE);
    file_size += (su2double)8*sizeof(passivedouble);

  }

  /*--- All ranks close the file after writing. ---*/

  MPI_File_close(&fhw);

#endif

  /*--- Compute and store the write time. ---*/

  StopTime = SU2_MPI::Wtime();

  UsedTime = StopTime-StartTime;

  /*--- Communicate the total file size for the restart ---*/

#ifdef HAVE_MPI
  su2double my_file_size = file_size;
  SU2_MPI::Allreduce(&my_file_size, &file_size, 1,
                     MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif

  /*--- Compute and store the bandwidth ---*/

  Bandwidth = file_size/(1.0e6)/UsedTime;
  config->SetRestart_Bandwidth_Agg(config->GetRestart_Bandwidth_Agg()+Bandwidth);

  if ((rank == MASTER_NODE) && (wrt_perf)) {
    cout << "Wrote " << file_size/1.0e6 << " MB to disk in ";
    cout << UsedTime << " s. (" << Bandwidth << " MB/s)." << endl;
  }

  /*--- Free temporary data buffer for writing the binary file. ---*/

  delete [] buf;

}

void COutputLegacy::WriteCSV_Slice(CConfig *config, CGeometry *geometry,
                             CSolver *FlowSolver, unsigned long iExtIter,
                             unsigned short val_iZone, unsigned short val_direction) {

  /*--- This routine is for exporting slices of field data for 2D cartesian
   grids. It assumes that the grid points lie on lines of constant x-
   or y-coordinates. It is a simple way to export slices or profiles on
   these meshes for use in verification and validation work. It will
   eventually be replaced by general routines for probing/slicing.  ---*/

  int DIRECTION = (int)val_direction;

  su2double coordMin, coordMax;
  coordMin = config->GetStations_Bounds(0);
  coordMax = config->GetStations_Bounds(1);

  unsigned short iVar;
  unsigned long iPoint, iVertex, Global_Index;
  char cstr[200];

  int rank = MASTER_NODE, iProcessor, nProcessor = SINGLE_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);
#endif

  bool isPeriodic;

  int *Local_Halo = new int[geometry->GetnPoint()];
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
    Local_Halo[iPoint] = !geometry->node[iPoint]->GetDomain();

  for (unsigned short iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) {

      /*--- Checking for less than or equal to the rank, because there may
       be some periodic halo nodes that send info to the same rank. ---*/

      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        isPeriodic = ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() > 0) &&
                      (geometry->vertex[iMarker][iVertex]->GetRotation_Type() % 2 == 1));
        if (isPeriodic) Local_Halo[iPoint] = false;
      }
    }
  }

  /*--- Sum total number of nodes that belong to the domain ---*/

  unsigned long nTotalPoint;
  unsigned long nLocalPoint = 0;
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
    if (Local_Halo[iPoint] == false)
      nLocalPoint++;

#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&nLocalPoint, &nTotalPoint, 1,
                     MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#else
  nTotalPoint = nLocalPoint;
#endif

  unsigned long *npoint_procs  = new unsigned long[nProcessor];
  unsigned long *nPoint_Linear = new unsigned long[nProcessor+1];

  unsigned long total_pt_accounted = 0;
  for (int ii = 0; ii < nProcessor; ii++) {
    npoint_procs[ii] = nTotalPoint/nProcessor;
    total_pt_accounted = total_pt_accounted + npoint_procs[ii];
  }

  /*--- Get the number of remainder points after the even division. ---*/

  unsigned long rem_points = nTotalPoint-total_pt_accounted;
  for (unsigned long ii = 0; ii < rem_points; ii++) {
    npoint_procs[ii]++;
  }

  /*--- Store the point offsets for each rank. ---*/

  nPoint_Linear[0] = 0;
  for (int ii = 1; ii < nProcessor; ii++) {
    nPoint_Linear[ii] = nPoint_Linear[ii-1] + npoint_procs[ii-1];
  }
  nPoint_Linear[nProcessor] = nTotalPoint;

  unsigned long Buffer_Send_nVertex[1], *Buffer_Recv_nVertex = NULL;
  unsigned long nLocalVertex_Surface = 0;
  unsigned long MaxLocalVertex_Surface = 0;

  /*--- Find the max number of vertices we will send from among all
   partitions and set up buffers. The master node will handle the
   writing of the CSV file after gathering all of the data. ---*/

  nLocalVertex_Surface = 0;
  for (iPoint = 0; iPoint < nParallel_Poin; iPoint++) {

    /*--- Global Index of the current point. (note outer loop over procs) ---*/

    Global_Index = iPoint + nPoint_Linear[rank];

    /*--- Only write original domain points, i.e., exclude any periodic
     or halo nodes, even if they are output in the viz. files. ---*/

    if (Global_Index < geometry->GetGlobal_nPointDomain()) {
      if ((Parallel_Data[DIRECTION][iPoint] > coordMin) &&
          (Parallel_Data[DIRECTION][iPoint] < coordMax)) {
        nLocalVertex_Surface++;
      }
    }
  }

  /*--- Communicate the number of local vertices on each partition
   to the master node ---*/

  Buffer_Send_nVertex[0] = nLocalVertex_Surface;
  if (rank == MASTER_NODE) Buffer_Recv_nVertex = new unsigned long [nProcessor];

#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&nLocalVertex_Surface, &MaxLocalVertex_Surface, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
  SU2_MPI::Gather(&Buffer_Send_nVertex, 1, MPI_UNSIGNED_LONG, Buffer_Recv_nVertex, 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
#else
  MaxLocalVertex_Surface = nLocalVertex_Surface;
  Buffer_Recv_nVertex[0] = Buffer_Send_nVertex[0];
#endif

  /*--- Send and Recv buffers ---*/

  su2double *Buffer_Send_Data = new su2double [MaxLocalVertex_Surface*nVar_Par];
  su2double *Buffer_Recv_Data = NULL;

  unsigned long *Buffer_Send_GlobalIndex = new unsigned long [MaxLocalVertex_Surface];
  unsigned long *Buffer_Recv_GlobalIndex = NULL;

  /*--- Prepare the receive buffers on the master node only. ---*/

  if (rank == MASTER_NODE) {
    Buffer_Recv_Data = new su2double [nProcessor*MaxLocalVertex_Surface*nVar_Par];
    Buffer_Recv_GlobalIndex  = new unsigned long [nProcessor*MaxLocalVertex_Surface];
  }

  /*--- Loop over all vertices in this partition and load the
   data of the specified type into the buffer to be sent to
   the master node. ---*/

  nLocalVertex_Surface = 0;
  for (iPoint = 0; iPoint < nParallel_Poin; iPoint++) {

    /*--- Global Index of the current point. (note outer loop over procs) ---*/

    Global_Index = iPoint + nPoint_Linear[rank];

    /*--- Only write original domain points, i.e., exclude any periodic
     or halo nodes, even if they are output in the viz. files. ---*/

    if (Global_Index < geometry->GetGlobal_nPointDomain()) {
      if ((Parallel_Data[DIRECTION][iPoint] > coordMin) &&
          (Parallel_Data[DIRECTION][iPoint] < coordMax)) {
        Buffer_Send_GlobalIndex[nLocalVertex_Surface] = Global_Index;
        for (iVar = 0; iVar < nVar_Par; iVar++) {
          Buffer_Send_Data[nLocalVertex_Surface*nVar_Par+iVar] = Parallel_Data[iVar][iPoint];
        }
        nLocalVertex_Surface++;
      }
    }
  }

  /*--- Send the information to the master node ---*/

#ifdef HAVE_MPI
  SU2_MPI::Gather(Buffer_Send_Data, MaxLocalVertex_Surface*nVar_Par, MPI_DOUBLE, Buffer_Recv_Data, MaxLocalVertex_Surface*nVar_Par, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  SU2_MPI::Gather(Buffer_Send_GlobalIndex, MaxLocalVertex_Surface, MPI_UNSIGNED_LONG, Buffer_Recv_GlobalIndex, MaxLocalVertex_Surface, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
#else
  for (iVertex = 0; iVertex < Buffer_Recv_nVertex[0]; iVertex++) {
    Buffer_Recv_GlobalIndex[iVertex] = Buffer_Send_GlobalIndex[iVertex];
    for (iVar = 0; iVar < nVar_Par; iVar++) {
      Buffer_Recv_Data[iVertex*nVar_Par+iVar] = Buffer_Send_Data[iVertex*nVar_Par+iVar];
    }
  }
#endif

  /*--- The master node unpacks the data and writes the surface CSV file ---*/

  if (rank == MASTER_NODE) {

    /*--- Write file name with extension if unsteady ---*/
    char buffer[50];
    string filename = "slice";
    if (DIRECTION == 0) {
      SPRINTF (buffer, "_vert.csv");
    } else if (DIRECTION == 1) {
      SPRINTF (buffer, "_hori.csv");
    }
    ofstream SurfFlow_file;

    /*--- Write file name with extension if unsteady ---*/
    strcpy (cstr, filename.c_str());
    strcat (cstr, buffer);
    SurfFlow_file.precision(15);
    SurfFlow_file.open(cstr, ios::out);

    /*--- Global index is first, then the rest of the data. We already have the names. ---*/

    SurfFlow_file << "\"Global_Index\",";
    for (iVar = 0; iVar < Variable_Names.size()-1; iVar++) {
      SurfFlow_file << "\"" << Variable_Names[iVar] << "\",";
    }
    SurfFlow_file << "\"" << Variable_Names[Variable_Names.size()-1] << "\"" << endl;

    /*--- Loop through all of the collected data and write each node's values ---*/

    unsigned long Total_Index;
    for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
      for (iVertex = 0; iVertex < Buffer_Recv_nVertex[iProcessor]; iVertex++) {

        /*--- Current index position and global index ---*/

        Total_Index  = iProcessor*MaxLocalVertex_Surface+iVertex;
        Global_Index = Buffer_Recv_GlobalIndex[Total_Index];

        /*--- Write the the data ---*/

        SurfFlow_file << scientific << Global_Index;
        Total_Index  = iProcessor*MaxLocalVertex_Surface*nVar_Par+iVertex*nVar_Par;
        for (iVar = 0; iVar < nVar_Par; iVar++) {
          SurfFlow_file << scientific << ", " << Buffer_Recv_Data[Total_Index+iVar];
        }
        SurfFlow_file << endl;

      }
    }

    /*--- Close the CSV file ---*/

    SurfFlow_file.close();

    /*--- Release the recv buffers on the master node ---*/

    delete [] Buffer_Recv_Data;
    delete [] Buffer_Recv_GlobalIndex;
    delete [] Buffer_Recv_nVertex;

  }

  /*--- Release the memory for the remaining buffers and exit ---*/

  delete [] Buffer_Send_Data;
  delete [] Buffer_Send_GlobalIndex;

}

void COutputLegacy::DeallocateConnectivity_Parallel(CConfig *config, CGeometry *geometry, bool surf_sol) {

  /*--- Deallocate memory for connectivity data on each processor. ---*/

  if (surf_sol) {
    if (Conn_BoundLine_Par != NULL) delete [] Conn_BoundLine_Par;
    if (Conn_BoundTria_Par != NULL) delete [] Conn_BoundTria_Par;
    if (Conn_BoundQuad_Par != NULL) delete [] Conn_BoundQuad_Par;
  }
  else {
    if (Conn_Tria_Par != NULL) delete [] Conn_Tria_Par;
    if (Conn_Quad_Par != NULL) delete [] Conn_Quad_Par;
    if (Conn_Tetr_Par != NULL) delete [] Conn_Tetr_Par;
    if (Conn_Hexa_Par != NULL) delete [] Conn_Hexa_Par;
    if (Conn_Pris_Par != NULL) delete [] Conn_Pris_Par;
    if (Conn_Pyra_Par != NULL) delete [] Conn_Pyra_Par;
  }

}

void COutputLegacy::DeallocateData_Parallel(CConfig *config, CGeometry *geometry) {

  /*--- Deallocate memory for solution data ---*/

  for (unsigned short iVar = 0; iVar < nVar_Par; iVar++) {
    if (Parallel_Data[iVar] != NULL) delete [] Parallel_Data[iVar];
  }
  if (Parallel_Data != NULL) delete [] Parallel_Data;

  /*--- Deallocate the structures holding the linear partitioning ---*/

  if (Local_Halo_Sort != NULL) delete [] Local_Halo_Sort;

  if (beg_node != NULL) delete [] beg_node;
  if (end_node != NULL) delete [] end_node;

  if (nPointLinear     != NULL) delete [] nPointLinear;
  if (nPointCumulative != NULL) delete [] nPointCumulative;

}

void COutputLegacy::DeallocateSurfaceData_Parallel(CConfig *config, CGeometry *geometry) {

  /*--- Deallocate memory for surface solution data ---*/

  for (unsigned short iVar = 0; iVar < nVar_Par; iVar++) {
    if (Parallel_Surf_Data[iVar] != NULL) delete [] Parallel_Surf_Data[iVar];
  }
  if (Parallel_Surf_Data != NULL) delete [] Parallel_Surf_Data;

}

void COutputLegacy::SpecialOutput_AnalyzeSurface(CSolver *solver, CGeometry *geometry, CConfig *config, bool output) {

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

  bool compressible   = config->GetKind_Regime() == COMPRESSIBLE;
  bool incompressible = config->GetKind_Regime() == INCOMPRESSIBLE;
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

        if (geometry->node[iPoint]->GetDomain()) {

          geometry->vertex[iMarker][iVertex]->GetNormal(Vector);

          if (axisymmetric) {
            if (geometry->node[iPoint]->GetCoord(1) != 0.0)
              AxiFactor = 2.0*PI_NUMBER*geometry->node[iPoint]->GetCoord(1);
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
            if (config->GetKind_DensityModel() == VARIABLE) {
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
    SU2_MPI::Allreduce(Surface_MassFlow_Local, Surface_MassFlow_Total, nMarker_Analyze, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(Surface_Mach_Local, Surface_Mach_Total, nMarker_Analyze, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(Surface_Temperature_Local, Surface_Temperature_Total, nMarker_Analyze, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(Surface_Density_Local, Surface_Density_Total, nMarker_Analyze, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(Surface_Enthalpy_Local, Surface_Enthalpy_Total, nMarker_Analyze, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(Surface_NormalVelocity_Local, Surface_NormalVelocity_Total, nMarker_Analyze, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(Surface_StreamVelocity2_Local, Surface_StreamVelocity2_Total, nMarker_Analyze, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(Surface_TransvVelocity2_Local, Surface_TransvVelocity2_Total, nMarker_Analyze, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(Surface_Pressure_Local, Surface_Pressure_Total, nMarker_Analyze, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(Surface_TotalTemperature_Local, Surface_TotalTemperature_Total, nMarker_Analyze, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(Surface_TotalPressure_Local, Surface_TotalPressure_Total, nMarker_Analyze, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(Surface_Area_Local, Surface_Area_Total, nMarker_Analyze, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(Surface_MassFlow_Abs_Local, Surface_MassFlow_Abs_Total, nMarker_Analyze, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
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

void COutputLegacy::MergeConnectivity_FEM(CConfig *config, CGeometry *geometry, unsigned short val_iZone) {

  /*--- Flags identifying the types of files to be written. ---*/

  bool Wrt_Vol = config->GetWrt_Vol_Sol();
  bool Wrt_Srf = config->GetWrt_Srf_Sol();

  /*--- Merge connectivity for each type of element (excluding halos). Note
   that we only need to merge the connectivity once, as it does not change
   during computation. Check whether the base file has been written. ---*/

  /*--- Merge volumetric grid. ---*/

  if (Wrt_Vol) {

    if ((rank == MASTER_NODE) && (size != SINGLE_NODE) && (nGlobal_Tria != 0))
      cout <<"Merging volumetric triangle grid connectivity." << endl;
    MergeVolumetricConnectivity_FEM(config, geometry, TRIANGLE    );

    if ((rank == MASTER_NODE) && (size != SINGLE_NODE) && (nGlobal_Quad != 0))
      cout <<"Merging volumetric quadrilateral grid connectivity." << endl;
    MergeVolumetricConnectivity_FEM(config, geometry, QUADRILATERAL   );

    if ((rank == MASTER_NODE) && (size != SINGLE_NODE) && (nGlobal_Tetr != 0))
      cout <<"Merging volumetric tetrahedron grid connectivity." << endl;
    MergeVolumetricConnectivity_FEM(config, geometry, TETRAHEDRON );

    if ((rank == MASTER_NODE) && (size != SINGLE_NODE) && (nGlobal_Hexa != 0))
      cout <<"Merging volumetric hexahedron grid connectivity." << endl;
    MergeVolumetricConnectivity_FEM(config, geometry, HEXAHEDRON  );

    if ((rank == MASTER_NODE) && (size != SINGLE_NODE) && (nGlobal_Pris != 0))
      cout <<"Merging volumetric prism grid connectivity." << endl;
    MergeVolumetricConnectivity_FEM(config, geometry, PRISM       );

    if ((rank == MASTER_NODE) && (size != SINGLE_NODE) && (nGlobal_Pyra != 0))
      cout <<"Merging volumetric pyramid grid connectivity." << endl;
    MergeVolumetricConnectivity_FEM(config, geometry, PYRAMID     );

  }

  /*--- Merge surface grid. ---*/

  if (Wrt_Srf) {

    if ((rank == MASTER_NODE) && (size != SINGLE_NODE) && (nGlobal_Line != 0))
      cout <<"Merging surface line grid connectivity." << endl;
    MergeSurfaceConnectivity_FEM(config, geometry, LINE);

    if ((rank == MASTER_NODE) && (size != SINGLE_NODE) && (nGlobal_BoundTria != 0))
      cout <<"Merging surface triangle grid connectivity." << endl;
    MergeSurfaceConnectivity_FEM(config, geometry, TRIANGLE);

    if ((rank == MASTER_NODE) && (size != SINGLE_NODE) && (nGlobal_BoundQuad != 0))
      cout <<"Merging surface quadrilateral grid connectivity." << endl;
    MergeSurfaceConnectivity_FEM(config, geometry, QUADRILATERAL);

  }

  /*--- Update total number of volume elements after merge. ---*/

  nGlobal_Elem = nGlobal_Tria + nGlobal_Quad + nGlobal_Tetr +
  nGlobal_Hexa + nGlobal_Pyra + nGlobal_Pris;

  /*--- Update total number of surface elements after merge. ---*/

  nSurf_Elem = nGlobal_Line + nGlobal_BoundTria + nGlobal_BoundQuad;

}

void COutputLegacy::MergeCoordinates_FEM(CConfig *config, CGeometry *geometry) {

  /*--- Local variables needed on all processors ---*/

  unsigned short iDim;
  unsigned long iPoint;

  /*--- Create an object of the class CMeshFEM_DG and retrieve the necessary
   geometrical information for the FEM DG solver. ---*/
  CMeshFEM_DG *DGGeometry = dynamic_cast<CMeshFEM_DG *>(geometry);

  unsigned short nDim = DGGeometry->GetnDim();

  unsigned long nVolElemOwned = DGGeometry->GetNVolElemOwned();
  CVolumeElementFEM *volElem = DGGeometry->GetVolElem();

  /*--- Create the map from the global DOF ID to the local index. ---*/
  vector<su2double> DOFsCoords;
  vector<unsigned long> globalID;

  /*--- Set the global ID and the coordinates of the local DOFs
        by looping over the owned volume elements. ---*/
  unsigned long nLocalPoint = 0;
  for(unsigned long l=0; l<nVolElemOwned; ++l) {
    for(unsigned short j=0; j<volElem[l].nDOFsSol; ++j) {

      const unsigned long globalIndex = volElem[l].offsetDOFsSolGlobal + j;
      globalID.push_back(globalIndex);

      const su2double *coor = volElem[l].coorSolDOFs.data() + j*nDim;
      for(unsigned short k=0; k<nDim; ++k) DOFsCoords.push_back(coor[k]);

      nLocalPoint++;
    }
  }

#ifndef HAVE_MPI

  /*--- In serial, the single process has access to all geometry, so simply
   load the coordinates into the data structure. ---*/

  // need to double check for halos in parallel

  nGlobal_Poin = nLocalPoint;
  nGlobal_Doma = nLocalPoint;

  /*--- Allocate the coordinates data structure. ---*/

  Coords = new su2double*[nDim];
  for (iDim = 0; iDim < nDim; iDim++) {
    Coords[iDim] = new su2double[nGlobal_Poin];
  }

  /*--- Loop over the mesh to collect the coords of the local points ---*/

  for (iPoint = 0; iPoint < nLocalPoint; iPoint++) {

    /*--- Check if the node belongs to the domain (i.e, not a halo node).
     Sort by the global index, even in serial there is a renumbering (e.g. RCM). ---*/

    /*--- Retrieve the current coordinates at this node. ---*/

    unsigned long iGlobal_Index = globalID[iPoint];

    for (iDim = 0; iDim < nDim; iDim++) {
      Coords[iDim][iGlobal_Index] = DOFsCoords[iPoint*nDim+iDim];

      /*--- If US system, the output should be in inches ---*/

      if ((config->GetSystemMeasurements() == US) && (config->GetKind_SU2() != SU2_DEF)) {
        Coords[iDim][iGlobal_Index] *= 12.0;
      }
    }
  }

#else

  /*--- MPI preprocessing ---*/
  int iProcessor, nProcessor = size;
  unsigned long jPoint;

  /*--- Local variables needed for merging the geometry with MPI. ---*/

  unsigned long Buffer_Send_nPoin[1], *Buffer_Recv_nPoin = NULL;
  unsigned long MaxLocalPoint = 0;
  unsigned long iGlobal_Index = 0, nBuffer_Scalar = 0;

  if (rank == MASTER_NODE) Buffer_Recv_nPoin = new unsigned long[nProcessor];

  /*--- Search all send/recv boundaries on this partition for any periodic
   nodes that were part of the original domain. We want to recover these
   for visualization purposes. ---*/

  // Need to double check halo layers

  Buffer_Send_nPoin[0] = nLocalPoint;

  /*--- Communicate the total number of nodes on this domain. ---*/

  SU2_MPI::Gather(&Buffer_Send_nPoin, 1, MPI_UNSIGNED_LONG,
                  Buffer_Recv_nPoin, 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&nLocalPoint, &MaxLocalPoint, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);

  if (rank == MASTER_NODE) {
    nGlobal_Doma = 0;
    for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
      nGlobal_Doma += Buffer_Recv_nPoin[iProcessor];
    }
  }
  nBuffer_Scalar = MaxLocalPoint;

  /*--- Send and Recv buffers. ---*/

  su2double *Buffer_Send_X = new su2double[MaxLocalPoint];
  su2double *Buffer_Recv_X = NULL;

  su2double *Buffer_Send_Y = new su2double[MaxLocalPoint];
  su2double *Buffer_Recv_Y = NULL;

  su2double *Buffer_Send_Z = NULL, *Buffer_Recv_Z = NULL;
  if (nDim == 3) Buffer_Send_Z = new su2double[MaxLocalPoint];

  unsigned long *Buffer_Send_GlobalIndex = new unsigned long[MaxLocalPoint];
  unsigned long *Buffer_Recv_GlobalIndex = NULL;

  /*--- Prepare the receive buffers in the master node only. ---*/

  if (rank == MASTER_NODE) {

    Buffer_Recv_X = new su2double[nProcessor*MaxLocalPoint];
    Buffer_Recv_Y = new su2double[nProcessor*MaxLocalPoint];
    if (nDim == 3) Buffer_Recv_Z = new su2double[nProcessor*MaxLocalPoint];
    Buffer_Recv_GlobalIndex = new unsigned long[nProcessor*MaxLocalPoint];

    /*--- Sum total number of nodes to be written and allocate arrays ---*/
    nGlobal_Poin = 0;
    for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
      nGlobal_Poin += Buffer_Recv_nPoin[iProcessor];
    }
    Coords = new su2double*[nDim];
    for (iDim = 0; iDim < nDim; iDim++) {
      Coords[iDim] = new su2double[nGlobal_Poin];
    }
  }

  /*--- Main communication routine. Loop over each coordinate and perform
   the MPI comm. Temporary 1-D buffers are used to send the coordinates at
   all nodes on each partition to the master node. These are then unpacked
   by the master and sorted by global index in one large n-dim. array. ---*/

  /*--- Loop over this partition to collect the coords of the local points. ---*/
  jPoint = 0;
  for (iPoint = 0; iPoint < nLocalPoint; iPoint++) {

    /*--- Check for halos and write only if requested ---*/
    //    if (!Local_Halo[iPoint] || Wrt_Halo) {

    /*--- Load local coords into the temporary send buffer. These were stored above ---*/
    Buffer_Send_X[jPoint] = DOFsCoords[iPoint*nDim+0];
    Buffer_Send_Y[jPoint] = DOFsCoords[iPoint*nDim+1];
    if (nDim == 3) Buffer_Send_Z[jPoint] = DOFsCoords[iPoint*nDim+2];

    /*--- If US system, the output should be in inches ---*/

    if ((config->GetSystemMeasurements() == US) && (config->GetKind_SU2() != SU2_DEF)) {
      Buffer_Send_X[jPoint] *= 12.0;
      Buffer_Send_Y[jPoint] *= 12.0;
      if (nDim == 3) Buffer_Send_Z[jPoint] *= 12.0;
    }

    /*--- Store the global index for this local node. ---*/
    Buffer_Send_GlobalIndex[jPoint] = globalID[iPoint];

    /*--- Increment jPoint as the counter. We need this because iPoint
     may include halo nodes that we skip over during this loop. ---*/
    jPoint++;
    //    }
  }

  /*--- Gather the coordinate data on the master node using MPI. ---*/

  SU2_MPI::Gather(Buffer_Send_X, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_X, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  SU2_MPI::Gather(Buffer_Send_Y, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Y, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  if (nDim == 3) {
    SU2_MPI::Gather(Buffer_Send_Z, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Z, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  }
  SU2_MPI::Gather(Buffer_Send_GlobalIndex, nBuffer_Scalar, MPI_UNSIGNED_LONG, Buffer_Recv_GlobalIndex, nBuffer_Scalar, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);

  /*--- The master node unpacks and sorts this variable by global index ---*/

  if (rank == MASTER_NODE) {
    jPoint = 0;
    for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
      for (iPoint = 0; iPoint < Buffer_Recv_nPoin[iProcessor]; iPoint++) {
        /*--- Get global index, then loop over each variable and store ---*/
        iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];
        if (iGlobal_Index >= nGlobal_Poin){
          cout << iGlobal_Index << " " << nGlobal_Poin << endl;
        }
        Coords[0][iGlobal_Index] = Buffer_Recv_X[jPoint];
        Coords[1][iGlobal_Index] = Buffer_Recv_Y[jPoint];
        if (nDim == 3) Coords[2][iGlobal_Index] = Buffer_Recv_Z[jPoint];
        jPoint++;
      }
      /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
      jPoint = (iProcessor+1)*nBuffer_Scalar;
    }
  }

  /*--- Immediately release the temporary data buffers. ---*/

  //  delete [] Local_Halo;
  delete [] Buffer_Send_X;
  delete [] Buffer_Send_Y;
  if (Buffer_Send_Z != NULL) delete [] Buffer_Send_Z;
  delete [] Buffer_Send_GlobalIndex;
  if (rank == MASTER_NODE) {
    delete [] Buffer_Recv_X;
    delete [] Buffer_Recv_Y;
    if (Buffer_Recv_Z != NULL)  delete [] Buffer_Recv_Z;
    delete [] Buffer_Recv_GlobalIndex;
    delete [] Buffer_Recv_nPoin;
  }

#endif

}

void COutputLegacy::MergeVolumetricConnectivity_FEM(CConfig *config, CGeometry *geometry, unsigned short Elem_Type) {

  int iProcessor;
  unsigned short NODES_PER_ELEMENT = 0;
  unsigned long iNode, jNode;
  unsigned long iElem = 0;
  unsigned long nLocalElem = 0, nElem_Total = 0;

  unsigned long jElem;

  unsigned long Buffer_Send_nElem[1], *Buffer_Recv_nElem = NULL;
  unsigned long nBuffer_Scalar = 0;
  unsigned long kNode = 0, kElem = 0;
  unsigned long MaxLocalElem = 0;

  bool Wrt_Halo = config->GetWrt_Halo();
  bool *Write_Elem = NULL;

  int *Conn_Elem = NULL;

  /*--- Create an object of the class CMeshFEM_DG and retrieve the necessary
   geometrical information for the FEM DG solver. ---*/
  CMeshFEM_DG *DGGeometry = dynamic_cast<CMeshFEM_DG *>(geometry);

  unsigned long nVolElemOwned = DGGeometry->GetNVolElemOwned();

  CVolumeElementFEM *volElem = DGGeometry->GetVolElem();

  const CFEMStandardElement *standardElementsSol = DGGeometry->GetStandardElementsSol();

  /* Define the vectors for the connectivity of the local linear subelements. */
  vector<unsigned long> volumeConn;

  /*--- Create the map from the global DOF ID to the local index. ---*/
  map<unsigned long, unsigned long> mapLocal2Global;
  unsigned long ii = 0;
  for(unsigned long i=0; i<nVolElemOwned; ++i) {
    for(unsigned short j=0; j<volElem[i].nDOFsSol; ++j, ++ii) {
      mapLocal2Global[ii] = volElem[i].offsetDOFsSolGlobal+j;
    }
  }

  /*--- Counter for keeping track of the number of this element locally. ---*/
  nLocalElem = 0;
  for(unsigned long i=0; i<nVolElemOwned; ++i) {

    /* Determine the necessary data from the corresponding standard elem,
     such as the number of linear sub elements, the number of DOFs per
     linear sub element and the corresponding local connectivity. */
    const unsigned short ind       = volElem[i].indStandardElement;
    const unsigned short VTK_Type1 = standardElementsSol[ind].GetVTK_Type1();
    const unsigned short VTK_Type2 = standardElementsSol[ind].GetVTK_Type2();

    /*--- Only store the linear sub elements if they are of
     the current type that we are merging. ---*/
    if (Elem_Type == VTK_Type1) {

      const unsigned short nSubElems       = standardElementsSol[ind].GetNSubElemsType1();
      const unsigned short nDOFsPerSubElem = standardElementsSol[ind].GetNDOFsPerSubElem(VTK_Type1);
      const unsigned short *connSubElems   = standardElementsSol[ind].GetSubConnType1();

      /* Loop over the number of subfaces and store the required data. */
      unsigned short kk = 0;
      for(unsigned short j=0; j<nSubElems; ++j) {

        /*--- Store the global index for the surface conn ---*/
        for(unsigned short k=0; k<nDOFsPerSubElem; ++k, ++kk)
          volumeConn.push_back(mapLocal2Global[volElem[i].offsetDOFsSolLocal+connSubElems[kk]]);

        nLocalElem++;
      }

    } else if (Elem_Type == VTK_Type2) {

      const unsigned short nSubElems       = standardElementsSol[ind].GetNSubElemsType2();
      const unsigned short nDOFsPerSubElem = standardElementsSol[ind].GetNDOFsPerSubElem(VTK_Type2);
      const unsigned short *connSubElems   = standardElementsSol[ind].GetSubConnType2();

      /* Loop over the number of subfaces and store the required data. */
      unsigned short kk = 0;
      for(unsigned short j=0; j<nSubElems; ++j) {

        /*--- Store the global index for the surface conn ---*/
        for(unsigned short k=0; k<nDOFsPerSubElem; ++k, ++kk)
          volumeConn.push_back(mapLocal2Global[volElem[i].offsetDOFsSolLocal+connSubElems[kk]]);

        nLocalElem++;
      }
    }
  }

  /*--- Store the local number of this element type and the number of nodes
   per this element type. In serial, this will be the total number of this
   element type in the entire mesh. In parallel, it is the number on only
   the current partition. ---*/

  switch (Elem_Type) {
    case TRIANGLE:
      NODES_PER_ELEMENT = N_POINTS_TRIANGLE;
      break;
    case QUADRILATERAL:
      NODES_PER_ELEMENT = N_POINTS_QUADRILATERAL;
      break;
    case TETRAHEDRON:
      NODES_PER_ELEMENT = N_POINTS_TETRAHEDRON;
      break;
    case HEXAHEDRON:
      NODES_PER_ELEMENT = N_POINTS_HEXAHEDRON;
      break;
    case PRISM:
      NODES_PER_ELEMENT = N_POINTS_PRISM;
      break;
    case PYRAMID:
      NODES_PER_ELEMENT = N_POINTS_PYRAMID;
      break;
    default:
      SU2_MPI::Error("Unrecognized element type", CURRENT_FUNCTION);
  }

  /*--- Find the max number of this element type among all
   partitions and set up buffers. ---*/

  Buffer_Send_nElem[0] = nLocalElem;
  if (rank == MASTER_NODE) Buffer_Recv_nElem = new unsigned long[size];

#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&nLocalElem, &MaxLocalElem, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
  SU2_MPI::Gather(&Buffer_Send_nElem, 1, MPI_UNSIGNED_LONG, Buffer_Recv_nElem, 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
#else
  MaxLocalElem = nLocalElem;
  Buffer_Recv_nElem[0] = Buffer_Send_nElem[0];
#endif

  nBuffer_Scalar = MaxLocalElem*NODES_PER_ELEMENT;

  /*--- Send and Recv buffers ---*/

  unsigned long *Buffer_Send_Elem = new unsigned long[nBuffer_Scalar];
  unsigned long *Buffer_Recv_Elem = NULL;

  unsigned short *Buffer_Send_Halo = new unsigned short[MaxLocalElem];
  unsigned short *Buffer_Recv_Halo = NULL;

  /*--- Prepare the receive buffers on the master node only. ---*/

  if (rank == MASTER_NODE) {
    Buffer_Recv_Elem = new unsigned long[size*nBuffer_Scalar];
    Buffer_Recv_Halo = new unsigned short[size*MaxLocalElem];
    if (MaxLocalElem > 0) Conn_Elem = new int[size*MaxLocalElem*NODES_PER_ELEMENT];
  }

  /*--- Loop over all elements in this partition and load the
   elements of the current type into the buffer to be sent to
   the master node. ---*/

  jNode = 0;
  for (iElem = 0; iElem < nLocalElem; iElem++) {

    /*--- Loop over all nodes in this element and load the
     connectivity into the send buffer. ---*/
    Buffer_Send_Halo[iElem] = false;
    for (iNode = 0; iNode < NODES_PER_ELEMENT; iNode++) {

      /*--- Store the global index values directly. ---*/

      Buffer_Send_Elem[jNode] = volumeConn[iElem*NODES_PER_ELEMENT+iNode];

      /*--- Increment jNode as the counter. We need this because iElem
       may include other elements that we skip over during this loop. ---*/

      jNode++;
    }
  }

  /*--- Gather the element connectivity information. ---*/

#ifdef HAVE_MPI
  SU2_MPI::Gather(Buffer_Send_Elem, nBuffer_Scalar, MPI_UNSIGNED_LONG, Buffer_Recv_Elem, nBuffer_Scalar, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
  SU2_MPI::Gather(Buffer_Send_Halo, MaxLocalElem, MPI_UNSIGNED_SHORT, Buffer_Recv_Halo, MaxLocalElem, MPI_UNSIGNED_SHORT, MASTER_NODE, MPI_COMM_WORLD);
#else
  for (unsigned long iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Elem[iPoint] = Buffer_Send_Elem[iPoint];
  for (unsigned long iPoint = 0; iPoint < MaxLocalElem; iPoint++) Buffer_Recv_Halo[iPoint] = Buffer_Send_Halo[iPoint];
#endif

  /*--- The master node unpacks and sorts the connectivity. ---*/

  if (rank == MASTER_NODE) {

    /*---  We need to remove any duplicate elements (halo cells) that
     exist on multiple partitions. Start by initializing all elements
     to the "write" state by using a boolean array. ---*/

    Write_Elem = new bool[size*MaxLocalElem];
    for (iElem = 0; iElem < size*MaxLocalElem; iElem++) {
      Write_Elem[iElem] = true;
    }

    /*--- Remove the rind layer from the solution only if requested ---*/

    if (!Wrt_Halo) {

      /*--- Loop for flagging duplicate elements so that they are not
       included in the final connectivity list. ---*/

      kElem = 0;
      for (iProcessor = 0; iProcessor < size; iProcessor++) {
        for (iElem = 0; iElem < Buffer_Recv_nElem[iProcessor]; iElem++) {

          /*--- Check if this element was marked as a halo. ---*/
          if (Buffer_Recv_Halo[kElem+iElem])
            Write_Elem[kElem+iElem] = false;

        }
        kElem = (iProcessor+1)*MaxLocalElem;
      }
    }

    /*--- Store the unique connectivity list for this element type. ---*/

    jNode = 0; kNode = 0; jElem = 0; nElem_Total = 0;
    for (iProcessor = 0; iProcessor < size; iProcessor++) {
      for (iElem = 0; iElem < Buffer_Recv_nElem[iProcessor]; iElem++) {

        /*--- Only write the elements that were flagged for it. ---*/
        if (Write_Elem[jElem+iElem]) {

          /*--- Increment total count for this element type ---*/
          nElem_Total++;

          /*--- Get global index, then loop over each variable and store.
           Note that we are adding one to the index value because CGNS/Tecplot
           use 1-based indexing.---*/

          for (iNode = 0; iNode < NODES_PER_ELEMENT; iNode++) {
            Conn_Elem[kNode] = (int)Buffer_Recv_Elem[jNode+iElem*NODES_PER_ELEMENT+iNode] + 1;
            kNode++;
          }
        }
      }
      /*--- Adjust jNode to index of next proc's data in the buffers. ---*/
      jElem = (iProcessor+1)*MaxLocalElem;
      jNode = (iProcessor+1)*nBuffer_Scalar;
    }
  }

  /*--- Immediately release the temporary buffers. ---*/
  delete [] Buffer_Send_Elem;
  delete [] Buffer_Send_Halo;
  if (rank == MASTER_NODE) {
    delete [] Buffer_Recv_nElem;
    delete [] Buffer_Recv_Elem;
    delete [] Buffer_Recv_Halo;
    delete [] Write_Elem;
  }

  /*--- Store the particular global element count in the class data,
   and set the class data pointer to the connectivity array. ---*/

  if (rank == MASTER_NODE) {
    switch (Elem_Type) {
      case TRIANGLE:
        nGlobal_Tria = nElem_Total;
        if (nGlobal_Tria > 0) Conn_Tria = Conn_Elem;
        break;
      case QUADRILATERAL:
        nGlobal_Quad = nElem_Total;
        if (nGlobal_Quad > 0) Conn_Quad = Conn_Elem;
        break;
      case TETRAHEDRON:
        nGlobal_Tetr = nElem_Total;
        if (nGlobal_Tetr > 0) Conn_Tetr = Conn_Elem;
        break;
      case HEXAHEDRON:
        nGlobal_Hexa = nElem_Total;
        if (nGlobal_Hexa > 0) Conn_Hexa = Conn_Elem;
        break;
      case PRISM:
        nGlobal_Pris = nElem_Total;
        if (nGlobal_Pris > 0) Conn_Pris = Conn_Elem;
        break;
      case PYRAMID:
        nGlobal_Pyra = nElem_Total;
        if (nGlobal_Pyra > 0) Conn_Pyra = Conn_Elem;
        break;
      default:
        SU2_MPI::Error("Unrecognized element type", CURRENT_FUNCTION);
    }
  }

}

void COutputLegacy::MergeSurfaceConnectivity_FEM(CConfig *config, CGeometry *geometry, unsigned short Elem_Type) {

  unsigned short NODES_PER_ELEMENT = 0;

  unsigned long iNode, jNode;
  unsigned long iElem = 0;
  unsigned long nLocalElem = 0, nElem_Total = 0;

  int iProcessor;
  unsigned long jElem;

  unsigned long Buffer_Send_nElem[1], *Buffer_Recv_nElem = NULL;
  unsigned long nBuffer_Scalar = 0;
  unsigned long kNode = 0, kElem = 0;
  unsigned long MaxLocalElem = 0;

  bool Wrt_Halo = config->GetWrt_Halo();
  bool *Write_Elem = NULL;

  int *Conn_Elem = NULL;

  /*--- Create an object of the class CMeshFEM_DG and retrieve the necessary
   geometrical information for the FEM DG solver. ---*/
  CMeshFEM_DG *DGGeometry = dynamic_cast<CMeshFEM_DG *>(geometry);

  unsigned long nVolElemOwned = DGGeometry->GetNVolElemOwned();

  CVolumeElementFEM *volElem = DGGeometry->GetVolElem();

  const CBoundaryFEM *boundaries = DGGeometry->GetBoundaries();

  const CFEMStandardBoundaryFace *standardBoundaryFacesSol = DGGeometry->GetStandardBoundaryFacesSol();

  /*--- Create the map from the global DOF ID to the local index. ---*/
  map<unsigned long, unsigned long> mapLocal2Global;
  unsigned long ii = 0;
  for(unsigned long i=0; i<nVolElemOwned; ++i) {
    for(unsigned short j=0; j<volElem[i].nDOFsSol; ++j, ++ii) {
      mapLocal2Global[ii] = volElem[i].offsetDOFsSolGlobal+j;
    }
  }

  /* Define the vectors for the connectivity of the local linear subelements. */
  vector<unsigned long>  surfaceConn;

  /*--- Counter for keeping track of the number of this element locally. ---*/
  nLocalElem = 0;

  /* Loop over the boundary markers. */
  for(unsigned short iMarker=0; iMarker < config->GetnMarker_All(); ++iMarker) {
    if( !boundaries[iMarker].periodicBoundary ) {

      /* Check for markers to be plotted. */
      if (config->GetMarker_All_Plotting(iMarker) == YES) {

        /* Loop over the surface elements of this marker. */
        const vector<CSurfaceElementFEM> &surfElem = boundaries[iMarker].surfElem;
        for(unsigned long i=0; i<surfElem.size(); ++i) {

          /* Determine the necessary data from the corresponding standard face,
           such as the number of linear subfaces, the number of DOFs per
           linear subface and the corresponding local connectivity. */
          const unsigned short ind           = surfElem[i].indStandardElement;
          const unsigned short VTK_Type      = standardBoundaryFacesSol[ind].GetVTK_Type();

          /*--- Only store the linear sub-elements if they are of
           the current type that we are merging. ---*/
          if (VTK_Type == Elem_Type) {

            const unsigned short nSubFaces     = standardBoundaryFacesSol[ind].GetNSubFaces();
            const unsigned short nDOFsPerFace  = standardBoundaryFacesSol[ind].GetNDOFsPerSubFace();
            const unsigned short *connSubFaces = standardBoundaryFacesSol[ind].GetSubFaceConn();

            /* Loop over the number of subfaces and store the required data. */
            unsigned short ii = 0;
            for(unsigned short j=0; j<nSubFaces; ++j) {

              /*--- Store the global index for the surface conn ---*/
              for(unsigned short k=0; k<nDOFsPerFace; ++k, ++ii)
                surfaceConn.push_back(mapLocal2Global[surfElem[i].DOFsSolFace[connSubFaces[ii]]]);

              nLocalElem++;
            }
          }
        }
      }
    }
  }

  switch (Elem_Type) {
    case LINE:
      NODES_PER_ELEMENT = N_POINTS_LINE;
      break;
    case TRIANGLE:
      NODES_PER_ELEMENT = N_POINTS_TRIANGLE;
      break;
    case QUADRILATERAL:
      NODES_PER_ELEMENT = N_POINTS_QUADRILATERAL;
      break;
    default:
      SU2_MPI::Error("Unrecognized element type", CURRENT_FUNCTION);
  }

  /*--- Find the max number of this element type among all
   partitions and set up buffers. ---*/

  Buffer_Send_nElem[0] = nLocalElem;
  if (rank == MASTER_NODE) Buffer_Recv_nElem = new unsigned long[size];

#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&nLocalElem, &MaxLocalElem, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
  SU2_MPI::Gather(&Buffer_Send_nElem, 1, MPI_UNSIGNED_LONG, Buffer_Recv_nElem, 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
#else
  MaxLocalElem = nLocalElem;
  Buffer_Recv_nElem[0] = Buffer_Send_nElem[0];
#endif

  nBuffer_Scalar = MaxLocalElem*NODES_PER_ELEMENT;

  /*--- Send and Recv buffers ---*/

  unsigned long *Buffer_Send_Elem = new unsigned long[nBuffer_Scalar];
  unsigned long *Buffer_Recv_Elem = NULL;

  unsigned short *Buffer_Send_Halo = new unsigned short[MaxLocalElem];
  unsigned short *Buffer_Recv_Halo = NULL;

  /*--- Prepare the receive buffers on the master node only. ---*/

  if (rank == MASTER_NODE) {
    Buffer_Recv_Elem = new unsigned long[size*nBuffer_Scalar];
    Buffer_Recv_Halo = new unsigned short[size*MaxLocalElem];
    if (MaxLocalElem > 0) Conn_Elem = new int[size*MaxLocalElem*NODES_PER_ELEMENT];
  }

  /*--- Loop over all elements in this partition and load the
   elements of the current type into the buffer to be sent to
   the master node. ---*/
  jNode = 0;
  for (iElem = 0; iElem < nLocalElem; iElem++) {

    /*--- Loop over all nodes in this element and load the
     connectivity into the send buffer. ---*/
    Buffer_Send_Halo[iElem] = false;
    for (iNode = 0; iNode < NODES_PER_ELEMENT; iNode++) {

      /*--- Store the global index values directly. ---*/
      Buffer_Send_Elem[jNode] = surfaceConn[iElem*NODES_PER_ELEMENT+iNode];

      /*--- Increment jNode as the counter for the buffer array. ---*/
      jNode++;
    }
  }

  /*--- Gather the element connectivity information. ---*/

#ifdef HAVE_MPI
  SU2_MPI::Gather(Buffer_Send_Elem, nBuffer_Scalar, MPI_UNSIGNED_LONG, Buffer_Recv_Elem, nBuffer_Scalar, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
  SU2_MPI::Gather(Buffer_Send_Halo, MaxLocalElem, MPI_UNSIGNED_SHORT, Buffer_Recv_Halo, MaxLocalElem, MPI_UNSIGNED_SHORT, MASTER_NODE, MPI_COMM_WORLD);
#else
  for (unsigned long iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Elem[iPoint] = Buffer_Send_Elem[iPoint];
  for (unsigned long iPoint = 0; iPoint < MaxLocalElem; iPoint++) Buffer_Recv_Halo[iPoint] = Buffer_Send_Halo[iPoint];
#endif

  /*--- The master node unpacks and sorts the connectivity. ---*/

  if (rank == MASTER_NODE) {

    /*---  We need to remove any duplicate elements (halo cells) that
     exist on multiple partitions. Start by initializing all elements
     to the "write" state by using a boolean array. ---*/

    Write_Elem = new bool[size*MaxLocalElem];
    for (iElem = 0; iElem < size*MaxLocalElem; iElem++) {
      Write_Elem[iElem] = true;
    }

    /*--- Remove the rind layer from the solution only if requested ---*/

    if (!Wrt_Halo) {

      /*--- Loop for flagging duplicate elements so that they are not
       included in the final connectivity list. ---*/

      kElem = 0;
      for (iProcessor = 0; iProcessor < size; iProcessor++) {
        for (iElem = 0; iElem < Buffer_Recv_nElem[iProcessor]; iElem++) {

          /*--- Check if this element was marked as a halo. ---*/
          if (Buffer_Recv_Halo[kElem+iElem])
            Write_Elem[kElem+iElem] = false;

        }
        kElem = (iProcessor+1)*MaxLocalElem;
      }
    }

    /*--- Store the unique connectivity list for this element type. ---*/

    jNode = 0; kNode = 0; jElem = 0; nElem_Total = 0;
    for (iProcessor = 0; iProcessor < size; iProcessor++) {
      for (iElem = 0; iElem < Buffer_Recv_nElem[iProcessor]; iElem++) {

        /*--- Only write the elements that were flagged for it. ---*/
        if (Write_Elem[jElem+iElem]) {

          /*--- Increment total count for this element type ---*/
          nElem_Total++;

          /*--- Get global index, then loop over each variable and store.
           Note that we are adding one to the index value because CGNS/Tecplot
           use 1-based indexing.---*/

          for (iNode = 0; iNode < NODES_PER_ELEMENT; iNode++) {
            Conn_Elem[kNode] = (int)Buffer_Recv_Elem[jNode+iElem*NODES_PER_ELEMENT+iNode] + 1;
            kNode++;
          }
        }
      }
      /*--- Adjust jNode to index of next proc's data in the buffers. ---*/
      jElem = (iProcessor+1)*MaxLocalElem;
      jNode = (iProcessor+1)*nBuffer_Scalar;
    }
  }

  /*--- Immediately release the temporary buffers. ---*/
  delete [] Buffer_Send_Elem;
  delete [] Buffer_Send_Halo;
  if (rank == MASTER_NODE) {
    delete [] Buffer_Recv_nElem;
    delete [] Buffer_Recv_Elem;
    delete [] Buffer_Recv_Halo;
    delete [] Write_Elem;
  }

  /*--- Store the particular global element count in the class data,
   and set the class data pointer to the connectivity array. ---*/

  if (rank == MASTER_NODE) {
    switch (Elem_Type) {
      case LINE:
        nGlobal_Line = nElem_Total;
        if (nGlobal_Line > 0) Conn_Line = Conn_Elem;
        break;
      case TRIANGLE:
        nGlobal_BoundTria = nElem_Total;
        if (nGlobal_BoundTria > 0) Conn_BoundTria = Conn_Elem;
        break;
      case QUADRILATERAL:
        nGlobal_BoundQuad = nElem_Total;
        if (nGlobal_BoundQuad > 0) Conn_BoundQuad = Conn_Elem;
        break;
      default:
        SU2_MPI::Error("Unrecognized element type", CURRENT_FUNCTION);
    }
  }
}

void COutputLegacy::MergeSolution_FEM(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned short val_iZone) {

  unsigned short Kind_Solver  = config->GetKind_Solver();
  unsigned short iVar = 0, FirstIndex = NONE, SecondIndex = NONE, ThirdIndex = NONE;
  unsigned short nVar_First = 0, nVar_Second = 0, nVar_Third = 0;

  unsigned long iPoint = 0, jPoint = 0;

  unsigned long Buffer_Send_nPoint[1], *Buffer_Recv_nPoint = NULL;
  unsigned long nLocalPoint = 0, MaxLocalPoint = 0;
  unsigned long iGlobal_Index = 0, nBuffer_Scalar = 0;

  int iProcessor;

  /*--- Prepare send buffers for the conservative variables. Need to
   find the total number of conservative variables and also the
   index for their particular solution container. ---*/

  switch (Kind_Solver) {
    case FEM_EULER : case FEM_NAVIER_STOKES: case FEM_LES: FirstIndex = FLOW_SOL; SecondIndex = NONE; ThirdIndex = NONE; break;
    default: SecondIndex = NONE; ThirdIndex = NONE; break;
  }

  nVar_First = solver[FirstIndex]->GetnVar();
  if (SecondIndex != NONE) nVar_Second = solver[SecondIndex]->GetnVar();
  if (ThirdIndex != NONE) nVar_Third = solver[ThirdIndex]->GetnVar();
  nVar_Consv = nVar_First + nVar_Second + nVar_Third;
  nVar_Total = nVar_Consv;

  /*--- Create an object of the class CMeshFEM_DG and retrieve the necessary
   geometrical information for the FEM DG solver. ---*/
  CMeshFEM_DG *DGGeometry = dynamic_cast<CMeshFEM_DG *>(geometry);

  unsigned long nVolElemOwned = DGGeometry->GetNVolElemOwned();
  CVolumeElementFEM *volElem = DGGeometry->GetVolElem();

  /*--- Create the map from the global DOF ID to the local index. ---*/
  map<unsigned long, unsigned long> mapLocal2Global;
  vector<su2double> DOFsSol;
  vector<unsigned long> globalID;

  /*--- Update the solution by looping over the owned volume elements. ---*/
  for(unsigned long l=0; l<nVolElemOwned; ++l) {

    /* Set the pointers for the residual and solution for this element. */
    const unsigned long offset  = nVar_Consv*volElem[l].offsetDOFsSolLocal;
    su2double *solDOFs          = solver[FirstIndex]->GetVecSolDOFs() + offset;

    /* Loop over the DOFs for this element and store the solution. */

    unsigned int i = 0;
    for(unsigned short j=0; j<volElem[l].nDOFsSol; ++j) {

      nLocalPoint++;
      const unsigned long globalIndex = volElem[l].offsetDOFsSolGlobal + j;
      globalID.push_back(globalIndex);

      for(unsigned short iVar=0; iVar<nVar_Consv; ++iVar, ++i) {
        DOFsSol.push_back(solDOFs[i]);
      }
    }
  }
  Buffer_Send_nPoint[0] = nLocalPoint;

  /*--- Each processor sends its local number of nodes to the master. ---*/

  if (rank == MASTER_NODE) Buffer_Recv_nPoint = new unsigned long[size];

#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&nLocalPoint, &MaxLocalPoint, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
  SU2_MPI::Gather(&Buffer_Send_nPoint, 1, MPI_UNSIGNED_LONG, Buffer_Recv_nPoint, 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
#else
  MaxLocalPoint = nLocalPoint;
  Buffer_Recv_nPoint[0] = Buffer_Send_nPoint[0];
#endif

  nBuffer_Scalar = MaxLocalPoint;

  /*--- Send and Recv buffers. ---*/

  su2double *Buffer_Send_Var = new su2double[MaxLocalPoint];
  su2double *Buffer_Recv_Var = NULL;

  unsigned long *Buffer_Send_GlobalIndex = new unsigned long[MaxLocalPoint];
  unsigned long *Buffer_Recv_GlobalIndex = NULL;

  /*--- Prepare the receive buffers in the master node only. ---*/

  if (rank == MASTER_NODE) {

    Buffer_Recv_Var = new su2double[size*MaxLocalPoint];
    Buffer_Recv_GlobalIndex = new unsigned long[size*MaxLocalPoint];

    /*--- Sum total number of nodes to be written and allocate arrays ---*/
    nGlobal_Poin = 0;
    for (iProcessor = 0; iProcessor < size; iProcessor++) {
      nGlobal_Poin += Buffer_Recv_nPoint[iProcessor];
    }
    Data = new su2double*[nVar_Total];
    for (iVar = 0; iVar < nVar_Total; iVar++) {
      Data[iVar] = new su2double[nGlobal_Poin];
    }
  }

  /*--- Main communication routine. Loop over each variable that has
   been requested by the user and perform the MPI comm. Temporary
   1-D buffers are used to send the solution for each variable at all
   nodes on each partition to the master node. These are then unpacked
   by the master and sorted by global index in one large n-dim. array. ---*/

  for (iVar = 0; iVar < nVar_Consv; iVar++) {

    /*--- Loop over this partition to collect the current variable ---*/

    jPoint = 0;
    for (iPoint = 0; iPoint < nLocalPoint; iPoint++) {

      /*--- Get this variable into the temporary send buffer. ---*/

      Buffer_Send_Var[jPoint] = DOFsSol[iPoint*nVar_Consv+iVar];

      /*--- Only send/recv the volumes & global indices during the first loop ---*/

      if (iVar == 0) {
        Buffer_Send_GlobalIndex[jPoint] = globalID[iPoint];
      }

      jPoint++;

    }

    /*--- Gather the data on the master node. ---*/

#ifdef HAVE_MPI
    SU2_MPI::Gather(Buffer_Send_Var, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Var, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
#else
    for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Var[iPoint] = Buffer_Send_Var[iPoint];
#endif

    if (iVar == 0) {
#ifdef HAVE_MPI
      SU2_MPI::Gather(Buffer_Send_GlobalIndex, nBuffer_Scalar, MPI_UNSIGNED_LONG, Buffer_Recv_GlobalIndex, nBuffer_Scalar, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
#else
      for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_GlobalIndex[iPoint] = Buffer_Send_GlobalIndex[iPoint];
#endif
    }

    /*--- The master node unpacks and sorts this variable by global index ---*/

    if (rank == MASTER_NODE) {
      jPoint = 0;
      for (iProcessor = 0; iProcessor < size; iProcessor++) {
        for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++) {

          /*--- Get global index, then loop over each variable and store ---*/

          iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];

          Data[iVar][iGlobal_Index] = Buffer_Recv_Var[jPoint];

          jPoint++;

        }
        /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
        jPoint = (iProcessor+1)*nBuffer_Scalar;
      }
    }

  }

  /*--- Immediately release the temporary buffers. ---*/

  delete [] Buffer_Send_Var;
  delete [] Buffer_Send_GlobalIndex;
  if (rank == MASTER_NODE) {
    delete [] Buffer_Recv_nPoint;
    delete [] Buffer_Recv_Var;
    delete [] Buffer_Recv_GlobalIndex;
  }

}

void COutputLegacy::MergeBaselineSolution_FEM(CConfig *config, CGeometry *geometry, CSolver *solver, unsigned short val_iZone) {

  /*--- Local variables needed on all processors ---*/
  unsigned short iVar;
  unsigned long iPoint = 0, jPoint = 0;

  unsigned long Buffer_Send_nPoint[1], *Buffer_Recv_nPoint = NULL;
  unsigned long nLocalPoint = 0, MaxLocalPoint = 0;
  unsigned long iGlobal_Index = 0, nBuffer_Scalar = 0;

  int iProcessor;

  /*--- We know the number of fields from reading the restart. ---*/

  nVar_Total = config->fields.size() - 1;

  /*--- Create an object of the class CMeshFEM_DG and retrieve the necessary
   geometrical information for the FEM DG solver. ---*/
  CMeshFEM_DG *DGGeometry = dynamic_cast<CMeshFEM_DG *>(geometry);

  unsigned long nVolElemOwned = DGGeometry->GetNVolElemOwned();
  CVolumeElementFEM *volElem = DGGeometry->GetVolElem();

  /*--- Create the map from the global DOF ID to the local index. ---*/
  vector<su2double> DOFsSol;
  vector<unsigned long> globalID;

  /*--- Update the solution by looping over the owned volume elements. ---*/
  for(unsigned long l=0; l<nVolElemOwned; ++l) {

    /* Set the pointers for the residual and solution for this element. */
    const unsigned long offset  = nVar_Total*volElem[l].offsetDOFsSolLocal;
    su2double *solDOFs          = solver->GetVecSolDOFs() + offset;

    /* Loop over the DOFs for this element and store the solution. */

    unsigned int i = 0;
    for(unsigned short j=0; j<volElem[l].nDOFsSol; ++j) {

      nLocalPoint++;
      const unsigned long globalIndex = volElem[l].offsetDOFsSolGlobal + j;
      globalID.push_back(globalIndex);

      for(unsigned short iVar=0; iVar<nVar_Total; ++iVar, ++i) {
        DOFsSol.push_back(solDOFs[i]);
      }
    }
  }
  Buffer_Send_nPoint[0] = nLocalPoint;

  /*--- Each processor sends its local number of nodes to the master. ---*/

  if (rank == MASTER_NODE) Buffer_Recv_nPoint = new unsigned long[size];

#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&nLocalPoint, &MaxLocalPoint, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
  SU2_MPI::Gather(&Buffer_Send_nPoint, 1, MPI_UNSIGNED_LONG, Buffer_Recv_nPoint, 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
#else
  MaxLocalPoint = nLocalPoint;
  Buffer_Recv_nPoint[0] = Buffer_Send_nPoint[0];
#endif

  nBuffer_Scalar = MaxLocalPoint;

  /*--- Send and Recv buffers. ---*/

  su2double *Buffer_Send_Var = new su2double[MaxLocalPoint];
  su2double *Buffer_Recv_Var = NULL;

  unsigned long *Buffer_Send_GlobalIndex = new unsigned long[MaxLocalPoint];
  unsigned long *Buffer_Recv_GlobalIndex = NULL;

  /*--- Prepare the receive buffers in the master node only. ---*/

  if (rank == MASTER_NODE) {

    Buffer_Recv_Var = new su2double[size*MaxLocalPoint];
    Buffer_Recv_GlobalIndex = new unsigned long[size*MaxLocalPoint];

    /*--- Sum total number of nodes to be written and allocate arrays ---*/
    nGlobal_Poin = 0;
    for (iProcessor = 0; iProcessor < size; iProcessor++) {
      nGlobal_Poin += Buffer_Recv_nPoint[iProcessor];
    }
    Data = new su2double*[nVar_Total];
    for (iVar = 0; iVar < nVar_Total; iVar++) {
      Data[iVar] = new su2double[nGlobal_Poin];
    }
  }

  /*--- Main communication routine. Loop over each variable that has
   been requested by the user and perform the MPI comm. Temporary
   1-D buffers are used to send the solution for each variable at all
   nodes on each partition to the master node. These are then unpacked
   by the master and sorted by global index in one large n-dim. array. ---*/

  for (iVar = 0; iVar < nVar_Total; iVar++) {

    /*--- Loop over this partition to collect the current variable ---*/

    jPoint = 0;
    for (iPoint = 0; iPoint < nLocalPoint; iPoint++) {

      /*--- Get this variable into the temporary send buffer. ---*/

      Buffer_Send_Var[jPoint] = DOFsSol[iPoint*nVar_Total+iVar];

      /*--- Only send/recv the volumes & global indices during the first loop ---*/

      if (iVar == 0) {
        Buffer_Send_GlobalIndex[jPoint] = globalID[iPoint];
      }

      jPoint++;

    }

    /*--- Gather the data on the master node. ---*/

#ifdef HAVE_MPI
    SU2_MPI::Gather(Buffer_Send_Var, nBuffer_Scalar, MPI_DOUBLE, Buffer_Recv_Var, nBuffer_Scalar, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
#else
    for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_Var[iPoint] = Buffer_Send_Var[iPoint];
#endif

    if (iVar == 0) {
#ifdef HAVE_MPI
      SU2_MPI::Gather(Buffer_Send_GlobalIndex, nBuffer_Scalar, MPI_UNSIGNED_LONG, Buffer_Recv_GlobalIndex, nBuffer_Scalar, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
#else
      for (iPoint = 0; iPoint < nBuffer_Scalar; iPoint++) Buffer_Recv_GlobalIndex[iPoint] = Buffer_Send_GlobalIndex[iPoint];
#endif
    }

    /*--- The master node unpacks and sorts this variable by global index ---*/

    if (rank == MASTER_NODE) {
      jPoint = 0;
      for (iProcessor = 0; iProcessor < size; iProcessor++) {
        for (iPoint = 0; iPoint < Buffer_Recv_nPoint[iProcessor]; iPoint++) {

          /*--- Get global index, then loop over each variable and store ---*/

          iGlobal_Index = Buffer_Recv_GlobalIndex[jPoint];

          Data[iVar][iGlobal_Index] = Buffer_Recv_Var[jPoint];

          jPoint++;

        }
        /*--- Adjust jPoint to index of next proc's data in the buffers. ---*/
        jPoint = (iProcessor+1)*nBuffer_Scalar;
      }
    }

  }

  /*--- Immediately release the temporary buffers. ---*/

  delete [] Buffer_Send_Var;
  delete [] Buffer_Send_GlobalIndex;
  if (rank == MASTER_NODE) {
    delete [] Buffer_Recv_nPoint;
    delete [] Buffer_Recv_Var;
    delete [] Buffer_Recv_GlobalIndex;
  }

}

void COutputLegacy::SetResult_Files_FEM(CSolver ****solver_container, CGeometry ***geometry, CConfig **config,
                                  unsigned long iExtIter, unsigned short val_nZone) {

  unsigned short iZone;

  for (iZone = 0; iZone < val_nZone; iZone++) {

    /*--- Check whether we are writing any output at all. We can
     disable all output to avoid serial bottleneck at scale. ---*/
    if (config[iZone]->GetWrt_Output()) {

      /*--- Flags identifying the types of files to be written. ---*/

      bool Wrt_Vol = config[iZone]->GetWrt_Vol_Sol();
      bool Wrt_Srf = config[iZone]->GetWrt_Srf_Sol();

      /*--- Get the file output format ---*/

      unsigned short FileFormat = config[iZone]->GetTabular_FileFormat();

      /*--- Merge the node coordinates and connectivity, if necessary. This
       is only performed if a volume solution file is requested, and it
       is active by default. ---*/

      if (Wrt_Vol || Wrt_Srf) {
        if (rank == MASTER_NODE) cout << endl << "Merging connectivities in the Master node." << endl;
        MergeConnectivity_FEM(config[iZone], geometry[iZone][MESH_0], iZone);
      }

      /*--- Merge coordinates of all grid nodes (excluding ghost points).
       The grid coordinates are always merged and included first in the
       restart files. ---*/

      if (rank == MASTER_NODE) cout << "Merging coordinates in the Master node." << endl;
      MergeCoordinates_FEM(config[iZone], geometry[iZone][MESH_0]);

      if ((rank == MASTER_NODE) && (Wrt_Vol || Wrt_Srf)) {
        if (FileFormat == TECPLOT_BINARY) {
          if (rank == MASTER_NODE) cout << "Writing Tecplot binary volume and surface mesh files." << endl;
//          SetTecplotBinary_DomainMesh(config[iZone], geometry[iZone][MESH_0], iZone);
//          SetTecplotBinary_SurfaceMesh(config[iZone], geometry[iZone][MESH_0], iZone);
        }
      }

      /*--- Merge the solution data needed for volume solutions and restarts ---*/

      if (rank == MASTER_NODE) cout << "Merging solution in the Master node." << endl;
      MergeSolution_FEM(config[iZone], geometry[iZone][MESH_0], solver_container[iZone][MESH_0], iZone);

      /*--- Write restart, or Tecplot files using the merged data.
       This data lives only on the master, and these routines are currently
       executed by the master proc alone (as if in serial). ---*/

      if (rank == MASTER_NODE) {

        /*--- Write a native restart file ---*/

        if (rank == MASTER_NODE) cout << "Writing SU2 native restart file." << endl;
        SetRestart(config[iZone], geometry[iZone][MESH_0], solver_container[iZone][MESH_0] , iZone);

        if (Wrt_Vol) {

          switch (FileFormat) {

            case TECPLOT:

              /*--- Write a Tecplot ASCII file ---*/

              if (rank == MASTER_NODE) cout << "Writing Tecplot ASCII file volume solution file." << endl;
//              SetTecplotASCII(config[iZone], geometry[iZone][MESH_0], solver_container[iZone][MESH_0], iZone, val_nZone, false);
              DeallocateConnectivity(config[iZone], geometry[iZone][MESH_0], false);
              break;


            case TECPLOT_BINARY:

              /*--- Write a Tecplot binary solution file ---*/

              if (rank == MASTER_NODE) cout << "Writing Tecplot binary volume solution file." << endl;
//              SetTecplotBinary_DomainSolution(config[iZone], geometry[iZone][MESH_0], iZone);
              DeallocateConnectivity(config[iZone], geometry[iZone][MESH_0], false);
              break;


            case PARAVIEW:

              /*--- Write a Paraview ASCII file ---*/

              if (rank == MASTER_NODE) cout << "Writing Paraview ASCII volume solution file." << endl;
//              SetParaview_ASCII(config[iZone], geometry[iZone][MESH_0], iZone, val_nZone, false);
              DeallocateConnectivity(config[iZone], geometry[iZone][MESH_0], false);
              break;

            case PARAVIEW_BINARY:

              /*--- Write a ParaView ASCII file instead for now. ---*/

              if (rank == MASTER_NODE) cout << "ParaView binary volume files not available in this mode." << endl;
              if (rank == MASTER_NODE) cout << " Writing ParaView ASCII volume solution file instead." << endl;
//              SetParaview_ASCII(config[iZone], geometry[iZone][MESH_0], iZone, val_nZone, false);
              DeallocateConnectivity(config[iZone], geometry[iZone][MESH_0], false);
              break;

            default:
              break;
          }

        }

        if (Wrt_Srf) {

          switch (FileFormat) {

            case TECPLOT:

              /*--- Write a Tecplot ASCII file ---*/

              if (rank == MASTER_NODE) cout << "Writing Tecplot ASCII surface solution file." << endl;
//              SetTecplotASCII(config[iZone], geometry[iZone][MESH_0], solver_container[iZone][MESH_0] , iZone, val_nZone, true);
              DeallocateConnectivity(config[iZone], geometry[iZone][MESH_0], true);
              break;

            case TECPLOT_BINARY:

              /*--- Write a Tecplot binary solution file ---*/

              if (rank == MASTER_NODE) cout << "Writing Tecplot binary surface solution file." << endl;
//              SetTecplotBinary_SurfaceSolution(config[iZone], geometry[iZone][MESH_0], iZone);
              DeallocateConnectivity(config[iZone], geometry[iZone][MESH_0], true);
              break;

            case PARAVIEW:

              /*--- Write a Paraview ASCII file ---*/

              if (rank == MASTER_NODE) cout << "Writing Paraview ASCII surface solution file." << endl;
//              SetParaview_ASCII(config[iZone], geometry[iZone][MESH_0], iZone, val_nZone, true);
              DeallocateConnectivity(config[iZone], geometry[iZone][MESH_0], true);
              break;

            case PARAVIEW_BINARY:

              /*--- Write a ParaView ASCII file instead for now. ---*/

              if (rank == MASTER_NODE) cout << "ParaView binary surface files not available in this mode." << endl;
              if (rank == MASTER_NODE) cout << " Writing ParaView ASCII surface solution file instead." << endl;
//              SetParaview_ASCII(config[iZone], geometry[iZone][MESH_0], iZone, val_nZone, true);
              DeallocateConnectivity(config[iZone], geometry[iZone][MESH_0], true);
              break;

            default:
              break;
          }

        }

        /*--- Release memory needed for merging the solution data. ---*/

        DeallocateCoordinates(config[iZone], geometry[iZone][MESH_0]);
        DeallocateSolution(config[iZone], geometry[iZone][MESH_0]);

      }

      /*--- Final broadcast (informing other procs that the base output
       file was written). ---*/

#ifdef HAVE_MPI
      SU2_MPI::Bcast(&wrote_base_file, 1, MPI_UNSIGNED_SHORT, MASTER_NODE, MPI_COMM_WORLD);
      SU2_MPI::Bcast(&wrote_surf_file, 1, MPI_UNSIGNED_SHORT, MASTER_NODE, MPI_COMM_WORLD);
#endif

    } else {
      if (rank == MASTER_NODE) cout << endl << "Restart and solution output disabled." << endl;
    }
  }
}

void COutputLegacy::SetBaselineResult_Files_FEM(CSolver ***solver, CGeometry ***geometry, CConfig **config,
                                          unsigned long iExtIter, unsigned short val_nZone) {

  unsigned short iZone, iInst, nInst;

  for (iZone = 0; iZone < val_nZone; iZone++) {

    nInst = config[iZone]->GetnTimeInstances();

    for (iInst = 0; iInst < nInst; iInst++) {

      config[iZone]->SetiInst(iInst);

      /*--- Flags identifying the types of files to be written. ---*/

      bool Wrt_Vol = config[iZone]->GetWrt_Vol_Sol();
      bool Wrt_Srf = config[iZone]->GetWrt_Srf_Sol();

      /*--- Get the file output format ---*/

      unsigned short FileFormat = config[iZone]->GetTabular_FileFormat();

      /*--- Merge the node coordinates and connectivity if necessary. This
       is only performed if a volume solution file is requested, and it
       is active by default. ---*/

      if ((Wrt_Vol || Wrt_Srf)) {
        if (rank == MASTER_NODE) cout << "Merging connectivities in the Master node." << endl;
        MergeConnectivity_FEM(config[iZone], geometry[iZone][iInst], iZone);
      }

      /*--- Merge the solution data needed for volume solutions and restarts ---*/

      if ((Wrt_Vol || Wrt_Srf)) {
        if (rank == MASTER_NODE) cout << "Merging solution in the Master node." << endl;
        MergeBaselineSolution_FEM(config[iZone], geometry[iZone][iInst], solver[iZone][iInst], iZone);
      }

      /*--- Write restart, Tecplot or Paraview files using the merged data.
       This data lives only on the master, and these routines are currently
       executed by the master proc alone (as if in serial). ---*/

      if (rank == MASTER_NODE) {

        if (Wrt_Vol) {

          switch (FileFormat) {

            case TECPLOT:

              /*--- Write a Tecplot ASCII file ---*/

              if (rank == MASTER_NODE) cout << "Writing Tecplot ASCII file (volume grid)." << endl;
//              SetTecplotASCII(config[iZone], geometry[iZone][iInst], &solver[iZone][iInst], iZone, val_nZone, false);
              DeallocateConnectivity(config[iZone], geometry[iZone][iInst], false);
              break;


            case TECPLOT_BINARY:

              /*--- Write a Tecplot binary solution file ---*/

              if (rank == MASTER_NODE) cout << "Writing Tecplot Binary file (volume grid)." << endl;
//              SetTecplotBinary_DomainMesh(config[iZone], geometry[iZone][iInst], iZone);
//              SetTecplotBinary_DomainSolution(config[iZone], geometry[iZone][iInst], iZone);
              break;


            case PARAVIEW:

              /*--- Write a Paraview ASCII file ---*/

              if (rank == MASTER_NODE) cout << "Writing Paraview ASCII file (volume grid)." << endl;
//              SetParaview_ASCII(config[iZone], geometry[iZone][iInst], iZone, val_nZone, false);
              DeallocateConnectivity(config[iZone], geometry[iZone][iInst], false);
              break;

            case PARAVIEW_BINARY:

              /*--- Write a ParaView ASCII file instead for now. ---*/

              if (rank == MASTER_NODE) cout << "ParaView binary (volume grid) files not available in this mode." << endl;
              if (rank == MASTER_NODE) cout << " Writing ParaView ASCII file instead." << endl;
//              SetParaview_ASCII(config[iZone], geometry[iZone][iInst], iZone, val_nZone, false);
              DeallocateConnectivity(config[iZone], geometry[iZone][iInst], false);
              break;

            default:
              break;
          }

        }

        if (Wrt_Srf) {

          switch (FileFormat) {

            case TECPLOT:

              /*--- Write a Tecplot ASCII file ---*/

              if (rank == MASTER_NODE) cout << "Writing Tecplot ASCII file (surface grid)." << endl;
//              SetTecplotASCII(config[iZone], geometry[iZone][iInst], &solver[iZone][iInst], iZone, val_nZone, true);
              DeallocateConnectivity(config[iZone], geometry[iZone][iInst], true);
              break;

            case TECPLOT_BINARY:

              /*--- Write a Tecplot binary solution file ---*/

              if (rank == MASTER_NODE) cout << "Writing Tecplot Binary file (surface grid)." << endl;
//              SetTecplotBinary_SurfaceMesh(config[iZone], geometry[iZone][iInst], iZone);
//              SetTecplotBinary_SurfaceSolution(config[iZone], geometry[iZone][iInst], iZone);
              break;

            case PARAVIEW:

              /*--- Write a Paraview ASCII file ---*/

              if (rank == MASTER_NODE) cout << "Writing Paraview ASCII file (surface grid)." << endl;
//              SetParaview_ASCII(config[iZone], geometry[iZone][iInst], iZone, val_nZone, true);
              DeallocateConnectivity(config[iZone], geometry[iZone][iInst], true);
              break;

            case PARAVIEW_BINARY:

              /*--- Write a ParaView ASCII file instead for now. ---*/

              if (rank == MASTER_NODE) cout << "ParaView binary (surface grid) files not available in this mode." << endl;
              if (rank == MASTER_NODE) cout << " Writing ParaView ASCII file instead." << endl;
//              SetParaview_ASCII(config[iZone], geometry[iZone][iInst], iZone, val_nZone, true);
              DeallocateConnectivity(config[iZone], geometry[iZone][iInst], true);
              break;

            default:
              break;
          }
        }

        if (FileFormat == TECPLOT_BINARY) {
          if (!wrote_base_file)
            DeallocateConnectivity(config[iZone], geometry[iZone][iInst], false);
          if (!wrote_surf_file)
            DeallocateConnectivity(config[iZone], geometry[iZone][iInst], wrote_surf_file);
        }

        if (Wrt_Vol || Wrt_Srf)
          DeallocateSolution(config[iZone], geometry[iZone][iInst]);
      }

      /*--- Final broadcast (informing other procs that the base output
       file was written). ---*/

#ifdef HAVE_MPI
      SU2_MPI::Bcast(&wrote_base_file, 1, MPI_UNSIGNED_SHORT, MASTER_NODE, MPI_COMM_WORLD);
#endif
    }
  }
}

void COutputLegacy::LoadLocalData_FEM(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned short val_iZone) {

  unsigned short iDim;
  unsigned short Kind_Solver = config->GetKind_Solver();
  unsigned short nDim = geometry->GetnDim();

  unsigned long iVar, jVar;
  unsigned long iPoint, jPoint, FirstIndex = NONE, SecondIndex = NONE;
  unsigned long nVar_First = 0, nVar_Second = 0, nVar_Consv_Par = 0;

  stringstream varname;

  /*--- Use a switch statement to decide how many solver containers we have
   in this zone for output. ---*/

  switch (Kind_Solver) {
    case FEM_EULER : case FEM_NAVIER_STOKES: case FEM_LES: FirstIndex = FLOW_SOL; SecondIndex = NONE; break;
    case FEM_RANS : FirstIndex = FLOW_SOL; SecondIndex = TURB_SOL; break;
    default: SecondIndex = NONE; break;
  }

  nVar_First = solver[FirstIndex]->GetnVar();
  if (SecondIndex != NONE) nVar_Second = solver[SecondIndex]->GetnVar();
  nVar_Consv_Par = nVar_First + nVar_Second;

  /*--------------------------------------------------------------------------*/
  /*--- Step 1: Register the variables that will be output. To register a  ---*/
  /*---         variable, two things are required. First, increment the    ---*/
  /*---         counter for the number of variables (nVar_Par), which      ---*/
  /*---         controls the size of the data structure allocation, i.e.,  ---*/
  /*---         the number of columns in an nPoint x nVar structure.       ---*/
  /*---         Second, add a name for the variable to the vector that     ---*/
  /*---         holds the string names.                                    ---*/
  /*--------------------------------------------------------------------------*/

  /*--- All output files first need the grid coordinates. ---*/

  nVar_Par  = 1; Variable_Names.push_back("x");
  nVar_Par += 1; Variable_Names.push_back("y");
  if (geometry->GetnDim() == 3) {
    nVar_Par += 1; Variable_Names.push_back("z");
  }

  /*--- At a mininum, the restarts and visualization files need the
   conservative variables, so these follow next. ---*/

  nVar_Par += nVar_Consv_Par;

  Variable_Names.push_back("Density");
  Variable_Names.push_back("X-Momentum");
  Variable_Names.push_back("Y-Momentum");
  if (geometry->GetnDim() == 3) Variable_Names.push_back("Z-Momentum");
  Variable_Names.push_back("Energy");

  /*--- Eventually, turbulence model from second container goes here. ---*/


  /*--- Add Pressure, Temperature, Cp, Mach. ---*/

  nVar_Par += 1;
  Variable_Names.push_back("Pressure");

  nVar_Par += 3;
  Variable_Names.push_back("Temperature");
  if (config->GetTabular_FileFormat() == TAB_CSV){
    Variable_Names.push_back("Pressure_Coefficient");
  } else {
    Variable_Names.push_back("C<sub>p</sub>");
  }
  Variable_Names.push_back("Mach");

  if ((Kind_Solver == FEM_NAVIER_STOKES) || (Kind_Solver == FEM_LES)){
    nVar_Par += 1;
    Variable_Names.push_back("Laminar_Viscosity");
  }
  if ((Kind_Solver == FEM_LES) && (config->GetKind_SGS_Model() != IMPLICIT_LES)){
    nVar_Par += 1;
    Variable_Names.push_back("Eddy_Viscosity");
  }

  if (solver[FLOW_SOL]->VerificationSolution) {
    if (solver[FLOW_SOL]->VerificationSolution->ExactSolutionKnown()) {
      nVar_Par += 2*nVar_Consv_Par;
      Variable_Names.push_back("Verification_Density");
      Variable_Names.push_back("Verification_Momentum_x");
      Variable_Names.push_back("Verification_Momentum_y");
      if (geometry->GetnDim() == 3) Variable_Names.push_back("Verification_Momentum_z");
      Variable_Names.push_back("Verification_Energy");
      Variable_Names.push_back("Error_Density");
      Variable_Names.push_back("Error_Momentum_x");
      Variable_Names.push_back("Error_Momentum_y");
      if (geometry->GetnDim() == 3) Variable_Names.push_back("Error_Momentum_z");
      Variable_Names.push_back("Error_Energy");
    }
  }

    /*--- New variables get registered here before the end of the loop. ---*/


  /*--- Create an object of the class CMeshFEM_DG and retrieve the necessary
   geometrical information for the FEM DG solver. ---*/

  CMeshFEM_DG *DGGeometry = dynamic_cast<CMeshFEM_DG *>(geometry);

  unsigned long nVolElemOwned = DGGeometry->GetNVolElemOwned();
  CVolumeElementFEM *volElem  = DGGeometry->GetVolElem();

  /*--- Get a pointer to the fluid model class from the DG-FEM solver
   so that we can access the states below. ---*/

  CFluidModel *DGFluidModel = solver[FLOW_SOL]->GetFluidModel();

  /*--- Allocate the local data structure now that we know how many
   variables are in the output. ---*/

  Local_Data = new su2double*[nLocalPoint_Sort];
  for (iPoint = 0; iPoint < nLocalPoint_Sort; iPoint++) {
    Local_Data[iPoint] = new su2double[nVar_Par];
  }

  /*--------------------------------------------------------------------------*/
  /*--- Step 2: Loop over all grid nodes and load up the desired data for  ---*/
  /*---         the restart and vizualization files. Note that we need to  ---*/
  /*---         increment the iVar variable after each variable load.      ---*/
  /*---         The idea is that we're filling up the columns of field     ---*/
  /*---         data for each iPoint (row) of the data structure.          ---*/
  /*---         This data will then be sorted, communicated, and written   ---*/
  /*---         to files automatically after this routine. Note that the   ---*/
  /*---         ordering of the data loading MUST match the order of the   ---*/
  /*---         variable registration above for the files to be correct.   ---*/
  /*--------------------------------------------------------------------------*/

  jPoint = 0;

  /*--- Access the solution by looping over the owned volume elements. ---*/

  for(unsigned long l=0; l<nVolElemOwned; ++l) {

    /* Set the pointers for the solution for this element. */

    const unsigned long offset = nVar_First*volElem[l].offsetDOFsSolLocal;
    su2double *solDOFs         = solver[FirstIndex]->GetVecSolDOFs() + offset;

    for(unsigned short j=0; j<volElem[l].nDOFsSol; ++j) {

      /*--- Restart the column index with each new point. ---*/

      iVar = 0;

      /*--- Get the conservative variables for this particular DOF. ---*/

      const su2double *U = solDOFs + j*nVar_First;

      /*--- Load the coordinate values of the solution DOFs. ---*/

      const su2double *coor = volElem[l].coorSolDOFs.data() + j*nDim;
      for(unsigned short k=0; k<nDim; ++k) {
        Local_Data[jPoint][iVar] = coor[k];
        iVar++;
      }

      /* Load the conservative variables. */

      for(jVar=0; jVar < nVar_First; ++jVar) {
        Local_Data[jPoint][iVar] = U[jVar];
        iVar++;
      }

      /*--- Prepare the primitive states. ---*/

      const su2double DensityInv = 1.0/U[0];
      su2double vel[3], Velocity2 = 0.0;
      for(iDim=0; iDim<nDim; ++iDim) {
        vel[iDim] = U[iDim+1]*DensityInv;
        Velocity2 += vel[iDim]*vel[iDim];
      }
      su2double StaticEnergy = U[nDim+1]*DensityInv - 0.5*Velocity2;
      DGFluidModel->SetTDState_rhoe(U[0], StaticEnergy);

      /*--- Load data for the pressure, temperature, Cp, and Mach variables. ---*/

      Local_Data[jPoint][iVar] = DGFluidModel->GetPressure(); iVar++;
      Local_Data[jPoint][iVar] = DGFluidModel->GetTemperature(); iVar++;
      Local_Data[jPoint][iVar] = DGFluidModel->GetCp(); iVar++;
      Local_Data[jPoint][iVar] = sqrt(Velocity2)/DGFluidModel->GetSoundSpeed(); iVar++;

      if ((Kind_Solver == FEM_NAVIER_STOKES) || (Kind_Solver == FEM_LES)){
        Local_Data[jPoint][iVar] = DGFluidModel->GetLaminarViscosity(); iVar++;
      }
      if ((Kind_Solver == FEM_LES) && (config->GetKind_SGS_Model() != IMPLICIT_LES)){
        // todo: Export Eddy instead of Laminar viscosity
        Local_Data[jPoint][iVar] = DGFluidModel->GetLaminarViscosity(); iVar++;
      }

      if (solver[FLOW_SOL]->VerificationSolution) {
        if (solver[FLOW_SOL]->VerificationSolution->ExactSolutionKnown()) {

          /*--- Get the physical time if necessary. ---*/
          su2double time = 0.0;
          if (config->GetTime_Marching()) time = config->GetPhysicalTime();

          /* Get the verification solution. */
          su2double mmsSol[5];
          solver[FLOW_SOL]->VerificationSolution->GetSolution(coor, time, mmsSol);
          for (jVar = 0; jVar < nVar_First; jVar++) {
            Local_Data[jPoint][iVar] = mmsSol[jVar];
            iVar++;
          }

          /* Get local error from the verification solution class. */
          su2double error[5];
          solver[FLOW_SOL]->VerificationSolution->GetLocalError(coor, time, U, error);
          for (jVar = 0; jVar < nVar_First; jVar++) {
            Local_Data[jPoint][iVar] = error[jVar];
            iVar++;
          }
        }
      }

      /*--- New variables can be loaded to the Local_Data structure here,
         assuming they were registered above correctly. ---*/

    }

    /*--- Increment the point counter. ---*/

    jPoint++;
  }


}

void COutputLegacy::PrepareOffsets(CConfig *config, CGeometry *geometry) {

  unsigned long iPoint;

  /*--- Bool to distinguish between the FVM and FEM solvers. ---*/

  const bool fem_solver = config->GetFEMSolver();

  /*--- Reset point sorting counters ---*/

  nGlobalPoint_Sort = 0;
  nLocalPoint_Sort  = 0;

  /*--- Prepare the offsets for the FV solver ---*/

  if (!fem_solver) {

    /*--- Search all send/recv boundaries on this partition for any periodic
     nodes that were part of the original domain. We want to recover these
     for visualization purposes. ---*/

    unsigned long iVertex;
    bool isPeriodic;

    Local_Halo_Sort = new int[geometry->GetnPoint()];
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
      Local_Halo_Sort[iPoint] = !geometry->node[iPoint]->GetDomain();

    for (unsigned short iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) {

        /*--- Checking for less than or equal to the rank, because there may
         be some periodic halo nodes that send info to the same rank. ---*/

        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          isPeriodic = ((geometry->vertex[iMarker][iVertex]->GetRotation_Type() > 0) &&
                        (geometry->vertex[iMarker][iVertex]->GetRotation_Type() % 2 == 1));
          if (isPeriodic) Local_Halo_Sort[iPoint] = false;
        }
      }
    }

    /*--- Sum total number of nodes that belong to the domain ---*/

    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
      if (Local_Halo_Sort[iPoint] == false)
        nLocalPoint_Sort++;

#ifdef HAVE_MPI
    SU2_MPI::Allreduce(&nLocalPoint_Sort, &nGlobalPoint_Sort, 1,
                       MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#else
    nGlobalPoint_Sort = nLocalPoint_Sort;
#endif

    /*--- Set a final variable for the number of points in the restart
     file. We do not write the periodic points for the FV solver, even if
     they show up in the viz. files. ---*/

    nPoint_Restart = geometry->GetGlobal_nPointDomain();

  } else {

    /*--- Create an object of the class CMeshFEM_DG and retrieve the necessary
     geometrical information for the FEM DG solver. ---*/

    CMeshFEM_DG *DGGeometry = dynamic_cast<CMeshFEM_DG *>(geometry);

    unsigned long nVolElemOwned = DGGeometry->GetNVolElemOwned();
    CVolumeElementFEM *volElem  = DGGeometry->GetVolElem();

    /*--- Update the solution by looping over the owned volume elements. ---*/

    for(unsigned long l=0; l<nVolElemOwned; ++l) {

      /* Count up the number of local points we have for allocating storage. */

      for(unsigned short j=0; j<volElem[l].nDOFsSol; ++j) {
        nLocalPoint_Sort++;
      }
    }

    Local_Halo_Sort = new int[nLocalPoint_Sort];
    for (iPoint = 0; iPoint < nLocalPoint_Sort; iPoint++)
      Local_Halo_Sort[iPoint] = false;

#ifdef HAVE_MPI
    SU2_MPI::Allreduce(&nLocalPoint_Sort, &nGlobalPoint_Sort, 1,
                       MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#else
    nGlobalPoint_Sort = nLocalPoint_Sort;
#endif

    /*--- Set a final variable for the number of points in the restart
     file. We do not write the periodic points for the FV solver, even if
     they show up in the viz. files. ---*/

    nPoint_Restart = nGlobalPoint_Sort;

  }

  /*--- Now that we know the actual number of points we need to output,
   compute the number of points that will be on each processor.
   This is a linear partitioning with the addition of a simple load
   balancing for any remainder points. ---*/

  beg_node = new unsigned long[size];
  end_node = new unsigned long[size];

  nPointLinear     = new unsigned long[size];
  nPointCumulative = new unsigned long[size+1];

  unsigned long total_points = 0;
  for (int ii = 0; ii < size; ii++) {
    nPointLinear[ii]   = nGlobalPoint_Sort/size;
    total_points += nPointLinear[ii];
  }

  /*--- Get the number of remainder points after the even division. ---*/

  unsigned long remainder = nGlobalPoint_Sort - total_points;
  for (unsigned long ii = 0; ii < remainder; ii++) {
    nPointLinear[ii]++;
  }

  /*--- Store the local number of nodes on each proc in the linear
   partitioning, the beginning/end index, and the linear partitioning
   within an array in cumulative storage format. ---*/

  beg_node[0] = 0;
  end_node[0] = beg_node[0] + nPointLinear[0];
  nPointCumulative[0] = 0;
  for (int ii = 1; ii < size; ii++) {
    beg_node[ii] = end_node[ii-1];
    end_node[ii] = beg_node[ii] + nPointLinear[ii];
    nPointCumulative[ii] = nPointCumulative[ii-1] + nPointLinear[ii-1];
  }
  nPointCumulative[size] = nGlobalPoint_Sort;

}

void COutputLegacy::SortConnectivity_FEM(CConfig *config, CGeometry *geometry, unsigned short val_iZone) {

  /*--- Flags identifying the types of files to be written. ---*/

  bool Wrt_Vol = config->GetWrt_Vol_Sol();
  bool Wrt_Srf = config->GetWrt_Srf_Sol();

  /*--- Sort connectivity for each type of element (excluding halos). Note
   In these routines, we sort the connectivity into a linear partitioning
   across all processors based on the global index of the grid nodes. ---*/

  /*--- Sort volumetric grid connectivity. ---*/

  if (Wrt_Vol) {

    if ((rank == MASTER_NODE) && (size != SINGLE_NODE))
      cout <<"Sorting volumetric grid connectivity." << endl;

    SortVolumetricConnectivity_FEM(config, geometry, TRIANGLE     );
    SortVolumetricConnectivity_FEM(config, geometry, QUADRILATERAL);
    SortVolumetricConnectivity_FEM(config, geometry, TETRAHEDRON  );
    SortVolumetricConnectivity_FEM(config, geometry, HEXAHEDRON   );
    SortVolumetricConnectivity_FEM(config, geometry, PRISM        );
    SortVolumetricConnectivity_FEM(config, geometry, PYRAMID      );

  }

  /*--- Sort surface grid connectivity. ---*/

  if (Wrt_Srf) {

    if ((rank == MASTER_NODE) && (size != SINGLE_NODE))
      cout <<"Sorting surface grid connectivity." << endl;

    SortSurfaceConnectivity_FEM(config, geometry, LINE         );
    SortSurfaceConnectivity_FEM(config, geometry, TRIANGLE     );
    SortSurfaceConnectivity_FEM(config, geometry, QUADRILATERAL);

  }

  /*--- Reduce the total number of cells we will be writing in the output files. ---*/

  unsigned long nTotal_Elem = nParallel_Tria + nParallel_Quad + nParallel_Tetr + nParallel_Hexa + nParallel_Pris + nParallel_Pyra;
  unsigned long nTotal_Surf_Elem = nParallel_Line + nParallel_BoundTria + nParallel_BoundQuad;
#ifndef HAVE_MPI
  nGlobal_Elem_Par = nTotal_Elem;
  nSurf_Elem_Par   = nTotal_Surf_Elem;
#else
  SU2_MPI::Allreduce(&nTotal_Elem, &nGlobal_Elem_Par, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&nTotal_Surf_Elem, &nSurf_Elem_Par, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#endif

}

void COutputLegacy::SortVolumetricConnectivity_FEM(CConfig *config, CGeometry *geometry, unsigned short Elem_Type) {

  /* Determine the number of nodes for this element type. */
  unsigned short NODES_PER_ELEMENT = 0;
  switch (Elem_Type) {
    case TRIANGLE:
      NODES_PER_ELEMENT = N_POINTS_TRIANGLE;
      break;
    case QUADRILATERAL:
      NODES_PER_ELEMENT = N_POINTS_QUADRILATERAL;
      break;
    case TETRAHEDRON:
      NODES_PER_ELEMENT = N_POINTS_TETRAHEDRON;
      break;
    case HEXAHEDRON:
      NODES_PER_ELEMENT = N_POINTS_HEXAHEDRON;
      break;
    case PRISM:
      NODES_PER_ELEMENT = N_POINTS_PRISM;
      break;
    case PYRAMID:
      NODES_PER_ELEMENT = N_POINTS_PYRAMID;
      break;
    default:
      SU2_MPI::Error("Unrecognized element type", CURRENT_FUNCTION);
  }

  /*--- Create an object of the class CMeshFEM_DG and retrieve the necessary
        geometrical information for the FEM DG solver. ---*/
  CMeshFEM_DG *DGGeometry = dynamic_cast<CMeshFEM_DG *>(geometry);

  unsigned long nVolElemOwned = DGGeometry->GetNVolElemOwned();
  CVolumeElementFEM *volElem  = DGGeometry->GetVolElem();

  const CFEMStandardElement *standardElementsSol = DGGeometry->GetStandardElementsSol();

  /*--- Determine the number of sub-elements on this rank. ---*/
  unsigned long nSubElem_Local = 0;
  for(unsigned long i=0; i<nVolElemOwned; ++i) {

    /* Determine the necessary data from the corresponding standard elem. */
    const unsigned short ind       = volElem[i].indStandardElement;
    const unsigned short VTK_Type1 = standardElementsSol[ind].GetVTK_Type1();
    const unsigned short VTK_Type2 = standardElementsSol[ind].GetVTK_Type2();

     /* Only store the linear sub elements if they are of
        the current type that we are storing. */
     if(Elem_Type == VTK_Type1) nSubElem_Local += standardElementsSol[ind].GetNSubElemsType1();
     if(Elem_Type == VTK_Type2) nSubElem_Local += standardElementsSol[ind].GetNSubElemsType2();
  }

  /* Allocate the memory to store the connectivity if the size is
     larger than zero. */
  int *Conn_SubElem = NULL;
  if(nSubElem_Local > 0) Conn_SubElem = new int[nSubElem_Local*NODES_PER_ELEMENT];

  /*--- Loop again over the local volume elements and store the global
        connectivities of the sub-elements. Note one is added to the
        index value, because visualization softwares typically use
        1-based indexing. ---*/
  unsigned long kNode = 0;
  for(unsigned long i=0; i<nVolElemOwned; ++i) {

    /* Determine the necessary data from the corresponding standard elem. */
    const unsigned short ind       = volElem[i].indStandardElement;
    const unsigned short VTK_Type1 = standardElementsSol[ind].GetVTK_Type1();
    const unsigned short VTK_Type2 = standardElementsSol[ind].GetVTK_Type2();

    /* Check if the first sub-element is of the required type. */
    if(Elem_Type == VTK_Type1) {

      /* Get the number of sub-elements and the local connectivity of
         the sub-elements. */
      const unsigned short nSubElems     = standardElementsSol[ind].GetNSubElemsType1();
      const unsigned short *connSubElems = standardElementsSol[ind].GetSubConnType1();

      /* Store the global connectivities. */
      const unsigned short kk = NODES_PER_ELEMENT*nSubElems;
      for(unsigned short k=0; k<kk; ++k, ++kNode)
        Conn_SubElem[kNode] = connSubElems[k] + volElem[i].offsetDOFsSolGlobal + 1;
    }

    /* Check if the second sub-element is of the required type. */
    if(Elem_Type == VTK_Type2) {

      /* Get the number of sub-elements and the local connectivity of
         the sub-elements. */
      const unsigned short nSubElems     = standardElementsSol[ind].GetNSubElemsType2();
      const unsigned short *connSubElems = standardElementsSol[ind].GetSubConnType2();

      /* Store the global connectivities. */
      const unsigned short kk = NODES_PER_ELEMENT*nSubElems;
      for(unsigned short k=0; k<kk; ++k, ++kNode)
        Conn_SubElem[kNode] = connSubElems[k] + volElem[i].offsetDOFsSolGlobal + 1;
    }
  }

  /*--- Store the particular global element count in the class data,
        and set the class data pointer to the connectivity array. ---*/
  switch (Elem_Type) {
    case TRIANGLE:
      nParallel_Tria = nSubElem_Local;
      Conn_Tria_Par = Conn_SubElem;
      break;
    case QUADRILATERAL:
      nParallel_Quad = nSubElem_Local;
      Conn_Quad_Par = Conn_SubElem;
      break;
    case TETRAHEDRON:
      nParallel_Tetr = nSubElem_Local;
      Conn_Tetr_Par = Conn_SubElem;
      break;
    case HEXAHEDRON:
      nParallel_Hexa = nSubElem_Local;
      Conn_Hexa_Par = Conn_SubElem;
      break;
    case PRISM:
      nParallel_Pris = nSubElem_Local;
      Conn_Pris_Par = Conn_SubElem;
      break;
    case PYRAMID:
      nParallel_Pyra = nSubElem_Local;
      Conn_Pyra_Par = Conn_SubElem;
      break;
    default:
      SU2_MPI::Error("Unrecognized element type", CURRENT_FUNCTION);
      break;
  }
}

void COutputLegacy::SortSurfaceConnectivity_FEM(CConfig *config, CGeometry *geometry, unsigned short Elem_Type) {

  /* Determine the number of nodes for this element type. */
  unsigned short NODES_PER_ELEMENT = 0;
  switch (Elem_Type) {
    case LINE:
      NODES_PER_ELEMENT = N_POINTS_LINE;
      break;
    case TRIANGLE:
      NODES_PER_ELEMENT = N_POINTS_TRIANGLE;
      break;
    case QUADRILATERAL:
      NODES_PER_ELEMENT = N_POINTS_QUADRILATERAL;
      break;
    default:
      SU2_MPI::Error("Unrecognized element type", CURRENT_FUNCTION);
  }

  /*--- Create an object of the class CMeshFEM_DG and retrieve the necessary
        geometrical information for the FEM DG solver. ---*/
  CMeshFEM_DG *DGGeometry = dynamic_cast<CMeshFEM_DG *>(geometry);

  unsigned long nVolElemOwned = DGGeometry->GetNVolElemOwned();
  CVolumeElementFEM *volElem  = DGGeometry->GetVolElem();

  const CBoundaryFEM *boundaries = DGGeometry->GetBoundaries();
  const CFEMStandardBoundaryFace *standardBoundaryFacesSol = DGGeometry->GetStandardBoundaryFacesSol();

  /*--- Create the map from the global DOF ID to the local index.
        Note one is added to the index value, because visualization
        softwares typically use 1-based indexing. ---*/
  vector<unsigned long> globalID;
  for(unsigned long l=0; l<nVolElemOwned; ++l) {
    for(unsigned short j=0; j<volElem[l].nDOFsSol; ++j) {
      const unsigned long globalIndex = volElem[l].offsetDOFsSolGlobal + j + 1;
      globalID.push_back(globalIndex);
    }
  }

  /*--- Determine the number of sub-elements on this rank by looping
        over the surface elements of the boundary markers that must
        be plotted. ---*/
  unsigned long nSubElem_Local = 0;
  for(unsigned short iMarker=0; iMarker < config->GetnMarker_All(); ++iMarker) {
    if( !boundaries[iMarker].periodicBoundary ) {
      if (config->GetMarker_All_Plotting(iMarker) == YES) {
        const vector<CSurfaceElementFEM> &surfElem = boundaries[iMarker].surfElem;
        for(unsigned long i=0; i<surfElem.size(); ++i) {
          const unsigned short ind      = surfElem[i].indStandardElement;
          const unsigned short VTK_Type = standardBoundaryFacesSol[ind].GetVTK_Type();
          if(Elem_Type == VTK_Type) nSubElem_Local += standardBoundaryFacesSol[ind].GetNSubFaces();
        }
      }
    }
  }

  /* Allocate the memory to store the connectivity if the size is
     larger than zero. */
  int *Conn_SubElem = NULL;
  if(nSubElem_Local > 0) Conn_SubElem = new int[nSubElem_Local*NODES_PER_ELEMENT];

  /*--- Repeat the loop over the surface elements of the boundary markers
        that must be plotted, but now store the connectivity. ---*/
  unsigned long kNode = 0;
  for(unsigned short iMarker=0; iMarker < config->GetnMarker_All(); ++iMarker) {
    if( !boundaries[iMarker].periodicBoundary ) {
      if (config->GetMarker_All_Plotting(iMarker) == YES) {
        const vector<CSurfaceElementFEM> &surfElem = boundaries[iMarker].surfElem;

        /* Loop over the surface elements of this boundary marker. */
        for(unsigned long i=0; i<surfElem.size(); ++i) {

          /* Check if this is the element type to be stored. */
          const unsigned short ind      = surfElem[i].indStandardElement;
          const unsigned short VTK_Type = standardBoundaryFacesSol[ind].GetVTK_Type();
          if(Elem_Type == VTK_Type) {

            /* Get the number of sub-elements and the local connectivity of
               the sub-elements. */
            const unsigned short nSubFaces     = standardBoundaryFacesSol[ind].GetNSubFaces();
            const unsigned short *connSubFaces = standardBoundaryFacesSol[ind].GetSubFaceConn();

            /* Store the global connectivities. */
            const unsigned short kk = NODES_PER_ELEMENT*nSubFaces;
            for(unsigned short k=0; k<kk; ++k, ++kNode)
              Conn_SubElem[kNode] = globalID[surfElem[i].DOFsSolFace[connSubFaces[k]]];
          }
        }
      }
    }
  }

  /*--- Store the particular global element count in the class data,
        and set the class data pointer to the connectivity array. ---*/
  switch (Elem_Type) {
    case LINE:
      nParallel_Line = nSubElem_Local;
      Conn_BoundLine_Par = Conn_SubElem;
      break;
    case TRIANGLE:
      nParallel_BoundTria = nSubElem_Local;
      Conn_BoundTria_Par = Conn_SubElem;
      break;
    case QUADRILATERAL:
      nParallel_BoundQuad = nSubElem_Local;
      Conn_BoundQuad_Par = Conn_SubElem;
      break;
    default:
      SU2_MPI::Error("Unrecognized element type", CURRENT_FUNCTION);
      break;
  }
}

void COutputLegacy::SortOutputData_FEM(CConfig *config, CGeometry *geometry) {

  unsigned long iProcessor;
  unsigned long iPoint, Global_Index;

  /* For convenience, set the total number of variables stored at each DOF. */

  int VARS_PER_POINT = nVar_Par;

#ifdef HAVE_MPI
  SU2_MPI::Request *send_req, *recv_req;
  SU2_MPI::Status status;
  int ind;
#endif

  /*--- Create an object of the class CMeshFEM_DG and retrieve the necessary
   geometrical information for the FEM DG solver. ---*/

  CMeshFEM_DG *DGGeometry = dynamic_cast<CMeshFEM_DG *>(geometry);

  unsigned long nVolElemOwned = DGGeometry->GetNVolElemOwned();
  CVolumeElementFEM *volElem  = DGGeometry->GetVolElem();

  /*--- Create the map from the global DOF ID to the local index. ---*/

  //map<unsigned long, unsigned long> mapLocal2Global;
  vector<unsigned long> globalID;

  /*--- Update the solution by looping over the owned volume elements. ---*/
  for(unsigned long l=0; l<nVolElemOwned; ++l) {

    /* Loop over the DOFs for this element and store the solution. */
    for(unsigned short j=0; j<volElem[l].nDOFsSol; ++j) {
      const unsigned long globalIndex = volElem[l].offsetDOFsSolGlobal + j;
      globalID.push_back(globalIndex);
    }
  }

  /*--- We start with the grid nodes distributed across all procs with
   no particular ordering assumed. We need to loop through our local partition
   and decide how many nodes we must send to each other rank in order to
   have all nodes sorted according to a linear partitioning of the grid
   nodes, i.e., rank 0 holds the first ~ nGlobalPoint()/nProcessors nodes.
   First, initialize a counter and flag. ---*/

  int *nPoint_Send = new int[size+1]; nPoint_Send[0] = 0;
  int *nPoint_Recv = new int[size+1]; nPoint_Recv[0] = 0;
  int *nPoint_Flag = new int[size];

  for (int ii=0; ii < size; ii++) {
    nPoint_Send[ii] = 0;
    nPoint_Recv[ii] = 0;
    nPoint_Flag[ii]= -1;
  }
  nPoint_Send[size] = 0; nPoint_Recv[size] = 0;

  for (iPoint = 0; iPoint < nLocalPoint_Sort; iPoint++ ) {

    /*--- Get the global index of the current point. ---*/

    Global_Index = globalID[iPoint];

    /*--- Search for the processor that owns this point ---*/

    iProcessor = Global_Index/nPointLinear[0];
    if (iProcessor >= (unsigned long)size)
      iProcessor = (unsigned long)size-1;
    if (Global_Index >= nPointCumulative[iProcessor])
      while(Global_Index >= nPointCumulative[iProcessor+1]) iProcessor++;
    else
      while(Global_Index <  nPointCumulative[iProcessor])   iProcessor--;

    /*--- If we have not visted this node yet, increment our
     number of elements that must be sent to a particular proc. ---*/

    if (nPoint_Flag[iProcessor] != (int)iPoint) {
      nPoint_Flag[iProcessor] = (int)iPoint;
      nPoint_Send[iProcessor+1]++;
    }

  }

  /*--- Communicate the number of nodes to be sent/recv'd amongst
   all processors. After this communication, each proc knows how
   many cells it will receive from each other processor. ---*/

#ifdef HAVE_MPI
  SU2_MPI::Alltoall(&(nPoint_Send[1]), 1, MPI_INT,
                    &(nPoint_Recv[1]), 1, MPI_INT, MPI_COMM_WORLD);
#else
  nPoint_Recv[1] = nPoint_Send[1];
#endif

  /*--- Prepare to send coordinates. First check how many
   messages we will be sending and receiving. Here we also put
   the counters into cumulative storage format to make the
   communications simpler. ---*/

  int nSends = 0, nRecvs = 0;
  for (int ii=0; ii < size; ii++) nPoint_Flag[ii] = -1;

  for (int ii = 0; ii < size; ii++) {
    if ((ii != rank) && (nPoint_Send[ii+1] > 0)) nSends++;
    if ((ii != rank) && (nPoint_Recv[ii+1] > 0)) nRecvs++;

    nPoint_Send[ii+1] += nPoint_Send[ii];
    nPoint_Recv[ii+1] += nPoint_Recv[ii];
  }

  /*--- Allocate memory to hold the connectivity that we are sending. ---*/

  su2double *connSend = NULL;
  connSend = new su2double[VARS_PER_POINT*nPoint_Send[size]];
  for (int ii = 0; ii < VARS_PER_POINT*nPoint_Send[size]; ii++)
    connSend[ii] = 0;

  /*--- Allocate arrays for sending the global ID. ---*/

  unsigned long *idSend = new unsigned long[nPoint_Send[size]];
  for (int ii = 0; ii < nPoint_Send[size]; ii++)
    idSend[ii] = 0;

  /*--- Create an index variable to keep track of our index
   positions as we load up the send buffer. ---*/

  unsigned long *index = new unsigned long[size];
  for (int ii=0; ii < size; ii++) index[ii] = VARS_PER_POINT*nPoint_Send[ii];

  unsigned long *idIndex = new unsigned long[size];
  for (int ii=0; ii < size; ii++) idIndex[ii] = nPoint_Send[ii];

  /*--- Loop through our elements and load the elems and their
   additional data that we will send to the other procs. ---*/

  for (iPoint = 0; iPoint < nLocalPoint_Sort; iPoint++) {

    /*--- Get the index of the current point. ---*/

    Global_Index = globalID[iPoint];

    /*--- Search for the processor that owns this point. ---*/

    iProcessor = Global_Index/nPointLinear[0];
    if (iProcessor >= (unsigned long)size)
      iProcessor = (unsigned long)size-1;
    if (Global_Index >= nPointCumulative[iProcessor])
      while(Global_Index >= nPointCumulative[iProcessor+1]) iProcessor++;
    else
      while(Global_Index <  nPointCumulative[iProcessor])   iProcessor--;

    /*--- Load data into the buffer for sending. ---*/

    if (nPoint_Flag[iProcessor] != (int)iPoint) {

      nPoint_Flag[iProcessor] = (int)iPoint;
      unsigned long nn = index[iProcessor];

      /*--- Load the data values. ---*/

      for (unsigned short kk = 0; kk < VARS_PER_POINT; kk++) {
        connSend[nn] = Local_Data[iPoint][kk]; nn++;
      }

      /*--- Load the global ID (minus offset) for sorting the
       points once they all reach the correct processor. ---*/

      nn = idIndex[iProcessor];
      idSend[nn] = Global_Index - beg_node[iProcessor];

      /*--- Increment the index by the message length ---*/

      index[iProcessor]  += VARS_PER_POINT;
      idIndex[iProcessor]++;

    }

  }

  /*--- Free memory after loading up the send buffer. ---*/

  delete [] index;
  delete [] idIndex;

  /*--- Allocate the memory that we need for receiving the conn
   values and then cue up the non-blocking receives. Note that
   we do not include our own rank in the communications. We will
   directly copy our own data later. ---*/

  su2double *connRecv = NULL;
  connRecv = new su2double[VARS_PER_POINT*nPoint_Recv[size]];
  for (int ii = 0; ii < VARS_PER_POINT*nPoint_Recv[size]; ii++)
    connRecv[ii] = 0;

  unsigned long *idRecv = new unsigned long[nPoint_Recv[size]];
  for (int ii = 0; ii < nPoint_Recv[size]; ii++)
    idRecv[ii] = 0;

#ifdef HAVE_MPI
  /*--- We need double the number of messages to send both the conn.
   and the global IDs. ---*/

  send_req = new SU2_MPI::Request[2*nSends];
  recv_req = new SU2_MPI::Request[2*nRecvs];

  unsigned long iMessage = 0;
  for (int ii=0; ii<size; ii++) {
    if ((ii != rank) && (nPoint_Recv[ii+1] > nPoint_Recv[ii])) {
      int ll     = VARS_PER_POINT*nPoint_Recv[ii];
      int kk     = nPoint_Recv[ii+1] - nPoint_Recv[ii];
      int count  = VARS_PER_POINT*kk;
      int source = ii;
      int tag    = ii + 1;
      SU2_MPI::Irecv(&(connRecv[ll]), count, MPI_DOUBLE, source, tag,
                     MPI_COMM_WORLD, &(recv_req[iMessage]));
      iMessage++;
    }
  }

  /*--- Launch the non-blocking sends of the connectivity. ---*/

  iMessage = 0;
  for (int ii=0; ii<size; ii++) {
    if ((ii != rank) && (nPoint_Send[ii+1] > nPoint_Send[ii])) {
      int ll = VARS_PER_POINT*nPoint_Send[ii];
      int kk = nPoint_Send[ii+1] - nPoint_Send[ii];
      int count  = VARS_PER_POINT*kk;
      int dest = ii;
      int tag    = rank + 1;
      SU2_MPI::Isend(&(connSend[ll]), count, MPI_DOUBLE, dest, tag,
                     MPI_COMM_WORLD, &(send_req[iMessage]));
      iMessage++;
    }
  }

  /*--- Repeat the process to communicate the global IDs. ---*/

  iMessage = 0;
  for (int ii=0; ii<size; ii++) {
    if ((ii != rank) && (nPoint_Recv[ii+1] > nPoint_Recv[ii])) {
      int ll     = nPoint_Recv[ii];
      int kk     = nPoint_Recv[ii+1] - nPoint_Recv[ii];
      int count  = kk;
      int source = ii;
      int tag    = ii + 1;
      SU2_MPI::Irecv(&(idRecv[ll]), count, MPI_UNSIGNED_LONG, source, tag,
                     MPI_COMM_WORLD, &(recv_req[iMessage+nRecvs]));
      iMessage++;
    }
  }

  /*--- Launch the non-blocking sends of the global IDs. ---*/

  iMessage = 0;
  for (int ii=0; ii<size; ii++) {
    if ((ii != rank) && (nPoint_Send[ii+1] > nPoint_Send[ii])) {
      int ll = nPoint_Send[ii];
      int kk = nPoint_Send[ii+1] - nPoint_Send[ii];
      int count  = kk;
      int dest   = ii;
      int tag    = rank + 1;
      SU2_MPI::Isend(&(idSend[ll]), count, MPI_UNSIGNED_LONG, dest, tag,
                     MPI_COMM_WORLD, &(send_req[iMessage+nSends]));
      iMessage++;
    }
  }
#endif

  /*--- Copy my own rank's data into the recv buffer directly. ---*/

  int mm = VARS_PER_POINT*nPoint_Recv[rank];
  int ll = VARS_PER_POINT*nPoint_Send[rank];
  int kk = VARS_PER_POINT*nPoint_Send[rank+1];

  for (int nn=ll; nn<kk; nn++, mm++) connRecv[mm] = connSend[nn];

  mm = nPoint_Recv[rank];
  ll = nPoint_Send[rank];
  kk = nPoint_Send[rank+1];

  for (int nn=ll; nn<kk; nn++, mm++) idRecv[mm] = idSend[nn];

  /*--- Wait for the non-blocking sends and recvs to complete. ---*/

#ifdef HAVE_MPI
  int number = 2*nSends;
  for (int ii = 0; ii < number; ii++)
    SU2_MPI::Waitany(number, send_req, &ind, &status);

  number = 2*nRecvs;
  for (int ii = 0; ii < number; ii++)
    SU2_MPI::Waitany(number, recv_req, &ind, &status);

  delete [] send_req;
  delete [] recv_req;
#endif

  /*--- Store the connectivity for this rank in the proper data
   structure before post-processing below. First, allocate the
   appropriate amount of memory for this section. ---*/

  Parallel_Data = new su2double*[VARS_PER_POINT];
  for (int jj = 0; jj < VARS_PER_POINT; jj++) {
    Parallel_Data[jj] = new su2double[nPoint_Recv[size]];
    for (int ii = 0; ii < nPoint_Recv[size]; ii++) {
      Parallel_Data[jj][idRecv[ii]] = connRecv[ii*VARS_PER_POINT+jj];
    }
  }

  /*--- Store the total number of local points my rank has for
   the current section after completing the communications. ---*/

  nParallel_Poin = nPoint_Recv[size];

  /*--- Reduce the total number of points we will write in the output files. ---*/

#ifndef HAVE_MPI
  nGlobal_Poin_Par = nParallel_Poin;
#else
  SU2_MPI::Allreduce(&nParallel_Poin, &nGlobal_Poin_Par, 1,
                     MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#endif

  /*--- Free temporary memory from communications ---*/

  delete [] connSend;
  delete [] connRecv;
  delete [] idSend;
  delete [] idRecv;
  delete [] nPoint_Recv;
  delete [] nPoint_Send;
  delete [] nPoint_Flag;

  for (iPoint = 0; iPoint < nLocalPoint_Sort; iPoint++)
    delete [] Local_Data[iPoint];
  delete [] Local_Data;

}

void COutputLegacy::SortOutputData_Surface_FEM(CConfig *config, CGeometry *geometry) {

  const int VARS_PER_POINT = nVar_Par;

  /*---------------------------------------------------*/
  /*--- Step 1: Determine the global DOF ID's of the   */
  /*---         locally stored surface connectivities. */
  /*---------------------------------------------------*/

  /* Loop over the surface connectivities and store global
     DOF ID's in a vector. Subtract 1, because the stored
     connectivities are 1 based. */
  std::vector<unsigned long> globalSurfaceDOFIDs;
  globalSurfaceDOFIDs.reserve(nParallel_Line*N_POINTS_LINE +
                              nParallel_BoundTria*N_POINTS_TRIANGLE +
                              nParallel_BoundQuad*N_POINTS_QUADRILATERAL);

  for(unsigned long i=0; i<(nParallel_Line*N_POINTS_LINE); ++i) {
    const unsigned long globalID = Conn_BoundLine_Par[i]-1;
    globalSurfaceDOFIDs.push_back(globalID);
  }

  for(unsigned long i=0; i<(nParallel_BoundTria*N_POINTS_TRIANGLE); ++i) {
    const unsigned long globalID = Conn_BoundTria_Par[i]-1;
    globalSurfaceDOFIDs.push_back(globalID);
  }

  for(unsigned long i=0; i<(nParallel_BoundQuad*N_POINTS_QUADRILATERAL); ++i) {
    const unsigned long globalID = Conn_BoundQuad_Par[i]-1;
    globalSurfaceDOFIDs.push_back(globalID);
  }

  /* Sort globalSurfaceDOFIDs in increasing order and remove the
     multiple entries. */
  sort(globalSurfaceDOFIDs.begin(), globalSurfaceDOFIDs.end());
  vector<unsigned long>::iterator lastEntry;
  lastEntry = unique(globalSurfaceDOFIDs.begin(), globalSurfaceDOFIDs.end());
  globalSurfaceDOFIDs.erase(lastEntry, globalSurfaceDOFIDs.end());

  /*-------------------------------------------------------------------*/
  /*--- Step 2: Communicate the information of globalSurfaceDOFIDs  ---*/
  /*---         to the ranks, which actually store this information ---*/
  /*---         in Parallel_Data.                                   ---*/
  /*-------------------------------------------------------------------*/

  /* Allocate the memory for the first index of the communication buffers. */
  vector<vector<unsigned long> > sendBuf(size, vector<unsigned long>(0));

  /* Loop over the entries of globalSurfaceDOFIDs and fill the
     communication buffers accordingly. */
  for(unsigned long i=0; i<globalSurfaceDOFIDs.size(); ++i) {

    /* Search for the processor that owns this point. */
    unsigned long iProcessor = globalSurfaceDOFIDs[i]/nPointLinear[0];
    if (iProcessor >= (unsigned long)size)
      iProcessor = (unsigned long)size-1;
    if (globalSurfaceDOFIDs[i] >= nPointCumulative[iProcessor])
      while(globalSurfaceDOFIDs[i] >= nPointCumulative[iProcessor+1]) ++iProcessor;
    else
      while(globalSurfaceDOFIDs[i] <  nPointCumulative[iProcessor])   --iProcessor;

    /* Store the global ID in the send buffer for iProcessor. */
    sendBuf[iProcessor].push_back(globalSurfaceDOFIDs[i]);
  }

  /* Determine the number of DOFs to be sent to each processor. */
  int nRankSend = 0;
  vector<unsigned long> nDOFSend(size);
  for(int i=0; i<size; ++i) {
    nDOFSend[i] = sendBuf[i].size();
    if(nDOFSend[i] && (i != rank)) ++nRankSend;
  }

  /* Communicate nDOFSend using Alltoall. */
  vector<unsigned long> nDOFRecv(size);

#ifdef HAVE_MPI
  SU2_MPI::Alltoall(nDOFSend.data(), 1, MPI_UNSIGNED_LONG,
                    nDOFRecv.data(), 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
#else
  nDOFRecv[rank] = nDOFSend[rank];
#endif

  /* Determine the number of messages this rank will receive. */
  int nRankRecv = 0;
  for(int i=0; i<size; ++i) {
    if(nDOFRecv[i] && (i != rank)) ++nRankRecv;
  }

  /* Define the receive buffers. */
  vector<vector<unsigned long> > recvBuf(size, vector<unsigned long>(0));

#ifdef HAVE_MPI
  /* Launch the non-blocking sends. Do not send anything to myself. */
  vector<SU2_MPI::Request> sendReq(nRankSend);
  nRankSend = 0;
  for(int i=0; i<size; ++i) {
    if(nDOFSend[i] && (i != rank)) {
      SU2_MPI::Isend(sendBuf[i].data(), nDOFSend[i], MPI_UNSIGNED_LONG,
                     i, rank, MPI_COMM_WORLD, &sendReq[nRankSend]);
      ++nRankSend;
    }
  }

  /* Launch the non-blocking receives. Do not receive anything from myself. */
  vector<SU2_MPI::Request> recvReq(nRankRecv);
  nRankRecv = 0;
  for(int i=0; i<size; ++i) {
    if(nDOFRecv[i] && (i != rank)) {
      recvBuf[i].resize(nDOFRecv[i]);
      SU2_MPI::Irecv(recvBuf[i].data(), nDOFRecv[i], MPI_UNSIGNED_LONG,
                     i, i, MPI_COMM_WORLD, &recvReq[nRankRecv]);
      ++nRankRecv;
    }
  }

#endif

  /* Copy my own data. */
  recvBuf[rank] = sendBuf[rank];

#ifdef HAVE_MPI
  /* Complete the non-blocking communication. */
  SU2_MPI::Waitall(nRankSend, sendReq.data(), MPI_STATUSES_IGNORE);
  SU2_MPI::Waitall(nRankRecv, recvReq.data(), MPI_STATUSES_IGNORE);
#endif

  /*-------------------------------------------------------------------*/
  /*--- Step 2: Determine the number of locally stored surface DOFs ---*/
  /*---         and store the data in Parallel_Surf_Data.           ---*/
  /*-------------------------------------------------------------------*/

  /* Store all received surface ID in one vector and sort it. */
  globalSurfaceDOFIDs.resize(0);
  for(int i=0; i<size; ++i) {
    if( nDOFRecv[i] )
      globalSurfaceDOFIDs.insert(globalSurfaceDOFIDs.end(),
                                  recvBuf[i].begin(), recvBuf[i].end());
  }

  sort(globalSurfaceDOFIDs.begin(), globalSurfaceDOFIDs.end());
  lastEntry = unique(globalSurfaceDOFIDs.begin(), globalSurfaceDOFIDs.end());
  globalSurfaceDOFIDs.erase(lastEntry, globalSurfaceDOFIDs.end());

  /* Allocate the memory for Parallel_Surf_Data. */
  nSurf_Poin_Par = globalSurfaceDOFIDs.size();

  Parallel_Surf_Data = new su2double*[VARS_PER_POINT];
  for(int jj=0; jj<VARS_PER_POINT; jj++)
    Parallel_Surf_Data[jj] = new su2double[nSurf_Poin_Par];

  /* Determine the local index of the global surface DOFs and
     copy the data into Parallel_Surf_Data. */
  for(unsigned long i=0; i<nSurf_Poin_Par; ++i) {
    const unsigned long ii = globalSurfaceDOFIDs[i] - nPointCumulative[rank];

    for(int jj=0; jj<VARS_PER_POINT; jj++)
      Parallel_Surf_Data[jj][i] = Parallel_Data[jj][ii];
  }

  /*--- Reduce the total number of surf points we have. This will be
        needed for writing the surface solution files later. ---*/
#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&nSurf_Poin_Par, &nGlobal_Surf_Poin, 1,
                     MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#else
  nGlobal_Surf_Poin = nSurf_Poin_Par;
#endif

  /*-------------------------------------------------------------------*/
  /*--- Step 3: Modify the surface connectivities, such that only   ---*/
  /*---         the data of surface DOFs needs to be written.       ---*/
  /*-------------------------------------------------------------------*/

  /* Determine the offset for my surface DOFs. */
  unsigned long offsetSurfaceDOFs = 0;
#ifdef HAVE_MPI
  vector<unsigned long> nSurfaceDOFsRanks(size);

  SU2_MPI::Allgather(&nSurf_Poin_Par, 1, MPI_UNSIGNED_LONG,
                     nSurfaceDOFsRanks.data(), 1, MPI_UNSIGNED_LONG,
                     MPI_COMM_WORLD);

  for(int i=0; i<rank; ++i) offsetSurfaceDOFs += nSurfaceDOFsRanks[i];
#endif

  /* Create the map from the global volume numbering to the new global
     surface numbering. */
  map<unsigned long, unsigned long> mapGlobalVol2Surf;
  for(unsigned long i=0; i<nSurf_Poin_Par; ++i)
    mapGlobalVol2Surf[globalSurfaceDOFIDs[i]] = offsetSurfaceDOFs + i;

  /* Fill the receive buffers with the modified global surface DOF numbers,
     such that this information can be returned to the original ranks. */
  for(int i=0; i<size; ++i) {
    for(unsigned long j=0; j<nDOFRecv[i]; ++j)
      recvBuf[i][j] = mapGlobalVol2Surf.find(recvBuf[i][j])->second;
  }

  /* Copy the original send buffers, because that information is
     still needed. */
  vector<vector<unsigned long> > originalSendBuf = sendBuf;

#ifdef HAVE_MPI
  /* Launch the non-blocking sends for the reverse communication, where the
     receive buffers must be used for sending. Do not send anything to myself. */
  nRankRecv = 0;
  for(int i=0; i<size; ++i) {
    if(nDOFRecv[i] && (i != rank)) {
      SU2_MPI::Isend(recvBuf[i].data(), nDOFRecv[i], MPI_UNSIGNED_LONG,
                     i, rank+1, MPI_COMM_WORLD, &recvReq[nRankRecv]);
      ++nRankRecv;
    }
  }

  /* Launch the non-blocking receives for the reverse communication, where the
     send buffers must be used for receiving.  Do not receive anything from myself. */
  nRankSend = 0;
  for(int i=0; i<size; ++i) {
    if(nDOFSend[i] && (i != rank)) {
      SU2_MPI::Irecv(sendBuf[i].data(), nDOFSend[i], MPI_UNSIGNED_LONG,
                     i, i+1, MPI_COMM_WORLD, &sendReq[nRankSend]);
      ++nRankSend;
    }
  }

#endif

  /* Copy my own data. */
  sendBuf[rank] = recvBuf[rank];

#ifdef HAVE_MPI
  /* Complete the non-blocking communication. */
  SU2_MPI::Waitall(nRankRecv, recvReq.data(), MPI_STATUSES_IGNORE);
  SU2_MPI::Waitall(nRankSend, sendReq.data(), MPI_STATUSES_IGNORE);
#endif

  /* Clear the map mapGlobalVol2Surf and fill it with the data
     needed to modify the surface connectivity. Note that 1 is added,
     because the visualization softwares typically use 1 based indexing. */
  mapGlobalVol2Surf.clear();

  for(int i=0; i<size; ++i) {
    for(unsigned long j=0; j<nDOFSend[i]; ++j)
      mapGlobalVol2Surf[originalSendBuf[i][j]+1] = sendBuf[i][j]+1;
  }

  /* Modify the locally stored surface connectivities. */
  for(unsigned long i=0; i<(nParallel_Line*N_POINTS_LINE); ++i)
    Conn_BoundLine_Par[i] = mapGlobalVol2Surf.find(Conn_BoundLine_Par[i])->second;

  for(unsigned long i=0; i<(nParallel_BoundTria*N_POINTS_TRIANGLE); ++i)
    Conn_BoundTria_Par[i] = mapGlobalVol2Surf.find(Conn_BoundTria_Par[i])->second;

  for(unsigned long i=0; i<(nParallel_BoundQuad*N_POINTS_QUADRILATERAL); ++i)
    Conn_BoundQuad_Par[i] = mapGlobalVol2Surf.find(Conn_BoundQuad_Par[i])->second;
}
