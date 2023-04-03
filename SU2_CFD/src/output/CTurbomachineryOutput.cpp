/*!
 * \file CNEMOCompOutput.cpp
 * \brief Main subroutines for turbomachinery output
 * \author F. Palacios, T. Economon, J. Kelly
 * \version 7.5.1 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2023, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../include/output/CTurbomachineryOutput.hpp"

#include "../../../Common/include/geometry/CGeometry.hpp"
#include "../../include/solvers/CBaselineSolver.hpp"
#include "../../include/fluid/CCoolProp.hpp"

CTurbomachineryOutput::CTurbomachineryOutput(CConfig *config) {

    rank = SU2_MPI::GetRank();
    size = SU2_MPI::GetSize();

    unsigned short iDim, iSpan, iMarker;

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

CTurbomachineryOutput::~CTurbomachineryOutput(void) {
    /* delete pointers initialized at construction*/
  /* Coords and Conn_*(Connectivity) have their own dealloc functions */
  /* Data is taken care of in DeallocateSolution function */

  /*--- Delete turboperformance pointers initiliazed at constrction  ---*/
  unsigned short iMarker, iSpan;
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