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

void CTurbomachineryOutput::ComputeTurboPerformance(CSolver *solver_container, CGeometry *geometry, CConfig *config) {
  /* I'm not sure this belongs here but it condenses turbo perf calculations into one file for now'*/

  CFluidModel *FluidModel;
  unsigned short nDim = geometry->GetnDim();
  unsigned short iMarkerTP, iSpan, iDim, iStage, iBlade;
  unsigned short nMarkerTP = config->GetnMarker_Turbomachinery();
  FluidModel = solver_container->GetFluidModel();
  su2double area, absVel2, soundSpeed, mach, tangVel, tangVel2, *relVel, relVel2;
  su2double relPressureIn, relPressureOut, enthalpyOutIs, relVelOutIs2;
  relVel = new su2double[nDim];
  su2double muLam, kine, omega, nu;
  bool turbulent = ((config->GetKind_Solver() == MAIN_SOLVER::RANS) || (config->GetKind_Solver() == MAIN_SOLVER::DISC_ADJ_RANS));
  bool menter_sst       = (config->GetKind_Turb_Model() == TURB_MODEL::SST);

  unsigned short nBladesRow, nStages;

  nBladesRow = config->GetnMarker_Turbomachinery();
  nStages    = SU2_TYPE::Int(nBladesRow/2);


  /*--- Compute BC imposed value for convergence monitoring ---*/
  for(iMarkerTP = 0; iMarkerTP < nMarkerTP; iMarkerTP++ ){
    for(iSpan = 0; iSpan < config->GetnSpan_iZones(iMarkerTP) + 1; iSpan++){
      if(config->GetRampOutletPressure() && config->GetInnerIter() > 0){
        PressureOut_BC[iMarkerTP][iSpan] = config->GetMonitorOutletPressure()/config->GetPressure_Ref();
      }
      FluidModel->SetTDState_PT(config->GetTotalPressureIn_BC(), config->GetTotalTemperatureIn_BC());
      TotalEnthalpyIn_BC[iMarkerTP][iSpan] = FluidModel->GetStaticEnergy()+ FluidModel->GetPressure()/FluidModel->GetDensity();
      EntropyIn_BC[iMarkerTP][iSpan]       = FluidModel->GetEntropy();
    }
  }

  /*--- Compute performance for each blade ---*/
  for(iMarkerTP = 0; iMarkerTP < nMarkerTP; iMarkerTP++ ){
    for(iSpan = 0; iSpan < config->GetnSpan_iZones(iMarkerTP) + 1; iSpan++){


      /*--- INFLOW ---*/
      /*--- Retrieve Inflow primitive quantities ---*/
      DensityIn[iMarkerTP][iSpan]          = solver_container->GetDensityIn(iMarkerTP, iSpan);
      PressureIn[iMarkerTP][iSpan]         = solver_container->GetPressureIn(iMarkerTP, iSpan);

      absVel2 = 0.0;

      for (iDim = 0; iDim < nDim; iDim++){
        TurboVelocityIn[iMarkerTP][iSpan][iDim]    = solver_container->GetTurboVelocityIn(iMarkerTP, iSpan)[iDim];
        absVel2   += TurboVelocityIn[iMarkerTP][iSpan][iDim]*TurboVelocityIn[iMarkerTP][iSpan][iDim];
      }
      TurboVelocityIn[iMarkerTP][iSpan][nDim] = sqrt(absVel2);

      TRadius[iMarkerTP][iSpan]  = geometry->GetTurboRadiusIn(iMarkerTP, iSpan);
      area                       = geometry->GetSpanAreaIn(iMarkerTP, iSpan);

      /*--- Compute static Inflow quantities ---*/
      FluidModel->SetTDState_Prho(PressureIn[iMarkerTP][iSpan], DensityIn[iMarkerTP][iSpan]);
      EntropyIn[iMarkerTP][iSpan]          = FluidModel->GetEntropy();
      MassFlowIn[iMarkerTP][iSpan]         = config->GetnBlades(iMarkerTP)*DensityIn[iMarkerTP][iSpan]*TurboVelocityIn[iMarkerTP][iSpan][0]*area;
      AbsFlowAngleIn[iMarkerTP][iSpan]     = atan(TurboVelocityIn[iMarkerTP][iSpan][1]/TurboVelocityIn[iMarkerTP][iSpan][0]);
      EnthalpyIn[iMarkerTP][iSpan]         = FluidModel->GetStaticEnergy() + PressureIn[iMarkerTP][iSpan]/DensityIn[iMarkerTP][iSpan];
      soundSpeed                           = FluidModel->GetSoundSpeed();


      /*--- Compute Total Inflow quantities ---*/
      TotalEnthalpyIn[iMarkerTP][iSpan]    = EnthalpyIn[iMarkerTP][iSpan] + 0.5*absVel2;
      FluidModel->SetTDState_hs(TotalEnthalpyIn[iMarkerTP][iSpan], EntropyIn[iMarkerTP][iSpan]);
      TotalPressureIn[iMarkerTP][iSpan]    = FluidModel->GetPressure();
      TotalTemperatureIn[iMarkerTP][iSpan] = FluidModel->GetTemperature();

      /*--- Retrieve Inflow relative quantities ---*/
      tangVel = geometry->GetTangGridVelIn(iMarkerTP, iSpan);
      tangVel2 = tangVel*tangVel;

      for (iDim = 0; iDim < nDim; iDim++){
        relVel[iDim] = TurboVelocityIn[iMarkerTP][iSpan][iDim];
      }
      relVel[1] -= tangVel;

      relVel2 = 0.0;
      for (iDim = 0; iDim < nDim; iDim++){
        relVel2 += relVel[iDim]*relVel[iDim];
      }

      /*--- Compute Total relative Inflow quantities ---*/
      RothalpyIn[iMarkerTP][iSpan]  = EnthalpyIn[iMarkerTP][iSpan] + 0.5*relVel2 - 0.5*tangVel2;
      FluidModel->SetTDState_hs(RothalpyIn[iMarkerTP][iSpan], EntropyIn[iMarkerTP][iSpan]);
      relPressureIn   = FluidModel->GetPressure();

      /*--- Compute kinematic relative Inflow quantities ---*/
      FlowAngleIn[iMarkerTP][iSpan]    = atan(relVel[1]/relVel[0]);
      mach          = 0.0;
      for (iDim = 0; iDim < nDim; iDim++){
        MachIn[iMarkerTP][iSpan][iDim] = relVel[iDim]/soundSpeed;
        mach = MachIn[iMarkerTP][iSpan][iDim]*MachIn[iMarkerTP][iSpan][iDim];
      }
      MachIn[iMarkerTP][iSpan][nDim]   = sqrt(mach);

      /*--- Compute Turbulent Inflow quantities ---*/
      if(turbulent){
        FluidModel->SetTDState_Prho(PressureIn[iMarkerTP][iSpan], DensityIn[iMarkerTP][iSpan]);
        muLam  = FluidModel->GetLaminarViscosity();
        if(menter_sst){
          kine   = solver_container->GetKineIn(iMarkerTP, iSpan);
          omega  = solver_container->GetOmegaIn(iMarkerTP, iSpan);
          TurbIntensityIn[iMarkerTP][iSpan]     =  sqrt(2.0/3.0*kine/absVel2);
          Turb2LamViscRatioIn[iMarkerTP][iSpan] = DensityIn[iMarkerTP][iSpan]*kine/(muLam*omega);
//          TurbIntensityIn[iMarkerTP][iSpan]     =  kine;
//          Turb2LamViscRatioIn[iMarkerTP][iSpan] = omega;
        }
        else{
          nu = solver_container->GetNuIn(iMarkerTP, iSpan);
          NuFactorIn[iMarkerTP][iSpan]          = nu*DensityIn[iMarkerTP][iSpan]/muLam;
          if (config->GetSAParsedOptions().bc) {
            NuFactorIn[iMarkerTP][iSpan]        = nu*DensityIn[iMarkerTP][iSpan]/muLam/0.005;
          }
        }
      }

      /*--- OUTFLOW ---*/
      /*--- Retrieve Outflow primitive quantities ---*/
      DensityOut[iMarkerTP][iSpan]         = solver_container->GetDensityOut(iMarkerTP, iSpan);
      PressureOut[iMarkerTP][iSpan]        = solver_container->GetPressureOut(iMarkerTP, iSpan);
      absVel2 = 0.0;

      for (iDim = 0; iDim < nDim; iDim++){
        TurboVelocityOut[iMarkerTP][iSpan][iDim]    = solver_container->GetTurboVelocityOut(iMarkerTP, iSpan)[iDim];
        absVel2   += TurboVelocityOut[iMarkerTP][iSpan][iDim]*TurboVelocityOut[iMarkerTP][iSpan][iDim];
      }
      TurboVelocityOut[iMarkerTP][iSpan][nDim] = sqrt(absVel2);


      for (iDim = 0; iDim < 3; iDim++){
      }
      area   = geometry->GetSpanAreaOut(iMarkerTP, iSpan);


      /*--- Compute all the Outflow quantities ---*/
      FluidModel->SetTDState_Prho(PressureOut[iMarkerTP][iSpan], DensityOut[iMarkerTP][iSpan]);
      EntropyOut[iMarkerTP][iSpan]          = FluidModel->GetEntropy();
      MassFlowOut[iMarkerTP][iSpan]         = config->GetnBlades(iMarkerTP)*DensityOut[iMarkerTP][iSpan]*TurboVelocityOut[iMarkerTP][iSpan][0]*area;
      AbsFlowAngleOut[iMarkerTP][iSpan]     = atan(TurboVelocityOut[iMarkerTP][iSpan][1]/TurboVelocityOut[iMarkerTP][iSpan][0]);
      EnthalpyOut[iMarkerTP][iSpan]         = FluidModel->GetStaticEnergy() + PressureOut[iMarkerTP][iSpan]/DensityOut[iMarkerTP][iSpan];
      soundSpeed                            = FluidModel->GetSoundSpeed();

      /*--- Compute Total Outflow quantities ---*/
      TotalEnthalpyOut[iMarkerTP][iSpan]    = EnthalpyOut[iMarkerTP][iSpan] + 0.5*absVel2;
      FluidModel->SetTDState_hs(TotalEnthalpyOut[iMarkerTP][iSpan], EntropyOut[iMarkerTP][iSpan]);
      TotalPressureOut[iMarkerTP][iSpan]    = FluidModel->GetPressure();
      TotalTemperatureOut[iMarkerTP][iSpan] = FluidModel->GetTemperature();

      /*--- Retrieve relative Outflow  quantities ---*/
      tangVel  = geometry->GetTangGridVelOut(iMarkerTP, iSpan);
      tangVel2 = tangVel*tangVel;

      for (iDim = 0; iDim < nDim; iDim++){
        relVel[iDim] = TurboVelocityOut[iMarkerTP][iSpan][iDim];
      }
      relVel[1] -= tangVel;

      relVel2 = 0.0;
      for (iDim = 0; iDim < nDim; iDim++){
        relVel2 += relVel[iDim]*relVel[iDim];
      }

      /*--- Compute Total relative Outflow quantities ---*/
      RothalpyOut[iMarkerTP][iSpan] = EnthalpyOut[iMarkerTP][iSpan] + 0.5*relVel2 - 0.5*tangVel2;
      FluidModel->SetTDState_hs(RothalpyOut[iMarkerTP][iSpan], EntropyOut[iMarkerTP][iSpan]);
      relPressureOut  = FluidModel->GetPressure();

      /*--- Compute isentropic Outflow quantities ---*/
      FluidModel->SetTDState_Ps(PressureOut[iMarkerTP][iSpan], EntropyIn[iMarkerTP][iSpan]);
      enthalpyOutIs   = FluidModel->GetStaticEnergy() + PressureOut[iMarkerTP][iSpan]/FluidModel->GetDensity();
      relVelOutIs2    = 2*(RothalpyOut[iMarkerTP][iSpan] - enthalpyOutIs) + tangVel2;


      /*--- Compute kinematic relative Outflow quantities ---*/
      FlowAngleOut[iMarkerTP][iSpan] = atan(relVel[1]/relVel[0]);
      mach   = 0.0;
      for (iDim = 0; iDim < nDim; iDim++){
        MachOut[iMarkerTP][iSpan][iDim] = relVel[iDim]/soundSpeed;
        mach = MachOut[iMarkerTP][iSpan][iDim]*MachOut[iMarkerTP][iSpan][iDim];
      }
      MachOut[iMarkerTP][iSpan][nDim]   = sqrt(mach);

      /*--- Compute Turbulent Outflow quantities ---*/
      if(turbulent){
        FluidModel->SetTDState_Prho(PressureOut[iMarkerTP][iSpan], DensityOut[iMarkerTP][iSpan]);
        muLam  = FluidModel->GetLaminarViscosity();
        if(menter_sst){
          kine   = solver_container->GetKineOut(iMarkerTP, iSpan);
          omega  = solver_container->GetOmegaOut(iMarkerTP, iSpan);
          TurbIntensityOut[iMarkerTP][iSpan]     =  sqrt(2.0/3.0*kine/absVel2);
          Turb2LamViscRatioOut[iMarkerTP][iSpan] = DensityOut[iMarkerTP][iSpan]*kine/(muLam*omega);
//          TurbIntensityOut[iMarkerTP][iSpan]     =  kine;
//          Turb2LamViscRatioOut[iMarkerTP][iSpan] = omega;
        }
        else{
          nu = solver_container->GetNuOut(iMarkerTP, iSpan);
          NuFactorOut[iMarkerTP][iSpan]          = nu*DensityOut[iMarkerTP][iSpan]/muLam;
          if (config->GetSAParsedOptions().bc) {
            NuFactorOut[iMarkerTP][iSpan]        = nu*DensityOut[iMarkerTP][iSpan]/muLam/0.005;
          }
        }
      }

      /*--- TURBO-PERFORMANCE---*/
      EntropyGen[iMarkerTP][iSpan]         = (EntropyOut[iMarkerTP][iSpan] - EntropyIn[iMarkerTP][iSpan])/abs(EntropyIn_BC[iMarkerTP][iSpan] + 1);
      EulerianWork[iMarkerTP][iSpan]       = TotalEnthalpyIn[iMarkerTP][iSpan] - TotalEnthalpyOut[iMarkerTP][iSpan];
      TotalPressureLoss[iMarkerTP][iSpan]  = (relPressureIn - relPressureOut)/(relPressureIn - PressureOut[iMarkerTP][iSpan]);
      KineticEnergyLoss[iMarkerTP][iSpan]  = 2*(EnthalpyOut[iMarkerTP][iSpan] - enthalpyOutIs)/relVelOutIs2;
      PressureRatio[iMarkerTP][iSpan]      = TotalPressureOut[iMarkerTP][iSpan]/TotalPressureIn[iMarkerTP][iSpan];
      EnthalpyOutIs[iMarkerTP][iSpan]      = (pow(TotalPressureOut[iMarkerTP][iSpan]/TotalPressureIn[iMarkerTP][iSpan], 0.4/1.4) - 1.0)/(TotalTemperatureOut[iMarkerTP][iSpan]/TotalTemperatureIn[iMarkerTP][iSpan] -1.0);
    }
  }

  if(nBladesRow > 1){
    /*--- Compute performance for each stage ---*/

    EulerianWork[nBladesRow + nStages][nSpanWiseSections]           = 0.0;
    /*---Comnpute performance for each stage---*/
    for(iStage = 0; iStage < nStages; iStage++ ){
      FluidModel->SetTDState_Ps(PressureOut[iStage*2 +1][config->GetnSpan_iZones(iStage*2 +1)], EntropyIn[iStage*2][config->GetnSpan_iZones(iStage*2)]);
      EnthalpyOutIs[nBladesRow + iStage][nSpanWiseSections]         = FluidModel->GetStaticEnergy() + PressureOut[iStage*2 +1][config->GetnSpan_iZones(iStage*2 +1)]/FluidModel->GetDensity();
      FluidModel->SetTDState_Prho(PressureOut[iStage*2 +1][config->GetnSpan_iZones(iStage*2 +1)], DensityOut[iStage*2 +1][config->GetnSpan_iZones(iStage*2 +1)]);
      absVel2 = 0.0;
      for (iDim = 0; iDim<nDim; iDim++)
        absVel2 += TurboVelocityOut[iStage*2 +1][config->GetnSpan_iZones(iStage*2 +1)][iDim]*TurboVelocityOut[iStage*2 +1][config->GetnSpan_iZones(iStage*2 +1)][iDim];
      TotalEnthalpyOutIs[nBladesRow + iStage][nSpanWiseSections]    = EnthalpyOutIs[nBladesRow + iStage][nSpanWiseSections] + 0.5*absVel2;

      TotalTotalEfficiency[nBladesRow + iStage][nSpanWiseSections]  = (TotalEnthalpyIn[iStage*2][config->GetnSpan_iZones(iStage*2)] - TotalEnthalpyOut[iStage*2 + 1][config->GetnSpan_iZones(iStage*2 +1)]);
      TotalTotalEfficiency[nBladesRow + iStage][nSpanWiseSections]  /= (TotalEnthalpyIn[iStage*2][config->GetnSpan_iZones(iStage*2)] - TotalEnthalpyOutIs[nBladesRow + iStage][nSpanWiseSections]);
      TotalStaticEfficiency[nBladesRow + iStage][nSpanWiseSections] = (TotalEnthalpyIn[iStage*2][config->GetnSpan_iZones(iStage*2)] - TotalEnthalpyOut[iStage*2 + 1][config->GetnSpan_iZones(iStage*2+1)]);
      TotalStaticEfficiency[nBladesRow + iStage][nSpanWiseSections] /= (TotalEnthalpyIn[iStage*2][config->GetnSpan_iZones(iStage*2)] - EnthalpyOutIs[nBladesRow + iStage][nSpanWiseSections]);
      PressureRatio[nBladesRow + iStage][nSpanWiseSections]         = (PressureRatio[iStage*2][config->GetnSpan_iZones(iStage*2)]*PressureOut[iStage*2][config->GetnSpan_iZones(iStage*2)]/PressureOut[iStage*2 + 1][config->GetnSpan_iZones(iStage*2+1)]);
      MassFlowIn[nBladesRow + iStage][nSpanWiseSections]            = MassFlowIn[iStage*2][config->GetnSpan_iZones(iStage*2)];
      MassFlowOut[nBladesRow + iStage][nSpanWiseSections]           = MassFlowOut[iStage*2 + 1][config->GetnSpan_iZones(iStage*2+1)];
      EntropyGen[nBladesRow + iStage][nSpanWiseSections]            = EntropyGen[iStage*2 + 1][config->GetnSpan_iZones(iStage*2 +1)] + EntropyGen[iStage*2][config->GetnSpan_iZones(iStage*2)];

    }

    /*---Compute turbo performance for full machine---*/
    FluidModel->SetTDState_Ps(PressureOut[nBladesRow-1][config->GetnSpan_iZones(nBladesRow-1)], EntropyIn[0][config->GetnSpan_iZones(0)]);
    EnthalpyOutIs[nBladesRow + nStages][nSpanWiseSections]          = FluidModel->GetStaticEnergy() + PressureOut[nBladesRow-1][config->GetnSpan_iZones(nBladesRow-1)]/FluidModel->GetDensity();
    FluidModel->SetTDState_Prho(PressureOut[nBladesRow-1][config->GetnSpan_iZones(nBladesRow-1)], DensityOut[nBladesRow-1][config->GetnSpan_iZones(nBladesRow-1)]);
    absVel2 = 0.0;
    for (iDim = 0; iDim<nDim;iDim++) absVel2 += TurboVelocityOut[nBladesRow-1][config->GetnSpan_iZones(nBladesRow-1)][iDim]*TurboVelocityOut[nBladesRow-1][config->GetnSpan_iZones(nBladesRow-1)][iDim];
    TotalEnthalpyOutIs[nBladesRow + nStages][nSpanWiseSections]     = EnthalpyOutIs[nBladesRow + nStages][nSpanWiseSections] + 0.5*absVel2;

    TotalTotalEfficiency[nBladesRow + nStages][nSpanWiseSections]   = (TotalEnthalpyIn[0][config->GetnSpan_iZones(0)] - TotalEnthalpyOut[nBladesRow-1][config->GetnSpan_iZones(nBladesRow-1)]);
    TotalTotalEfficiency[nBladesRow + nStages][nSpanWiseSections]  /= (TotalEnthalpyIn[0][config->GetnSpan_iZones(0)] - TotalEnthalpyOutIs[nBladesRow + nStages][nSpanWiseSections]);
    TotalStaticEfficiency[nBladesRow +nStages][nSpanWiseSections]   = (TotalEnthalpyIn[0][config->GetnSpan_iZones(0)] - TotalEnthalpyOut[nBladesRow-1][config->GetnSpan_iZones(nBladesRow-1)]);
    TotalStaticEfficiency[nBladesRow +nStages][nSpanWiseSections]  /= (TotalEnthalpyIn[0][config->GetnSpan_iZones(0)] - EnthalpyOutIs[nBladesRow + nStages][nSpanWiseSections]);
    PressureRatio[nBladesRow + nStages][nSpanWiseSections]          = PressureRatio[0][config->GetnSpan_iZones(0)]*PressureOut[0][config->GetnSpan_iZones(0)]/PressureOut[nBladesRow-1][config->GetnSpan_iZones(nBladesRow-1)];
    MassFlowIn[nBladesRow + nStages][nSpanWiseSections]             = MassFlowIn[0][config->GetnSpan_iZones(0)];
    MassFlowOut[nBladesRow + nStages][nSpanWiseSections]            = MassFlowOut[nBladesRow-1][config->GetnSpan_iZones(nBladesRow-1)];

    EntropyGen[nBladesRow + nStages][nSpanWiseSections]             = 0.0;
    for(iBlade = 0; iBlade < nBladesRow; iBlade++ ){
      EntropyGen[nBladesRow + nStages][nSpanWiseSections]          += EntropyGen[iBlade][config->GetnSpan_iZones(iBlade)];
    }
  }

  delete [] relVel;
}