/*!
 * \file output_physics.cpp
 * \brief Main subroutines to compute physical output quantities such as CL, CD, entropy generation, mass flow, ecc... .
 * \author S. Vitale
 * \version 5.0.0 "Raven"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
 *                      Dr. Thomas D. Economon (economon@stanford.edu).
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *                 Prof. Edwin van der Weide's group at the University of Twente.
 *                 Prof. Vincent Terrapon's group at the University of Liege.
 *
 * Copyright (C) 2012-2017 SU2, the open-source CFD code.
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


void COutput::ComputeTurboPerformance(CSolver *solver_container, CGeometry *geometry, CConfig *config) {

  CFluidModel *FluidModel;
  unsigned short nDim = geometry->GetnDim();
  unsigned short iTimeInstance = config->GetiZone();
  unsigned short iMarkerTP, iSpan, iDim, iStage, iBlade;
  unsigned short nMarkerTP = config->GetnMarker_Turbomachinery();
  FluidModel = solver_container->GetFluidModel();
  su2double area, absVel2, soundSpeed, mach, tangVel, tangVel2, *relVel, relVel2;
  su2double relPressureIn, relPressureOut, enthalpyOutIs, relVelOutIs2;
  relVel = new su2double[nDim];

  unsigned short nBladesRow, nStages;

  nBladesRow = config->GetnMarker_Turbomachinery();
  nStages    = SU2_TYPE::Int(nBladesRow/2);


  /*--- Compute BC imposed value for convergence monitoring ---*/
  for(iMarkerTP = 0; iMarkerTP < nMarkerTP; iMarkerTP++ ){
    for(iSpan = 0; iSpan < config->GetnSpan_iZones(iMarkerTP) + 1; iSpan++){
      if(config->GetRampOutletPressure() && config->GetExtIter() > 0){
        PressureOut_BC[nMarkerTurboPerf*iTimeInstance + iMarkerTP][iSpan] = config->GetMonitorOutletPressure()/config->GetPressure_Ref();
      }
      FluidModel->SetTDState_PT(config->GetTotalPressureIn_BC(), config->GetTotalTemperatureIn_BC());
      TotalEnthalpyIn_BC[nMarkerTurboPerf*iTimeInstance + iMarkerTP][iSpan] = FluidModel->GetStaticEnergy()+ FluidModel->GetPressure()/FluidModel->GetDensity();
      EntropyIn_BC[nMarkerTurboPerf*iTimeInstance + iMarkerTP][iSpan] = FluidModel->GetEntropy();
    }
  }

  /*--- Compute performance for each blade ---*/
  for(iMarkerTP = 0; iMarkerTP < nMarkerTP; iMarkerTP++ ){
    for(iSpan = 0; iSpan < config->GetnSpan_iZones(iMarkerTP) + 1; iSpan++){


      /*--- INFLOW ---*/
      /*--- Retrieve Inflow primitive quantities ---*/
      DensityIn[nMarkerTurboPerf*iTimeInstance + iMarkerTP][iSpan]          = solver_container->GetDensityIn(iMarkerTP, iSpan);
      PressureIn[nMarkerTurboPerf*iTimeInstance + iMarkerTP][iSpan]         = solver_container->GetPressureIn(iMarkerTP, iSpan);

      absVel2 = 0.0;

      for (iDim = 0; iDim < nDim; iDim++){
        TurboVelocityIn[nMarkerTurboPerf*iTimeInstance + iMarkerTP][iSpan][iDim] = solver_container->GetTurboVelocityIn(iMarkerTP, iSpan)[iDim];
        absVel2   += TurboVelocityIn[nMarkerTurboPerf*iTimeInstance + iMarkerTP][iSpan][iDim]*TurboVelocityIn[nMarkerTurboPerf*iTimeInstance + iMarkerTP][iSpan][iDim];
      }
      TurboVelocityIn[nMarkerTurboPerf*iTimeInstance + iMarkerTP][iSpan][nDim] = sqrt(absVel2);

      TRadius[nMarkerTurboPerf*iTimeInstance + iMarkerTP][iSpan] = geometry->GetTurboRadiusIn(iMarkerTP, iSpan);
      area = geometry->GetSpanAreaIn(iMarkerTP, iSpan);

      /*--- Compute static Inflow quantities ---*/
      FluidModel->SetTDState_Prho(PressureIn[nMarkerTurboPerf*iTimeInstance + iMarkerTP][iSpan], DensityIn[nMarkerTurboPerf*iTimeInstance + iMarkerTP][iSpan]);
      EntropyIn[nMarkerTurboPerf*iTimeInstance + iMarkerTP][iSpan]  = FluidModel->GetEntropy();
      MassFlowIn[nMarkerTurboPerf*iTimeInstance + iMarkerTP][iSpan] =  config->GetnBlades(iTimeInstance)*DensityIn[nMarkerTurboPerf*iTimeInstance + iMarkerTP][iSpan];
      MassFlowIn[nMarkerTurboPerf*iTimeInstance + iMarkerTP][iSpan] *= TurboVelocityIn[nMarkerTurboPerf*iTimeInstance + iMarkerTP][iSpan][0]*area;
      AbsFlowAngleIn[nMarkerTurboPerf*iTimeInstance + iMarkerTP][iSpan]     = atan(TurboVelocityIn[nMarkerTurboPerf*iTimeInstance + iMarkerTP][iSpan][1]/TurboVelocityIn[nMarkerTurboPerf*iTimeInstance + iMarkerTP][iSpan][0]);
      EnthalpyIn[nMarkerTurboPerf*iTimeInstance + iMarkerTP][iSpan]         = FluidModel->GetStaticEnergy() + PressureIn[nMarkerTurboPerf*iTimeInstance + iMarkerTP][iSpan]/DensityIn[nMarkerTurboPerf*iTimeInstance + iMarkerTP][iSpan];
      soundSpeed                           = FluidModel->GetSoundSpeed();


      /*--- Compute Total Inflow quantities ---*/
      TotalEnthalpyIn[nMarkerTurboPerf*iTimeInstance + iMarkerTP][iSpan]    = EnthalpyIn[nMarkerTurboPerf*iTimeInstance + iMarkerTP][iSpan] + 0.5*absVel2;
      FluidModel->SetTDState_hs(TotalEnthalpyIn[nMarkerTurboPerf*iTimeInstance + iMarkerTP][iSpan], EntropyIn[nMarkerTurboPerf*iTimeInstance + iMarkerTP][iSpan]);
      TotalPressureIn[nMarkerTurboPerf*iTimeInstance + iMarkerTP][iSpan]    = FluidModel->GetPressure();
      TotalTemperatureIn[nMarkerTurboPerf*iTimeInstance + iMarkerTP][iSpan] = FluidModel->GetTemperature();

      /*--- Retrieve Inflow relative quantities ---*/
      tangVel = geometry->GetTangGridVelIn(iMarkerTP, iSpan);
      tangVel2 = tangVel*tangVel;

      for (iDim = 0; iDim < nDim; iDim++){
        relVel[iDim] = TurboVelocityIn[nMarkerTurboPerf*iTimeInstance + iMarkerTP][iSpan][iDim];
      }
      relVel[1] -= tangVel;

      relVel2 = 0.0;
      for (iDim = 0; iDim < nDim; iDim++){
        relVel2 += relVel[iDim]*relVel[iDim];
      }

      /*--- Compute Total relative Inflow quantities ---*/
      RothalpyIn[nMarkerTurboPerf*iTimeInstance + iMarkerTP][iSpan]  = EnthalpyIn[nMarkerTurboPerf*iTimeInstance + iMarkerTP][iSpan] + 0.5*relVel2 - 0.5*tangVel2;
      FluidModel->SetTDState_hs(RothalpyIn[nMarkerTurboPerf*iTimeInstance + iMarkerTP][iSpan], EntropyIn[nMarkerTurboPerf*iTimeInstance + iMarkerTP][iSpan]);
      relPressureIn   = FluidModel->GetPressure();

      /*--- Compute kinematic relative Inflow quantities ---*/
      FlowAngleIn[nMarkerTurboPerf*iTimeInstance + iMarkerTP][iSpan]   = atan(relVel[1]/relVel[0]);
      mach          = 0.0;
      for (iDim = 0; iDim < nDim; iDim++){
        MachIn[nMarkerTurboPerf*iTimeInstance + iMarkerTP][iSpan][iDim]       = relVel[iDim]/soundSpeed;
        mach = MachIn[nMarkerTurboPerf*iTimeInstance + iMarkerTP][iSpan][iDim]*MachIn[nMarkerTurboPerf*iTimeInstance + iMarkerTP][iSpan][iDim];
      }
      MachIn[nMarkerTurboPerf*iTimeInstance + iMarkerTP][iSpan][nDim]            = sqrt(mach);



      /*--- OUTFLOW ---*/
      /*--- Retrieve Outflow primitive quantities ---*/
      DensityOut[nMarkerTurboPerf*iTimeInstance + iMarkerTP][iSpan]         = solver_container->GetDensityOut(iMarkerTP, iSpan);
      PressureOut[nMarkerTurboPerf*iTimeInstance + iMarkerTP][iSpan]         = solver_container->GetPressureOut(iMarkerTP, iSpan);
      absVel2 = 0.0;

      for (iDim = 0; iDim < nDim; iDim++){
        TurboVelocityOut[nMarkerTurboPerf*iTimeInstance + iMarkerTP][iSpan][iDim]    = solver_container->GetTurboVelocityOut(iMarkerTP, iSpan)[iDim];
        absVel2   += TurboVelocityOut[nMarkerTurboPerf*iTimeInstance + iMarkerTP][iSpan][iDim]*TurboVelocityOut[nMarkerTurboPerf*iTimeInstance + iMarkerTP][iSpan][iDim];
      }
      TurboVelocityOut[nMarkerTurboPerf*iTimeInstance + iMarkerTP][iSpan][nDim] = sqrt(absVel2);


      for (iDim = 0; iDim < 3; iDim++){
      }
      //			TRadius[iMarkerTP][iSpan]            = geometry->GetTurboRadiusIn(iMarkerTP, iSpan);
      area																 = geometry->GetSpanAreaOut(iMarkerTP, iSpan);


      /*--- Compute all the Outflow quantities ---*/
      FluidModel->SetTDState_Prho(PressureOut[nMarkerTurboPerf*iTimeInstance + iMarkerTP][iSpan], DensityOut[nMarkerTurboPerf*iTimeInstance + iMarkerTP][iSpan]);
      EntropyOut[nMarkerTurboPerf*iTimeInstance + iMarkerTP][iSpan]				 = FluidModel->GetEntropy();
      MassFlowOut[nMarkerTurboPerf*iTimeInstance + iMarkerTP][iSpan]  = config->GetnBlades(iTimeInstance)*DensityOut[nMarkerTurboPerf*iTimeInstance + iMarkerTP][iSpan];
      MassFlowOut[nMarkerTurboPerf*iTimeInstance + iMarkerTP][iSpan] *= TurboVelocityOut[nMarkerTurboPerf*iTimeInstance + iMarkerTP][iSpan][0]*area;
      AbsFlowAngleOut[nMarkerTurboPerf*iTimeInstance + iMarkerTP][iSpan]     = atan(TurboVelocityOut[nMarkerTurboPerf*iTimeInstance + iMarkerTP][iSpan][1]/TurboVelocityOut[nMarkerTurboPerf*iTimeInstance + iMarkerTP][iSpan][0]);
      EnthalpyOut[nMarkerTurboPerf*iTimeInstance + iMarkerTP][iSpan]         = FluidModel->GetStaticEnergy() + PressureOut[nMarkerTurboPerf*iTimeInstance + iMarkerTP][iSpan]/DensityOut[nMarkerTurboPerf*iTimeInstance + iMarkerTP][iSpan];
      soundSpeed                           = FluidModel->GetSoundSpeed();

      /*--- Compute Total Outflow quantities ---*/
      TotalEnthalpyOut[nMarkerTurboPerf*iTimeInstance + iMarkerTP][iSpan]    = EnthalpyOut[nMarkerTurboPerf*iTimeInstance + iMarkerTP][iSpan] + 0.5*absVel2;
      FluidModel->SetTDState_hs(TotalEnthalpyOut[nMarkerTurboPerf*iTimeInstance + iMarkerTP][iSpan], EntropyOut[nMarkerTurboPerf*iTimeInstance + iMarkerTP][iSpan]);
      TotalPressureOut[nMarkerTurboPerf*iTimeInstance + iMarkerTP][iSpan]    = FluidModel->GetPressure();
      TotalTemperatureOut[nMarkerTurboPerf*iTimeInstance + iMarkerTP][iSpan] = FluidModel->GetTemperature();

      /*--- Retrieve relative Outflow  quantities ---*/
      tangVel = geometry->GetTangGridVelOut(iMarkerTP, iSpan);
      tangVel2 = tangVel*tangVel;

      for (iDim = 0; iDim < nDim; iDim++){
        relVel[iDim] = TurboVelocityOut[nMarkerTurboPerf*iTimeInstance + iMarkerTP][iSpan][iDim];
      }
      relVel[1] -= tangVel;

      relVel2 = 0.0;
      for (iDim = 0; iDim < nDim; iDim++){
        relVel2 += relVel[iDim]*relVel[iDim];
      }

      /*--- Compute Total relative Outflow quantities ---*/
      RothalpyOut[nMarkerTurboPerf*iTimeInstance + iMarkerTP][iSpan] = EnthalpyOut[nMarkerTurboPerf*iTimeInstance + iMarkerTP][iSpan] + 0.5*relVel2 - 0.5*tangVel2;
      FluidModel->SetTDState_hs(RothalpyOut[nMarkerTurboPerf*iTimeInstance + iMarkerTP][iSpan], EntropyOut[nMarkerTurboPerf*iTimeInstance + iMarkerTP][iSpan]);
      relPressureOut    = FluidModel->GetPressure();

      /*--- Compute isentropic Outflow quantities ---*/
      FluidModel->SetTDState_Ps(PressureOut[nMarkerTurboPerf*iTimeInstance + iMarkerTP][iSpan], EntropyIn[nMarkerTurboPerf*iTimeInstance + iMarkerTP][iSpan]);
      enthalpyOutIs = FluidModel->GetStaticEnergy() + PressureOut[nMarkerTurboPerf*iTimeInstance + iMarkerTP][iSpan]/FluidModel->GetDensity();
      relVelOutIs2  = 2*(RothalpyOut[nMarkerTurboPerf*iTimeInstance + iMarkerTP][iSpan] - enthalpyOutIs) + tangVel2;


      /*--- Compute kinematic relative Outflow quantities ---*/
      FlowAngleOut[nMarkerTurboPerf*iTimeInstance + iMarkerTP][iSpan]   = atan(relVel[1]/relVel[0]);
      mach          = 0.0;
      for (iDim = 0; iDim < nDim; iDim++){
        MachOut[nMarkerTurboPerf*iTimeInstance + iMarkerTP][iSpan][iDim]       = relVel[iDim]/soundSpeed;
        mach = MachOut[nMarkerTurboPerf*iTimeInstance + iMarkerTP][iSpan][iDim]*MachOut[nMarkerTurboPerf*iTimeInstance + iMarkerTP][iSpan][iDim];
      }
      MachOut[nMarkerTurboPerf*iTimeInstance + iMarkerTP][iSpan][nDim]            = sqrt(mach);



      /*--- TURBO-PERFORMANCE---*/
      EntropyGen[nMarkerTurboPerf*iTimeInstance + iMarkerTP][iSpan] = (EntropyOut[nMarkerTurboPerf*iTimeInstance + iMarkerTP][iSpan] -
          EntropyIn[nMarkerTurboPerf*iTimeInstance + iMarkerTP][iSpan])/abs(EntropyIn_BC[nMarkerTurboPerf*iTimeInstance + iMarkerTP][iSpan] + 1);
      EulerianWork[nMarkerTurboPerf*iTimeInstance + iMarkerTP][iSpan]       = TotalEnthalpyIn[nMarkerTurboPerf*iTimeInstance + iMarkerTP][iSpan] - TotalEnthalpyOut[nMarkerTurboPerf*iTimeInstance + iMarkerTP][iSpan];
      TotalPressureLoss[nMarkerTurboPerf*iTimeInstance + iMarkerTP][iSpan]  = (relPressureIn - relPressureOut)/(relPressureIn - PressureOut[nMarkerTurboPerf*iTimeInstance + iMarkerTP][iSpan]);
      KineticEnergyLoss[nMarkerTurboPerf*iTimeInstance + iMarkerTP][iSpan]  = 2*(EnthalpyOut[nMarkerTurboPerf*iTimeInstance + iMarkerTP][iSpan] - enthalpyOutIs)/relVelOutIs2;
      PressureRatio[nMarkerTurboPerf*iTimeInstance + iMarkerTP][iSpan]      = TotalPressureOut[nMarkerTurboPerf*iTimeInstance + iMarkerTP][iSpan]/TotalPressureIn[nMarkerTurboPerf*iTimeInstance + iMarkerTP][iSpan];

    }
  }

  if(nBladesRow > 1){
    /*--- Compute performance for each stage ---*/

    EulerianWork[nMarkerTurboPerf*iTimeInstance +  nBladesRow + nStages][nSpanWiseSections]           = 0.0;
    /*---Comnpute performance for each stage---*/
    for(iStage = 0; iStage < nStages; iStage++ ){
      FluidModel->SetTDState_Ps(PressureOut[nMarkerTurboPerf*iTimeInstance + iStage*2 +1][config->GetnSpan_iZones(iStage*2 +1)], EntropyIn[nMarkerTurboPerf*iTimeInstance + iStage*2][config->GetnSpan_iZones(iStage*2)]);
      EnthalpyOutIs[nMarkerTurboPerf*iTimeInstance + nBladesRow + iStage][nSpanWiseSections]         = FluidModel->GetStaticEnergy() + PressureOut[nMarkerTurboPerf*iTimeInstance + iStage*2 +1][config->GetnSpan_iZones(iStage*2 +1)]/FluidModel->GetDensity();
      FluidModel->SetTDState_Prho(PressureOut[nMarkerTurboPerf*iTimeInstance + iStage*2 +1][config->GetnSpan_iZones(iStage*2 +1)], DensityOut[nMarkerTurboPerf*iTimeInstance + iStage*2 +1][config->GetnSpan_iZones(iStage*2 +1)]);
      absVel2 = 0.0;
      for (iDim = 0; iDim<nDim; iDim++)
        absVel2 += TurboVelocityOut[nMarkerTurboPerf*iTimeInstance + iStage*2 +1][config->GetnSpan_iZones(iStage*2 +1)][iDim]*TurboVelocityOut[nMarkerTurboPerf*iTimeInstance + iStage*2 +1][config->GetnSpan_iZones(iStage*2 +1)][iDim];
      TotalEnthalpyOutIs[nMarkerTurboPerf*iTimeInstance + nBladesRow + iStage][nSpanWiseSections]    = EnthalpyOutIs[nMarkerTurboPerf*iTimeInstance + nBladesRow + iStage][nSpanWiseSections] + 0.5*absVel2;

      TotalTotalEfficiency[nMarkerTurboPerf*iTimeInstance + nBladesRow + iStage][nSpanWiseSections]  = (TotalEnthalpyIn[nMarkerTurboPerf*iTimeInstance + iStage*2][config->GetnSpan_iZones(iStage*2)] - TotalEnthalpyOut[nMarkerTurboPerf*iTimeInstance + iStage*2 + 1][config->GetnSpan_iZones(iStage*2 +1)]);
      TotalTotalEfficiency[nMarkerTurboPerf*iTimeInstance + nBladesRow + iStage][nSpanWiseSections]  /= (TotalEnthalpyIn[nMarkerTurboPerf*iTimeInstance + iStage*2][config->GetnSpan_iZones(iStage*2)] - TotalEnthalpyOutIs[nMarkerTurboPerf*iTimeInstance + nBladesRow + iStage][nSpanWiseSections]);
      TotalStaticEfficiency[nMarkerTurboPerf*iTimeInstance + nBladesRow + iStage][nSpanWiseSections] = (TotalEnthalpyIn[nMarkerTurboPerf*iTimeInstance + iStage*2][config->GetnSpan_iZones(iStage*2)] - TotalEnthalpyOut[nMarkerTurboPerf*iTimeInstance + iStage*2 + 1][config->GetnSpan_iZones(iStage*2+1)]);
      TotalStaticEfficiency[nMarkerTurboPerf*iTimeInstance + nBladesRow + iStage][nSpanWiseSections] /= (TotalEnthalpyIn[nMarkerTurboPerf*iTimeInstance + iStage*2][config->GetnSpan_iZones(iStage*2)] - EnthalpyOutIs[nMarkerTurboPerf*iTimeInstance + nBladesRow + iStage][nSpanWiseSections]);
      PressureRatio[nMarkerTurboPerf*iTimeInstance + nBladesRow + iStage][nSpanWiseSections]         = (PressureRatio[nMarkerTurboPerf*iTimeInstance + iStage*2][config->GetnSpan_iZones(iStage*2)]*PressureOut[nMarkerTurboPerf*iTimeInstance + iStage*2][config->GetnSpan_iZones(iStage*2)]/PressureOut[nMarkerTurboPerf*iTimeInstance + iStage*2 + 1][config->GetnSpan_iZones(iStage*2+1)]);
      MassFlowIn[nMarkerTurboPerf*iTimeInstance + nBladesRow + iStage][nSpanWiseSections]            = MassFlowIn[nMarkerTurboPerf*iTimeInstance + iStage*2][config->GetnSpan_iZones(iStage*2)];
      MassFlowOut[nMarkerTurboPerf*iTimeInstance + nBladesRow + iStage][nSpanWiseSections]           = MassFlowOut[nMarkerTurboPerf*iTimeInstance + iStage*2 + 1][config->GetnSpan_iZones(iStage*2+1)];
      EntropyGen[nMarkerTurboPerf*iTimeInstance + nBladesRow + iStage][nSpanWiseSections]            = EntropyGen[nMarkerTurboPerf*iTimeInstance + iStage*2 + 1][config->GetnSpan_iZones(iStage*2 +1)] + EntropyGen[nMarkerTurboPerf*iTimeInstance + iStage*2][config->GetnSpan_iZones(iStage*2)];

    }

    /*---Compute turbo performance for full machine---*/
    FluidModel->SetTDState_Ps(PressureOut[nMarkerTurboPerf*iTimeInstance + nBladesRow-1][config->GetnSpan_iZones(nBladesRow-1)], EntropyIn[nMarkerTurboPerf*iTimeInstance][config->GetnSpan_iZones(0)]);
    EnthalpyOutIs[nMarkerTurboPerf*iTimeInstance + nBladesRow + nStages][nSpanWiseSections]          = FluidModel->GetStaticEnergy() + PressureOut[nMarkerTurboPerf*iTimeInstance + nBladesRow-1][config->GetnSpan_iZones(nBladesRow-1)]/FluidModel->GetDensity();
    FluidModel->SetTDState_Prho(PressureOut[nMarkerTurboPerf*iTimeInstance + nBladesRow-1][config->GetnSpan_iZones(nBladesRow-1)], DensityOut[nMarkerTurboPerf*iTimeInstance + nBladesRow-1][config->GetnSpan_iZones(nBladesRow-1)]);
    absVel2 = 0.0;
    for (iDim = 0; iDim<nDim;iDim++) absVel2 += TurboVelocityOut[nMarkerTurboPerf*iTimeInstance + nBladesRow-1][config->GetnSpan_iZones(nBladesRow-1)][iDim]*TurboVelocityOut[nMarkerTurboPerf*iTimeInstance + nBladesRow-1][config->GetnSpan_iZones(nBladesRow-1)][iDim];
    TotalEnthalpyOutIs[nMarkerTurboPerf*iTimeInstance + nBladesRow + nStages][nSpanWiseSections]     = EnthalpyOutIs[nMarkerTurboPerf*iTimeInstance + nBladesRow + nStages][nSpanWiseSections] + 0.5*absVel2;

    TotalTotalEfficiency[nMarkerTurboPerf*iTimeInstance + nBladesRow + nStages][nSpanWiseSections]   = (TotalEnthalpyIn[nMarkerTurboPerf*iTimeInstance][config->GetnSpan_iZones(0)] - TotalEnthalpyOut[nMarkerTurboPerf*iTimeInstance + nBladesRow-1][config->GetnSpan_iZones(nBladesRow-1)]);
    TotalTotalEfficiency[nMarkerTurboPerf*iTimeInstance + nBladesRow + nStages][nSpanWiseSections]  /= (TotalEnthalpyIn[nMarkerTurboPerf*iTimeInstance][config->GetnSpan_iZones(0)] - TotalEnthalpyOutIs[nMarkerTurboPerf*iTimeInstance + nBladesRow + nStages][nSpanWiseSections]);
    TotalStaticEfficiency[nMarkerTurboPerf*iTimeInstance + nBladesRow +nStages][nSpanWiseSections]   = (TotalEnthalpyIn[nMarkerTurboPerf*iTimeInstance][config->GetnSpan_iZones(0)] - TotalEnthalpyOut[nMarkerTurboPerf*iTimeInstance + nBladesRow-1][config->GetnSpan_iZones(nBladesRow-1)]);
    TotalStaticEfficiency[nMarkerTurboPerf*iTimeInstance + nBladesRow +nStages][nSpanWiseSections]  /= (TotalEnthalpyIn[nMarkerTurboPerf*iTimeInstance][config->GetnSpan_iZones(0)] - EnthalpyOutIs[nMarkerTurboPerf*iTimeInstance + nBladesRow + nStages][nSpanWiseSections]);
    PressureRatio[nMarkerTurboPerf*iTimeInstance + nBladesRow + nStages][nSpanWiseSections]          = PressureRatio[nMarkerTurboPerf*iTimeInstance][config->GetnSpan_iZones(0)]*PressureOut[nMarkerTurboPerf*iTimeInstance][config->GetnSpan_iZones(0)]/PressureOut[nMarkerTurboPerf*iTimeInstance + nBladesRow-1][config->GetnSpan_iZones(nBladesRow-1)];
    MassFlowIn[nMarkerTurboPerf*iTimeInstance + nBladesRow + nStages][nSpanWiseSections]             = MassFlowIn[nMarkerTurboPerf*iTimeInstance][config->GetnSpan_iZones(0)];
    MassFlowOut[nMarkerTurboPerf*iTimeInstance + nBladesRow + nStages][nSpanWiseSections]            = MassFlowOut[nMarkerTurboPerf*iTimeInstance + nBladesRow-1][config->GetnSpan_iZones(nBladesRow-1)];

    EntropyGen[nMarkerTurboPerf*iTimeInstance + nBladesRow + nStages][nSpanWiseSections]     = 0.0;
    for(iBlade = 0; iBlade < nBladesRow; iBlade++ ){
      EntropyGen[nMarkerTurboPerf*iTimeInstance + nBladesRow + nStages][nSpanWiseSections]  += EntropyGen[iBlade][config->GetnSpan_iZones(iBlade)];
    }
  }

  delete [] relVel;
}

void COutput::ComputeAvgTurboPerformance_HB(CConfig *config, unsigned short nTimeInstances, unsigned short iGeomZone) {
  unsigned short iTimeInstance, nStages;
  unsigned short nBladesRow = config->GetnMarker_Turbomachinery();
  nStages    = SU2_TYPE::Int(nBladesRow/2);

  EntropyGenAverage_HB     [iGeomZone] = 0.;
  TotPressureLossAverage_HB[iGeomZone] = 0.;
  TotalTotalEfficiencyAverage_HB = 0.;
  TotalStaticEfficiencyAverage_HB= 0.;

  for (iTimeInstance = 0; iTimeInstance < nTimeInstances; iTimeInstance++ ){
    EntropyGenAverage_HB[iGeomZone]      += EntropyGen[iTimeInstance * nMarkerTurboPerf + iGeomZone][nSpanWiseSections];
    TotPressureLossAverage_HB[iGeomZone] += TotalPressureLoss[iTimeInstance * nMarkerTurboPerf + iGeomZone][nSpanWiseSections];
  }

  EntropyGenAverage_HB[iGeomZone] /= nTimeInstances;
  TotPressureLossAverage_HB[iGeomZone] /= nTimeInstances;
//  cout << EntropyGenAverage_HB[iGeomZone] << endl;
//  cout << TotPressureLossAverage_HB[iGeomZone] << endl;

  if (nBladesRow > 1){
    for (iTimeInstance = 0; iTimeInstance < nTimeInstances; iTimeInstance++ ){
      TotalTotalEfficiencyAverage_HB += TotalTotalEfficiency[nMarkerTurboPerf*iTimeInstance + nBladesRow + nStages][nSpanWiseSections]  ;
      TotalStaticEfficiencyAverage_HB += TotalStaticEfficiency[nMarkerTurboPerf*iTimeInstance + nBladesRow + nStages][nSpanWiseSections] ;
    }

    TotalTotalEfficiencyAverage_HB  /= nTimeInstances;
    TotalStaticEfficiencyAverage_HB /= nTimeInstances;


  }

}

su2double COutput::GetEntropyGen(unsigned short iMarkerTP, unsigned short iSpan){return EntropyGen[iMarkerTP][iSpan];}




