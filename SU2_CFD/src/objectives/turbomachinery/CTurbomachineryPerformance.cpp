/*!
 * \file CTurbomachineryPerformance.cpp
 * \brief Source of the Turbomachinery Performance class.
 * \author S. Vitale
 * \version 7.0.1 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2019, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../../include/objectives/turbomachinery/CTurbomachineryPerformance.hpp"

CTurbomachineryPrimitiveState::CTurbomachineryPrimitiveState() {
  Density = Pressure = TangVelocity = 0.0;
}
CTurbomachineryPrimitiveState::CTurbomachineryPrimitiveState(vector<su2double> TurboPrimitive,
                                                             unsigned short nDim,
                                                             su2double tangVel) : Density(TurboPrimitive[0]), Pressure(TurboPrimitive[1]), TangVelocity(tangVel) {
  // Velocity.assign(TurboPrimitive+2, TurboPrimitive + nDim+2 );
  Velocity = {TurboPrimitive.begin()+2, TurboPrimitive.end()};
}

CTurbomachineryCombinedPrimitiveStates::CTurbomachineryCombinedPrimitiveStates(const CTurbomachineryPrimitiveState& inletPrimitiveState, 
                                        const CTurbomachineryPrimitiveState&  outletPrimitiveState) : InletPrimitiveState(inletPrimitiveState), OutletPrimitiveState(outletPrimitiveState){}

CTurbomachineryState::CTurbomachineryState(){
  Density = Pressure = Entropy = Enthalpy = Temperature = TotalTemperature = TotalPressure = TotalEnthalpy = 0.0;
  AbsFlowAngle = FlowAngle = MassFlow = Rothalpy = TotalRelPressure = 0.0;
  Area = Radius = 0.0;
}

CTurbomachineryState::CTurbomachineryState(unsigned short nDim, su2double area, su2double radius) : CTurbomachineryState() {
  generate_n(back_inserter(Velocity),nDim, [] {return 0.0;});
  generate_n(back_inserter(RelVelocity),nDim, [] {return 0.0;});
  generate_n(back_inserter(Mach),nDim, [] {return 0.0;});
  generate_n(back_inserter(RelMach),nDim, [] {return 0.0;});
  Area = area;
  Radius = radius;
}

void CTurbomachineryState::ComputeState(CFluidModel& fluidModel, const CTurbomachineryPrimitiveState &primitiveState) {

  /*--- Assign new primitive values ---*/
  Density = primitiveState.GetDensity();
  Pressure = primitiveState.GetPressure();
  auto velocity = primitiveState.GetVelocity();
  Velocity.assign(velocity.begin(), velocity.end());
  auto tangVel = primitiveState.GetTangVelocity();

  /*--- Compute static TD quantities ---*/
  fluidModel.SetTDState_Prho(Pressure, Density);
  Entropy = fluidModel.GetEntropy();
  Enthalpy = fluidModel.GetStaticEnergy() + Pressure / Density;
  su2double soundSpeed = fluidModel.GetSoundSpeed();

  /*--- Compute total TD quantities ---*/
  TotalEnthalpy = Enthalpy + 0.5 * GetVelocityValue() * GetVelocityValue();
  fluidModel.SetTDState_hs(TotalEnthalpy, Entropy);
  TotalPressure = fluidModel.GetPressure();
  TotalTemperature = fluidModel.GetTemperature();

  /*--- Compute absolute kinematic quantities---*/
  MassFlow = Density * Velocity[0] * Area;
  AbsFlowAngle = atan(Velocity[1] / Velocity[0]);
  Mach.assign(Velocity.begin(), Velocity.end());
  std::for_each(Mach.begin(), Mach.end(), [&](su2double &el) { el /= soundSpeed; });

  /*--- Compute relative kinematic quantities ---*/
  auto tangVel2 = tangVel * tangVel;
  RelVelocity.assign(Velocity.begin(), Velocity.end());
  RelVelocity[1] -= tangVel;
  auto relVel2 = GetRelVelocityValue();
  FlowAngle = atan(RelVelocity[1] / RelVelocity[0]);
  RelMach.assign(RelVelocity.begin(), RelVelocity.end());
  std::for_each(RelMach.begin(), RelMach.end(), [&](su2double &el) { el /= soundSpeed; });

  /*--- Compute total relative TD quantities ---*/
  Rothalpy = Enthalpy + 0.5 * relVel2 - 0.5 * tangVel2;
  fluidModel.SetTDState_hs(Rothalpy, Entropy);
  TotalRelPressure = fluidModel.GetPressure();

  /*--- Compute isentropic quantities ---*/
  fluidModel.SetTDState_Ps(Pressure, Entropy);
  
}

CTurbomachineryBladePerformance::CTurbomachineryBladePerformance(CFluidModel& fluidModel,
                                                                 unsigned short nDim,
                                                                 su2double areaIn,
                                                                 su2double radiusIn,
                                                                 su2double areaOut,
                                                                 su2double radiusOut) : FluidModel(fluidModel) {
   InletState               = CTurbomachineryState(nDim, areaIn, radiusIn);
   OutletState              = CTurbomachineryState(nDim, areaOut,radiusOut);
}


CTurbineBladePerformance::CTurbineBladePerformance(CFluidModel& fluidModel,
                                                  unsigned short nDim,
                                                  su2double areaIn,
                                                  su2double radiusIn,
                                                  su2double areaOut,
                                                  su2double radiusOut) : CTurbomachineryBladePerformance(fluidModel, nDim, areaIn, radiusIn, areaOut, radiusOut){}

void CTurbineBladePerformance::ComputePerformance(const CTurbomachineryCombinedPrimitiveStates &primitives) {

  /*--- Compute Inlet and Outlet state ---*/
  InletState.ComputeState(FluidModel, primitives.GetInletPrimitiveState());
  OutletState.ComputeState(FluidModel, primitives.GetOutletPrimitiveState());

  /*--- Compute isentropic Outflow quantities ---*/
  FluidModel.SetTDState_Ps(OutletState.GetPressure(), InletState.GetEntropy());
  auto enthalpyOutIs = FluidModel.GetStaticEnergy() + OutletState.GetPressure() / FluidModel.GetDensity();
  auto tangVel = primitives.GetOutletPrimitiveState().GetTangVelocity();
  auto relVelOutIs2 = 2 * (OutletState.GetRothalpy() - enthalpyOutIs) + tangVel * tangVel;

  /*--- Compute performance ---*/
  EntropyGen = (OutletState.GetEntropy() - InletState.GetEntropy()); // / abs(InletState.GetEntropy() + 1);
  EulerianWork = InletState.GetTotalEnthalpy() - OutletState.GetTotalEnthalpy();
  TotalPressureLoss = (InletState.GetTotalRelPressure() - OutletState.GetTotalRelPressure()) /
                      (InletState.GetTotalRelPressure() - OutletState.GetPressure());
  KineticEnergyLoss = 2 * (OutletState.GetEnthalpy() - enthalpyOutIs) / relVelOutIs2;
}

CCompressorBladePerformance::CCompressorBladePerformance(CFluidModel& fluidModel,
                                                        unsigned short nDim,
                                                        su2double areaIn,
                                                        su2double radiusIn,
                                                        su2double areaOut,
                                                        su2double radiusOut) : CTurbomachineryBladePerformance(fluidModel, nDim, areaIn, radiusIn, areaOut, radiusOut){}

void CCompressorBladePerformance::ComputePerformance(const CTurbomachineryCombinedPrimitiveStates &primitives) {
  /*--- Compute Inlet and Outlet state ---*/
  InletState.ComputeState(FluidModel, primitives.GetInletPrimitiveState());
  OutletState.ComputeState(FluidModel, primitives.GetOutletPrimitiveState());

  /*--- Compute isentropic Outflow quantities ---*/
  FluidModel.SetTDState_Ps(OutletState.GetPressure(), InletState.GetEntropy());
  auto enthalpyOutIs = FluidModel.GetStaticEnergy() + OutletState.GetPressure() / FluidModel.GetDensity();
  auto tangVel = primitives.GetOutletPrimitiveState().GetTangVelocity();
  auto relVelOutIs2 = 2 * (OutletState.GetRothalpy() - enthalpyOutIs) + tangVel * tangVel;

  /*--- Compute performance ---*/
  EntropyGen = (OutletState.GetEntropy() - InletState.GetEntropy()); // / abs(InletState.GetEntropy() + 1);
  EulerianWork = InletState.GetTotalEnthalpy() - OutletState.GetTotalEnthalpy();
  TotalPressureLoss = (InletState.GetTotalRelPressure() - OutletState.GetTotalRelPressure()) /
                      (InletState.GetTotalRelPressure() - OutletState.GetPressure());
  KineticEnergyLoss = 2 * (OutletState.GetEnthalpy() - enthalpyOutIs) / relVelOutIs2;
}

CPropellorBladePerformance::CPropellorBladePerformance(CFluidModel& fluidModel,
                                                        unsigned short nDim,
                                                        su2double areaIn,
                                                        su2double radiusIn,
                                                        su2double areaOut,
                                                        su2double radiusOut) : CTurbomachineryBladePerformance(fluidModel, nDim, areaIn, radiusIn, areaOut, radiusOut){}

void CPropellorBladePerformance::ComputePerformance(const CTurbomachineryCombinedPrimitiveStates &primitives) {

}

CTurbomachineryPerformance::CTurbomachineryPerformance(const CConfig& config,
                                                       const CGeometry& geometry,
                                                       CFluidModel& fluidModel) {
  unsigned short nBladesRow = config.GetnMarker_Turbomachinery();
  unsigned short nDim = geometry.GetnDim();
  unsigned short nStages = SU2_TYPE::Int(nBladesRow / 2);

  for (unsigned short iBladeRow = 0; iBladeRow < nBladesRow; iBladeRow++) {
    vector <shared_ptr<CTurbomachineryBladePerformance>> bladeSpanPerformances;
    auto nSpan = config.GetnSpanWiseSections();
    for (unsigned short iSpan = 0; iSpan < nSpan + 1; iSpan++) {

      auto areaIn = geometry.GetSpanAreaIn(iBladeRow, iSpan);
      auto areaOut = geometry.GetSpanAreaOut(iBladeRow, iSpan);
      auto radiusIn = geometry.GetTurboRadiusIn(iBladeRow, iSpan);
      auto radiusOut = geometry.GetTurboRadiusOut(iBladeRow, iSpan);

      // TODO: I have a feeling this should not be in such a for loop, to be discussed with Salvo (Nitish)
      SU2_OMP_PARALLEL
      {
        const int thread = omp_get_thread_num();

      /* Switch between the Turbomachinery Performance Kind */
      // TODO: This needs to be fixed
      switch (config.GetKind_TurboPerf(iBladeRow)) {
        
        case TURBINE:
          bladeSpanPerformances.push_back(make_shared<CTurbineBladePerformance>(fluidModel, nDim, areaIn, radiusIn, areaOut,
                                                                      radiusOut));
          break;

        case COMPRESSOR:
          bladeSpanPerformances.push_back(make_shared<CCompressorBladePerformance>(fluidModel, nDim, areaIn, radiusIn, areaOut,
                                                                      radiusOut));
          break;

        case PROPELLOR:
          bladeSpanPerformances.push_back(make_shared<CPropellorBladePerformance>(fluidModel, nDim, areaIn, radiusIn, areaOut,
                                                                      radiusOut));
          break;

        default:
          bladeSpanPerformances.push_back(make_shared<CTurbineBladePerformance>(fluidModel, nDim, areaIn, radiusIn, areaOut,
                                                                      radiusOut));
          break;
        }
      }
    }
    BladesPerformances.push_back(bladeSpanPerformances);
  }
}

void CTurbomachineryPerformance::ComputeTurbomachineryPerformance(
  vector <vector<CTurbomachineryCombinedPrimitiveStates>> const bladesPrimitives) {
  boost::for_each(BladesPerformances, bladesPrimitives, boost::bind(&ComputePerBlade, _1, _2));
}

void CTurbomachineryPerformance::ComputePerBlade(
  vector <shared_ptr<CTurbomachineryBladePerformance>> const bladePerformances,
  vector <CTurbomachineryCombinedPrimitiveStates> const bladePrimitives) {
  boost::for_each(bladePerformances, bladePrimitives, boost::bind(&ComputePerSpan, _1, _2));
}

void CTurbomachineryPerformance::ComputePerSpan(shared_ptr <CTurbomachineryBladePerformance> const spanPerformances,
                                                const CTurbomachineryCombinedPrimitiveStates &spanPrimitives) {
  spanPerformances->ComputePerformance(spanPrimitives);
}

CTurbomachineryStagePerformance::CTurbomachineryStagePerformance(CFluidModel& fluid) : fluidModel(fluid) {
  su2double TotalStaticEfficiency=0.0, 
            TotalTotalEfficiency=0.0, 
            KineticEnergyLoss=0.0, 
            TotalPressureLoss=0.0, 
            EntropyGen=0.0, 
            PressureRatio=0.0, 
            EulerianWork=0.0;
  // CFluidModel fluidModel(fluid);
}

void CTurbomachineryStagePerformance::ComputePerformanceStage(CTurbomachineryState InState, CTurbomachineryState OutState, CConfig* config) {

  switch (config->GetKind_TurboPerf(ZONE_0)) {
  case TURBINE:
    /* code */
    ComputeTurbineStagePerformance(InState, OutState);
    break;

  default:
    break;
  }
}

void CTurbomachineryStagePerformance::ComputeTurbineStagePerformance(CTurbomachineryState InState, CTurbomachineryState OutState) {

    /*--- Compute isentropic Outflow quantities ---*/
  fluidModel.SetTDState_Ps(OutState.GetPressure(), InState.GetEntropy());
  auto enthalpyOutIs = fluidModel.GetStaticEnergy() + OutState.GetPressure() / fluidModel.GetDensity();
  auto tangVel = 0.0; //primitives.GetOutletPrimitiveState().GetTangVelocity();
  auto relVelOutIs2 = 2 * (OutState.GetRothalpy() - enthalpyOutIs) + tangVel * tangVel;
  su2double IseTotalEnthalpyOut = enthalpyOutIs + 0.5*OutState.GetVelocityValue()*OutState.GetVelocityValue();

  /*--- Compute performance ---*/
  EntropyGen = (OutState.GetEntropy() - InState.GetEntropy())/InState.GetEntropy() * 100;
  EulerianWork = InState.GetTotalEnthalpy() - OutState.GetTotalEnthalpy();
  TotalPressureLoss = (InState.GetTotalRelPressure() - OutState.GetTotalRelPressure()) /
                      (InState.GetTotalRelPressure() - OutState.GetPressure());
  KineticEnergyLoss = 2 * (OutState.GetEnthalpy() - enthalpyOutIs) / relVelOutIs2;

  TotalStaticEfficiency = EulerianWork/(InState.GetEnthalpy()-enthalpyOutIs) * 100;
  TotalTotalEfficiency =  EulerianWork/(InState.GetEnthalpy()-IseTotalEnthalpyOut) * 100;
  PressureRatio = InState.GetTotalPressure()/OutState.GetTotalPressure();
}

void CTurbomachineryStagePerformance::ComputeCompressorStagePerformance(CTurbomachineryState InState, CTurbomachineryState OutState) {

  /*--- Compute isentropic Outflow quantities ---*/
  fluidModel.SetTDState_Ps(OutState.GetPressure(), InState.GetEntropy());
  auto enthalpyOutIs = fluidModel.GetStaticEnergy() + OutState.GetPressure() / fluidModel.GetDensity();
  auto tangVel = 0.0; //primitives.GetOutletPrimitiveState().GetTangVelocity();
  auto relVelOutIs2 = 2 * (OutState.GetRothalpy() - enthalpyOutIs) + tangVel * tangVel;
  su2double IseTotalEnthalpyOut = enthalpyOutIs + 0.5*OutState.GetVelocityValue()*OutState.GetVelocityValue();

  /*--- Compute performance ---*/
  EntropyGen = (OutState.GetEntropy() - InState.GetEntropy()) / InState.GetEntropy() * 100;
  EulerianWork = OutState.GetTotalEnthalpy() - InState.GetTotalEnthalpy();
  TotalPressureLoss = (InState.GetTotalRelPressure() - OutState.GetTotalRelPressure()) /
                      (InState.GetTotalRelPressure() - OutState.GetPressure());
                      
  KineticEnergyLoss = 2 * (OutState.GetEnthalpy() - enthalpyOutIs) / relVelOutIs2;

  /*--- Efficiency Computation ---*/
  TotalStaticEfficiency = EulerianWork/(InState.GetEnthalpy()-enthalpyOutIs) * 100;
  TotalTotalEfficiency =  EulerianWork/(InState.GetEnthalpy()-IseTotalEnthalpyOut) * 100;

  /*--- Pressure Ratio ---*/
  PressureRatio = OutState.GetTotalPressure()/InState.GetTotalPressure();

}

CTurbomachineryStagePerformance::~CTurbomachineryStagePerformance() {}