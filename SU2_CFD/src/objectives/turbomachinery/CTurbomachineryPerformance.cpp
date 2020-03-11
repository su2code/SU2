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

CTurbomachineryPrimitiveState::CTurbomachineryPrimitiveState(){
  Density = Pressure = TangVelocity = 0.0;
}
CTurbomachineryPrimitiveState::CTurbomachineryPrimitiveState(su2double density, su2double pressure, su2double *velocity, unsigned short nDim, su2double tangVel) : 
                              Density(density), Pressure(pressure), TangVelocity(tangVel){
  Velocity.assign(velocity, velocity + nDim );
}

CTurbomachineryPrimitiveState::~CTurbomachineryPrimitiveState(){}

CTurbomachineryCombinedPrimitiveStates::CTurbomachineryCombinedPrimitiveStates(const CTurbomachineryPrimitiveState& inletPrimitiveState, 
                                        const CTurbomachineryPrimitiveState&  outletPrimitiveState) : InletPrimitiveState(inletPrimitiveState), OutletPrimitiveState(outletPrimitiveState){}

CTurbomachineryCombinedPrimitiveStates::~CTurbomachineryCombinedPrimitiveStates(){}

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

CTurbomachineryState::~CTurbomachineryState(){}


void CTurbomachineryState::ComputeState(CFluidModel *const fluidModel, const CTurbomachineryPrimitiveState &primitiveState) {

  /*--- Assign new primitive values ---*/
  Density = primitiveState.GetDensity();
  Pressure = primitiveState.GetPressure();
  auto velocity = primitiveState.GetVelocity();
  Velocity.assign(velocity.begin(), velocity.end());
  auto tangVel = primitiveState.GetTangVelocity();

  /*--- Compute static TD quantities ---*/
  fluidModel->SetTDState_Prho(Pressure, Density);
  Entropy = fluidModel->GetEntropy();
  Enthalpy = fluidModel->GetStaticEnergy() + Pressure / Density;
  su2double soundSpeed = fluidModel->GetSoundSpeed();

  /*--- Compute total TD quantities ---*/
  TotalEnthalpy = Enthalpy + 0.5 * GetVelocityValue() * GetVelocityValue();
  fluidModel->SetTDState_hs(TotalEnthalpy, Entropy);
  TotalPressure = fluidModel->GetPressure();
  TotalTemperature = fluidModel->GetTemperature();

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
  fluidModel->SetTDState_hs(Rothalpy, Entropy);
  TotalRelPressure = fluidModel->GetPressure();

}

CTurbomachineryBladePerformance::CTurbomachineryBladePerformance(CFluidModel* const fluidModel, unsigned short nDim, su2double areaIn, su2double radiusIn, su2double areaOut, su2double radiusOut) : FluidModel(fluidModel){
   InletState               = CTurbomachineryState(nDim, areaIn, radiusIn);
   OutletState              = CTurbomachineryState(nDim, areaOut,radiusOut);
}

CTurbomachineryBladePerformance::~CTurbomachineryBladePerformance(){}

CTurbineBladePerformance::CTurbineBladePerformance(CFluidModel* const fluidModel, unsigned short nDim, su2double areaIn, su2double radiusIn, su2double areaOut, su2double radiusOut) : CTurbomachineryBladePerformance(fluidModel, nDim, areaIn, radiusIn, areaOut, radiusOut){}

void CTurbineBladePerformance::ComputePerformance(const CTurbomachineryCombinedPrimitiveStates &primitives) {

  /*--- Compute Inlet and Outlet state ---*/
  InletState.ComputeState(FluidModel, primitives.GetInletPrimitiveState());
  OutletState.ComputeState(FluidModel, primitives.GetOutletPrimitiveState());

  /*--- Compute isentropic Outflow quantities ---*/
  FluidModel->SetTDState_Ps(OutletState.GetPressure(), InletState.GetEntropy());
  auto enthalpyOutIs = FluidModel->GetStaticEnergy() + OutletState.GetPressure() / FluidModel->GetDensity();
  auto tangVel = primitives.GetOutletPrimitiveState().GetTangVelocity();
  auto relVelOutIs2 = 2 * (OutletState.GetRothalpy() - enthalpyOutIs) + tangVel * tangVel;

  /*--- Compute performance ---*/
  EntropyGen = (OutletState.GetEntropy() - InletState.GetEntropy()) / abs(InletState.GetEntropy() + 1);
  EulerianWork = InletState.GetTotalEnthalpy() - OutletState.GetTotalEnthalpy();
  TotalPressureLoss = (InletState.GetTotalRelPressure() - OutletState.GetTotalRelPressure()) /
                      (InletState.GetTotalRelPressure() - OutletState.GetPressure());
  KineticEnergyLoss = 2 * (OutletState.GetEnthalpy() - enthalpyOutIs) / relVelOutIs2;
}

CTurbineBladePerformance::~CTurbineBladePerformance() {}
// CCompressorBladePerformance::CCompressorBladePerformance(unsigned short nDim, su2double area, su2double radius) : CTurbomachineryBladePerformance(nDim, area, radius){}

CTurbomachineryPerformance::CTurbomachineryPerformance(CConfig *const config, CGeometry *const geometry,
                                                       CFluidModel *const fluidModel) {
  unsigned short nBladesRow = config->GetnMarker_Turbomachinery();
  unsigned short nDim = geometry->GetnDim();
  unsigned short nStages = SU2_TYPE::Int(nBladesRow / 2);

  for (unsigned short iBladeRow = 0; iBladeRow < nBladesRow; iBladeRow++) {
    vector <shared_ptr<CTurbomachineryBladePerformance>> bladeSpanPerformances;
    auto nSpan = config->GetnSpan_iZones(iBladeRow);
    for (unsigned short iSpan = 0; iSpan < nSpan + 1; iSpan++) {

      auto areaIn = geometry->GetSpanAreaIn(iBladeRow, iSpan);
      auto areaOut = geometry->GetSpanAreaOut(iBladeRow, iSpan);
      auto radiusIn = geometry->GetTurboRadiusIn(iBladeRow, iSpan);
      auto radiusOut = geometry->GetTurboRadiusOut(iBladeRow, iSpan);
      auto bladePerformance = make_shared<CTurbineBladePerformance>(fluidModel, nDim, areaIn, radiusIn, areaOut,
                                                                    radiusOut);;
      bladeSpanPerformances.push_back(bladePerformance);
    }
    BladesPerformances.push_back(bladeSpanPerformances);
  }
}

CTurbomachineryPerformance::~CTurbomachineryPerformance() {}

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