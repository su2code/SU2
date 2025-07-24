/*!
 * \file CTurboOutput.hpp
 * \brief Headers of the Turbomachinery Performance class.
 * \author S. Vitale, N. Anand
 * \version 7.5.1 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2025, SU2 Contributors (cf. AUTHORS.md)
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

#pragma once

#include <vector>
#include <memory>
#include <algorithm>
#include <iterator>
#include <numeric>
#include "../../../Common/include/geometry/CGeometry.hpp"
#include "../../../Common/include/CConfig.hpp"
#include "../fluid/CFluidModel.hpp"

/*!
 * \brief Class containing the required primitive variables for initiating a turboperformance calculation
 * \author S. Vitale, N. Anand, J. Kelly
 * \ingroup Output
 */
class CTurbomachineryPrimitiveState {
 private:
  su2double Density, Pressure, TangVelocity;
  vector<su2double> Velocity;

 public:
  CTurbomachineryPrimitiveState();

  CTurbomachineryPrimitiveState(vector<su2double> TurboPrimitives, unsigned short nDim, su2double tangVel);

  const su2double& GetDensity() const { return Density; }

  const su2double& GetPressure() const { return Pressure; }

  const su2double& GetTangVelocity() const { return TangVelocity; }

  const std::vector<su2double>& GetVelocity() const { return Velocity; }
};

/*!
 * \brief Class containing the combined primitive inlet and outlet states for a given blade
 * \author S. Vitale, N. Anand, J. Kelly
 * \ingroup Output
 */
class CTurbomachineryCombinedPrimitiveStates {
 private:
  CTurbomachineryPrimitiveState InletPrimitiveState;
  CTurbomachineryPrimitiveState OutletPrimitiveState;

 public:
  CTurbomachineryCombinedPrimitiveStates(const CTurbomachineryPrimitiveState& inletPrimitiveState,
                                         const CTurbomachineryPrimitiveState& outletPrimitiveState);

  const CTurbomachineryPrimitiveState& GetInletPrimitiveState() const& { return InletPrimitiveState; }

  const CTurbomachineryPrimitiveState& GetOutletPrimitiveState() const& { return OutletPrimitiveState; }
};

/*!
 * \brief Class containing state information for a turbomachine
 * \author S. Vitale, N. Anand, J. Kelly
 * \ingroup Output
 */
class CTurbomachineryState {
 private:
  su2double Density, Pressure, Entropy, Enthalpy, Temperature, TotalTemperature, TotalPressure, TotalEnthalpy;
  su2double AbsFlowAngle, FlowAngle, MassFlow, Rothalpy, TangVelocity, TotalRelPressure;
  vector<su2double> Velocity, RelVelocity, Mach, RelMach;
  su2double Area, Radius;

 public:
  CTurbomachineryState();

  CTurbomachineryState(unsigned short nDim, su2double area, su2double radius);

  void ComputeState(CFluidModel& fluidModel, const CTurbomachineryPrimitiveState& primitiveState);

  const su2double& GetDensity() const { return Density; }

  const su2double& GetPressure() const { return Pressure; }

  const su2double& GetEntropy() const { return Entropy; }

  const su2double& GetEnthalpy() const { return Enthalpy; }

  const su2double& GetTemperature() const { return Temperature; }

  const su2double& GetTotalTemperature() const { return TotalTemperature; }

  const su2double& GetTotalPressure() const { return TotalPressure; }

  const su2double& GetTotalRelPressure() const { return TotalRelPressure; }

  const su2double& GetTotalEnthalpy() const { return TotalEnthalpy; }

  const su2double& GetAbsFlowAngle() const { return AbsFlowAngle; }

  const su2double& GetFlowAngle() const { return FlowAngle; }

  const su2double& GetMassFlow() const { return MassFlow; }

  const su2double& GetRothalpy() const { return Rothalpy; }

  const su2double& GetTangVelocity() const { return TangVelocity; }

  const vector<su2double>& GetVelocity() const { return Velocity; }

  const vector<su2double>& GetMach() const { return Mach; }

  su2double GetVelocityValue() const { return Norm(Velocity); }

  su2double GetMachValue() const { return Norm(Mach); }

  su2double GetRelVelocityValue() const { return Norm(RelVelocity); }

  su2double GetRelMachValue() const { return Norm(RelMach); }

  su2double Norm(const vector<su2double>& u) const {
    su2double accum = 0.;
    for (auto i = 0u; i < u.size(); ++i) {
      accum += u[i] * u[i];
    }
    return sqrt(accum);
  }
};

/*!
 * \brief Class containing additional states and performance calculation routines for blades in different turbomachines
 * \author S. Vitale, N. Anand, J. Kelly
 * \ingroup Output
 */
class CTurbomachineryBladePerformance {
 protected:
  CTurbomachineryState InletState;
  CTurbomachineryState OutletState;
  su2double KineticEnergyLoss, TotalPressureLoss, EntropyGen, PressureRatio, EulerianWork;
  CFluidModel& FluidModel;

 public:
  CTurbomachineryBladePerformance(CFluidModel& fluidModel, unsigned short nDim, su2double areaIn, su2double radiusIn, su2double areaOut, su2double radiusOut);

  virtual ~CTurbomachineryBladePerformance() = default;

  virtual void ComputePerformance(const CTurbomachineryCombinedPrimitiveStates& primitives){};

  const CTurbomachineryState& GetInletState() { return InletState; }

  const CTurbomachineryState& GetOutletState() { return OutletState; }

  const su2double& GetKineticEnergyLoss() const { return KineticEnergyLoss; }

  const su2double& GetTotalPressureLoss() const { return TotalPressureLoss; }

  const su2double& GetEntropyGen() const { return EntropyGen; }

  const su2double& GetPressureRatio() const { return PressureRatio; }

  const su2double& GetEulerianWork() const { return EulerianWork; }
};

class CTurbineBladePerformance : public CTurbomachineryBladePerformance {
 public:
  CTurbineBladePerformance(CFluidModel& fluidModel, unsigned short nDim, su2double areaIn, su2double radiusIn, su2double areaOut, su2double radiusOut);

  void ComputePerformance(const CTurbomachineryCombinedPrimitiveStates& primitives) override;
};

class CCompressorBladePerformance : public CTurbomachineryBladePerformance {
 public:
  CCompressorBladePerformance(CFluidModel& fluidModel, unsigned short nDim, su2double areaIn, su2double radiusIn, su2double areaOut, su2double radiusOut);

  void ComputePerformance(const CTurbomachineryCombinedPrimitiveStates& primitives) override;
};

class CPropellorBladePerformance : public CTurbomachineryBladePerformance {
 public:
  CPropellorBladePerformance(CFluidModel& fluidModel, unsigned short nDim, su2double areaIn, su2double radiusIn, su2double areaOut, su2double radiusOut);

  void ComputePerformance(const CTurbomachineryCombinedPrimitiveStates& primitives) override;
};

/*!
 * \brief Class for computng full stage performance
 * \author S. Vitale, N. Anand, J. Kelly
 * \ingroup Output
 */
class CTurbomachineryStagePerformance {
 protected:
  su2double TotalStaticEfficiency, TotalTotalEfficiency, NormEntropyGen, TotalStaticPressureRatio, TotalTotalPressureRatio, EulerianWork, TotalPressureLoss, KineticEnergyLoss;
  CFluidModel& fluidModel;

 public:
  CTurbomachineryStagePerformance(CFluidModel& fluid);

  virtual ~CTurbomachineryStagePerformance() = default;

  virtual void ComputePerformanceStage(const CTurbomachineryState& InState, const CTurbomachineryState& OutState, const CConfig* config);

  virtual void ComputeTurbineStagePerformance(const CTurbomachineryState& InState, const CTurbomachineryState& OutState);

  virtual void ComputeCompressorStagePerformance(const CTurbomachineryState& InState, const CTurbomachineryState& OutState);

  su2double GetTotalStaticEfficiency() const { return TotalStaticEfficiency; }

  su2double GetTotalTotalEfficiency() const { return TotalTotalEfficiency; }

  su2double GetEulerianWork() const { return EulerianWork; }

  su2double GetNormEntropyGen() const { return NormEntropyGen; }

  su2double GetTotalStaticPressureRatio() const { return TotalStaticPressureRatio; }

  su2double GetTotalTotalPressureRatio() const { return TotalTotalPressureRatio; }

  su2double GetTotalPressureLoss() const { return TotalPressureLoss; }

  su2double GetKineticEnergyLoss() const { return KineticEnergyLoss; }
};

/*!
 * \brief Class for handling the calculation of turboperformance variables across a blade, span and full machine
 * \author S. Vitale, N. Anand, J. Kelly
 * \ingroup Output
 */
class CTurboOutput {
 private:
  vector<vector<shared_ptr<CTurbomachineryBladePerformance>>> BladesPerformances;

  static void ComputePerBlade(vector<shared_ptr<CTurbomachineryBladePerformance>> const bladePerformances, vector<CTurbomachineryCombinedPrimitiveStates> const bladePrimitives);

  static void ComputePerSpan(shared_ptr<CTurbomachineryBladePerformance> const spanPerformances, const CTurbomachineryCombinedPrimitiveStates& spanPrimitives);
  
 public:
  CTurboOutput(CConfig** config, const CGeometry& geometry, CFluidModel& fluidModel);

  const vector<vector<shared_ptr<CTurbomachineryBladePerformance>>>& GetBladesPerformances() const { return BladesPerformances; }

  void ComputeTurbomachineryPerformance(vector<vector<CTurbomachineryCombinedPrimitiveStates>> const primitives);
};