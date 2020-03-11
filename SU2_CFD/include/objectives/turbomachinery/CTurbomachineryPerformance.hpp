/*!
 * \file CTurbomachineryPerformance.hpp
 * \brief Headers of the Turbomachinery Performance class.
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


#pragma once

#include <vector>
#include <memory>
#include <algorithm>
#include <iterator>
#include <numeric>
#include <boost/bind.hpp>
#include <boost/range/algorithm_ext/for_each.hpp>
#include "../../../../Common/include/mpi_structure.hpp"
#include "../../../../Common/include/geometry/CGeometry.hpp"
#include "../../../../Common/include/CConfig.hpp"
#include "../../fluid_model.hpp"

using namespace std;

class CTurbomachineryPrimitiveState {
  private:
  su2double Density, Pressure, TangVelocity;
  vector <su2double> Velocity;
  public:
  CTurbomachineryPrimitiveState();

  CTurbomachineryPrimitiveState(su2double density, su2double pressure, su2double *velocity, unsigned short nDim,
                                su2double tangVel);

  ~CTurbomachineryPrimitiveState();

  su2double GetDensity() const { return Density; }

  su2double GetPressure() const { return Pressure; }

  su2double GetTangVelocity() const { return TangVelocity; }

  vector <su2double> GetVelocity() const { return Velocity; }

};

class CTurbomachineryCombinedPrimitiveStates {
  private:
  CTurbomachineryPrimitiveState InletPrimitiveState;
  CTurbomachineryPrimitiveState OutletPrimitiveState;
  public:
  CTurbomachineryCombinedPrimitiveStates();

  CTurbomachineryCombinedPrimitiveStates(const CTurbomachineryPrimitiveState &inletPrimitiveState,
                                         const CTurbomachineryPrimitiveState &outletPrimitiveState);

  ~CTurbomachineryCombinedPrimitiveStates();

  CTurbomachineryPrimitiveState GetInletPrimitiveState() const { return InletPrimitiveState; }

  CTurbomachineryPrimitiveState GetOutletPrimitiveState() const { return OutletPrimitiveState; }
};

class CTurbomachineryState {
  private:
  su2double Density, Pressure, Entropy, Enthalpy, Temperature, TotalTemperature, TotalPressure, TotalEnthalpy;
  su2double AbsFlowAngle, FlowAngle, MassFlow, Rothalpy, TotalRelPressure;
  vector <su2double> Velocity, RelVelocity, Mach, RelMach;
  su2double Area, Radius;


  public:
  CTurbomachineryState();

  CTurbomachineryState(unsigned short nDim, su2double area, su2double radius);

  ~CTurbomachineryState();

  void ComputeState(CFluidModel *const fluidModel, const CTurbomachineryPrimitiveState &primitiveState);

  su2double GetDensity() const { return Density; }

  su2double GetPressure() const { return Pressure; }

  su2double GetEntropy() const { return Entropy; }

  su2double GetEnthalpy() const { return Enthalpy; }

  su2double GetTemperature() const { return Temperature; }

  su2double GetTotalTemerature() const { return TotalTemperature; }

  su2double GetTotalPressure() const { return TotalPressure; }

  su2double GetTotalRelPressure() const { return TotalRelPressure; }

  su2double GetTotalEnthalpy() const { return TotalEnthalpy; }

  su2double GetAbsFlowAngle() const { return AbsFlowAngle; }

  su2double GetFlowAngle() const { return FlowAngle; }

  su2double GetMassFlow() const { return MassFlow; }

  su2double GetRothalpy() const { return Rothalpy; }

  vector <su2double> GetVelocity() const { return Velocity; }

  vector <su2double> GetMach() const { return Mach; }

  su2double GetVelocityValue() const {
    return sqrt(inner_product(Velocity.begin(), Velocity.end(), Velocity.begin(), 0.0));
  }

  su2double GetMachValue() const { return sqrt(inner_product(Mach.begin(), Mach.end(), Mach.begin(), 0.0)); }

  su2double GetRelVelocityValue() const {
    return sqrt(inner_product(RelVelocity.begin(), RelVelocity.end(), RelVelocity.begin(), 0.0));
  }

  su2double GetRelMachValue() const {
    return sqrt(inner_product(RelMach.begin(), RelMach.end(), RelMach.begin(), 0.0));
  }
};


class CTurbomachineryBladePerformance {
  protected:
  CTurbomachineryState InletState;
  CTurbomachineryState OutletState;
  su2double KineticEnergyLoss, TotalPressureLoss, EntropyGen, PressureRatio, EulerianWork;
  CFluidModel *FluidModel;

  public:
  CTurbomachineryBladePerformance(CFluidModel *const fluidModel, unsigned short nDim, su2double areaIn,
                                  su2double radiusIn, su2double areaOut, su2double radiusOut);

  virtual ~CTurbomachineryBladePerformance();

  virtual void ComputePerformance(const CTurbomachineryCombinedPrimitiveStates &primitives) {};

  CTurbomachineryState GetInletState() const { return InletState; }

  CTurbomachineryState GetOutletState() const { return OutletState; }

  su2double GetKineticEnergyLoss() const { return KineticEnergyLoss; }

  su2double GetTotalPressureLoss() const { return TotalPressureLoss; }

  su2double GetEntropyGen() const { return EntropyGen; }

  su2double GetPressureRatio() const { return PressureRatio; }

  su2double GetEulerianWork() const { return EulerianWork; }

};

class CTurbineBladePerformance : public CTurbomachineryBladePerformance {

  public:
  CTurbineBladePerformance(CFluidModel *const fluidModel, unsigned short nDim, su2double areaIn, su2double radiusIn,
                           su2double areaOut, su2double radiusOut);

  ~CTurbineBladePerformance();

  void ComputePerformance(const CTurbomachineryCombinedPrimitiveStates &primitives) override;

};

class CCompressorBladePerformance : public CTurbomachineryBladePerformance {

  public:
  CCompressorBladePerformance(CFluidModel *const fluidModel, unsigned short nDim, su2double areaIn, su2double radiusIn,
                              su2double areaOut, su2double radiusOut);

  ~CCompressorBladePerformance();

  void ComputePerformance(const CTurbomachineryCombinedPrimitiveStates &primitives) override;

};
// class CTurbomachineryStagePerformance :  CTurbomachineryBladePerformance {
//    protected: 
//       su2double TotalStaticEfficiency, TotalTotalEfficiency, KineticEnergyLoss, TotalPressureLoss, EntropyGen, PressureRatio, EulerianWork;
//    public:
//       CTurbomachineryStagePerformance();
//       ~CTurbomachineryStagePerformance();
//       virtual void ComputePerformance(CFluidModel* fluidModel, su2double tangVel);
//       su2double GetTotalStaticEfficiency() const { return TotalStaticEfficiency; }
//       su2double GetTotalTotalEfficiency() const { return TotalTotalEfficiency; }
// };

// class CTurbineStagePerformance : CTurbomachineryStagePerformance {

//    public:
//       CTurbineStagePerformance();
//       ~CTurbineStagePerformance();
//       void ComputePerformance(CFluidModel* fluidModel, su2double tangVel) override;

// };

// class CCompressorStagePerformance : CTurbomachineryStagePerformance {

//    public:
//       CCompressorStagePerformance();
//       ~CCompressorStagePerformance();
//       void ComputePerformance(CFluidModel* fluidModel, su2double tangVel) override;

// };

class CTurbomachineryPerformance {
  private:
  vector <vector<shared_ptr < CTurbomachineryBladePerformance>>>
  BladesPerformances;

  static void ComputePerBlade(vector <shared_ptr<CTurbomachineryBladePerformance>> const bladePerformances,
                              vector <CTurbomachineryCombinedPrimitiveStates> const bladePrimitives);

  static void ComputePerSpan(shared_ptr <CTurbomachineryBladePerformance> const spanPerformances,
                             const CTurbomachineryCombinedPrimitiveStates &spanPrimitives);
  // vector<shared_ptr<CTurbomachineryStagePerformance>> StagePerformances;
  // shared_ptr<CTurbomachineryStagePerformance> MachinePerformances;
  public:
  CTurbomachineryPerformance(CConfig *const config, CGeometry *const geometry, CFluidModel *const fluidModel);

  ~CTurbomachineryPerformance();

  vector <vector<shared_ptr < CTurbomachineryBladePerformance>>>

  GetBladesPerformances() const { return BladesPerformances; }

  // vector<shared_ptr<CTurbomachineryStagePerformance>> GetStagePerformances() const { return StagePerformances; }
  // shared_ptr<CTurbomachineryStagePerformance> GetMachinePerformances() const { return MachinePerformances; }
  void ComputeTurbomachineryPerformance(vector <vector<CTurbomachineryCombinedPrimitiveStates>> const primitives);
};