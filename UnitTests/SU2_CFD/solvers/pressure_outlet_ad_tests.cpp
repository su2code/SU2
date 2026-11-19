/*!
 * \file pressure_outlet_ad_tests.cpp
 * \brief Differentiation of subsonic, sonic-limited and supersonic pressure outlets.
 * \version 8.5.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md)
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

#include <memory>
#include <tuple>
#include <utility>
#include "catch.hpp"
#include "../../../SU2_CFD/include/fluid/CFluidModel.hpp"
#include "../../../Common/include/toolboxes/geometry_toolbox.hpp"
#include "../../../SU2_CFD/include/fluid/CIdealGas.hpp"
#include "../../../SU2_CFD/include/fluid/CPengRobinson.hpp"
#include "../../../SU2_CFD/include/fluid/CVanDerWaalsGas.hpp"

TEST_CASE("Pressure outlet state derivatives match finite differences", "[Pressure-outlet-AD]") {
  // Check the analytic ideal-gas path and integrated Peng-Robinson/van der Waals paths with the same inputs.
  for (auto kind : {0u, 1u, 2u}) {
    std::unique_ptr<CFluidModel> gas;
    if (kind == 0) gas.reset(new CIdealGas(1.4, 35.152));
    if (kind == 1) gas.reset(new CPengRobinson(1.0125, 35.152, 1415200, 564.1, .529));
    if (kind == 2) gas.reset(new CVanDerWaalsGas(1.0125, 35.152, 1415200, 564.1));
    const double rho = 230, p = 1500112;
    gas->SetTDState_Prho(p, rho);
    const double e = SU2_TYPE::GetValue(gas->GetStaticEnergy()), a = SU2_TYPE::GetValue(gas->GetSoundSpeed());
    const su2double normal[2] = {1, 0};
    // Cover reverse subsonic flow, ordinary outflow, near-sonic flow and supersonic extrapolation;
    // the lowest pressure ratio makes subsonic outward expansion stop at the moving sonic limit.
    for (const auto mach : {-.3, .3, .99, 1.4})
      for (const auto ratio : {1., .95, .001}) {
        const double values[6] = {
            rho, rho * (mach * a + 20), rho * 5, rho * (e + .5 * (pow(mach * a + 20, 2) + 25)), ratio * p, 20};
        // The first four inputs are conservative variables; the last two are pressure and grid speed.
        auto evaluate = [&](su2double* input, su2double* output) {
          return gas->ComputePressureOutletState(2, input, normal, input[5], input[4], output, 0, kind == 0 ? 1.4 : 0.);
        };
        for (auto column = 0u; column < 6; ++column) {
          const double scale[6] = {rho, rho * a, rho * a, rho * a * a, p, a};
          const double delta = 2e-5 * fmax(fabs(values[column]), scale[column]);
          su2double plus[6], minus[6], fp[4], fm[4];
          for (auto k = 0u; k < 6; ++k) plus[k] = minus[k] = values[k];
          plus[column] += delta;
          minus[column] -= delta;
          REQUIRE(evaluate(plus, fp) == CFluidModel::PressureOutletStatus::Valid);
          REQUIRE(evaluate(minus, fm) == CFluidModel::PressureOutletStatus::Valid);
          // Finite differences rerun the complete map, so the reference includes movement of the sonic root.
          for (auto row = 0u; row < 4; ++row) {
            const double expected = SU2_TYPE::GetValue((fp[row] - fm[row]) / (2 * delta));
            su2double input[6], output[4];
            for (auto k = 0u; k < 6; ++k) input[k] = values[k];
#ifdef CODI_REVERSE_TYPE
            // Record one scalar input/output pair per tape to keep sensitivities independent.
            AD::Reset();
            AD::StartRecording();
            AD::RegisterInput(input[column]);
#else
            // Forward mode propagates a unit seed for the selected input column.
            SU2_TYPE::SetDerivative(input[column], 1);
#endif
            REQUIRE(evaluate(input, output) == CFluidModel::PressureOutletStatus::Valid);
#ifdef CODI_REVERSE_TYPE
            AD::RegisterOutput(output[row]);
            AD::StopRecording();
            SU2_TYPE::SetDerivative(output[row], 1);
            AD::ComputeAdjoint();
            const double derivative = SU2_TYPE::GetDerivative(input[column]);
            AD::Reset();
#else
            const double derivative = SU2_TYPE::GetDerivative(output[row]);
#endif
            INFO("gas=" << kind << " Mach=" << mach << " pressure ratio=" << ratio << " row=" << row
                        << " column=" << column);
            CHECK(derivative == Approx(expected).epsilon(3e-4).margin(2e-4));
          }
        }
      }
  }
}

namespace {
// Compare the derivative of every output with central finite differences of the complete map, one input at a time.
template <class Map>
void CheckMapDerivatives(const Map& evaluate, const double* values, const double* scale) {
  for (auto column = 0u; column < 6; ++column) {
    const double delta = 2e-5 * fmax(fabs(values[column]), scale[column]);
    su2double plus[6], minus[6], fp[4], fm[4];
    for (auto k = 0u; k < 6; ++k) plus[k] = minus[k] = values[k];
    plus[column] += delta;
    minus[column] -= delta;
    REQUIRE(evaluate(plus, fp) == CFluidModel::PressureOutletStatus::Valid);
    REQUIRE(evaluate(minus, fm) == CFluidModel::PressureOutletStatus::Valid);
    for (auto row = 0u; row < 4; ++row) {
      const double expected = SU2_TYPE::GetValue((fp[row] - fm[row]) / (2 * delta));
      su2double input[6], output[4];
      for (auto k = 0u; k < 6; ++k) input[k] = values[k];
#ifdef CODI_REVERSE_TYPE
      AD::Reset();
      AD::StartRecording();
      AD::RegisterInput(input[column]);
#else
      SU2_TYPE::SetDerivative(input[column], 1);
#endif
      REQUIRE(evaluate(input, output) == CFluidModel::PressureOutletStatus::Valid);
#ifdef CODI_REVERSE_TYPE
      AD::RegisterOutput(output[row]);
      AD::StopRecording();
      SU2_TYPE::SetDerivative(output[row], 1);
      AD::ComputeAdjoint();
      const double derivative = SU2_TYPE::GetDerivative(input[column]);
      AD::Reset();
#else
      const double derivative = SU2_TYPE::GetDerivative(output[row]);
#endif
      INFO("row=" << row << " column=" << column);
      // The margin also covers the finite-difference noise of the adaptive EOS integration (about 1e-10 relative).
      const double noise = 1e-9 * fmax(fabs(SU2_TYPE::GetValue(fp[row])), fabs(SU2_TYPE::GetValue(fm[row]))) / delta;
      CHECK(derivative == Approx(expected).epsilon(3e-4).margin(2e-4 + noise));
    }
  }
}
}  // namespace

TEST_CASE("Outlet backflow state derivatives match finite differences", "[Pressure-outlet-AD]") {
  for (auto kind : {0u, 1u, 2u}) {
    std::unique_ptr<CFluidModel> gas;
    if (kind == 0) gas.reset(new CIdealGas(1.4, 35.152));
    if (kind == 1) gas.reset(new CPengRobinson(1.0125, 35.152, 1415200, 564.1, .529));
    if (kind == 2) gas.reset(new CVanDerWaalsGas(1.0125, 35.152, 1415200, 564.1));
    const double rho = 230, p = 1500112;
    gas->SetTDState_Prho(p, rho);
    const double e = SU2_TYPE::GetValue(gas->GetStaticEnergy()), a = SU2_TYPE::GetValue(gas->GetSoundSpeed());
    const su2double temperature = SU2_TYPE::GetValue(gas->GetTemperature()), normal[2] = {1, 0};
    const double gamma = kind == 0 ? 1.4 : 0;
    // Inputs: four conservative variables, then reservoir total pressure (or imposed pressure) and grid speed.
    const double scale[6] = {rho, rho * a, rho * a, rho * a * a, p, a};
    auto reservoir = [&](su2double* input, su2double* output) {
      return gas->ComputeReservoirInflowState(2, input, normal, input[5], input[4], temperature, nullptr, output, 0,
                                              gamma);
    };
    auto sonic = [&](su2double* input, su2double* output) {
      return gas->ComputeSonicBackflowState(2, input, normal, input[5], input[4], output, 0, gamma);
    };
    // Subsonic reservoir inflow (the root carries the derivatives), the sonic reservoir limit and the reservoir at
    // rest (both independent of the interior), then the sonic fallback without data.
    // Each case gives the normal Mach number relative to the grid, the pressure ratio and the grid speed.
    for (const auto& [mach, ratio, grid] :
         {std::tuple<double, double, double>{-.3, 1.05, -20}, {-.05, 1.05, -20}, {-1.3, 1.05, -20}, {-.02, .9, 0}}) {
      CAPTURE(kind, mach, ratio, grid);
      const double u = mach * a + grid;
      const double values[6] = {rho, rho * u, rho * 5, rho * (e + .5 * (u * u + 25)), ratio * p, grid};
      CheckMapDerivatives(reservoir, values, scale);
    }
    // Inflow with the interior tangential velocity: explicit total pressure, and zero total pressure, where input[4]
    // is the imposed pressure that the reservoir has at zero normal velocity.
    auto tangentialReservoir = [&](su2double* input, su2double* output, bool zeroTotal) {
      const su2double tangential[2] = {0, input[2] / input[0]};
      return gas->ComputeReservoirInflowState(2, input, normal, input[5], zeroTotal ? su2double(0) : input[4],
                                              temperature, nullptr, output, 0, gamma, tangential, input[4]);
    };
    for (const auto zeroTotal : {false, true})
      for (const auto mach : {-.3, -.05}) {
        CAPTURE(kind, zeroTotal, mach);
        const double grid = -20, u = mach * a + grid, ratio = zeroTotal ? 1. : 1.05;
        const double values[6] = {rho, rho * u, rho * .3 * a, rho * (e + .5 * (u * u + .09 * a * a)), ratio * p, grid};
        CheckMapDerivatives([&](su2double* in, su2double* out) { return tangentialReservoir(in, out, zeroTotal); },
                            values, scale);
      }
    for (const auto ratio : {.95, 1.1}) {
      CAPTURE(kind, ratio);
      const double grid = 20, u = -1.4 * a + grid;
      const double values[6] = {rho, rho * u, rho * 5, rho * (e + .5 * (u * u + 25)), ratio * p, grid};
      CheckMapDerivatives(sonic, values, scale);
    }
  }
}
