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

#include <tuple>
#include "catch.hpp"
#include "pressure_outlet_test_utils.hpp"

namespace {
using Status = CEulerSolver::PressureOutletStatus;

// Compare the derivative of every output with central finite differences of the complete map, one input at a time.
// Inputs: four conservative variables, then the imposed pressure and the grid speed along the normal.
// The margin can include the relative noise of a map (e.g. of nested EOS integrations) divided by the step.
template <class Map>
void CheckMapDerivatives(const Map& evaluate, const double* values, const double* scale, double relativeNoise = 0) {
  for (auto column = 0u; column < 6; ++column) {
    const double delta = 2e-5 * fmax(fabs(values[column]), scale[column]);
    su2double plus[6], minus[6], fp[4], fm[4];
    for (auto k = 0u; k < 6; ++k) plus[k] = minus[k] = values[k];
    plus[column] += delta;
    minus[column] -= delta;
    REQUIRE(evaluate(plus, fp) == Status::Valid);
    REQUIRE(evaluate(minus, fm) == Status::Valid);
    // Finite differences rerun the complete map, so the reference includes the movement of the sonic point.
    for (auto row = 0u; row < 4; ++row) {
      const double expected = SU2_TYPE::GetValue((fp[row] - fm[row]) / (2 * delta));
      su2double input[6], output[4];
      for (auto k = 0u; k < 6; ++k) input[k] = values[k];
#ifdef CODI_REVERSE_TYPE
      // One scalar input/output pair per tape keeps the sensitivities independent.
      AD::Reset();
      AD::StartRecording();
      AD::RegisterInput(input[column]);
#else
      SU2_TYPE::SetDerivative(input[column], 1);
#endif
      REQUIRE(evaluate(input, output) == Status::Valid);
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
      const double noise =
          relativeNoise * fmax(fabs(SU2_TYPE::GetValue(fp[row])), fabs(SU2_TYPE::GetValue(fm[row]))) / delta;
      CHECK(derivative == Approx(expected).epsilon(3e-4).margin(2e-4 + noise));
    }
  }
}
}  // namespace

TEST_CASE("Pressure outlet state derivatives match finite differences", "[Pressure-outlet-AD]") {
  // The closed form of the ideal gas and the integration for Peng-Robinson and van der Waals, with the same inputs.
  for (auto kind : {0u, 1u, 2u}) {
    auto gas = MakeOutletTestGas(kind);
    const double rho = 230, p = 1500112;
    gas->SetTDState_Prho(p, rho);
    const double e = SU2_TYPE::GetValue(gas->GetStaticEnergy()), a = SU2_TYPE::GetValue(gas->GetSoundSpeed());
    const su2double normal[2] = {1, 0};
    const double scale[6] = {rho, rho * a, rho * a, rho * a * a, p, a};
    auto outlet = [&](su2double* input, su2double* output) {
      return CEulerSolver::ComputePressureOutletState(gas.get(), 2, input, normal, input[5], input[4], 0, kind == 0,
                                                      1.4, output);
    };
    auto sonic = [&](su2double* input, su2double* output) {
      return CEulerSolver::ComputeSonicBackflowState(gas.get(), 2, input, normal, input[5], input[4], 0, kind == 0, 1.4,
                                                     output);
    };
    // Subsonic reverse flow, outflow, near-sonic flow and supersonic outflow; the lowest pressure ratio stops the
    // expansion at the sonic point.
    for (const auto mach : {-.3, .3, .99, 1.4})
      for (const auto ratio : {1., .995, .95, .001}) {
        CAPTURE(kind, mach, ratio);
        const double u = mach * a + 20;
        const double values[6] = {rho, rho * u, rho * 5, rho * (e + .5 * (u * u + 25)), ratio * p, 20};
        CheckMapDerivatives(outlet, values, scale);
      }
    // Sonic backflow fallback of steady runs, with compression and expansion.
    for (const auto ratio : {.95, 1.1}) {
      CAPTURE(kind, ratio);
      const double u = -1.4 * a + 20;
      const double values[6] = {rho, rho * u, rho * 5, rho * (e + .5 * (u * u + 25)), ratio * p, 20};
      CheckMapDerivatives(sonic, values, scale);
    }
  }
}

TEST_CASE("Reservoir backflow state derivatives match finite differences", "[Pressure-outlet-AD]") {
  for (auto kind : {0u, 1u, 2u}) {
    auto gas = MakeOutletTestGas(kind);
    const double rho = 230, p = 1500112;
    gas->SetTDState_Prho(p, rho);
    const double e = SU2_TYPE::GetValue(gas->GetStaticEnergy()), a = SU2_TYPE::GetValue(gas->GetSoundSpeed());
    const double temperature = SU2_TYPE::GetValue(gas->GetTemperature());
    const su2double normal[2] = {1, 0};
    const double scale[6] = {rho, rho * a, rho * a, rho * a * a, p, a};
    // input[4] is the reservoir total pressure, or with a zero total pressure the imposed pressure. The solver
    // precomputes the total state of a positive total pressure (SetReservoirTotalState).
    auto reservoirMap = [&](su2double* input, su2double* output, bool zeroTotal, bool tangential, bool precomputed) {
      CEulerSolver::PressureOutletReservoir reservoir;
      reservoir.totalPressure = zeroTotal ? su2double(0) : input[4];
      reservoir.totalTemperature = temperature;
      reservoir.interiorTangential = tangential;
      if (precomputed) REQUIRE(CEulerSolver::SetReservoirTotalState(gas.get(), kind == 0, 1.4, reservoir));
      return CEulerSolver::ComputeReservoirInflowState(gas.get(), 2, input, normal, input[5], input[4], 0, kind == 0,
                                                       1.4, reservoir, output);
    };
    // Subsonic inflow (the root carries the derivatives), the sonic limit and the reservoir at rest (both independent
    // of the interior): normal Mach number relative to the grid, pressure ratio and grid speed.
    for (const auto precomputed : {false, true})
      for (const auto& [mach, ratio, grid] : {std::tuple<double, double, double>{-.3, 1.05, -20},
                                              {-.05, 1.05, -20},
                                              {-1.3, 1.05, -20},
                                              {-1.3, 1.05, 0},
                                              {-.02, .9, 0}}) {
        CAPTURE(kind, precomputed, mach, ratio, grid);
        const double u = mach * a + grid;
        const double values[6] = {rho, rho * u, rho * 5, rho * (e + .5 * (u * u + 25)), ratio * p, grid};
        CheckMapDerivatives(
            [&](su2double* in, su2double* out) { return reservoirMap(in, out, false, false, precomputed); }, values,
            scale, 1e-9);
      }
    // Inflow with the interior tangential velocity, with a positive and a zero total pressure.
    for (const auto zeroTotal : {false, true})
      for (const auto mach : {-.3, -.05}) {
        CAPTURE(kind, zeroTotal, mach);
        const double grid = -20, u = mach * a + grid, ratio = zeroTotal ? 1. : 1.05;
        const double values[6] = {rho, rho * u, rho * .3 * a, rho * (e + .5 * (u * u + .09 * a * a)), ratio * p, grid};
        CheckMapDerivatives(
            [&](su2double* in, su2double* out) { return reservoirMap(in, out, zeroTotal, true, false); }, values, scale,
            1e-9);
      }
  }
}
