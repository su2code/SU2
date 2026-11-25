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

#include "catch.hpp"
#include "pressure_outlet_test_utils.hpp"

namespace {
using Status = CEulerSolver::PressureOutletStatus;

// Compare the derivative of every output with central finite differences of the complete map, one input at a time.
// Inputs: four conservative variables, then the imposed pressure and the grid speed along the normal.
template <class Map>
void CheckMapDerivatives(const Map& evaluate, const double* values, const double* scale) {
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
      CHECK(derivative == Approx(expected).epsilon(3e-4).margin(2e-4));
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
