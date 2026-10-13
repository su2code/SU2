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
#include "../../../SU2_CFD/include/solvers/characteristic_pressure_outlet.hpp"
#include "../../../SU2_CFD/include/fluid/CIdealGas.hpp"
#include "../../../SU2_CFD/include/fluid/CPengRobinson.hpp"
#include "../../../SU2_CFD/include/fluid/CVanDerWaalsGas.hpp"

TEST_CASE("Pressure outlet state derivatives match finite differences", "[Pressure-outlet-AD]") {
  for (auto kind : {0u, 1u, 2u}) {
    std::unique_ptr<CFluidModel> gas;
    if (kind == 0) gas.reset(new CIdealGas(1.4, 35.152));
    if (kind == 1) gas.reset(new CPengRobinson(1.0125, 35.152, 1415200, 564.1, .529));
    if (kind == 2) gas.reset(new CVanDerWaalsGas(1.0125, 35.152, 1415200, 564.1));
    const double rho = 230, p = 1500112;
    gas->SetTDState_Prho(p, rho);
    const double e = SU2_TYPE::GetValue(gas->GetStaticEnergy()), a = SU2_TYPE::GetValue(gas->GetSoundSpeed());
    const su2double normal[2] = {1, 0};
    for (const auto mach : {.3, .99, 1.4})
      for (const auto ratio : {.95, .001}) {
        const double values[6] = {
            rho, rho * (mach * a + 20), rho * 5, rho * (e + .5 * (pow(mach * a + 20, 2) + 25)), ratio * p, 20};
        auto evaluate = [&](su2double* input, su2double* output) {
          return PressureOutlet::State(*gas, 2, input, normal, input[5], input[4], output, 0, kind == 0 ? 1.4 : 0.);
        };
        for (auto column = 0u; column < 6; ++column) {
          const double scale[6] = {rho, rho * a, rho * a, rho * a * a, p, a};
          const double delta = 2e-5 * fmax(fabs(values[column]), scale[column]);
          su2double plus[6], minus[6], fp[4], fm[4];
          for (auto k = 0u; k < 6; ++k) plus[k] = minus[k] = values[k];
          plus[column] += delta;
          minus[column] -= delta;
          REQUIRE(evaluate(plus, fp) == PressureOutlet::Status::Valid);
          REQUIRE(evaluate(minus, fm) == PressureOutlet::Status::Valid);
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
            REQUIRE(evaluate(input, output) == PressureOutlet::Status::Valid);
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
