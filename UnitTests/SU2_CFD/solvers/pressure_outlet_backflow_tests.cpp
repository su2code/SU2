/*!
 * \file pressure_outlet_backflow_tests.cpp
 * \brief Backflow states of static-pressure outlets: reservoir inflow and the sonic fallback.
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
#include "../../../SU2_CFD/include/fluid/CIdealGas.hpp"
#include "../../../SU2_CFD/include/fluid/CPengRobinson.hpp"
#include "../../../SU2_CFD/include/fluid/CVanDerWaalsGas.hpp"
#include "../../../Common/include/toolboxes/geometry_toolbox.hpp"
#include "pressure_outlet_test_utils.hpp"

namespace {
using Status = CFluidModel::PressureOutletStatus;

// Thermodynamic state of a conservative boundary state: pressure, entropy, enthalpy, sound speed.
struct BoundaryState {
  su2double pressure, entropy, enthalpy, sound, speed, velocity[3] = {};
};
BoundaryState Evaluate(CFluidModel& gas, unsigned short dim, const su2double* u, su2double k = 0) {
  BoundaryState b;
  for (auto d = 0u; d < dim; ++d) b.velocity[d] = u[d + 1] / u[0];
  b.speed = GeometryToolbox::Norm(dim, b.velocity);
  gas.SetTDState_rhoe(u[0], u[dim + 1] / u[0] - .5 * b.speed * b.speed - k);
  b.pressure = gas.GetPressure();
  b.entropy = gas.GetEntropy();
  b.enthalpy = gas.GetStaticEnergy() + b.pressure / u[0];
  b.sound = gas.GetSoundSpeed();
  return b;
}

double Value(su2double x) { return SU2_TYPE::GetValue(x); }
}  // namespace

TEST_CASE("Reservoir backflow keeps the reservoir state and the interior acoustic invariant",
          "[Pressure-outlet][Outlet-backflow]") {
  for (auto dim : {2u, 3u})
    for (auto kind : {0u, 1u, 2u}) {
      auto gas = MakeOutletTestGas(kind);
      const su2double p = 1500112, rho = 230, normal[3] = {1, 0, 0};
      gas->SetTDState_Prho(p, rho);
      const su2double e = gas->GetStaticEnergy(), a = gas->GetSoundSpeed(), T = gas->GetTemperature();
      // Reservoir 5% above the interior pressure at the interior temperature.
      const su2double totalPressure = 1.05 * p, totalTemperature = T;
      gas->SetTDState_PT(totalPressure, totalTemperature);
      const su2double totalEnthalpy = gas->GetStaticEnergy() + gas->GetPressure() / gas->GetDensity();
      const su2double totalEntropy = gas->GetEntropy();
      const su2double oblique[3] = {-.8, .6, 0}, zero[3] = {};
      // The reservoir is at rest in the lab frame; a grid moving into the domain keeps the backflow inward.
      for (const su2double grid : {0., -20.})
        for (const auto mach : {-.3, -.05})
          for (const su2double* direction : {zero, oblique}) {
            CAPTURE(dim, kind, grid, mach, direction == oblique);
            su2double velocity[3] = {mach * a + grid, 5, -2}, state[5] = {rho}, out[5];
            for (auto d = 0u; d < dim; ++d) state[d + 1] = rho * velocity[d];
            state[dim + 1] = rho * (e + .5 * GeometryToolbox::SquaredNorm(dim, velocity));
            REQUIRE(gas->ComputeReservoirInflowState(dim, state, normal, grid,
                                                     MakeReservoir(totalPressure, totalTemperature, direction), out, 0,
                                                     kind == 0 ? 1.4 : 0) == Status::Valid);
            const auto b = Evaluate(*gas, dim, out);
            // Reservoir total enthalpy and entropy.
            // The reference uses SetTDState_PT directly; the van der Waals (p,T) inversion is only accurate to ~1e-6.
            CHECK(Value(b.enthalpy + .5 * b.speed * b.speed) ==
                  Approx(Value(totalEnthalpy)).epsilon(kind == 2 ? 1e-6 : 1e-9));
            CHECK(Value(b.entropy) == Approx(Value(totalEntropy)).margin(1e-5));
            // Velocity along the given direction (the inward normal for a zero direction).
            su2double unit[3] = {-1, 0, 0};
            if (direction == oblique)
              for (auto d = 0u; d < dim; ++d) unit[d] = oblique[d] / GeometryToolbox::Norm(dim, oblique);
            CHECK(Value(GeometryToolbox::DotProduct(dim, b.velocity, unit)) == Approx(Value(b.speed)).epsilon(1e-10));
            // Subsonic inflow relative to the grid.
            const su2double q = b.velocity[0] - grid;
            CHECK(q < 0);
            CHECK(-q < b.sound);
            // The interior outgoing invariant: the pressure-outlet map at the same pressure gives the same speed.
            su2double reference[5];
            REQUIRE(gas->ComputePressureOutletState(dim, state, normal, grid, b.pressure, reference, 0,
                                                    kind == 0 ? 1.4 : 0) == Status::Valid);
            CHECK(Value(q) == Approx(Value(reference[1] / reference[0] - grid)).margin(1e-7 * Value(a)));
            if (kind == 0) {
              // Gamma=1.4: q + 5*a is constant along the interior isentrope.
              const su2double interiorSound = a * pow(b.pressure / p, 1 / 7.);
              CHECK(Value(q + 5 * interiorSound) == Approx(Value(mach * a + 5 * a)).epsilon(1e-10));
            }
          }
    }
}

TEST_CASE("Reservoir backflow is limited to the sonic state of the reservoir", "[Pressure-outlet][Outlet-backflow]") {
  for (auto kind : {0u, 1u, 2u}) {
    auto gas = MakeOutletTestGas(kind);
    const su2double p = 1500112, rho = 230, normal[2] = {1, 0};
    gas->SetTDState_Prho(p, rho);
    const su2double e = gas->GetStaticEnergy(), a = gas->GetSoundSpeed(), T = gas->GetTemperature();
    su2double sonic[4] = {};
    // A strong reservoir drives supersonic backflow; fully incoming interior flow carries no interior information.
    for (const auto& [mach, ratio] : {std::pair<double, double>{-.95, 3.}, {-1.3, 1.05}, {-1.6, 1.05}}) {
      CAPTURE(kind, mach, ratio);
      // For the dense real gases the strong reservoir still gives subsonic inflow; the other cases enter fully.
      if (kind != 0 && mach == -.95) continue;
      const su2double state[4] = {rho, rho * mach * a, rho * 5, rho * (e + .5 * (mach * mach * a * a + 25))};
      su2double out[4];
      REQUIRE(gas->ComputeReservoirInflowState(2, state, normal, 0, MakeReservoir(ratio * p, T, nullptr), out, 0,
                                               kind == 0 ? 1.4 : 0) == Status::Valid);
      const auto b = Evaluate(*gas, 2, out);
      CHECK(Value(b.speed) == Approx(Value(b.sound)).epsilon(1e-8));
      CHECK(Value(b.velocity[0]) < 0);
      CHECK(fabs(Value(b.velocity[1])) < 1e-12 * Value(a));
      if (kind == 0) {
        // Ideal gas: p* = p_t (2/(gamma+1))^(gamma/(gamma-1)).
        CHECK(Value(b.pressure) == Approx(Value(ratio * p * pow(2 / 2.4, 3.5))).epsilon(1e-10));
      }
      // The sonic state depends on the reservoir only.
      if (mach == -1.6)
        for (auto k = 0u; k < 4; ++k) CHECK(Value(out[k]) == Approx(Value(sonic[k])).epsilon(1e-12));
      if (mach == -1.3)
        for (auto k = 0u; k < 4; ++k) sonic[k] = out[k];
    }
  }
}

TEST_CASE("Reservoir backflow stops when the interior wave cannot draw flow in", "[Pressure-outlet][Outlet-backflow]") {
  for (auto kind : {0u, 1u, 2u}) {
    auto gas = MakeOutletTestGas(kind);
    const su2double p = 1500112, rho = 230, normal[2] = {1, 0};
    gas->SetTDState_Prho(p, rho);
    const su2double e = gas->GetStaticEnergy(), a = gas->GetSoundSpeed(), T = gas->GetTemperature();
    // The reservoir pressure is 10% below the interior, so the weak reverse flow would turn outward.
    const su2double state[4] = {rho, rho * -.02 * a, 0, rho * (e + .5 * .0004 * a * a)};
    su2double out[4];
    REQUIRE(gas->ComputeReservoirInflowState(2, state, normal, 0, MakeReservoir(.9 * p, T, nullptr), out, 0,
                                             kind == 0 ? 1.4 : 0) == Status::Valid);
    const auto b = Evaluate(*gas, 2, out);
    CHECK(Value(b.speed) == 0);
    CHECK(Value(b.pressure) == Approx(Value(.9 * p)).epsilon(1e-8));
    CHECK(Value(gas->GetTemperature()) == Approx(Value(T)).epsilon(1e-8));
  }
}

TEST_CASE("Reservoir backflow keeps frozen SST energy out of the EOS", "[Pressure-outlet][Outlet-backflow]") {
  for (auto kind : {0u, 1u}) {
    auto gas = MakeOutletTestGas(kind);
    const su2double p = 1500112, rho = 230, normal[3] = {0, 0, 1}, k = 17.5;
    gas->SetTDState_Prho(p, rho);
    const su2double e = gas->GetStaticEnergy(), a = gas->GetSoundSpeed(), T = gas->GetTemperature();
    su2double state[5] = {rho, 0, 0, -.2 * rho * a, rho * (e + .02 * a * a)}, plain[5], withK[5];
    REQUIRE(gas->ComputeReservoirInflowState(3, state, normal, 0, MakeReservoir(1.05 * p, T, nullptr), plain, 0,
                                             kind == 0 ? 1.4 : 0) == Status::Valid);
    state[4] += rho * k;
    REQUIRE(gas->ComputeReservoirInflowState(3, state, normal, 0, MakeReservoir(1.05 * p, T, nullptr), withK, k,
                                             kind == 0 ? 1.4 : 0) == Status::Valid);
    for (auto v = 0u; v < 4; ++v) CHECK(Value(withK[v]) == Approx(Value(plain[v])).epsilon(1e-9).margin(1e-9));
    CHECK(Value(withK[4]) == Approx(Value(plain[4] + k * plain[0])).epsilon(1e-9));
  }
}

TEST_CASE("Reservoir backflow rejects a direction that leaves the domain", "[Pressure-outlet][Outlet-backflow]") {
  CIdealGas gas(1.4, 287);
  const su2double normal[2] = {1, 0}, outward[3] = {1, .2, 0};
  const su2double state[4] = {1.2, -30 * 1.2, 0, 1e5 / .4 + .5 * 1.2 * 900}, rest[4] = {1.2, 0, 0, 1e5 / .4};
  su2double out[4];
  CHECK(gas.ComputeReservoirInflowState(2, state, normal, 0, MakeReservoir(1.01e5, 300, outward), out, 0, 1.4) ==
        Status::InvalidState);
  CHECK(gas.ComputeReservoirInflowState(2, rest, normal, 0, MakeReservoir(-1, 300, nullptr), out, 0, 1.4) ==
        Status::InvalidState);
}

TEST_CASE("Backflow fallback keeps interior entropy, tangential velocity and the imposed pressure",
          "[Pressure-outlet][Outlet-backflow]") {
  for (auto dim : {2u, 3u})
    for (auto kind : {0u, 1u, 2u}) {
      auto gas = MakeOutletTestGas(kind);
      const su2double p = 1500112, rho = 230, normal[3] = {1, 0, 0};
      gas->SetTDState_Prho(p, rho);
      const su2double e = gas->GetStaticEnergy(), a = gas->GetSoundSpeed(), s = gas->GetEntropy();
      for (const su2double grid : {0., 20.})
        for (const auto ratio : {.95, 1.1}) {
          CAPTURE(dim, kind, grid, ratio);
          su2double velocity[3] = {-1.4 * a + grid, 5, -2}, state[5] = {rho}, out[5];
          for (auto d = 0u; d < dim; ++d) state[d + 1] = rho * velocity[d];
          state[dim + 1] = rho * (e + .5 * GeometryToolbox::SquaredNorm(dim, velocity));
          su2double unused[5];
          REQUIRE(gas->ComputePressureOutletState(dim, state, normal, grid, ratio * p, unused, 0,
                                                  kind == 0 ? 1.4 : 0) == Status::Backflow);
          REQUIRE(gas->ComputeSonicBackflowState(dim, state, normal, grid, ratio * p, out, 0, kind == 0 ? 1.4 : 0) ==
                  Status::Valid);
          const auto b = Evaluate(*gas, dim, out);
          CHECK(Value(b.pressure) == Approx(Value(ratio * p)).epsilon(1e-7));
          CHECK(Value(b.entropy) == Approx(Value(s)).margin(1e-5));
          // Inflow, never faster than sonic relative to the grid. For the ideal gas the interior invariant would
          // give supersonic inflow here, so the sonic limit is active; dense gases can stay subsonic.
          const su2double q = b.velocity[0] - grid;
          CHECK(q < 0);
          CHECK(Value(-q) <= Value(b.sound) * (1 + 1e-9));
          if (kind == 0) CHECK(Value(q) == Approx(Value(-b.sound)).epsilon(1e-9));
          for (auto d = 1u; d < dim; ++d) CHECK(Value(b.velocity[d]) == Approx(Value(velocity[d])).epsilon(1e-12));
        }
    }
}

TEST_CASE("Backflow fallback is continuous with subsonic backflow", "[Pressure-outlet][Outlet-backflow]") {
  // Bisect the interior normal speed at which the pressure-outlet map switches from subsonic backflow to backflow
  // through every characteristic. With compression the boundary speed reaches -a; with expansion the interior speed
  // reaches -a first. On the fully incoming side the fallback must match the subsonic state on the other side.
  for (auto kind : {0u, 1u, 2u})
    for (const auto ratio : {1.05, .95}) {
      CAPTURE(kind, ratio);
      auto gas = MakeOutletTestGas(kind);
      const su2double p = 1500112, rho = 230, normal[2] = {1, 0};
      gas->SetTDState_Prho(p, rho);
      const su2double e = gas->GetStaticEnergy(), a = gas->GetSoundSpeed(), target = ratio * p;
      const su2double gamma = kind == 0 ? 1.4 : 0;
      auto state = [&](double mach, su2double* u) {
        u[0] = rho;
        u[1] = rho * mach * a;
        u[2] = rho * 3;
        u[3] = rho * (e + .5 * (mach * mach * a * a + 9));
      };
      double lower = -1.5, upper = -.5;
      su2double u[4], subsonic[4], fallback[4];
      for (auto iteration = 0u; iteration < 50; ++iteration) {
        const double middle = .5 * (lower + upper);
        state(middle, u);
        if (gas->ComputePressureOutletState(2, u, normal, 0, target, subsonic, 0, gamma) == Status::Backflow)
          lower = middle;
        else
          upper = middle;
      }
      state(upper, u);
      REQUIRE(gas->ComputePressureOutletState(2, u, normal, 0, target, subsonic, 0, gamma) == Status::Valid);
      state(lower, u);
      REQUIRE(gas->ComputeSonicBackflowState(2, u, normal, 0, target, fallback, 0, gamma) == Status::Valid);
      for (auto k = 0u; k < 4; ++k)
        CHECK(Value(fallback[k]) == Approx(Value(subsonic[k])).epsilon(1e-6).margin(1e-6 * Value(rho * a)));
    }
}

TEST_CASE("Outlet pressure and normal velocity are continuous at flow reversal for a plenum at the outlet pressure",
          "[Pressure-outlet][Outlet-backflow]") {
  // With the reservoir total pressure equal to the imposed static pressure, the boundary pressure and normal velocity
  // (hence the physical flux) are continuous where the outlet state switches between outflow and reservoir inflow;
  // density, temperature and tangential velocity jump, as across a contact discontinuity.
  for (auto kind : {0u, 1u, 2u}) {
    CAPTURE(kind);
    auto gas = MakeOutletTestGas(kind);
    const su2double p = 1500112, rho = 230, normal[2] = {1, 0};
    gas->SetTDState_Prho(p, rho);
    const su2double e = gas->GetStaticEnergy(), a = gas->GetSoundSpeed(), T = gas->GetTemperature();
    const su2double reservoir[5] = {T, 0, 0, 0, 0};
    for (const auto mach : {1e-7, -1e-7}) {
      CAPTURE(mach);
      const su2double state[4] = {rho, rho * mach * a, rho * 40, rho * (e + .5 * (mach * mach * a * a + 1600))};
      su2double out[4];
      REQUIRE(ExpectedOutletState(*gas, 2, state, normal, 0, p, out, 0, kind == 0 ? 1.4 : 0, reservoir, true) ==
              Status::Valid);
      const auto b = Evaluate(*gas, 2, out);
      CHECK(Value(b.pressure) == Approx(Value(p)).epsilon(1e-6));
      CHECK(fabs(Value(b.velocity[0])) < 1e-5 * Value(a));
    }
  }
}

TEST_CASE("Reservoir backflow with the interior tangential velocity keeps total enthalpy and entropy",
          "[Pressure-outlet][Outlet-backflow]") {
  for (auto dim : {2u, 3u})
    for (auto kind : {0u, 1u, 2u}) {
      auto gas = MakeOutletTestGas(kind);
      const su2double p = 1500112, rho = 230, normal[3] = {1, 0, 0};
      gas->SetTDState_Prho(p, rho);
      const su2double e = gas->GetStaticEnergy(), a = gas->GetSoundSpeed(), T = gas->GetTemperature();
      const su2double totalPressure = 1.05 * p;
      gas->SetTDState_PT(totalPressure, T);
      const su2double totalEnthalpy = gas->GetStaticEnergy() + gas->GetPressure() / gas->GetDensity();
      const su2double totalEntropy = gas->GetEntropy();
      for (const auto mach : {-.3, -.1})
        for (const auto swirl : {.05, .2}) {
          CAPTURE(dim, kind, mach, swirl);
          su2double velocity[3] = {mach * a, swirl * a, -.5 * swirl * a}, state[5] = {rho}, out[5], tangential[3] = {};
          for (auto d = 0u; d < dim; ++d) state[d + 1] = rho * velocity[d];
          state[dim + 1] = rho * (e + .5 * GeometryToolbox::SquaredNorm(dim, velocity));
          InteriorTangentialVelocity(dim, state, normal, tangential);
          REQUIRE(gas->ComputeReservoirInflowState(dim, state, normal, 0,
                                                   MakeReservoir(totalPressure, T, nullptr, tangential), out, 0,
                                                   kind == 0 ? 1.4 : 0) == Status::Valid);
          const auto b = Evaluate(*gas, dim, out);
          // The van der Waals (p,T) reference is accurate only to about 1e-6.
          CHECK(Value(b.enthalpy + .5 * b.speed * b.speed) ==
                Approx(Value(totalEnthalpy)).epsilon(kind == 2 ? 1e-6 : 1e-9));
          CHECK(Value(b.entropy) == Approx(Value(totalEntropy)).margin(1e-5));
          for (auto d = 1u; d < dim; ++d) CHECK(Value(b.velocity[d]) == Approx(Value(velocity[d])).epsilon(1e-10));
          // Subsonic normal inflow that meets the interior outgoing invariant.
          CHECK(b.velocity[0] < 0);
          CHECK(-b.velocity[0] < b.sound);
          su2double reference[5];
          REQUIRE(gas->ComputePressureOutletState(dim, state, normal, 0, b.pressure, reference, 0,
                                                  kind == 0 ? 1.4 : 0) == Status::Valid);
          CHECK(Value(b.velocity[0]) == Approx(Value(reference[1] / reference[0])).margin(1e-7 * Value(a)));
        }
    }
}

TEST_CASE("Interior tangential backflow is continuous at flow reversal for a zero total pressure",
          "[Pressure-outlet][Outlet-backflow]") {
  // With a zero total pressure the reservoir has the imposed pressure where the inflow has no normal velocity, so
  // pressure and the whole velocity are continuous at the reversal; only density and temperature jump.
  for (auto kind : {0u, 1u, 2u}) {
    CAPTURE(kind);
    auto gas = MakeOutletTestGas(kind);
    const su2double p = 1500112, rho = 230, normal[2] = {1, 0};
    gas->SetTDState_Prho(p, rho);
    const su2double e = gas->GetStaticEnergy(), a = gas->GetSoundSpeed(), T = gas->GetTemperature();
    const su2double reservoir[5] = {T, 0, 0, 0, 0};
    for (const auto mach : {1e-7, -1e-7}) {
      CAPTURE(mach);
      const su2double state[4] = {rho, rho * mach * a, rho * .3 * a, rho * (e + .5 * (mach * mach + .09) * a * a)};
      su2double out[4];
      REQUIRE(ExpectedOutletState(*gas, 2, state, normal, 0, p, out, 0, kind == 0 ? 1.4 : 0, reservoir, true, true) ==
              Status::Valid);
      const auto b = Evaluate(*gas, 2, out);
      CHECK(Value(b.pressure) == Approx(Value(p)).epsilon(1e-6));
      CHECK(fabs(Value(b.velocity[0])) < 1e-5 * Value(a));
      CHECK(Value(b.velocity[1]) == Approx(Value(.3 * a)).epsilon(1e-10));
    }
  }
}

TEST_CASE("Interior tangential backflow limits the normal speed to sonic", "[Pressure-outlet][Outlet-backflow]") {
  for (auto kind : {0u, 1u, 2u}) {
    CAPTURE(kind);
    auto gas = MakeOutletTestGas(kind);
    const su2double p = 1500112, rho = 230, normal[2] = {1, 0};
    gas->SetTDState_Prho(p, rho);
    const su2double e = gas->GetStaticEnergy(), a = gas->GetSoundSpeed(), T = gas->GetTemperature();
    const su2double state[4] = {rho, -1.3 * rho * a, .3 * rho * a, rho * (e + .5 * (1.69 + .09) * a * a)};
    su2double out[4], tangential[2];
    InteriorTangentialVelocity(2, state, normal, tangential);
    REQUIRE(gas->ComputeReservoirInflowState(2, state, normal, 0, MakeReservoir(1.05 * p, T, nullptr, tangential), out,
                                             0, kind == 0 ? 1.4 : 0) == Status::Valid);
    const auto b = Evaluate(*gas, 2, out);
    CHECK(Value(-b.velocity[0]) == Approx(Value(b.sound)).epsilon(1e-8));
    CHECK(Value(b.velocity[1]) == Approx(Value(.3 * a)).epsilon(1e-10));
  }
}

TEST_CASE("Reservoir backflow limits the normal speed relative to the grid to sonic",
          "[Pressure-outlet][Outlet-backflow]") {
  const su2double oblique[2] = {-.8, .6};
  for (auto kind : {0u, 1u, 2u})
    for (const su2double grid : {0., 20., -20.})
      for (const auto* direction : {static_cast<const su2double*>(nullptr), oblique}) {
        CAPTURE(kind, grid, direction != nullptr);
        auto gas = MakeOutletTestGas(kind);
        const su2double p = 1500112, rho = 230, normal[2] = {1, 0};
        gas->SetTDState_Prho(p, rho);
        const su2double e = gas->GetStaticEnergy(), a = gas->GetSoundSpeed(), T = gas->GetTemperature();
        // The interior flow enters through every characteristic, so the sonic state of the reservoir applies.
        const su2double u = -1.3 * a + grid;
        const su2double state[4] = {rho, rho * u, 0, rho * (e + .5 * u * u)};
        su2double out[4];
        REQUIRE(gas->ComputeReservoirInflowState(2, state, normal, grid, MakeReservoir(1.05 * p, T, direction), out, 0,
                                                 kind == 0 ? 1.4 : 0) == Status::Valid);
        const auto b = Evaluate(*gas, 2, out);
        CHECK(Value(b.velocity[0] - grid) == Approx(Value(-b.sound)).epsilon(1e-8));
        if (direction != nullptr) CHECK(Value(b.velocity[1] / b.velocity[0]) == Approx(-.75).epsilon(1e-10));
      }
}

TEST_CASE("Reservoir backflow takes the rest state when the grid outruns sound", "[Pressure-outlet][Outlet-backflow]") {
  const su2double oblique[2] = {-.8, .6};
  for (auto kind : {0u, 1u, 2u})
    for (const su2double gridFactor : {1.2, 3.})
      for (const auto* direction : {static_cast<const su2double*>(nullptr), oblique}) {
        CAPTURE(kind, gridFactor, direction != nullptr);
        auto gas = MakeOutletTestGas(kind);
        const su2double p = 1500112, rho = 230, normal[2] = {1, 0};
        gas->SetTDState_Prho(p, rho);
        const su2double e = gas->GetStaticEnergy(), a = gas->GetSoundSpeed(), T = gas->GetTemperature();
        gas->SetTDState_PT(1.05 * p, T);
        // The boundary moves outward faster than sound in the reservoir at rest, so no sonic inflow exists.
        const su2double grid = gridFactor * gas->GetSoundSpeed();
        for (const su2double mach : {-1.3, -.3}) {
          const su2double u = mach * a + grid;
          const su2double state[4] = {rho, rho * u, 0, rho * (e + .5 * u * u)};
          su2double out[4];
          REQUIRE(gas->ComputeReservoirInflowState(2, state, normal, grid, MakeReservoir(1.05 * p, T, direction), out,
                                                   0, kind == 0 ? 1.4 : 0) == Status::Valid);
          const auto b = Evaluate(*gas, 2, out);
          CHECK(Value(b.speed) == 0);
          CHECK(Value(b.pressure) == Approx(Value(1.05 * p)).epsilon(1e-8));
          CHECK(Value(gas->GetTemperature()) == Approx(Value(T)).epsilon(1e-6));
        }
      }
}
