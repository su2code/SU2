/*!
 * \file pressure_outlet_backflow_tests.cpp
 * \brief Inflow state of static-pressure outlets from a backflow reservoir (MARKER_OUTLET_BACKFLOW).
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

#include <utility>
#include "catch.hpp"
#include "pressure_outlet_test_utils.hpp"

namespace {
using Status = CEulerSolver::PressureOutletStatus;

// Reservoir with a total pressure and temperature, a direction (nullptr: inward normal) and the inflow type.
CEulerSolver::PressureOutletReservoir Reservoir(su2double totalPressure, su2double totalTemperature,
                                                const su2double* direction = nullptr, bool interiorTangential = false) {
  CEulerSolver::PressureOutletReservoir reservoir;
  reservoir.totalPressure = totalPressure;
  reservoir.totalTemperature = totalTemperature;
  if (direction != nullptr)
    for (auto d = 0u; d < 3; ++d) reservoir.direction[d] = direction[d];
  reservoir.interiorTangential = interiorTangential;
  return reservoir;
}

// Reservoir inflow state of the test gas (closed form for the ideal gas); the imposed pressure matters only for a
// zero total pressure.
Status Inflow(CFluidModel& gas, unsigned short kind, unsigned short dim, const su2double* state,
              const su2double* normal, su2double grid, const CEulerSolver::PressureOutletReservoir& reservoir,
              su2double* out, su2double turbulentK = 0, su2double pressure = OUTLET_TEST_P) {
  return CEulerSolver::ComputeReservoirInflowState(&gas, dim, state, normal, grid, pressure, turbulentK, kind == 0, 1.4,
                                                   reservoir, out);
}

// Thermodynamic state of a conservative boundary state.
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
      const su2double p = OUTLET_TEST_P, rho = OUTLET_TEST_RHO, normal[3] = {1, 0, 0};
      gas->SetTDState_Prho(p, rho);
      const su2double e = gas->GetStaticEnergy(), a = gas->GetSoundSpeed(), T = gas->GetTemperature();
      // Reservoir 5% above the interior pressure at the interior temperature.
      const su2double totalPressure = 1.05 * p;
      gas->SetTDState_PT(totalPressure, T);
      const su2double totalEnthalpy = gas->GetStaticEnergy() + gas->GetPressure() / gas->GetDensity();
      const su2double totalEntropy = gas->GetEntropy();
      const su2double oblique[3] = {-.8, .6, 0};
      // The reservoir is at rest in the lab frame; a grid moving into the domain keeps the backflow inward.
      for (const su2double grid : {0., -20.})
        for (const auto mach : {-.3, -.05})
          for (const su2double* direction : {static_cast<const su2double*>(nullptr), oblique}) {
            CAPTURE(dim, kind, grid, mach, direction == oblique);
            su2double velocity[3] = {mach * a + grid, 5, -2}, state[5] = {rho}, out[5];
            for (auto d = 0u; d < dim; ++d) state[d + 1] = rho * velocity[d];
            state[dim + 1] = rho * (e + .5 * GeometryToolbox::SquaredNorm(dim, velocity));
            REQUIRE(Inflow(*gas, kind, dim, state, normal, grid, Reservoir(totalPressure, T, direction), out) ==
                    Status::Valid);
            const auto b = Evaluate(*gas, dim, out);
            // Reservoir total enthalpy and entropy; the van der Waals (p,T) reference is accurate to about 1e-6.
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
            REQUIRE(CEulerSolver::ComputePressureOutletState(gas.get(), dim, state, normal, grid, b.pressure, 0,
                                                             kind == 0, 1.4, reference) == Status::Valid);
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
    const su2double p = OUTLET_TEST_P, rho = OUTLET_TEST_RHO, normal[2] = {1, 0};
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
      // With the precomputed total state of the solver as well.
      auto reservoir = Reservoir(ratio * p, T);
      for (const auto precomputed : {false, true}) {
        if (precomputed) REQUIRE(CEulerSolver::SetReservoirTotalState(gas.get(), kind == 0, 1.4, reservoir));
        REQUIRE(Inflow(*gas, kind, 2, state, normal, 0, reservoir, out) == Status::Valid);
        const auto b = Evaluate(*gas, 2, out);
        CHECK(Value(b.speed) == Approx(Value(b.sound)).epsilon(1e-8));
        CHECK(Value(b.velocity[0]) < 0);
        CHECK(fabs(Value(b.velocity[1])) < 1e-12 * Value(a));
        // Ideal gas: p* = p_t (2/(gamma+1))^(gamma/(gamma-1)).
        if (kind == 0) CHECK(Value(b.pressure) == Approx(Value(ratio * p * pow(2 / 2.4, 3.5))).epsilon(1e-10));
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
    const su2double p = OUTLET_TEST_P, rho = OUTLET_TEST_RHO, normal[2] = {1, 0};
    gas->SetTDState_Prho(p, rho);
    const su2double e = gas->GetStaticEnergy(), a = gas->GetSoundSpeed(), T = gas->GetTemperature();
    // The reservoir pressure is 10% below the interior, so the weak reverse flow would turn outward.
    const su2double state[4] = {rho, rho * -.02 * a, 0, rho * (e + .5 * .0004 * a * a)};
    su2double out[4];
    REQUIRE(Inflow(*gas, kind, 2, state, normal, 0, Reservoir(.9 * p, T), out) == Status::Valid);
    const auto b = Evaluate(*gas, 2, out);
    CHECK(Value(b.speed) == 0);
    CHECK(Value(b.pressure) == Approx(Value(.9 * p)).epsilon(1e-8));
    CHECK(Value(gas->GetTemperature()) == Approx(Value(T)).epsilon(1e-8));
  }
}

TEST_CASE("Reservoir backflow keeps the SST energy out of the EOS", "[Pressure-outlet][Outlet-backflow]") {
  for (auto kind : {0u, 1u}) {
    auto gas = MakeOutletTestGas(kind);
    const su2double p = OUTLET_TEST_P, rho = OUTLET_TEST_RHO, normal[3] = {0, 0, 1}, k = 17.5;
    gas->SetTDState_Prho(p, rho);
    const su2double e = gas->GetStaticEnergy(), a = gas->GetSoundSpeed(), T = gas->GetTemperature();
    su2double state[5] = {rho, 0, 0, -.2 * rho * a, rho * (e + .02 * a * a)}, plain[5], withK[5];
    REQUIRE(Inflow(*gas, kind, 3, state, normal, 0, Reservoir(1.05 * p, T), plain) == Status::Valid);
    state[4] += rho * k;
    REQUIRE(Inflow(*gas, kind, 3, state, normal, 0, Reservoir(1.05 * p, T), withK, k) == Status::Valid);
    for (auto v = 0u; v < 4; ++v) CHECK(Value(withK[v]) == Approx(Value(plain[v])).epsilon(1e-9).margin(1e-9));
    CHECK(Value(withK[4]) == Approx(Value(plain[4] + k * plain[0])).epsilon(1e-9));
  }
}

TEST_CASE("Reservoir backflow rejects invalid reservoir data", "[Pressure-outlet][Outlet-backflow]") {
  CIdealGas gas(1.4, 287);
  const su2double normal[2] = {1, 0}, outward[3] = {1, .2, 0};
  const su2double state[4] = {1.2, -30 * 1.2, 0, 1e5 / .4 + .5 * 1.2 * 900};
  su2double out[4];
  // A direction that leaves the domain, a negative total pressure, a zero total pressure with an invalid imposed
  // pressure, and a direction together with the interior tangential velocity.
  CHECK(Inflow(gas, 0, 2, state, normal, 0, Reservoir(1.01e5, 300, outward), out) == Status::InvalidState);
  CHECK(Inflow(gas, 0, 2, state, normal, 0, Reservoir(-1, 300), out) == Status::InvalidState);
  CHECK(Inflow(gas, 0, 2, state, normal, 0, Reservoir(0, 300), out, 0, 0) == Status::InvalidState);
  CHECK(Inflow(gas, 0, 2, state, normal, 0, Reservoir(1.01e5, 300, outward, true), out) == Status::InvalidState);
}

TEST_CASE("Pressure and normal velocity are continuous at flow reversal for a zero total pressure",
          "[Pressure-outlet][Outlet-backflow]") {
  // With a zero total pressure the reservoir has the imposed pressure at zero normal velocity, so the boundary
  // pressure and normal velocity (hence the physical flux) are continuous where the outlet state switches between
  // outflow and reservoir inflow; density and temperature jump, as across a contact discontinuity. With the interior
  // tangential velocity the whole velocity is continuous.
  for (auto kind : {0u, 1u, 2u})
    for (const auto tangential : {false, true}) {
      CAPTURE(kind, tangential);
      auto gas = MakeOutletTestGas(kind);
      const su2double p = OUTLET_TEST_P, rho = OUTLET_TEST_RHO, normal[2] = {1, 0};
      gas->SetTDState_Prho(p, rho);
      const su2double e = gas->GetStaticEnergy(), a = gas->GetSoundSpeed(), T = gas->GetTemperature();
      const auto reservoir = Reservoir(0, T, nullptr, tangential);
      for (const auto mach : {1e-7, -1e-7}) {
        CAPTURE(mach);
        const su2double state[4] = {rho, rho * mach * a, rho * .3 * a, rho * (e + .5 * (mach * mach + .09) * a * a)};
        su2double out[4];
        CHECK_FALSE(ExpectedSteadyOutletState(*gas, 2, state, normal, 0, p, 0, kind == 0, 1.4, out, &reservoir));
        const auto b = Evaluate(*gas, 2, out);
        CHECK(Value(b.pressure) == Approx(Value(p)).epsilon(1e-6));
        CHECK(fabs(Value(b.velocity[0])) < 1e-5 * Value(a));
        if (tangential || mach > 0) CHECK(Value(b.velocity[1]) == Approx(Value(.3 * a)).epsilon(1e-10));
      }
    }
}

TEST_CASE("Reservoir backflow with the interior tangential velocity keeps total enthalpy and entropy",
          "[Pressure-outlet][Outlet-backflow]") {
  for (auto dim : {2u, 3u})
    for (auto kind : {0u, 1u, 2u}) {
      auto gas = MakeOutletTestGas(kind);
      const su2double p = OUTLET_TEST_P, rho = OUTLET_TEST_RHO, normal[3] = {1, 0, 0};
      gas->SetTDState_Prho(p, rho);
      const su2double e = gas->GetStaticEnergy(), a = gas->GetSoundSpeed(), T = gas->GetTemperature();
      const su2double totalPressure = 1.05 * p;
      gas->SetTDState_PT(totalPressure, T);
      const su2double totalEnthalpy = gas->GetStaticEnergy() + gas->GetPressure() / gas->GetDensity();
      const su2double totalEntropy = gas->GetEntropy();
      for (const auto mach : {-.3, -.1})
        for (const auto swirl : {.05, .2}) {
          CAPTURE(dim, kind, mach, swirl);
          su2double velocity[3] = {mach * a, swirl * a, -.5 * swirl * a}, state[5] = {rho}, out[5];
          for (auto d = 0u; d < dim; ++d) state[d + 1] = rho * velocity[d];
          state[dim + 1] = rho * (e + .5 * GeometryToolbox::SquaredNorm(dim, velocity));
          REQUIRE(Inflow(*gas, kind, dim, state, normal, 0, Reservoir(totalPressure, T, nullptr, true), out) ==
                  Status::Valid);
          const auto b = Evaluate(*gas, dim, out);
          CHECK(Value(b.enthalpy + .5 * b.speed * b.speed) ==
                Approx(Value(totalEnthalpy)).epsilon(kind == 2 ? 1e-6 : 1e-9));
          CHECK(Value(b.entropy) == Approx(Value(totalEntropy)).margin(1e-5));
          for (auto d = 1u; d < dim; ++d) CHECK(Value(b.velocity[d]) == Approx(Value(velocity[d])).epsilon(1e-10));
          // Subsonic normal inflow that meets the interior outgoing invariant.
          CHECK(b.velocity[0] < 0);
          CHECK(-b.velocity[0] < b.sound);
          su2double reference[5];
          REQUIRE(CEulerSolver::ComputePressureOutletState(gas.get(), dim, state, normal, 0, b.pressure, 0, kind == 0,
                                                           1.4, reference) == Status::Valid);
          CHECK(Value(b.velocity[0]) == Approx(Value(reference[1] / reference[0])).margin(1e-7 * Value(a)));
        }
    }
}

TEST_CASE("Interior tangential backflow limits the normal speed to sonic", "[Pressure-outlet][Outlet-backflow]") {
  for (auto kind : {0u, 1u, 2u}) {
    CAPTURE(kind);
    auto gas = MakeOutletTestGas(kind);
    const su2double p = OUTLET_TEST_P, rho = OUTLET_TEST_RHO, normal[2] = {1, 0};
    gas->SetTDState_Prho(p, rho);
    const su2double e = gas->GetStaticEnergy(), a = gas->GetSoundSpeed(), T = gas->GetTemperature();
    const su2double state[4] = {rho, -1.3 * rho * a, .3 * rho * a, rho * (e + .5 * (1.69 + .09) * a * a)};
    su2double out[4];
    REQUIRE(Inflow(*gas, kind, 2, state, normal, 0, Reservoir(1.05 * p, T, nullptr, true), out) == Status::Valid);
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
        const su2double p = OUTLET_TEST_P, rho = OUTLET_TEST_RHO, normal[2] = {1, 0};
        gas->SetTDState_Prho(p, rho);
        const su2double e = gas->GetStaticEnergy(), a = gas->GetSoundSpeed(), T = gas->GetTemperature();
        // The interior flow enters through every characteristic, so the sonic state of the reservoir applies.
        const su2double u = -1.3 * a + grid;
        const su2double state[4] = {rho, rho * u, 0, rho * (e + .5 * u * u)};
        su2double out[4];
        REQUIRE(Inflow(*gas, kind, 2, state, normal, grid, Reservoir(1.05 * p, T, direction), out) == Status::Valid);
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
        const su2double p = OUTLET_TEST_P, rho = OUTLET_TEST_RHO, normal[2] = {1, 0};
        gas->SetTDState_Prho(p, rho);
        const su2double e = gas->GetStaticEnergy(), a = gas->GetSoundSpeed(), T = gas->GetTemperature();
        gas->SetTDState_PT(1.05 * p, T);
        // The boundary moves outward faster than sound in the reservoir at rest, so no sonic inflow exists.
        const su2double grid = gridFactor * gas->GetSoundSpeed();
        for (const su2double mach : {-1.3, -.3}) {
          const su2double u = mach * a + grid;
          const su2double state[4] = {rho, rho * u, 0, rho * (e + .5 * u * u)};
          su2double out[4];
          REQUIRE(Inflow(*gas, kind, 2, state, normal, grid, Reservoir(1.05 * p, T, direction), out) == Status::Valid);
          const auto b = Evaluate(*gas, 2, out);
          CHECK(Value(b.speed) == 0);
          CHECK(Value(b.pressure) == Approx(Value(1.05 * p)).epsilon(1e-8));
          CHECK(Value(gas->GetTemperature()) == Approx(Value(T)).epsilon(1e-6));
        }
      }
}
