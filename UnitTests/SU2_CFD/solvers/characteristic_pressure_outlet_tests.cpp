/*!
 * \file characteristic_pressure_outlet_tests.cpp
 * \brief Static-pressure state of the outlets and the BC_Riemann outlet.
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
#include "../../../SU2_CFD/include/numerics/flow/convection/hllc.hpp"
#include "../../../SU2_CFD/include/fluid/CDataDrivenFluid.hpp"

namespace {
// BC_Riemann static-pressure outlet on a moving grid, at the start-up state of an MDM nozzle. A low back pressure
// leaves the Peng-Robinson validity range under a linear characteristic jump, but not along the isentrope.
void RiemannPressureOutlet(su2double pressure, bool backflow = false) {
  UnitQuadTestCase test;
  std::ostringstream options;
  options << OutletTestOptions(2, 1, 0, true, "5,129,5")
          << "GRID_MOVEMENT= STEADY_TRANSLATION\nTRANSLATION_RATE= 20,0,0\nMARKER_EULER= (x_minus,y_minus,y_plus)\n"
          << "MARKER_MONITORING= (x_minus)\nMARKER_RIEMANN= (x_plus,STATIC_PRESSURE," << pressure << ",0,0,0,0)\n";
  test.config_options = options.str();
  test.InitConfig();
  test.InitGeometry();
  test.InitSolver();
  auto* solver = test.solver[FLOW_SOL];
  auto gas = MakeOutletTestGas(1);
  su2double state[4] = {139.72158881490299, 11271.637336535776, -2338.5053086132339, 221781581.74179047};
  if (backflow) {
    // Inflow faster than sound with the same thermal state.
    const auto oldVelocity = state[1] / state[0];
    state[1] = state[0] * -200;
    state[3] += .5 * state[0] * (40000 - oldVelocity * oldVelocity);
  }
  const su2double grid[2] = {20, 0};
  SetOutletTestState(test, *gas, state, grid, 0);
  const auto marker = OutletTestMarker(*test.config, "x_plus");
  SU2_OMP_PARALLEL {
    CUpwGeneralHLLC_Flow numerics(2, 4, test.config.get());
    solver->BC_Riemann(test.geometry.get(), test.solver, &numerics, nullptr, test.config.get(), marker);
  }
  CUpwGeneralHLLC_Flow numerics(2, 4, test.config.get());
  unsigned long fallbacks = 0;
  for (auto vertex = 0ul; vertex < test.geometry->nVertex[marker]; ++vertex) {
    const auto point = test.geometry->vertex[marker][vertex]->GetNode();
    if (!test.geometry->nodes->GetDomain(point)) continue;
    su2double normal[2], unit[2], boundary[4], flux[4];
    OutwardNormal(*test.geometry, marker, vertex, 2, normal, unit);
    if (ExpectedSteadyOutletState(*gas, 2, state, unit, GeometryToolbox::DotProduct(2, grid, unit), pressure, 0, false,
                                  0, boundary))
      ++fallbacks;
    PhysicalOutletFlux(*gas, numerics, 2, boundary, normal, grid, flux);
    CheckOutletPrimitives(solver->GetCharacPrimVar(marker, vertex), *gas, 2);
    for (auto k = 0u; k < 4; ++k)
      CHECK(SU2_TYPE::GetValue(solver->LinSysRes(point, k)) == Approx(SU2_TYPE::GetValue(flux[k])).epsilon(1e-8));
    CheckRiemannOutletJacobian(solver->Jacobian.GetBlock(point, point), *gas, numerics, 2, state, boundary, unit,
                               normal, grid);
  }
  // Full backflow and an invalid pressure use the steady fallback states, which the monitor counts.
  CHECK((fallbacks > 0) == (backflow || pressure <= 0));
  CHECK(CountOutletFallbacks(solver, test.config.get()) == fallbacks);
}
}  // namespace

TEST_CASE("Riemann pressure outlet follows the isentrope from the start-up state", "[Pressure-outlet]") {
  for (const auto pressure : {1400000., 100000., 1000.}) RiemannPressureOutlet(pressure);
}

TEST_CASE("Riemann pressure outlet uses the steady fallback states", "[Pressure-outlet]") {
  RiemannPressureOutlet(0);
  RiemannPressureOutlet(1400000, true);
}

TEST_CASE("Nonlinear pressure outflow preserves the acoustic invariant and frame", "[Pressure-outlet]") {
  for (auto dim : {2u, 3u}) {
    for (auto kind : {0u, 1u, 2u}) {
      auto gas = MakeOutletTestGas(kind);
      const su2double p = 1500112, rho = 230, normal[3] = {1, 0, 0}, boost[3] = {20, -7, 3};
      gas->SetTDState_Prho(p, rho);
      const auto e = gas->GetStaticEnergy(), sound = gas->GetSoundSpeed(), entropy = gas->GetEntropy();
      // The smallest ratio is 1000/1500112, reproducing the startup expansion past the sonic limit.
      for (const auto ratio : {1., .995, .9, 1.05, .0006666168959384}) {
        su2double out[2][5] = {};
        // Boost both flow and grid: the relative wave speeds and reconstructed thermodynamics must agree.
        for (auto frame = 0u; frame < 2; ++frame) {
          su2double velocity[3] = {.3 * sound, 5, -2}, state[5] = {rho};
          for (auto d = 0u; d < dim; ++d) {
            if (frame) velocity[d] += boost[d];
            state[d + 1] = rho * velocity[d];
          }
          state[dim + 1] = rho * (e + .5 * GeometryToolbox::SquaredNorm(dim, velocity));
          REQUIRE(CEulerSolver::ComputePressureOutletState(gas.get(), dim, state, normal, frame ? boost[0] : 0,
                                                           ratio * p, 0, false, 0,
                                                           out[frame]) == CEulerSolver::PressureOutletStatus::Valid);
          const auto density = out[frame][0];
          su2double v[3] = {};
          for (auto d = 0u; d < dim; ++d) v[d] = out[frame][d + 1] / density;
          gas->SetTDState_rhoe(density, out[frame][dim + 1] / density - .5 * GeometryToolbox::SquaredNorm(dim, v));
          const auto q = v[0] - (frame ? boost[0] : 0);
          CHECK(density > 0);
          CHECK(gas->GetTemperature() > 0);
          CHECK(gas->GetSoundSpeed2() > 0);
          CHECK(SU2_TYPE::GetValue(gas->GetEntropy()) == Approx(SU2_TYPE::GetValue(entropy)).margin(1e-5));
          if (ratio > .1) {
            CHECK(SU2_TYPE::GetValue(gas->GetPressure()) == Approx(SU2_TYPE::GetValue(ratio * p)).epsilon(1e-7));
            CHECK(q < gas->GetSoundSpeed());
          } else {
            CHECK(gas->GetPressure() > ratio * p);
            CHECK(SU2_TYPE::GetValue(q) == Approx(SU2_TYPE::GetValue(gas->GetSoundSpeed())).epsilon(1e-7));
          }
          if (kind == 0) {
            // Gamma=1.4 gives 2/(gamma-1)=5, so q+5*a provides an analytic integration check.
            const su2double aRatio = ratio > .1 ? pow(ratio, .4 / 2.8) : (.4 * .3 + 2) / 2.4;
            CHECK(SU2_TYPE::GetValue(density) == Approx(SU2_TYPE::GetValue(rho * pow(aRatio, 5))).epsilon(1e-7));
            CHECK(SU2_TYPE::GetValue(q + 5 * gas->GetSoundSpeed()) ==
                  Approx(SU2_TYPE::GetValue(5.3 * sound)).epsilon(1e-7));
          }
        }
        CHECK(SU2_TYPE::GetValue(out[1][0]) == Approx(SU2_TYPE::GetValue(out[0][0])).epsilon(1e-9));
        // In the boosted frame, rho*E gains boost dot momentum + rho*|boost|^2/2.
        su2double energyBoost = 0, boost2 = 0;
        for (auto d = 0u; d < dim; ++d) {
          CHECK(SU2_TYPE::GetValue(out[1][d + 1]) ==
                Approx(SU2_TYPE::GetValue(out[0][d + 1] + boost[d] * out[0][0])).epsilon(1e-8));
          energyBoost += boost[d] * out[0][d + 1];
          boost2 += boost[d] * boost[d];
        }
        CHECK(SU2_TYPE::GetValue(out[1][dim + 1]) ==
              Approx(SU2_TYPE::GetValue(out[0][dim + 1] + energyBoost + .5 * boost2 * out[0][0])).epsilon(1e-8));
      }
      su2double state[5] = {rho}, unchanged[5] = {};
      state[1] = 2 * sound * rho;
      state[dim + 1] = rho * (e + 2 * sound * sound);
      REQUIRE(CEulerSolver::ComputePressureOutletState(gas.get(), dim, state, normal, 0, 1000, 0, false, 0,
                                                       unchanged) == CEulerSolver::PressureOutletStatus::Valid);
      for (auto v = 0u; v < dim + 2; ++v) CHECK(unchanged[v] == state[v]);
      // Frozen SST k must remain in total energy without entering the EOS.
      state[1] = .3 * sound * rho;
      state[dim + 1] = rho * (e + .5 * .3 * .3 * sound * sound);
      su2double plain[5], withK[5];
      REQUIRE(CEulerSolver::ComputePressureOutletState(gas.get(), dim, state, normal, 0, .9 * p, 0, false, 0, plain) ==
              CEulerSolver::PressureOutletStatus::Valid);
      state[dim + 1] += 17.5 * rho;
      REQUIRE(CEulerSolver::ComputePressureOutletState(gas.get(), dim, state, normal, 0, .9 * p, 17.5, false, 0,
                                                       withK) == CEulerSolver::PressureOutletStatus::Valid);
      for (auto v = 0u; v < dim + 1; ++v) CHECK(withK[v] == Approx(SU2_TYPE::GetValue(plain[v])).epsilon(1e-8));
      CHECK(withK[dim + 1] == Approx(SU2_TYPE::GetValue(plain[dim + 1] + 17.5 * plain[0])).epsilon(1e-8));
      // Reversing subsonic flow retains the pressure-only closure; fully incoming flow needs inlet data.
      state[1] = -state[1];
      CHECK(CEulerSolver::ComputePressureOutletState(gas.get(), dim, state, normal, 0, p, 17.5, false, 0, unchanged) ==
            CEulerSolver::PressureOutletStatus::Valid);
      state[1] = -2 * sound * rho;
      state[dim + 1] = rho * (e + 2 * sound * sound + 17.5);
      CHECK(CEulerSolver::ComputePressureOutletState(gas.get(), dim, state, normal, 0, p, 17.5, false, 0, unchanged) ==
            CEulerSolver::PressureOutletStatus::Backflow);
    }
  }
}

namespace {
// Counts the EOS evaluations of a wrapped fluid model.
class CountingFluid final : public CFluidModel {
  CFluidModel& gas;

 public:
  unsigned long calls = 0;
  explicit CountingFluid(CFluidModel& wrapped) : gas(wrapped) {}
  void SetTDState_rhoe(su2double rho, su2double e) override {
    ++calls;
    gas.SetTDState_rhoe(rho, e);
    Density = gas.GetDensity();
    StaticEnergy = gas.GetStaticEnergy();
    Pressure = gas.GetPressure();
    Temperature = gas.GetTemperature();
    SoundSpeed2 = gas.GetSoundSpeed2();
    Entropy = gas.GetEntropy();
  }
};
}  // namespace

TEST_CASE("Small pressure changes take a single integration step", "[Pressure-outlet]") {
  // A converged outlet imposes almost the interior pressure: one RK4 step (4 stages and a final check), plus the
  // interior and boundary states.
  auto gas = MakeOutletTestGas(1);
  CountingFluid counting(*gas);
  const su2double p = 1500112, rho = 230, normal[2] = {1, 0};
  gas->SetTDState_Prho(p, rho);
  const su2double e = gas->GetStaticEnergy(), a = gas->GetSoundSpeed(),
                  state[4] = {rho, rho * .3 * a, 0, rho * (e + .045 * a * a)};
  for (const auto ratio : {1.001, .999, .995}) {
    CAPTURE(ratio);
    su2double boundary[4];
    counting.calls = 0;
    REQUIRE(CEulerSolver::ComputePressureOutletState(&counting, 2, state, normal, 0, ratio * p, 0, false, 0,
                                                     boundary) == CEulerSolver::PressureOutletStatus::Valid);
    CHECK(counting.calls == 7);
    gas->SetTDState_rhoe(boundary[0], boundary[3] / boundary[0] - .5 * pow(boundary[1] / boundary[0], 2));
    CHECK(SU2_TYPE::GetValue(gas->GetPressure()) == Approx(SU2_TYPE::GetValue(ratio * p)).epsilon(1e-9));
  }
}

TEST_CASE("Sonic backflow fallback keeps the interior entropy, tangential velocity and imposed pressure",
          "[Pressure-outlet]") {
  using Status = CEulerSolver::PressureOutletStatus;
  for (auto dim : {2u, 3u})
    for (auto kind : {0u, 1u, 2u}) {
      auto gas = MakeOutletTestGas(kind);
      const su2double p = 1500112, rho = 230, normal[3] = {1, 0, 0};
      gas->SetTDState_Prho(p, rho);
      const su2double e = gas->GetStaticEnergy(), a = gas->GetSoundSpeed(), s = gas->GetEntropy();
      for (const su2double grid : {0., 20.})
        for (const auto ratio : {.95, 1.1}) {
          CAPTURE(dim, kind, grid, ratio);
          su2double velocity[3] = {-1.4 * a + grid, 5, -2}, state[5] = {rho}, out[5], unused[5];
          for (auto d = 0u; d < dim; ++d) state[d + 1] = rho * velocity[d];
          state[dim + 1] = rho * (e + .5 * GeometryToolbox::SquaredNorm(dim, velocity));
          REQUIRE(CEulerSolver::ComputePressureOutletState(gas.get(), dim, state, normal, grid, ratio * p, 0, kind == 0,
                                                           1.4, unused) == Status::Backflow);
          REQUIRE(CEulerSolver::ComputeSonicBackflowState(gas.get(), dim, state, normal, grid, ratio * p, 0, kind == 0,
                                                          1.4, out) == Status::Valid);
          su2double v[3] = {};
          for (auto d = 0u; d < dim; ++d) v[d] = out[d + 1] / out[0];
          gas->SetTDState_rhoe(out[0], out[dim + 1] / out[0] - .5 * GeometryToolbox::SquaredNorm(dim, v));
          CHECK(SU2_TYPE::GetValue(gas->GetPressure()) == Approx(SU2_TYPE::GetValue(ratio * p)).epsilon(1e-7));
          CHECK(SU2_TYPE::GetValue(gas->GetEntropy()) == Approx(SU2_TYPE::GetValue(s)).margin(1e-5));
          // Inflow relative to the grid, never faster than sound. The ideal gas reaches the sonic limit here.
          const su2double q = v[0] - grid, sound = gas->GetSoundSpeed();
          CHECK(q < 0);
          CHECK(SU2_TYPE::GetValue(-q) <= SU2_TYPE::GetValue(sound) * (1 + 1e-9));
          if (kind == 0) CHECK(SU2_TYPE::GetValue(q) == Approx(SU2_TYPE::GetValue(-sound)).epsilon(1e-9));
          for (auto d = 1u; d < dim; ++d)
            CHECK(SU2_TYPE::GetValue(v[d]) == Approx(SU2_TYPE::GetValue(velocity[d])).epsilon(1e-12));
        }
    }
}

TEST_CASE("Sonic backflow fallback is continuous with subsonic reverse flow", "[Pressure-outlet]") {
  // Bisect the interior normal speed where the outlet state switches from subsonic reverse flow to flow that enters
  // through every characteristic; the fallback on one side must match the outlet state on the other side.
  using Status = CEulerSolver::PressureOutletStatus;
  for (auto kind : {0u, 1u, 2u})
    for (const auto ratio : {1.05, .95}) {
      CAPTURE(kind, ratio);
      auto gas = MakeOutletTestGas(kind);
      const su2double p = 1500112, rho = 230, normal[2] = {1, 0};
      gas->SetTDState_Prho(p, rho);
      const su2double e = gas->GetStaticEnergy(), a = gas->GetSoundSpeed(), target = ratio * p;
      auto interior = [&](double mach, su2double* u) {
        u[0] = rho;
        u[1] = rho * mach * a;
        u[2] = rho * 3;
        u[3] = rho * (e + .5 * (mach * mach * a * a + 9));
      };
      double lower = -1.5, upper = -.5;
      su2double u[4], subsonic[4], fallback[4];
      for (auto iteration = 0u; iteration < 50; ++iteration) {
        const double middle = .5 * (lower + upper);
        interior(middle, u);
        if (CEulerSolver::ComputePressureOutletState(gas.get(), 2, u, normal, 0, target, 0, kind == 0, 1.4, subsonic) ==
            Status::Backflow)
          lower = middle;
        else
          upper = middle;
      }
      interior(upper, u);
      REQUIRE(CEulerSolver::ComputePressureOutletState(gas.get(), 2, u, normal, 0, target, 0, kind == 0, 1.4,
                                                       subsonic) == Status::Valid);
      interior(lower, u);
      REQUIRE(CEulerSolver::ComputeSonicBackflowState(gas.get(), 2, u, normal, 0, target, 0, kind == 0, 1.4,
                                                      fallback) == Status::Valid);
      for (auto k = 0u; k < 4; ++k)
        CHECK(SU2_TYPE::GetValue(fallback[k]) ==
              Approx(SU2_TYPE::GetValue(subsonic[k])).epsilon(1e-6).margin(1e-6 * SU2_TYPE::GetValue(rho * a)));
    }
}

#ifdef HAVE_MLPCPP
TEST_CASE("Pressure outlet follows the learned EOS without pressure-density inversion", "[Pressure-outlet]") {
  // MLP_nitrogen_outlet.mlp maps (density, thermal energy) to entropy; its derivatives define the EOS.
  // This fixture checks that outlet integration can use that native input pair without a (p,rho) solve.
  std::stringstream options;
  options << "SOLVER= EULER\nCONV_NUM_METHOD_FLOW= ROE\nFLUID_MODEL= DATADRIVEN_FLUID\nINTERPOLATION_METHOD= "
             "MLP\nUSE_PINN= YES\n"
             "FILENAMES_INTERPOLATOR= (src/SU2/UnitTests/SU2_CFD/fluid/MLP_nitrogen_outlet.mlp)\n";
  CConfig config(options, SU2_COMPONENT::SU2_CFD, false);
  CDataDrivenFluid gas(&config, false);
  gas.SetTDState_rhoe(16, 130000);
  const su2double pressure = gas.GetPressure(), sound = gas.GetSoundSpeed(), entropy = gas.GetEntropy();
  const su2double normal[2] = {1, 0};
  for (const auto mach : {-.3, 0., .3}) {
    su2double state[4] = {16, 16 * mach * sound, 0, 16 * (130000 + .5 * mach * mach * sound * sound)}, boundary[4];
    REQUIRE(CEulerSolver::ComputePressureOutletState(&gas, 2, state, normal, 0, .98 * pressure, 0, false, 0,
                                                     boundary) == CEulerSolver::PressureOutletStatus::Valid);
    const su2double velocity = boundary[1] / boundary[0];
    gas.SetTDState_rhoe(boundary[0], boundary[3] / boundary[0] - .5 * velocity * velocity);
    CHECK(SU2_TYPE::GetValue(gas.GetPressure()) == Approx(SU2_TYPE::GetValue(.98 * pressure)).epsilon(1e-7));
    CHECK(SU2_TYPE::GetValue(gas.GetEntropy()) == Approx(SU2_TYPE::GetValue(entropy)).margin(1e-5));
    CHECK(gas.GetTemperature() > 0);
    CHECK(gas.GetSoundSpeed2() > 0);
  }
}
#endif
