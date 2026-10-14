/*!
 * \file characteristic_pressure_outlet_tests.cpp
 * \brief Regression for characteristic pressure boundary admissibility.
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
#include "../../UnitQuadTestCase.hpp"
#include "../../../SU2_CFD/include/numerics/flow/convection/hllc.hpp"
#include "../../../SU2_CFD/include/solvers/characteristic_pressure_outlet.hpp"
#include "../../../SU2_CFD/include/fluid/CVanDerWaalsGas.hpp"
#include "../../../SU2_CFD/include/fluid/CPengRobinson.hpp"
#include "../../../SU2_CFD/include/fluid/CDataDrivenFluid.hpp"

namespace {
void CharacteristicPressureFixture(su2double pressure, bool backflow = false) {
  UnitQuadTestCase test;
  std::ostringstream options;
  options << "SOLVER= EULER\nMESH_FORMAT= RECTANGLE\nMESH_BOX_SIZE= 5,129,5\n"
          << "REF_DIMENSIONALIZATION= DIMENSIONAL\nMACH_NUMBER= 0.1\n"
          << "FREESTREAM_OPTION= DENSITY_FS\nFREESTREAM_DENSITY= 230\nFREESTREAM_PRESSURE= 1500112\n"
          << "FLUID_MODEL= PR_GAS\nGAMMA_VALUE= 1.0125\nGAS_CONSTANT= 35.152\n"
          << "CRITICAL_PRESSURE= 1415200\nCRITICAL_TEMPERATURE= 564.1\nACENTRIC_FACTOR= 0.529\n"
          << "CONV_NUM_METHOD_FLOW= HLLC\nTIME_DISCRE_FLOW= EULER_IMPLICIT\n"
          << "GRID_MOVEMENT= STEADY_TRANSLATION\nTRANSLATION_RATE= 20,0,0\n"
          << "MARKER_EULER= (x_minus,y_minus,y_plus)\nMARKER_MONITORING= (x_minus)\n"
          << "MARKER_RIEMANN= (x_plus,STATIC_PRESSURE," << pressure << ",0,0,0,0)\n";
  test.config_options = options.str();
  test.InitConfig();
  test.InitGeometry();
  test.InitSolver();
  auto* solver = test.solver[FLOW_SOL];
  auto* nodes = static_cast<CEulerVariable*>(solver->GetNodes());
  CPengRobinson gas(1.0125, 35.152, 1415200, 564.1, 0.529);
  // Captured conservative startup state: low back pressure previously produced an invalid EOS state.
  su2double state[4] = {139.72158881490299, 11271.637336535776, -2338.5053086132339, 221781581.74179047};
  if (backflow) {
    const auto oldVelocity = state[1] / state[0];
    state[1] = state[0] * -200;
    // Change only normal velocity; adjust total energy to keep the original thermal state.
    state[3] += .5 * state[0] * (40000 - oldVelocity * oldVelocity);
  }
  const su2double grid[2] = {20, 0};
  for (auto p = 0ul; p < test.geometry->GetnPoint(); ++p) {
    nodes->SetSolution(p, state);
    nodes->SetSolution_Old(p, state);
    REQUIRE(nodes->SetPrimVar(p, &gas));
    nodes->SetSecondaryVar(p, &gas);
    test.geometry->nodes->SetGridVel(p, grid);
    for (auto k = 0u; k < 4; ++k) solver->LinSysRes(p, k) = 0;
    su2double zero[16] = {};
    solver->Jacobian.SetBlock(p, p, zero);
  }
  CUpwGeneralHLLC_Flow numerics(2, 4, test.config.get());
  for (auto marker = 0u; marker < test.config->GetnMarker_All(); ++marker)
    if (test.config->GetMarker_All_TagBound(marker) == "x_plus") SU2_OMP_PARALLEL {
        CUpwGeneralHLLC_Flow threadNumerics(2, 4, test.config.get());
        solver->BC_Riemann(test.geometry.get(), test.solver, &threadNumerics, nullptr, test.config.get(), marker);
      }
  for (auto p = 0ul; p < test.geometry->GetnPoint(); ++p)
    for (auto k = 0u; k < 4; ++k) CHECK(std::isfinite(SU2_TYPE::GetValue(solver->LinSysRes(p, k))));
  for (auto marker = 0u; marker < test.config->GetnMarker_All(); ++marker) {
    if (test.config->GetMarker_All_TagBound(marker) != "x_plus") continue;
    for (auto v = 0ul; v < test.geometry->nVertex[marker]; ++v) {
      const auto point = test.geometry->vertex[marker][v]->GetNode();
      if (!test.geometry->nodes->GetDomain(point)) continue;
      su2double normal[2], unit[2];
      test.geometry->vertex[marker][v]->GetNormal(normal);
      const auto area = GeometryToolbox::Norm(2, normal);
      for (auto d = 0u; d < 2; ++d) unit[d] = (normal[d] = -normal[d]) / area;
      // BC_Riemann uses the physical boundary flux, including transport by the moving grid.
      su2double boundary[4];
      REQUIRE(PressureOutlet::State(gas, 2, state, unit, GeometryToolbox::DotProduct(2, grid, unit), pressure,
                                    boundary) == PressureOutlet::Status::Valid);
      auto flux = [&](const su2double* boundary, su2double* output) {
        su2double velocity[2];
        for (auto d = 0u; d < 2; ++d) velocity[d] = boundary[d + 1] / boundary[0];
        const auto energy = boundary[3] / boundary[0];
        gas.SetTDState_rhoe(boundary[0], energy - .5 * GeometryToolbox::SquaredNorm(2, velocity));
        const auto p = gas.GetPressure(), h = energy + p / boundary[0];
        numerics.GetInviscidProjFlux(&boundary[0], velocity, &p, &h, normal, output);
        const auto speed = GeometryToolbox::DotProduct(2, grid, normal);
        for (auto k = 0u; k < 4; ++k) output[k] -= speed * boundary[k];
      };
      su2double reference[4];
      flux(boundary, reference);
      const auto* primitive = solver->GetCharacPrimVar(marker, v);
      CHECK(SU2_TYPE::GetValue(primitive[0]) == Approx(SU2_TYPE::GetValue(gas.GetTemperature())).epsilon(1e-8));
      CHECK(SU2_TYPE::GetValue(primitive[3]) == Approx(SU2_TYPE::GetValue(gas.GetPressure())).epsilon(1e-8));
      CHECK(SU2_TYPE::GetValue(primitive[4]) == Approx(SU2_TYPE::GetValue(gas.GetDensity())).epsilon(1e-8));
      CHECK(SU2_TYPE::GetValue(primitive[6]) == Approx(SU2_TYPE::GetValue(gas.GetSoundSpeed())).epsilon(1e-8));
      for (auto k = 0u; k < 4; ++k)
        CHECK(SU2_TYPE::GetValue(solver->LinSysRes(point, k)) ==
              Approx(SU2_TYPE::GetValue(reference[k])).epsilon(1e-8));
      for (auto column = 0u; column < 4; ++column) {
        // Freeze the wave basis and differentiate the physical boundary flux along outgoing increments.
        su2double plus[4], minus[4], fPlus[4], fMinus[4], increment[4];
        OutgoingOutletIncrement(gas, 2, state, unit, GeometryToolbox::DotProduct(2, grid, unit), 0, column, increment);
        const auto delta = 2e-5 * fabs(state[column]);
        for (auto k = 0u; k < 4; ++k) {
          plus[k] = boundary[k] + delta * increment[k];
          minus[k] = boundary[k] - delta * increment[k];
        }
        flux(plus, fPlus);
        flux(minus, fMinus);
        for (auto row = 0u; row < 4; ++row) {
          const auto expected = SU2_TYPE::GetValue((fPlus[row] - fMinus[row]) / (2 * delta));
          CHECK(solver->Jacobian.GetBlock(point, point)[row * 4 + column] ==
                Approx(expected).epsilon(3e-4).margin(1e-5));
        }
      }
    }
  }
}
}  // namespace

TEST_CASE("Characteristic pressure boundary reconstructs the captured startup states", "[Pressure-outlet]") {
  for (const auto pressure : {1400000., 100000., 1000.}) CharacteristicPressureFixture(pressure);
}

TEST_CASE("Nonlinear pressure outflow preserves the acoustic invariant and frame", "[Pressure-outlet]") {
  for (auto dim : {2u, 3u}) {
    for (auto kind : {0u, 1u, 2u}) {
      std::unique_ptr<CFluidModel> gas;
      if (kind == 0) gas.reset(new CIdealGas(1.4, 35.152));
      if (kind == 1) gas.reset(new CPengRobinson(1.0125, 35.152, 1415200, 564.1, .529));
      if (kind == 2) gas.reset(new CVanDerWaalsGas(1.0125, 35.152, 1415200, 564.1));
      const su2double p = 1500112, rho = 230, normal[3] = {1, 0, 0}, boost[3] = {20, -7, 3};
      gas->SetTDState_Prho(p, rho);
      const auto e = gas->GetStaticEnergy(), sound = gas->GetSoundSpeed(), entropy = gas->GetEntropy();
      // The smallest ratio is 1000/1500112, reproducing the startup expansion past the sonic limit.
      for (const auto ratio : {1., .9, 1.05, .0006666168959384}) {
        su2double out[2][5] = {};
        // Boost both flow and grid: the relative wave speeds and reconstructed thermodynamics must agree.
        for (auto frame = 0u; frame < 2; ++frame) {
          su2double velocity[3] = {.3 * sound, 5, -2}, state[5] = {rho};
          for (auto d = 0u; d < dim; ++d) {
            if (frame) velocity[d] += boost[d];
            state[d + 1] = rho * velocity[d];
          }
          state[dim + 1] = rho * (e + .5 * GeometryToolbox::SquaredNorm(dim, velocity));
          REQUIRE(PressureOutlet::State(*gas, dim, state, normal, frame ? boost[0] : 0, ratio * p, out[frame]) ==
                  PressureOutlet::Status::Valid);
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
      REQUIRE(PressureOutlet::State(*gas, dim, state, normal, 0, 1000, unchanged) == PressureOutlet::Status::Valid);
      for (auto v = 0u; v < dim + 2; ++v) CHECK(unchanged[v] == state[v]);
      // Frozen SST k must remain in total energy without entering the EOS.
      state[1] = .3 * sound * rho;
      state[dim + 1] = rho * (e + .5 * .3 * .3 * sound * sound);
      su2double plain[5], withK[5];
      REQUIRE(PressureOutlet::State(*gas, dim, state, normal, 0, .9 * p, plain) == PressureOutlet::Status::Valid);
      state[dim + 1] += 17.5 * rho;
      REQUIRE(PressureOutlet::State(*gas, dim, state, normal, 0, .9 * p, withK, 17.5) == PressureOutlet::Status::Valid);
      for (auto v = 0u; v < dim + 1; ++v) CHECK(withK[v] == Approx(SU2_TYPE::GetValue(plain[v])).epsilon(1e-8));
      CHECK(withK[dim + 1] == Approx(SU2_TYPE::GetValue(plain[dim + 1] + 17.5 * plain[0])).epsilon(1e-8));
      // Reversing subsonic flow retains the pressure-only closure; fully incoming flow needs inlet data.
      state[1] = -state[1];
      CHECK(PressureOutlet::State(*gas, dim, state, normal, 0, p, unchanged, 17.5) == PressureOutlet::Status::Valid);
      state[1] = -2 * sound * rho;
      state[dim + 1] = rho * (e + 2 * sound * sound + 17.5);
      CHECK(PressureOutlet::State(*gas, dim, state, normal, 0, p, unchanged, 17.5) == PressureOutlet::Status::Backflow);
    }
  }
}

// Hidden Catch cases terminate the process; the Python checker runs them separately to inspect the error.
TEST_CASE("Characteristic pressure boundary rejects nonpositive pressure", "[.][Boundary-invalid-pressure]") {
  CharacteristicPressureFixture(0);
}
TEST_CASE("Characteristic pressure boundary requires data for backflow", "[.][Boundary-backflow]") {
  CharacteristicPressureFixture(1400000, true);
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
    REQUIRE(PressureOutlet::State(gas, 2, state, normal, 0, .98 * pressure, boundary) == PressureOutlet::Status::Valid);
    const su2double velocity = boundary[1] / boundary[0];
    gas.SetTDState_rhoe(boundary[0], boundary[3] / boundary[0] - .5 * velocity * velocity);
    CHECK(SU2_TYPE::GetValue(gas.GetPressure()) == Approx(SU2_TYPE::GetValue(.98 * pressure)).epsilon(1e-7));
    CHECK(SU2_TYPE::GetValue(gas.GetEntropy()) == Approx(SU2_TYPE::GetValue(entropy)).margin(1e-5));
    CHECK(gas.GetTemperature() > 0);
    CHECK(gas.GetSoundSpeed2() > 0);
  }
}
#endif
