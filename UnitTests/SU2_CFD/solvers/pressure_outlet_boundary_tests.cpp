/*!
 * \file pressure_outlet_boundary_tests.cpp
 * \brief Pressure outlet classification and boundary Jacobian regression.
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
// BC_Outlet with an ideal gas, normal Mach number relative to the grid and pressure ratio to the interior.
void IdealOutlet(unsigned short dim, double normalMach, double pressureRatio, bool moving, bool threaded = false,
                 double turbulentK = 0) {
  UnitQuadTestCase test;
  std::ostringstream options;
  options << std::setprecision(17) << OutletTestOptions(dim, 0, turbulentK, true, threaded ? "5,17,17" : "5,5,5")
          << "MARKER_EULER= (x_minus,"
          << "y_minus,y_plus" << (dim == 3 ? ",z_minus,z_plus" : "") << ")\nMARKER_MONITORING= (x_minus)\n"
          << "MARKER_OUTLET= (x_plus," << pressureRatio * OUTLET_TEST_P << ")\n";
  if (moving) options << "GRID_MOVEMENT= STEADY_TRANSLATION\nTRANSLATION_RATE= 500,240,0\n";
  test.config_options = options.str();
  test.InitConfig();
  test.InitGeometry();
  test.InitSolver();
  auto* solver = test.solver[FLOW_SOL];
  auto gas = MakeOutletTestGas(0);
  const double p = OUTLET_TEST_P, a = sqrt(1.4 * p / OUTLET_TEST_RHO);
  const su2double grid[3] = {moving ? 500. : 0., moving ? 240. : 0., 0.};
  // A tangential Mach number of 1.5 makes the total speed supersonic even when the normal one is subsonic; the grid
  // velocity separates the lab-frame speed from the normal speed relative to the grid.
  const su2double velocity[3] = {normalMach * a + grid[0], 1.5 * a + grid[1], 0};
  su2double state[5];
  OutletTestState(*gas, dim, velocity, turbulentK, state);
  SetOutletTestState(test, *gas, state, moving ? grid : nullptr, turbulentK);
  const auto marker = OutletTestMarker(*test.config, "x_plus");
  SU2_OMP_PARALLEL {
    OutletFluxProbe probe(dim, test.config.get());
    solver->BC_Outlet(test.geometry.get(), test.solver, &probe, nullptr, test.config.get(), marker);
  }
  // Independent ideal-gas reference: gamma=1.4 gives q+5a=constant and rho proportional to a^5.
  auto reference = [&](const su2double* u, su2double* out) {
    const su2double density = u[0], q = u[1] / density - grid[0];
    su2double v[3] = {};
    for (auto d = 0u; d < dim; ++d) v[d] = u[d + 1] / density;
    const su2double pressure = .4 * (u[dim + 1] - density * (.5 * GeometryToolbox::SquaredNorm(dim, v) + turbulentK));
    const su2double sound = sqrt(1.4 * pressure / density);
    if (q >= sound) {
      for (auto k = 0u; k < dim + 2u; ++k) out[k] = u[k];
      return;
    }
    const su2double boundarySound =
        fmax(sound * pow(pressureRatio * p / pressure, .4 / 2.8), (.4 * q + 2 * sound) / 2.4);
    out[0] = density * pow(boundarySound / sound, 5);
    v[0] += 5 * (sound - boundarySound);
    for (auto d = 0u; d < dim; ++d) out[d + 1] = out[0] * v[d];
    out[dim + 1] =
        out[0] * (boundarySound * boundarySound / .56 + .5 * GeometryToolbox::SquaredNorm(dim, v) + turbulentK);
  };
  // Full backflow and an invalid pressure use the steady fallback states.
  su2double expected[5];
  const su2double outward[3] = {1, 0, 0};
  const bool fallback = normalMach <= -1 || pressureRatio <= 0;
  if (fallback)
    CHECK(ExpectedSteadyOutletState(*gas, dim, state, outward, grid[0], pressureRatio * p, turbulentK, true, 1.4,
                                    expected));
  else
    reference(state, expected);
  unsigned long domainVertices = 0;
  for (auto vertex = 0ul; vertex < test.geometry->nVertex[marker]; ++vertex) {
    const auto point = test.geometry->vertex[marker][vertex]->GetNode();
    if (!test.geometry->nodes->GetDomain(point)) continue;
    ++domainVertices;
    for (auto row = 0u; row < dim + 2u; ++row)
      CHECK(SU2_TYPE::GetValue(solver->LinSysRes(point, row)) ==
            Approx(SU2_TYPE::GetValue(expected[row])).epsilon(1e-8).margin(1e-8));
    // Only the flux derivative with respect to the interior state (unit here) enters the implicit matrix.
    for (auto entry = 0u; entry < (dim + 2u) * (dim + 2u); ++entry)
      CHECK(solver->Jacobian.GetBlock(point, point)[entry] == (entry % (dim + 3) == 0 ? 1. : 0.));
  }
  CHECK(CountOutletFallbacks(solver, test.config.get()) == (fallback ? domainVertices : 0));
}
}  // namespace

TEST_CASE("Ordinary pressure outlet uses relative normal characteristics and the interior flux Jacobian",
          "[Pressure-outlet]") {
  for (auto dim : {2u, 3u})
    for (const auto moving : {false, true})
      for (const auto normalMach : {-.3, 0., .3, .99, 1.4})
        for (const auto ratio : {1., .95, 1.02, .001}) IdealOutlet(dim, normalMach, ratio, moving);
}

TEST_CASE("Ordinary pressure outlet uses the steady fallback states", "[Pressure-outlet]") {
  IdealOutlet(2, .3, 0, false);
  for (const auto ratio : {.95, 1.1}) IdealOutlet(2, -1.4, ratio, true);
}

TEST_CASE("Pressure outlet reconstruction is thread local", "[Pressure-outlet]") {
  IdealOutlet(3, .3, .95, true, true);
}

TEST_CASE("Ordinary pressure outlet keeps the SST energy out of the EOS", "[Pressure-outlet]") {
  for (auto dim : {2u, 3u})
    for (auto mach : {-.3, .3, 1.4}) IdealOutlet(dim, mach, .95, true, false, 1500);
}
