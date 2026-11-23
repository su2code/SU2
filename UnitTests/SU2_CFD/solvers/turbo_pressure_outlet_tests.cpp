/*!
 * \file turbo_pressure_outlet_tests.cpp
 * \brief Turbomachinery pressure-outlet reconstruction and ALE regression.
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

namespace {
// Supplies a known spanwise pressure profile, independently of the radial-equilibrium update elsewhere.
class TurboOutletSolver final : public CEulerSolver {
 public:
  TurboOutletSolver(CGeometry* geometry, CConfig* config, double pressure, bool radial)
      : CEulerSolver(geometry, config, MESH_0, config->GetViscous()) {
    RadialEquilibriumPressure.resize(nMarker, config->GetnSpanWiseSections());
    ExtAverageTurboVelocity.resize(nDim);
    SetOutletPressure(pressure, radial);
  }
  void SetOutletPressure(double pressure, bool radial) {
    for (auto marker = 0u; marker < nMarker; ++marker)
      for (auto span = 0u; span < RadialEquilibriumPressure.cols(); ++span)
        RadialEquilibriumPressure[marker][span] = pressure * (radial ? .98 + .02 * span : 1.);
  }
};

// Fixed SST k, without a turbulence model.
class FrozenOutletTurbulence final : public CSolver {
  CVariable variables;
  CVariable* GetBaseClassPointerToNodes() override { return &variables; }

 public:
  FrozenOutletTurbulence(unsigned long points, const CConfig* config) : variables(points, 2, config) {
    SetBaseClassPointerToNodes();
  }
};

// Build the mesh, turbo vertices and solver once, then check every listed Mach number and outlet pressure ratio.
void TurboPressureOutlet(unsigned short dim, unsigned short kind, std::initializer_list<double> machs,
                         std::initializer_list<double> ratios, bool moving, bool radial, bool implicit = true,
                         bool totalConditions = false, bool threaded = false, double turbulentK = 0) {
  CAPTURE(dim, kind, moving, radial, implicit, totalConditions, threaded, turbulentK);
  // Axial turbo normals point along x in 2D and along z in 3D.
  const unsigned short normalDim = dim == 2 ? 0 : 2;
  const std::string inlet = dim == 2 ? "x_minus" : "z_minus", outlet = dim == 2 ? "x_plus" : "z_plus";
  const double p = OUTLET_TEST_P, firstRatio = *ratios.begin();
  UnitQuadTestCase test;
  std::ostringstream options;
  options << OutletTestOptions(dim, kind, turbulentK, implicit, threaded ? "5,129,5" : "5,5,5")
          << "MESH_BOX_OFFSET= 1,0,0\nMARKER_EULER= " << (dim == 2 ? "(y_minus,y_plus)" : "(x_minus,x_plus)")
          << "\nMARKER_MONITORING= (" << inlet << ")\nMARKER_TURBOMACHINERY= (" << inlet << ',' << outlet << ")\n"
          << "TURBOMACHINERY_KIND= AXIAL\nTURBO_PERF_KIND= TURBINE\nSPANWISE_KIND= EQUISPACED\n"
          << "NUM_SPANWISE_SECTIONS= 3\nMARKER_RIEMANN= (" << inlet << ",TOTAL_CONDITIONS_PT,1600000,600,1,0,0,"
          << outlet << ','
          << (totalConditions ? "TOTAL_CONDITIONS_PT"
              : radial        ? "RADIAL_EQUILIBRIUM"
                              : "STATIC_PRESSURE")
          << ',' << firstRatio * p << (totalConditions ? ",600,1,0,0)\n" : ",0,0,0,0)\n");
  if (moving) options << "GRID_MOVEMENT= STEADY_TRANSLATION\nTRANSLATION_RATE= 20,-7,30\n";
  if (dim == 3) options << "MARKER_PERIODIC= (y_minus,y_plus,0,0,0,0,0,0,0,1,0)\n";
  test.config_options = options.str();
  test.InitConfig();
  test.InitGeometry();
  auto* config = test.config.get();
  std::ostringstream quiet;
  std::cout.rdbuf(quiet.rdbuf());
  if (dim == 3) {
    test.geometry->MatchPeriodic(config, 1);
    test.geometry->PreprocessPeriodicComms(test.geometry.get(), config);
  }
  for (auto flag : {INFLOW, OUTFLOW}) test.geometry->ComputeNSpan(config, 0, flag, true);
  for (auto flag : {INFLOW, OUTFLOW}) test.geometry->SetTurboVertex(config, 0, flag, true);
  config->SetnSpanMaxAllZones(config->GetnSpanWiseSections());
  test.solver = new CSolver*[MAX_SOLS]();
  auto* solver = new TurboOutletSolver(test.geometry.get(), config, firstRatio * p, radial);
  test.solver[FLOW_SOL] = solver;
  std::cout.rdbuf(test.orig_buf);
  std::unique_ptr<CSolver> turbulence;
  if (turbulentK) {
    turbulence.reset(new FrozenOutletTurbulence(test.geometry->GetnPoint(), config));
    test.solver[TURB_SOL] = turbulence.get();
  }
  auto gas = MakeOutletTestGas(kind);
  gas->SetTDState_Prho(p, OUTLET_TEST_RHO);
  const su2double sound = gas->GetSoundSpeed();
  const su2double grid[3] = {moving ? 20. : 0., moving ? -7. : 0., moving ? 30. : 0.};
  const auto marker = OutletTestMarker(*config, outlet);
  const auto spans = test.geometry->GetnSpanWiseSections(OUTFLOW);
  CUpwGeneralHLLC_Flow numerics(dim, dim + 2, config);
  for (const auto ratio : ratios) {
    // Only static-pressure and radial-equilibrium outlets read these pressures.
    config->SetPressureOut_BC(ratio * p / config->GetPressure_Ref());
    solver->SetOutletPressure(ratio * p, radial);
    auto target = [&](unsigned short span) { return ratio * p * (radial ? .98 + .02 * span : 1.); };
    for (const auto mach : machs) {
      CAPTURE(mach, ratio);
      su2double velocity[3] = {3, -2, 1}, state[5];
      velocity[normalDim] = mach * sound;
      for (auto d = 0u; d < dim; ++d) velocity[d] += grid[d];
      OutletTestState(*gas, dim, velocity, turbulentK, state);
      SetOutletTestState(test, *gas, state, moving ? grid : nullptr, turbulentK);
      SU2_OMP_PARALLEL {
        CUpwGeneralHLLC_Flow threadNumerics(dim, dim + 2, config);
        OutletFluxProbe diffusion(dim, config, true);
        solver->BC_TurboRiemann(test.geometry.get(), test.solver, &threadNumerics, &diffusion, config, marker);
      }
      unsigned long fallbacks = 0;
      for (auto span = 0u; span < spans; ++span)
        for (auto vertex = 0ul; vertex < test.geometry->GetnVertexSpan(marker, span); ++vertex) {
          const auto* turbo = test.geometry->turbovertex[marker][span][vertex];
          const auto point = turbo->GetNode(), oldVertex = turbo->GetOldVertex();
          su2double normal[3], unused[3], unit[3], boundary[5], flux[5];
          OutwardNormal(*test.geometry, marker, oldVertex, dim, normal, unused);
          turbo->GetNormal(unit);
          if (ExpectedSteadyOutletState(*gas, dim, state, unit, GeometryToolbox::DotProduct(dim, grid, unit),
                                        target(span), 0, kind == 0, 1.4, boundary))
            ++fallbacks;
          // The moving-grid transport term appears once, in the explicit and the implicit path.
          PhysicalOutletFlux(*gas, numerics, dim, boundary, normal, moving ? grid : nullptr, flux);
          for (auto row = 0u; row < dim + 2u; ++row)
            CHECK(SU2_TYPE::GetValue(solver->LinSysRes(point, row)) ==
                  Approx(SU2_TYPE::GetValue(flux[row])).epsilon(1e-8).margin(1e-7));
          const auto* primitive = solver->GetCharacPrimVar(marker, oldVertex);
          CheckOutletPrimitives(primitive, *gas, dim);
          if (mach == .3 && ratio > .1)
            CHECK(SU2_TYPE::GetValue(primitive[dim + 1]) == Approx(target(span)).epsilon(1e-8));
          if (implicit)
            CheckRiemannOutletJacobian(solver->Jacobian.GetBlock(point, point), *gas, numerics, dim, state, boundary,
                                       unit, normal, moving ? grid : nullptr);
        }
      // Full backflow and an invalid pressure use the steady fallback states, which the monitor counts.
      CHECK((fallbacks > 0) == (mach <= -1 || ratio <= 0));
      CHECK(CountOutletFallbacks(solver, config) == fallbacks);
    }
  }
}
}  // namespace

TEST_CASE("TurboRiemann pressure outlets preserve admissibility and the ALE flux", "[Pressure-outlet][Turbo-outlet]") {
  for (auto dim : {2u, 3u})
    for (auto kind : {0u, 1u, 2u})
      for (const auto radial : {false, true})
        for (const auto moving : {false, true})
          TurboPressureOutlet(dim, kind, {-.3, 0., .3, .99, 1.4}, {.95, .0006666168959384}, moving, radial);
}

TEST_CASE("Explicit TurboRiemann pressure outlet uses one grid flux correction", "[Pressure-outlet][Turbo-outlet]") {
  TurboPressureOutlet(2, 0, {1.4}, {.95}, true, false, false);
}

TEST_CASE("TurboRiemann real-gas startup pressure remains admissible", "[Pressure-outlet][Turbo-startup]") {
  // This ratio prescribes 1000 Pa from an interior pressure of 1500112 Pa.
  TurboPressureOutlet(2, 1, {.3}, {.0006666168959384}, false, false);
}

TEST_CASE("Implicit TurboRiemann pressure outlet uses one grid flux correction", "[Pressure-outlet][Turbo-moving]") {
  TurboPressureOutlet(2, 0, {1.4}, {.95}, true, false);
}

TEST_CASE("TurboRiemann pressure reconstruction is thread local", "[Pressure-outlet][Turbo-outlet]") {
  // Use enough outlet vertices to enter the OpenMP work-sharing loop with separate EOS instances.
  TurboPressureOutlet(2, 1, {.3}, {.95}, true, false, true, false, true);
}

TEST_CASE("Other TurboRiemann modes retain one ALE correction", "[Pressure-outlet][Turbo-outlet]") {
  for (const auto implicit : {false, true}) TurboPressureOutlet(2, 0, {1.4}, {.95}, true, false, implicit, true);
}

TEST_CASE("TurboRiemann pressure outlets use the steady fallback states", "[Pressure-outlet][Turbo-outlet]") {
  TurboPressureOutlet(2, 1, {.3}, {0}, false, false);
  TurboPressureOutlet(3, 1, {-1.4}, {.95}, true, true);
}

TEST_CASE("TurboRiemann real-gas radial startup remains admissible", "[Pressure-outlet][Turbo-radial-startup]") {
  TurboPressureOutlet(3, 1, {.3}, {.0006666168959384}, false, true);
}

TEST_CASE("TurboRiemann pressure outlet treats the SST energy like the other Riemann kinds",
          "[Pressure-outlet][Turbo-outlet]") {
  for (auto kind : {0u, 1u}) TurboPressureOutlet(3, kind, {.3}, {.95}, true, true, true, false, false, 13);
}
