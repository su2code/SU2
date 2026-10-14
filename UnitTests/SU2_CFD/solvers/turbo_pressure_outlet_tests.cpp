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

#include <iomanip>
#include "catch.hpp"
#include "pressure_outlet_test_utils.hpp"
#include "../../UnitQuadTestCase.hpp"
#include "../../../SU2_CFD/include/variables/CNSVariable.hpp"
#include "../../../SU2_CFD/include/fluid/CPengRobinson.hpp"
#include "../../../SU2_CFD/include/fluid/CVanDerWaalsGas.hpp"
#include "../../../SU2_CFD/include/numerics/flow/convection/hllc.hpp"
#include "../../../SU2_CFD/include/fluid/CFluidModel.hpp"
#include "../../../Common/include/toolboxes/geometry_toolbox.hpp"

namespace {
// Supply a known spanwise pressure profile to test how the boundary consumes it, independently
// of the radial-equilibrium profile update elsewhere in the solver.
class TurboOutletSolver final : public CEulerSolver {
 public:
  TurboOutletSolver(CGeometry* geometry, CConfig* config, double pressure, bool radial)
      : CEulerSolver(geometry, config, MESH_0, config->GetViscous()) {
    const auto spans = config->GetnSpanWiseSections();
    RadialEquilibriumPressure.resize(nMarker, spans);
    ExtAverageTurboVelocity.resize(nDim);
    for (auto marker = 0u; marker < nMarker; ++marker)
      for (auto span = 0u; span < spans; ++span)
        RadialEquilibriumPressure[marker][span] = pressure * (radial ? .98 + .02 * span : 1.);
  }
};

// Exercise the viscous boundary path with zero diffusion so the checks measure only convective fluxes.
class ZeroOutletDiffusion final : public CNumerics {
  su2double residual[5] = {}, matrix[5][5] = {};
  const su2double* rows[5];

 public:
  ZeroOutletDiffusion(unsigned short dim, const CConfig* config) : CNumerics(dim, dim + 2, config) {
    for (auto row = 0u; row < nVar; ++row) rows[row] = matrix[row];
  }
  ResidualType<> ComputeResidual(const CConfig*) override { return {residual, rows, rows}; }
};

// Supply fixed SST k to verify thermal-energy bookkeeping without advancing a turbulence model.
class FrozenOutletTurbulence final : public CSolver {
  CVariable variables;
  CVariable* GetBaseClassPointerToNodes() override { return &variables; }

 public:
  FrozenOutletTurbulence(unsigned long points, const CConfig* config, double k) : variables(points, 2, config) {
    SetBaseClassPointerToNodes();
    for (auto point = 0ul; point < points; ++point) variables.SetSolution(point, 0, k);
  }
};

void TurboPressureOutlet(unsigned short dim, unsigned short kind, double mach, double ratio, bool moving, bool radial,
                         bool implicit = true, bool totalConditions = false, bool threaded = false,
                         double turbulentK = 0) {
  CAPTURE(dim, kind, mach, ratio, moving, radial, implicit, totalConditions, threaded, turbulentK);
  // Axial turbo normals point along x in the 2D fixture and z in the 3D fixture.
  const unsigned short normalDim = dim == 2 ? 0 : 2;
  const std::string inlet = dim == 2 ? "x_minus" : "z_minus", outlet = dim == 2 ? "x_plus" : "z_plus";
  const double p = 1500112, rho = 230;
  UnitQuadTestCase test;
  std::ostringstream options;
  options
      << std::setprecision(17) << "SOLVER= " << (turbulentK ? "RANS" : "EULER")
      << "\nMESH_FORMAT= " << (dim == 2 ? "RECTANGLE" : "BOX")
      << "\nMESH_BOX_SIZE= " << (threaded ? "5,129,5" : "5,5,5")
      << "\nMESH_BOX_OFFSET= 1,0,0\nREF_DIMENSIONALIZATION= DIMENSIONAL\n"
      << "MACH_NUMBER= 0.3\nFREESTREAM_OPTION= DENSITY_FS\nFREESTREAM_DENSITY= 230\n"
      << "FREESTREAM_PRESSURE= 1500112\nFLUID_MODEL= "
      << (kind == 0   ? "IDEAL_GAS"
          : kind == 1 ? "PR_GAS"
                      : "VW_GAS")
      << "\nGAMMA_VALUE= " << (kind == 0 ? "1.4" : "1.0125")
      << "\nGAS_CONSTANT= 35.152\nCRITICAL_PRESSURE= 1415200\nCRITICAL_TEMPERATURE= 564.1\nACENTRIC_FACTOR= 0.529\n"
      << "CONV_NUM_METHOD_FLOW= HLLC\nTIME_DISCRE_FLOW= " << (implicit ? "EULER_IMPLICIT" : "RUNGE-KUTTA_EXPLICIT")
      << "\nMARKER_EULER= " << (dim == 2 ? "(y_minus,y_plus)" : "(x_minus,x_plus)") << "\nMARKER_MONITORING= (" << inlet
      << ")\nMARKER_TURBOMACHINERY= (" << inlet << ',' << outlet << ")\n"
      << "TURBOMACHINERY_KIND= AXIAL\nTURBO_PERF_KIND= TURBINE\nSPANWISE_KIND= EQUISPACED\nNUM_SPANWISE_SECTIONS= 3\n"
      << "MARKER_RIEMANN= (" << inlet << ",TOTAL_CONDITIONS_PT,1600000,600,1,0,0," << outlet << ','
      << (totalConditions ? "TOTAL_CONDITIONS_PT"
          : radial        ? "RADIAL_EQUILIBRIUM"
                          : "STATIC_PRESSURE")
      << ',' << ratio * p << (totalConditions ? ",600,1,0,0)\n" : ",0,0,0,0)\n");
  if (moving) options << "GRID_MOVEMENT= STEADY_TRANSLATION\nTRANSLATION_RATE= 20,-7,30\n";
  if (dim == 3) options << "MARKER_PERIODIC= (y_minus,y_plus,0,0,0,0,0,0,0,1,0)\n";
  if (turbulentK)
    options << "KIND_TURB_MODEL= SST\nINIT_OPTION= TD_CONDITIONS\nREYNOLDS_NUMBER= 1E6\n"
               "VISCOSITY_MODEL= CONSTANT_VISCOSITY\nCONDUCTIVITY_MODEL= CONSTANT_CONDUCTIVITY\n";
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
  auto* solver = new TurboOutletSolver(test.geometry.get(), config, ratio * p, radial);
  test.solver[FLOW_SOL] = solver;
  std::cout.rdbuf(test.orig_buf);
  std::unique_ptr<CFluidModel> gas;
  if (kind == 0) gas.reset(new CIdealGas(1.4, 35.152));
  if (kind == 1) gas.reset(new CPengRobinson(1.0125, 35.152, 1415200, 564.1, .529));
  if (kind == 2) gas.reset(new CVanDerWaalsGas(1.0125, 35.152, 1415200, 564.1));
  if (turbulentK) {
    gas->SetLaminarViscosityModel(config);
    gas->SetThermalConductivityModel(config);
  }
  gas->SetTDState_Prho(p, rho);
  const su2double sound = gas->GetSoundSpeed(), e = gas->GetStaticEnergy();
  const su2double grid[3] = {moving ? 20. : 0., moving ? -7. : 0., moving ? 30. : 0.};
  su2double state[5] = {rho}, velocity[3] = {3, -2, 1};
  for (auto d = 0u; d < dim; ++d) {
    velocity[d] = (d == normalDim ? mach * sound : velocity[d]) + grid[d];
    state[d + 1] = rho * velocity[d];
  }
  state[dim + 1] = rho * (e + .5 * GeometryToolbox::SquaredNorm(dim, velocity) + turbulentK);
  std::unique_ptr<CSolver> turbulence;
  if (turbulentK) {
    turbulence.reset(new FrozenOutletTurbulence(test.geometry->GetnPoint(), config, turbulentK));
    test.solver[TURB_SOL] = turbulence.get();
  }
  auto* nodes = static_cast<CEulerVariable*>(solver->GetNodes());
  for (auto point = 0ul; point < test.geometry->GetnPoint(); ++point) {
    nodes->SetSolution(point, state);
    if (turbulentK) {
      auto* viscousNodes = static_cast<CNSVariable*>(nodes);
      REQUIRE(viscousNodes->SetPrimVar(point, 0, turbulentK, gas.get()));
      viscousNodes->SetSecondaryVar(point, gas.get());
    } else {
      REQUIRE(nodes->SetPrimVar(point, gas.get()));
      nodes->SetSecondaryVar(point, gas.get());
    }
    if (moving) test.geometry->nodes->SetGridVel(point, grid);
    for (auto k = 0u; k < dim + 2; ++k) solver->LinSysRes(point, k) = 0;
    su2double zero[25] = {};
    solver->Jacobian.SetBlock(point, point, zero);
  }
  for (auto marker = 0u; marker < config->GetnMarker_All(); ++marker) {
    if (config->GetMarker_All_TagBound(marker) != outlet) continue;
    const auto spans = test.geometry->GetnSpanWiseSections(OUTFLOW);
    auto target = [&](unsigned short span) { return ratio * p * (radial ? .98 + .02 * span : 1.); };
    SU2_OMP_PARALLEL {
      CUpwGeneralHLLC_Flow numerics(dim, dim + 2, config);
      ZeroOutletDiffusion diffusion(dim, config);
      solver->BC_TurboRiemann(test.geometry.get(), test.solver, &numerics, &diffusion, config, marker);
    }
    CHECK(solver->GetFluidModel()->GetPressure() > 0);
    CUpwGeneralHLLC_Flow numerics(dim, dim + 2, config);
    for (auto span = 0u; span < spans; ++span)
      for (auto vertex = 0ul; vertex < test.geometry->GetnVertexSpan(marker, span); ++vertex) {
        const auto* turbo = test.geometry->turbovertex[marker][span][vertex];
        const auto point = turbo->GetNode(), oldVertex = turbo->GetOldVertex();
        su2double normal[3], unit[3];
        test.geometry->vertex[marker][oldVertex]->GetNormal(normal);
        turbo->GetNormal(unit);
        for (auto d = 0u; d < dim; ++d) normal[d] = -normal[d];
        su2double boundary[5];
        REQUIRE(gas->ComputePressureOutletState(dim, state, unit, GeometryToolbox::DotProduct(dim, grid, unit),
                                                target(span), boundary, turbulentK,
                                                kind == 0 ? 1.4 : 0) == CFluidModel::PressureOutletStatus::Valid);
        auto flux = [&](const su2double* boundary, su2double* result) {
          su2double v[3];
          for (auto d = 0u; d < dim; ++d) v[d] = boundary[d + 1] / boundary[0];
          const su2double energy = boundary[dim + 1] / boundary[0];
          gas->SetTDState_rhoe(boundary[0], energy - .5 * GeometryToolbox::SquaredNorm(dim, v) - turbulentK);
          const su2double pressure = gas->GetPressure(), h = energy + pressure / boundary[0];
          numerics.GetInviscidProjFlux(&boundary[0], v, &pressure, &h, normal, result);
          // The moving-grid transport term occurs once, in both explicit and implicit boundary paths.
          const su2double speed = GeometryToolbox::DotProduct(dim, grid, normal);
          for (auto k = 0u; k < dim + 2; ++k) result[k] -= speed * boundary[k];
        };
        su2double expected[5];
        flux(boundary, expected);
        for (auto row = 0u; row < dim + 2; ++row) {
          CHECK(std::isfinite(SU2_TYPE::GetValue(solver->LinSysRes(point, row))));
          CHECK(SU2_TYPE::GetValue(solver->LinSysRes(point, row)) ==
                Approx(SU2_TYPE::GetValue(expected[row])).epsilon(1e-8).margin(1e-7));
        }
        const auto* primitive = solver->GetCharacPrimVar(marker, oldVertex);
        CHECK(SU2_TYPE::GetValue(primitive[0]) == Approx(SU2_TYPE::GetValue(gas->GetTemperature())).epsilon(1e-8));
        CHECK(SU2_TYPE::GetValue(primitive[dim + 1]) == Approx(SU2_TYPE::GetValue(gas->GetPressure())).epsilon(1e-8));
        CHECK(SU2_TYPE::GetValue(primitive[dim + 2]) == Approx(SU2_TYPE::GetValue(gas->GetDensity())).epsilon(1e-8));
        CHECK(SU2_TYPE::GetValue(primitive[dim + 4]) == Approx(SU2_TYPE::GetValue(gas->GetSoundSpeed())).epsilon(1e-8));
        if (mach == .3 && ratio > .1)
          CHECK(SU2_TYPE::GetValue(primitive[dim + 1]) == Approx(target(span)).epsilon(1e-8));
        if (!implicit) continue;
        // The implicit matrix uses a frozen characteristic projection, with local span pressure and SST k.
        for (auto column = 0u; column < dim + 2; ++column) {
          su2double plus[5], minus[5], fp[5], fm[5], increment[5];
          OutgoingOutletIncrement(*gas, dim, state, unit, GeometryToolbox::DotProduct(dim, grid, unit), turbulentK,
                                  column, increment);
          const double delta =
              2e-5 * fmax(fabs(SU2_TYPE::GetValue(state[column])), column == 0 ? rho : rho * SU2_TYPE::GetValue(sound));
          for (auto k = 0u; k < dim + 2; ++k) {
            plus[k] = boundary[k] + delta * increment[k];
            minus[k] = boundary[k] - delta * increment[k];
          }
          flux(plus, fp);
          flux(minus, fm);
          for (auto row = 0u; row < dim + 2; ++row)
            CHECK(solver->Jacobian.GetBlock(point, point)[row * (dim + 2) + column] ==
                  Approx(SU2_TYPE::GetValue((fp[row] - fm[row]) / (2 * delta))).epsilon(5e-4).margin(2e-5));
        }
      }
  }
}
}  // namespace

TEST_CASE("TurboRiemann pressure outlets preserve admissibility and the ALE flux", "[Pressure-outlet][Turbo-outlet]") {
  for (auto dim : {2u, 3u})
    for (auto kind : {0u, 1u, 2u})
      for (const auto radial : {false, true})
        for (const auto moving : {false, true})
          for (const auto mach : {-.3, 0., .3, .99, 1.4})
            for (const auto ratio : {.95, .0006666168959384})
              TurboPressureOutlet(dim, kind, mach, ratio, moving, radial);
}

TEST_CASE("Explicit TurboRiemann pressure outlet uses one grid flux correction", "[Pressure-outlet][Turbo-outlet]") {
  TurboPressureOutlet(2, 0, 1.4, .95, true, false, false);
}

TEST_CASE("TurboRiemann real-gas startup pressure remains admissible", "[Pressure-outlet][Turbo-startup]") {
  // This ratio prescribes 1000 Pa from an interior pressure of 1500112 Pa.
  TurboPressureOutlet(2, 1, .3, .0006666168959384, false, false);
}

TEST_CASE("Implicit TurboRiemann pressure outlet uses one grid flux correction", "[Pressure-outlet][Turbo-moving]") {
  TurboPressureOutlet(2, 0, 1.4, .95, true, false);
}

TEST_CASE("TurboRiemann pressure reconstruction is thread local", "[Pressure-outlet][Turbo-outlet]") {
  // Use enough outlet vertices to enter the OpenMP work-sharing loop with separate EOS instances.
  TurboPressureOutlet(2, 1, .3, .95, true, false, true, false, true);
}

TEST_CASE("Other TurboRiemann modes retain one ALE correction", "[Pressure-outlet][Turbo-outlet]") {
  for (const auto implicit : {false, true}) TurboPressureOutlet(2, 0, 1.4, .95, true, false, implicit, true);
}

TEST_CASE("TurboRiemann pressure outlet rejects nonpositive pressure", "[.][Turbo-invalid-pressure]") {
  TurboPressureOutlet(2, 1, .3, 0, false, false);
}

TEST_CASE("TurboRiemann radial pressure outlet requires data for backflow", "[.][Turbo-backflow]") {
  TurboPressureOutlet(3, 1, -1.4, .95, true, true);
}

TEST_CASE("TurboRiemann real-gas radial startup remains admissible", "[Pressure-outlet][Turbo-radial-startup]") {
  TurboPressureOutlet(3, 1, .3, .0006666168959384, false, true);
}

TEST_CASE("TurboRiemann pressure outlet separates frozen turbulence energy from the EOS",
          "[Pressure-outlet][Turbo-outlet]") {
  for (auto kind : {0u, 1u}) TurboPressureOutlet(3, kind, .3, .95, true, true, true, false, false, 13);
}
