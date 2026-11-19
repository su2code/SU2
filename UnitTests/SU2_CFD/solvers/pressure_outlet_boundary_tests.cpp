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
#include "../../UnitQuadTestCase.hpp"
#include "../../../SU2_CFD/include/fluid/CIdealGas.hpp"
#include "../../../SU2_CFD/include/variables/CNSVariable.hpp"

namespace {
// Unit flux matrices expose I + the outgoing-wave projection in the assembled implicit block.
class OutletStateProbe final : public CNumerics {
  su2double state[5] = {}, identity[5][5] = {};
  const su2double* rows[5];

 public:
  OutletStateProbe(unsigned short dim, const CConfig* config) : CNumerics(dim, dim + 2, config) {
    for (auto k = 0u; k < nVar; ++k) {
      rows[k] = identity[k];
      identity[k][k] = 1;
    }
  }
  ResidualType<> ComputeResidual(const CConfig*) override {
    const su2double rho = V_j[nDim + 2], p = V_j[nDim + 1];
    state[0] = rho;
    for (auto d = 0u; d < nDim; ++d) state[d + 1] = rho * V_j[d + 1];
    state[nDim + 1] = rho * V_j[nDim + 3] - p;
    return {state, rows, rows};
  }
};

// Backflow modes: with or without MARKER_OUTLET_BACKFLOW reservoir data, in steady or unsteady runs.
enum class Backflow { None, Reservoir, Unsteady, UnsteadyReservoir, ReservoirAtOutletPressure, ReservoirTangential };

void IdealOutlet(unsigned short dim, double normalMach, double pressureRatio, bool moving, bool threaded = false,
                 double turbulentK = 0, Backflow backflow = Backflow::None) {
  CAPTURE(dim, normalMach, pressureRatio, moving, threaded, turbulentK, static_cast<int>(backflow));
  UnitQuadTestCase test;
  std::ostringstream options;
  options << "SOLVER= " << (turbulentK ? "RANS" : "EULER") << "\nMESH_FORMAT= " << (dim == 2 ? "RECTANGLE" : "BOX")
          << "\nMESH_BOX_SIZE= " << (threaded ? "5,17,17" : "5,5,5") << "\nREF_DIMENSIONALIZATION= DIMENSIONAL\n"
          << "MACH_NUMBER= 0.3\nFLUID_MODEL= IDEAL_GAS\nGAMMA_VALUE= 1.4\nGAS_CONSTANT= 287\n"
          << "FREESTREAM_OPTION= DENSITY_FS\nFREESTREAM_DENSITY= 1.2\nFREESTREAM_PRESSURE= 100000\n"
          << "CONV_NUM_METHOD_FLOW= HLLC\nTIME_DISCRE_FLOW= EULER_IMPLICIT\n"
          << "MARKER_EULER= (x_minus,y_minus,y_plus" << (dim == 3 ? ",z_minus,z_plus" : "")
          << ")\nMARKER_MONITORING= (x_minus)\nMARKER_OUTLET= (x_plus," << pressureRatio * 100000 << ")\n";
  if (moving) options << "GRID_MOVEMENT= STEADY_TRANSLATION\nTRANSLATION_RATE= 500,240,0\n";
  const bool reservoir = backflow == Backflow::Reservoir || backflow == Backflow::UnsteadyReservoir ||
                         backflow == Backflow::ReservoirAtOutletPressure || backflow == Backflow::ReservoirTangential;
  // The tangential mode uses a zero total pressure: the imposed pressure at zero normal velocity.
  const bool tangentialMode = backflow == Backflow::ReservoirTangential;
  // A zero total pressure stands for the imposed outlet pressure.
  const double totalPressure = backflow == Backflow::ReservoirAtOutletPressure ? pressureRatio * 100000 : 101000;
  if (reservoir)
    options << "MARKER_OUTLET_BACKFLOW= (x_plus, 300, "
            << (backflow == Backflow::ReservoirAtOutletPressure || tangentialMode ? 0. : totalPressure)
            << ", 0, 0, 0)\n";
  if (tangentialMode) options << "OUTLET_BACKFLOW_TYPE= (INTERIOR_TANGENTIAL)\n";
  if (backflow == Backflow::Unsteady || backflow == Backflow::UnsteadyReservoir)
    options << "TIME_DOMAIN= YES\nTIME_MARCHING= DUAL_TIME_STEPPING-2ND_ORDER\nTIME_STEP= 1e-5\n";
  if (turbulentK)
    options << "KIND_TURB_MODEL= SST\nINIT_OPTION= TD_CONDITIONS\nREYNOLDS_NUMBER= 1E6\n"
               "VISCOSITY_MODEL= CONSTANT_VISCOSITY\nCONDUCTIVITY_MODEL= CONSTANT_CONDUCTIVITY\n";
  test.config_options = options.str();
  test.InitConfig();
  test.InitGeometry();
  test.InitSolver();
  auto* solver = test.solver[FLOW_SOL];
  auto* nodes = static_cast<CEulerVariable*>(solver->GetNodes());
  CIdealGas gas(1.4, 287);
  if (turbulentK) {
    gas.SetLaminarViscosityModel(test.config.get());
    gas.SetThermalConductivityModel(test.config.get());
  }
  const double rho = 1.2, p = 100000, a = sqrt(1.4 * p / rho);
  const su2double grid[3] = {moving ? 500. : 0., moving ? 240. : 0., 0.};
  // Tangential Mach 1.5 makes total speed supersonic even when normalMach is subsonic;
  // the grid boost also distinguishes lab-frame speed from the relative normal wave speed.
  su2double state[5] = {rho}, velocity[3] = {normalMach * a + grid[0], 1.5 * a + grid[1], 0};
  for (auto d = 0u; d < dim; ++d) state[d + 1] = rho * velocity[d];
  state[dim + 1] = p / .4 + rho * (.5 * GeometryToolbox::SquaredNorm(dim, velocity) + turbulentK);
  for (auto point = 0ul; point < test.geometry->GetnPoint(); ++point) {
    nodes->SetSolution(point, state);
    if (turbulentK) {
      test.solver[TURB_SOL]->GetNodes()->SetSolution(point, 0, turbulentK);
      REQUIRE(static_cast<CNSVariable*>(nodes)->SetPrimVar(point, 0, turbulentK, &gas));
    } else {
      REQUIRE(nodes->SetPrimVar(point, &gas));
    }
    if (moving) test.geometry->nodes->SetGridVel(point, grid);
    for (auto k = 0u; k < dim + 2; ++k) solver->LinSysRes(point, k) = 0;
    su2double zero[25] = {};
    solver->Jacobian.SetBlock(point, point, zero);
  }
  for (auto marker = 0u; marker < test.config->GetnMarker_All(); ++marker) {
    if (test.config->GetMarker_All_TagBound(marker) != "x_plus") continue;
    SU2_OMP_PARALLEL {
      OutletStateProbe probe(dim, test.config.get());
      solver->BC_Outlet(test.geometry.get(), test.solver, &probe, nullptr, test.config.get(), marker);
    }
    // Independent ideal-gas reference: gamma=1.4 gives q+5*a=constant and rho proportional to a^5.
    auto reference = [&](const su2double* u, su2double* out) {
      const su2double density = u[0], q = u[1] / density - grid[0];
      su2double v[3] = {};
      for (auto d = 0u; d < dim; ++d) v[d] = u[d + 1] / density;
      const su2double pressure = .4 * (u[dim + 1] - density * (.5 * GeometryToolbox::SquaredNorm(dim, v) + turbulentK));
      const su2double sound = sqrt(1.4 * pressure / density);
      if (q >= sound) {
        for (auto k = 0u; k < dim + 2; ++k) out[k] = u[k];
        return;
      }
      const su2double sonic = (.4 * q + 2 * sound) / 2.4;
      const su2double target = sound * pow(pressureRatio * p / pressure, .4 / 2.8);
      const su2double boundarySound = fmax(target, sonic);
      out[0] = density * pow(boundarySound / sound, 5);
      v[0] += 5 * (sound - boundarySound);
      for (auto d = 0u; d < dim; ++d) out[d + 1] = out[0] * v[d];
      out[dim + 1] =
          out[0] * (boundarySound * boundarySound / .56 + .5 * GeometryToolbox::SquaredNorm(dim, v) + turbulentK);
    };
    su2double expected[5] = {};
    if (normalMach > -1) reference(state, expected);
    // Backflow states come from the fluid-model maps, which pressure_outlet_backflow_tests.cpp checks.
    const su2double outward[3] = {1, 0, 0};
    const su2double reservoirData[5] = {300, tangentialMode ? 0. : totalPressure, 0, 0, 0};
    const bool steady = backflow != Backflow::Unsteady && backflow != Backflow::UnsteadyReservoir;
    const bool limited = backflow == Backflow::None && normalMach <= -1;
    // (Unsteady backflow without data is a hidden error case: the solver must stop with its own message.)
    if ((normalMach <= -1 && (reservoir || steady)) || (reservoir && expected[1] / expected[0] < grid[0]))
      REQUIRE(ExpectedOutletState(gas, dim, state, outward, grid[0], pressureRatio * p, expected, turbulentK, 1.4,
                                  reservoir ? reservoirData : nullptr, steady,
                                  tangentialMode) == CFluidModel::PressureOutletStatus::Valid);
    unsigned long domainVertices = 0;
    for (auto vertex = 0ul; vertex < test.geometry->nVertex[marker]; ++vertex) {
      const auto point = test.geometry->vertex[marker][vertex]->GetNode();
      if (!test.geometry->nodes->GetDomain(point)) continue;
      ++domainVertices;
      for (auto row = 0u; row < dim + 2; ++row) {
        CHECK(SU2_TYPE::GetValue(solver->LinSysRes(point, row)) ==
              Approx(SU2_TYPE::GetValue(expected[row])).epsilon(1e-8).margin(1e-8));
      }
      const su2double normal[3] = {1, 0, 0};
      for (auto column = 0u; column < dim + 2; ++column) {
        su2double increment[5];
        OutgoingOutletIncrement(gas, dim, state, normal, grid[0], turbulentK, column, increment);
        for (auto row = 0u; row < dim + 2; ++row) {
          CHECK(solver->Jacobian.GetBlock(point, point)[row * (dim + 2) + column] ==
                Approx(SU2_TYPE::GetValue(increment[row]) + (row == column ? 1. : 0.)).epsilon(1e-8).margin(1e-8));
        }
      }
    }
    // Steady runs count the vertices limited to sonic backflow; the monitor call also resets the count.
    std::ostringstream quiet;
    std::cout.rdbuf(quiet.rdbuf());
    auto* euler = static_cast<CEulerSolver*>(solver);
    euler->MonitorOutletBackflow(test.config.get());
    const auto counted = euler->GetOutletBackflowClamped();
    euler->MonitorOutletBackflow(test.config.get());
    const auto recounted = euler->GetOutletBackflowClamped();
    std::cout.rdbuf(test.orig_buf);
    CHECK(counted == (limited ? domainVertices : 0));
    CHECK(recounted == 0);
    if (limited) CHECK(quiet.str().find("limited to sonic inflow") != std::string::npos);
  }
}
}  // namespace

TEST_CASE("Ordinary pressure outlet uses relative normal characteristics and an outgoing-wave matrix",
          "[Pressure-outlet]") {
  for (auto dim : {2u, 3u})
    for (const auto moving : {false, true})
      for (const auto normalMach : {-.3, 0., .3, .99, 1.4})
        for (const auto ratio : {1., .95, 1.02, .001}) IdealOutlet(dim, normalMach, ratio, moving);
}

TEST_CASE("Ordinary pressure outlet rejects nonpositive pressure", "[.][Ordinary-invalid-pressure]") {
  IdealOutlet(2, .3, 0, false);
}
TEST_CASE("Ordinary pressure outlet requires data for unsteady backflow", "[.][Ordinary-backflow]") {
  IdealOutlet(2, -1.4, .95, false, false, 0, Backflow::Unsteady);
}

TEST_CASE("Ordinary pressure outlet limits steady backflow without data to sonic inflow", "[Pressure-outlet]") {
  for (auto dim : {2u, 3u})
    for (const auto moving : {false, true}) IdealOutlet(dim, -1.4, .95, moving);
}

TEST_CASE("Ordinary pressure outlet takes backflow from the MARKER_OUTLET_BACKFLOW reservoir",
          "[Pressure-outlet][Outlet-backflow]") {
  for (auto dim : {2u, 3u})
    for (const auto moving : {false, true})
      for (const auto normalMach : {-1.4, -.3, 0., .3})
        for (const auto ratio : {1., .95, 1.02})
          IdealOutlet(dim, normalMach, ratio, moving, false, 0, Backflow::Reservoir);
  IdealOutlet(3, -.3, 1.02, true, false, 1500, Backflow::Reservoir);
}

TEST_CASE("A zero backflow total pressure means the imposed outlet pressure", "[Pressure-outlet][Outlet-backflow]") {
  for (const auto normalMach : {-1.4, -.3, 0.})
    for (const auto ratio : {.95, 1.02})
      IdealOutlet(2, normalMach, ratio, false, false, 0, Backflow::ReservoirAtOutletPressure);
}

TEST_CASE("Ordinary outlet backflow can keep the interior tangential velocity", "[Pressure-outlet][Outlet-backflow]") {
  for (auto dim : {2u, 3u})
    for (const auto moving : {false, true})
      for (const auto normalMach : {-1.4, -.3, 0., .3})
        for (const auto ratio : {.95, 1.02})
          IdealOutlet(dim, normalMach, ratio, moving, false, 0, Backflow::ReservoirTangential);
  IdealOutlet(3, -.3, 1.02, true, false, 1500, Backflow::ReservoirTangential);
}

TEST_CASE("Unsteady outlet with a reservoir accepts backflow", "[Pressure-outlet][Outlet-backflow]") {
  IdealOutlet(2, -1.4, .95, false, false, 0, Backflow::UnsteadyReservoir);
  IdealOutlet(2, -.3, 1.02, false, false, 0, Backflow::UnsteadyReservoir);
}

TEST_CASE("Pressure outlet reconstruction is thread local", "[Pressure-outlet]") {
  IdealOutlet(3, .3, .95, true, true);
}

TEST_CASE("Ordinary outlet projection keeps SST energy in the characteristic basis", "[Pressure-outlet]") {
  for (auto dim : {2u, 3u})
    for (auto mach : {-.3, .3, 1.4}) IdealOutlet(dim, mach, .95, true, false, 1500);
}
