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
#include "../../UnitQuadTestCase.hpp"
#include "../../../SU2_CFD/include/fluid/CIdealGas.hpp"

namespace {
// A boundary-state probe isolates the chain rule from approximate numerical-flux Jacobians.
class OutletStateProbe final : public CNumerics {
  su2double state[5] = {}, zero[5][5] = {}, identity[5][5] = {};
  const su2double *left[5], *right[5];

 public:
  OutletStateProbe(unsigned short dim, const CConfig* config) : CNumerics(dim, dim + 2, config) {
    for (auto k = 0u; k < nVar; ++k) {
      left[k] = zero[k];
      right[k] = identity[k];
      identity[k][k] = 1;
    }
  }
  ResidualType<> ComputeResidual(const CConfig*) override {
    const su2double rho = V_j[nDim + 2], p = V_j[nDim + 1];
    state[0] = rho;
    for (auto d = 0u; d < nDim; ++d) state[d + 1] = rho * V_j[d + 1];
    state[nDim + 1] = rho * V_j[nDim + 3] - p;
    return {state, left, right};
  }
};

void IdealOutlet(unsigned short dim, double normalMach, double pressureRatio, bool moving, bool threaded = false) {
  UnitQuadTestCase test;
  std::ostringstream options;
  options << "SOLVER= EULER\nMESH_FORMAT= " << (dim == 2 ? "RECTANGLE" : "BOX")
          << "\nMESH_BOX_SIZE= " << (threaded ? "5,17,17" : "5,5,5") << "\nREF_DIMENSIONALIZATION= DIMENSIONAL\n"
          << "MACH_NUMBER= 0.3\nFLUID_MODEL= IDEAL_GAS\nGAMMA_VALUE= 1.4\nGAS_CONSTANT= 287\n"
          << "FREESTREAM_OPTION= DENSITY_FS\nFREESTREAM_DENSITY= 1.2\nFREESTREAM_PRESSURE= 100000\n"
          << "CONV_NUM_METHOD_FLOW= HLLC\nTIME_DISCRE_FLOW= EULER_IMPLICIT\n"
          << "MARKER_EULER= (x_minus,y_minus,y_plus" << (dim == 3 ? ",z_minus,z_plus" : "")
          << ")\nMARKER_MONITORING= (x_minus)\nMARKER_OUTLET= (x_plus," << pressureRatio * 100000 << ")\n";
  if (moving) options << "GRID_MOVEMENT= STEADY_TRANSLATION\nTRANSLATION_RATE= 500,240,0\n";
  test.config_options = options.str();
  test.InitConfig();
  test.InitGeometry();
  test.InitSolver();
  auto* solver = test.solver[FLOW_SOL];
  auto* nodes = static_cast<CEulerVariable*>(solver->GetNodes());
  CIdealGas gas(1.4, 287);
  const double rho = 1.2, p = 100000, a = sqrt(1.4 * p / rho);
  const su2double grid[3] = {moving ? 500. : 0., moving ? 240. : 0., 0.};
  su2double state[5] = {rho}, velocity[3] = {normalMach * a + grid[0], 1.5 * a + grid[1], 0};
  for (auto d = 0u; d < dim; ++d) state[d + 1] = rho * velocity[d];
  state[dim + 1] = p / .4 + .5 * rho * GeometryToolbox::SquaredNorm(dim, velocity);
  for (auto point = 0ul; point < test.geometry->GetnPoint(); ++point) {
    nodes->SetSolution(point, state);
    REQUIRE(nodes->SetPrimVar(point, &gas));
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
    // Closed-form ideal-gas invariant; tangential speed does not determine characteristic direction.
    auto reference = [&](const su2double* u, su2double* out) {
      const su2double density = u[0], q = u[1] / density - grid[0];
      su2double v[3] = {};
      for (auto d = 0u; d < dim; ++d) v[d] = u[d + 1] / density;
      const su2double pressure = .4 * (u[dim + 1] - .5 * density * GeometryToolbox::SquaredNorm(dim, v));
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
      out[dim + 1] = out[0] * (boundarySound * boundarySound / .56 + .5 * GeometryToolbox::SquaredNorm(dim, v));
    };
    su2double expected[5];
    reference(state, expected);
    for (auto vertex = 0ul; vertex < test.geometry->nVertex[marker]; ++vertex) {
      const auto point = test.geometry->vertex[marker][vertex]->GetNode();
      if (!test.geometry->nodes->GetDomain(point)) continue;
      for (auto row = 0u; row < dim + 2; ++row) {
        CHECK(SU2_TYPE::GetValue(solver->LinSysRes(point, row)) ==
              Approx(SU2_TYPE::GetValue(expected[row])).epsilon(1e-8).margin(1e-8));
      }
      for (auto column = 0u; column < dim + 2; ++column) {
        su2double plus[5], minus[5], fp[5], fm[5];
        const double delta = 2e-5 * fmax(fabs(SU2_TYPE::GetValue(state[column])), column == 0 ? rho : rho * a);
        for (auto k = 0u; k < dim + 2; ++k) plus[k] = minus[k] = state[k];
        plus[column] += delta;
        minus[column] -= delta;
        reference(plus, fp);
        reference(minus, fm);
        for (auto row = 0u; row < dim + 2; ++row) {
          CHECK(solver->Jacobian.GetBlock(point, point)[row * (dim + 2) + column] ==
                Approx(SU2_TYPE::GetValue((fp[row] - fm[row]) / (2 * delta))).epsilon(2e-4).margin(2e-5));
        }
      }
    }
  }
}
}  // namespace

TEST_CASE("Ordinary pressure outlet uses relative normal characteristics and the state chain rule",
          "[Pressure-outlet]") {
  for (auto dim : {2u, 3u})
    for (const auto moving : {false, true})
      for (const auto normalMach : {-.3, 0., .3, .99, 1.4})
        for (const auto ratio : {1., .95, 1.02, .001}) IdealOutlet(dim, normalMach, ratio, moving);
}

TEST_CASE("Ordinary pressure outlet rejects nonpositive pressure", "[.][Ordinary-invalid-pressure]") {
  IdealOutlet(2, .3, 0, false);
}
TEST_CASE("Ordinary pressure outlet requires data for backflow", "[.][Ordinary-backflow]") {
  IdealOutlet(2, -1.4, .95, false);
}

TEST_CASE("Pressure outlet reconstruction is thread local", "[Pressure-outlet]") {
  IdealOutlet(3, .3, .95, true, true);
}
