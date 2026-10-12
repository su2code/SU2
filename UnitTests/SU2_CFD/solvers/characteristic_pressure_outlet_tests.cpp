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
#include "../../UnitQuadTestCase.hpp"
#include "../../../SU2_CFD/include/numerics/flow/convection/hllc.hpp"
#include "../../../SU2_CFD/include/fluid/CPengRobinson.hpp"

namespace {
void CharacteristicPressureFixture(su2double pressure) {
  UnitQuadTestCase test;
  std::ostringstream options;
  options << "SOLVER= EULER\nMESH_FORMAT= RECTANGLE\nMESH_BOX_SIZE= 5,5,5\n"
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
  const su2double state[4] = {139.72158881490299, 11271.637336535776, -2338.5053086132339, 221781581.74179047};
  const su2double grid[2] = {20, 0};
  for (auto p = 0ul; p < test.geometry->GetnPoint(); ++p) {
    nodes->SetSolution(p, state);
    nodes->SetSolution_Old(p, state);
    REQUIRE(nodes->SetPrimVar(p, &gas));
    nodes->SetSecondaryVar(p, &gas);
    test.geometry->nodes->SetGridVel(p, grid);
    for (auto k = 0u; k < 4; ++k) solver->LinSysRes(p, k) = 0;
  }
  CUpwGeneralHLLC_Flow numerics(2, 4, test.config.get());
  for (auto marker = 0u; marker < test.config->GetnMarker_All(); ++marker)
    if (test.config->GetMarker_All_TagBound(marker) == "x_plus")
      solver->BC_Riemann(test.geometry.get(), test.solver, &numerics, nullptr, test.config.get(), marker);
  for (auto p = 0ul; p < test.geometry->GetnPoint(); ++p)
    for (auto k = 0u; k < 4; ++k) CHECK(std::isfinite(SU2_TYPE::GetValue(solver->LinSysRes(p, k))));
}
}  // namespace

TEST_CASE("Characteristic pressure boundary accepts an admissible state", "[Boundary-valid-state]") {
  CharacteristicPressureFixture(1400000);
}

// These fixtures abort through SU2_MPI::Error and are checked by a subprocess test.
TEST_CASE("Characteristic pressure boundary rejects negative density", "[.][Boundary-negative-density]") {
  CharacteristicPressureFixture(1000);
}

TEST_CASE("Characteristic pressure boundary rejects invalid thermodynamics", "[.][Boundary-invalid-thermo]") {
  CharacteristicPressureFixture(100000);
}
