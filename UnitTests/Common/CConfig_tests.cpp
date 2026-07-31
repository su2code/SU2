/*!
 * \file CConfig_tests.cpp
 * \brief Unit tests for configuration-level time semantics.
 * \author SU2 Contributors
 * \version 8.5.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2026, SU2 Contributors
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
#include <sstream>
#include "../../Common/include/CConfig.hpp"
#include "../../Common/include/grid_movement/CVolumetricMovement.hpp"
#include "../UnitQuadTestCase.hpp"

namespace {

std::unique_ptr<CConfig> MakeTimeConfig(const char* time_marching, const char* math_problem = "DIRECT") {
  std::stringstream options;
  options << "SOLVER= EULER\n"
          << "MATH_PROBLEM= " << math_problem << "\n"
          << "TIME_DOMAIN= YES\n"
          << "TIME_MARCHING= " << time_marching << "\n"
          << "TIME_STEP= 0.5\n"
          << "TIME_ITER= 10\n";
  return std::unique_ptr<CConfig>(new CConfig(options, SU2_COMPONENT::SU2_CFD, false));
}

}  // namespace

TEST_CASE("Physical time distinguishes a time-iteration label from its target state", "[CConfig][time]") {
  constexpr su2double delta_time = 0.25;

  SECTION("dual-time labels denote the state at the end of the step") {
    auto config = MakeTimeConfig("DUAL_TIME_STEPPING-1ST_ORDER");
    config->SetDelta_UnstTimeND(delta_time);

    CHECK(config->GetPhysicalTime(0) == Approx(delta_time));
    CHECK(config->GetPhysicalTime(7) == Approx(8.0 * delta_time));

    config->SetTimeIter(0);
    config->SetPhysicalTime(config->GetPhysicalTime(0));
    config->SetGlobalParam(config->GetKind_Solver(), RUNTIME_FLOW_SYS);

    CHECK(config->GetTimeIter() == 0);
    CHECK(config->GetCurrent_UnstTimeND() == Approx(config->GetPhysicalTime()));
    CHECK(config->GetPhysicalTimeAtEndOfTimeStep() == Approx(delta_time));
  }

  SECTION("second-order dual time has the same target-time convention") {
    auto config = MakeTimeConfig("DUAL_TIME_STEPPING-2ND_ORDER");
    config->SetDelta_UnstTimeND(delta_time);

    CHECK(config->GetPhysicalTime(0) == Approx(delta_time));
  }

  SECTION("stage-based time stepping starts at the initial time") {
    auto config = MakeTimeConfig("TIME_STEPPING");
    config->SetDelta_UnstTimeND(delta_time);

    CHECK(config->GetPhysicalTime(0) == Approx(0.0));
    CHECK(config->GetPhysicalTime(7) == Approx(7.0 * delta_time));

    config->SetPhysicalTime(config->GetPhysicalTime(0));
    CHECK(config->GetPhysicalTimeAtEndOfTimeStep() == Approx(delta_time));
  }

  SECTION("continuous adjoint time retains its reverse-time interpretation") {
    auto config = MakeTimeConfig("DUAL_TIME_STEPPING-1ST_ORDER", "CONTINUOUS_ADJOINT");
    config->SetDelta_UnstTimeND(delta_time);

    CHECK(config->GetPhysicalTime(0) == Approx(0.0));
  }
}

TEST_CASE("Rigid motion evaluates the first dual-time target state", "[CConfig][time][grid movement]") {
  UnitQuadTestCase test_case;
  test_case.AddOption("TIME_DOMAIN= YES");
  test_case.AddOption("TIME_MARCHING= DUAL_TIME_STEPPING-1ST_ORDER");
  test_case.AddOption("TIME_STEP= 0.25");
  test_case.AddOption("TIME_ITER= 1");
  test_case.AddOption("GRID_MOVEMENT= RIGID_MOTION");
  test_case.AddOption("PLUNGING_OMEGA= 0.0 1.0 0.0");
  test_case.AddOption("PLUNGING_AMPL= 0.0 0.1 0.0");
  test_case.InitConfig();

  constexpr su2double delta_time = 0.25;
  test_case.config->SetDelta_UnstTimeND(delta_time);
  test_case.config->SetLength_Ref(1.0);
  test_case.config->SetOmega_Ref(1.0);
  test_case.config->SetPhysicalTime(test_case.config->GetPhysicalTime(0));
  test_case.InitGeometry();

  const auto initial_y = test_case.geometry->nodes->GetCoord(0, 1);
  const auto amplitude = test_case.config->GetPlunging_Ampl(1) / test_case.config->GetLength_Ref();
  const auto omega = test_case.config->GetPlunging_Omega(1) / test_case.config->GetOmega_Ref();
  const auto expected_delta = -amplitude * sin(omega * delta_time);

  CVolumetricMovement movement(test_case.geometry.get());
  movement.Rigid_Plunging(test_case.geometry.get(), test_case.config.get(), ZONE_0, 0);

  CHECK(test_case.geometry->nodes->GetCoord(0, 1) - initial_y == Approx(expected_delta));
  CHECK(expected_delta != Approx(0.0));
}
