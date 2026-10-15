/*!
 * \file under_relaxation_ad_tests.cpp
 * \brief Tests density ownership and derivatives in compressible under-relaxation.
 * \author SU2 Contributors
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

TEST_CASE("Compressible under-relaxation owns its AD density values", "[UnderRelaxationAD]") {
  UnitQuadTestCase test;
  test.config_options =
      "SOLVER= EULER\nMESH_FORMAT= RECTANGLE\nMESH_BOX_SIZE= 5,5,5\n"
      "MARKER_EULER= (x_minus,x_plus,y_minus,y_plus)\n"
      "TIME_DISCRE_FLOW= EULER_IMPLICIT\nMAX_UPDATE_FLOW= 0.01\n";
  test.InitConfig();
  test.InitGeometry();
  test.InitSolver();
  auto* solver = static_cast<CEulerSolver*>(test.solver[FLOW_SOL]);
  auto* nodes = solver->GetNodes();
  for (auto p = 0ul; p < test.geometry->GetnPoint(); ++p)
    for (auto v = 0u; v < 4; ++v) solver->LinSysSol(p, v) = 0.0;

  for (const bool densityLimited : {true, false}) {
    su2double rho = 2.0;
#ifdef CODI_REVERSE_TYPE
    AD::StartRecording();
    AD::RegisterInput(rho);
#else
    SU2_TYPE::SetDerivative(rho, 1.0);
#endif
    nodes->SetSolution(0, 0, rho);
    nodes->SetSolution(0, 1, 3.0);
    nodes->SetSolution(0, 2, 0.0);
    nodes->SetSolution(0, 3, 10.0);
    solver->LinSysSol(0, 0) = densityLimited ? 0.1 : 0.0;
    solver->LinSysSol(0, 1) = densityLimited ? 0.2 : 0.0;
    solver->LinSysSol(0, 3) = densityLimited ? 0.4 : 4.0;
    solver->ComputeUnderRelaxationFactor(test.config.get());
    su2double factor = nodes->GetUnderRelaxation(0);
#ifdef CODI_REVERSE_TYPE
    AD::RegisterOutput(factor);
    AD::StopRecording();
    SU2_TYPE::SetDerivative(factor, 1.0);
    AD::ComputeAdjoint();
    const auto derivative = SU2_TYPE::GetDerivative(rho);
#else
    const auto derivative = SU2_TYPE::GetDerivative(factor);
#endif
    CHECK(SU2_TYPE::GetValue(factor) == Approx(densityLimited ? 0.2 : 0.019375));
    CHECK(derivative == Approx(densityLimited ? 0.1 : 0.0028125));
#ifdef CODI_REVERSE_TYPE
    AD::Reset();
#endif
  }
}
