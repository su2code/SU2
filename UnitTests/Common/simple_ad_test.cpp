/*!
 * \file simple_ad_test.cpp
 * \brief Tests the algorithmic differentiation in SU2. Aside from testing
 *        basic functionality, this also serves as a regression test
 *        to make sure that AD works within unit testing.
 * \author C. Pederson
 * \version 8.0.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2023, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../Common/include/basic_types/datatype_structure.hpp"

su2double func(const su2double& x) { return x * x * x; }

/*---
 * This test case is based off of Tutorial 2 in the CoDiPack
 * documentation:
 *    https://www.scicomp.uni-kl.de/codi/dc/d7e/Tutorial2.html
 * SU2 wrapper functions have been substituted for the CoDiPack calls.
 * ---*/
TEST_CASE("Simple AD Test", "[AD tests]") {
  su2double x = 4.0;

  AD::StartRecording();
  AD::RegisterInput(x);

  su2double y = func(x);

  AD::RegisterOutput(y);
  AD::StopRecording();
  SU2_TYPE::SetDerivative(y, 1.0);
  AD::ComputeAdjoint();

  CHECK(SU2_TYPE::GetValue(y) == Approx(64));
  CHECK(SU2_TYPE::GetDerivative(x) == Approx(48));
}
