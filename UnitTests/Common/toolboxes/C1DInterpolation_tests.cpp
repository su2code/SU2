/*!
 * \file C1DInterpolation_tests.cpp
 * \brief Unit tests for splines and what not.
 * \author P. Gomes
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
#include <sstream>
#include <iomanip>
#include <vector>
#include "../../../Common/include/toolboxes/C1DInterpolation.hpp"

su2double myPoly(su2double x) { return 1 + x * (-1 + x * (-1 + x)); }

TEST_CASE("C1DInterpolation", "[Toolboxes]") {
  std::vector<su2double> x{{-1.5, -1.0, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 2.0}}, y;

  for (auto v : x) y.push_back(myPoly(v));

  /*--- piece-wise linear ---*/
  CLinearInterpolation L(x, y);

  /*--- natural spline ---*/
  CCubicSpline S0(x, y);

  /*--- analytical end conditions ---*/
  CCubicSpline S1(x, y, CCubicSpline::SECOND, -11.0, CCubicSpline::FIRST, 7.0);

  CCubicSpline S2(x, y, CCubicSpline::FIRST, 8.75, CCubicSpline::SECOND, 10.0);

  /*--- at the knots ---*/
  for (auto v : x) {
    auto ref = myPoly(v);
    CHECK(S0(v) == ref);
    CHECK(S1(v) == ref);
    CHECK(S2(v) == ref);
  }

  /*--- away from the knots ---*/
  for (auto& v : x) v = std::min(v + 0.1, 2.0);

  for (auto v : x) {
    auto ref = myPoly(v);
    CHECK(S0(v) == Approx(ref).epsilon(0.1));
    CHECK(S1(v) == Approx(ref));
    CHECK(S2(v) == Approx(ref));
  }

  /*--- Checks that intervals are mapped correctly ---*/
  for (size_t i = 1; i < x.size() - 2; ++i) {
    CHECK(L(x[i]) == Approx(0.5 * (y[i] + y[i + 1])));
  }
}
