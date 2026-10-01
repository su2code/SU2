/*!
 * \file limiters.cpp
 * \brief Unit tests for the limiter functions.
 * \author A. Rausa
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
#include "../../SU2_CFD/include/solvers/CSolver.hpp"
#include "../../SU2_CFD/include/limiters/computeLimiters.hpp"

/*--- The limiters are evaluated with the projection and the difference to the neighbor maximum (both >= 0)
 * and with the projection and the difference to the neighbor minimum (both <= 0). A limiter function must
 * then give the same value for (proj, delta) and (-proj, -delta). ---*/

TEST_CASE("Nishikawa limiters are symmetric", "[Limiters]") {
  using Helpers = LimiterHelpers<>;
  const su2double eps = 1e-3;

  for (const su2double proj : {0.1, 0.5, 1.0, 2.0}) {
    for (const su2double delta : {0.1, 0.5, 1.0, 1.9}) {
      CHECK(Helpers::r3Function(-proj, -delta, eps) == Approx(Helpers::r3Function(proj, delta, eps)));
      CHECK(Helpers::r4Function(-proj, -delta, eps) == Approx(Helpers::r4Function(proj, delta, eps)));
      CHECK(Helpers::r5Function(-proj, -delta, eps) == Approx(Helpers::r5Function(proj, delta, eps)));
    }
  }

  /*--- Nishikawa (AIAA 2022-1374): with a = |delta| = 1 and b = |proj| = 1, S4 = 6 and R4 = 7/8. ---*/
  CHECK(Helpers::r4Function(1.0, 1.0, 0.0) == Approx(0.875));
  CHECK(Helpers::r4Function(-1.0, -1.0, 0.0) == Approx(0.875));
}
