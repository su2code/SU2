/*!
 * \file random_toolbox_tests.cpp
 * \brief Unit tests for the random toolbox.
 * \author A. Passariello
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
#include <cstdint>
#include "../../Common/include/toolboxes/random_toolbox.hpp"

TEST_CASE("HashToUniform", "[Toolboxes]") {
  /*--- Basic checks: no zero, no >1 ---*/
  uint64_t simple_values[] = {0ULL, 1ULL, 123ULL, UINT64_MAX};
  for (auto v : simple_values) {
    auto u = RandomToolbox::HashToUniform(v);
    CHECK(u > 0.0);
    CHECK(u <= 1.0);
  }

  /*--- Worst-case mantissa: top 53 bits set to 1 ---*/
  uint64_t x = (uint64_t(-1) >> 11) << 11;  // force top 53 bits = all ones
  x |= 0x7FF;                               // set lower bits to avoid masking issues

  double u = RandomToolbox::HashToUniform(x);

  CHECK(u > 0.0);
  CHECK(u <= 1.0);

  /*--- Check exact theoretical min and max ---*/
  constexpr double inv53 = 1.0 / 9007199254740992.0;

  CHECK(RandomToolbox::HashToUniform(0ULL) == Approx(inv53));
  CHECK(RandomToolbox::HashToUniform(UINT64_MAX) == Approx(1.0));
}
