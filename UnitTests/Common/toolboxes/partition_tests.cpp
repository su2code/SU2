/*!
 * \file partition_tests.cpp
 * \brief Unit tests for the halo block message buffers and the partition imbalance measure.
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
#include "../../../Common/include/toolboxes/halo_block_buffer.hpp"
#include "../../../Common/include/toolboxes/partition_load.hpp"

TEST_CASE("Halo block buffer keeps neighboring large indices apart", "[Toolboxes]") {
  /*--- 2^24 and 2^24 + 1 are the same float. ---*/
  constexpr unsigned long big = 16777216;
  const float blockA[2] = {1.f, 2.f}, blockB[2] = {3.f, 4.f};

  CHaloBlockBuffer<float> buf;
  const auto pos = buf.BeginPoint();
  buf.AddBlock(pos, big, blockA, 2);
  buf.AddBlock(pos, big + 1, blockB, 2);
  const auto pos2 = buf.BeginPoint();
  buf.AddBlock(pos2, 7, blockB, 2);

  CHaloBlockReader<float> reader(buf.index.data(), buf.value.data(), 2);
  unsigned long col;
  REQUIRE(reader.NextPoint() == 2);
  const auto* a = reader.NextBlock(col);
  CHECK(col == big);
  CHECK(a[0] == 1.f);
  const auto* b = reader.NextBlock(col);
  CHECK(col == big + 1);
  CHECK(b[0] == 3.f);
  REQUIRE(reader.NextPoint() == 1);
  reader.NextBlock(col);
  CHECK(col == 7);
}

TEST_CASE("Partition imbalance flags moving a dominant line", "[Toolboxes]") {
  CHECK(PartitionImbalance({255, 255, 255, 255}) == Approx(1.0));
  CHECK(PartitionImbalance({930, 30, 30, 30}) == Approx(3.64706).epsilon(1e-4));
  CHECK(PartitionImbalance({0, 0}) == Approx(1.0));
}
