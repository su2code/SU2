/*!
 * \file ndflattener_tests.cpp
 * \brief Unit tests for NdFlattener template classes.
 * \author M. Aehle
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
#include "../../Common/include/toolboxes/ndflattener.hpp"

TEST_CASE("NdFlattener Test", "[NdFlattener]") {
  int rank_;
  SU2_MPI::Comm_rank(SU2_MPI::GetComm(), &rank_);
  int size_;
  SU2_MPI::Comm_size(SU2_MPI::GetComm(), &size_);
  const size_t rank = rank_, size = size_;

  /*-- Provide non-flat array --*/
  su2double** A = new su2double*[2];
  A[0] = new su2double[2];
  A[0][0] = 0.0;
  A[0][1] = 1.0;
  A[1] = new su2double[3 + rank];
  for (size_t i = 0; i < 3 + rank; i++) A[1][i] = 2.0 + rank + i;

  /*-- Accessor --*/
  auto f = std::make_pair((size_t)2, [rank, A](int i) {
    return std::make_pair((size_t)(i == 0 ? 2 : (3 + rank)), [rank, A, i](int j) { return A[i][j]; });
  });

  /*-- Read into flattening structure --*/
  NdFlattener<2> nd2;
  nd2.initialize_or_refresh(f);

  /*-- Modify A -> this should not alter nd2 at this point --*/
  A[0][0] = 0.5;

  /*-- Check structure --*/
  REQUIRE(nd2.size() == 2);
  REQUIRE(nd2[0][0] == 0.0);
  REQUIRE(nd2[0][1] == 1.0);
  REQUIRE(nd2[1].size() == 3 + rank);
  for (size_t i = 0; i < 3 + rank; i++) {
    REQUIRE(nd2[1][i] == 2.0 + rank + i);
  }

  /*-- Modify structure. --*/
  nd2[0][0] = 0.7;
  nd2[0].data()[1] = 1.7;

  /*-- gather flattening structures of all processes --*/
  NdFlattener<3> nd3(Nd_MPI_Environment(), nd2);

  /*-- Check gathered structure, non-const look-up. --*/
  REQUIRE(nd3.size() == size);
  for (size_t r = 0; r < size; r++) {
    REQUIRE(nd3[r].size() == 2);
    REQUIRE(nd3[r][0][0] == 0.7);
    REQUIRE(nd3[r][0][1] == 1.7);
    REQUIRE(nd3[r][1].size() == 3 + r);
    for (size_t i = 0; i < 3 + r; i++) {
      REQUIRE(nd3[r][1][i] == 2.0 + r + i);
    }
  }

  /*-- Check gathered structure, const look-up. --*/
  const NdFlattener<3>& nd3_const = nd3;
  REQUIRE(nd3_const.size() == size);
  for (size_t r = 0; r < size; r++) {
    REQUIRE(nd3_const[r].size() == 2);
    REQUIRE(nd3_const[r][0][0] == 0.7);
    REQUIRE(nd3_const[r][0][1] == 1.7);
    REQUIRE(nd3_const[r][1].size() == 3 + r);
    for (size_t i = 0; i < 3 + r; i++) {
      REQUIRE(nd3_const[r][1][i] == 2.0 + r + i);
    }
  }

  /*-- Reread modified A and check again. --*/
  nd2.initialize_or_refresh(f);
  nd3.refresh(Nd_MPI_Environment(), nd2);
  REQUIRE(nd3.size() == size);
  for (size_t r = 0; r < size; r++) {
    REQUIRE(nd3[r].size() == 2);
    REQUIRE(nd3[r][0][0] == 0.5);
    REQUIRE(nd3[r][0][1] == 1.0);
    REQUIRE(nd3[r][1].size() == 3 + r);
    for (size_t i = 0; i < 3 + r; i++) {
      REQUIRE(nd3[r][1][i] == 2.0 + r + i);
    }
  }

  /*-- Check 1D functionality, communicating unsigned long data. --*/
  NdFlattener<1, unsigned long> a1(
      std::make_pair((size_t)(3 + rank), [rank, A](int i) { return (unsigned long)(2 + rank + i); }));
  const NdFlattener<1, unsigned long>& a1_const = a1;
  REQUIRE(a1.size() == 3 + rank);
  for (size_t i = 0; i < 3 + rank; i++) {
    REQUIRE(a1[i] == 2 + rank + i);
    REQUIRE(a1_const[i] == 2 + rank + i);
  }
  a1[0] = 1;
  REQUIRE(a1_const.data()[0] == 1);
  a1.data()[0] = 2 + rank;
  REQUIRE(a1_const[0] == 2 + rank);
  const NdFlattener<2, unsigned long> a2_const(Nd_MPI_Environment(MPI_UNSIGNED_LONG), a1);
  REQUIRE(a2_const.size() == size);
  for (size_t r = 0; r < size; r++) {
    REQUIRE(a2_const[r].size() == 3 + r);
    for (size_t i = 0; i < 3 + r; i++) {
      REQUIRE(a2_const[r][i] == 2 + r + i);
    }
  }

  /*-- free stuff --*/
  delete[] A[0];
  delete[] A[1];
  delete[] A;
}
