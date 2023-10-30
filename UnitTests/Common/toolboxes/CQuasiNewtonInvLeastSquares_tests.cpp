/*!
 * \file CQuasiNewtonInvLeastSquares_tests.cpp
 * \brief Unit tests for the CQuasiNewtonInvLeastSquares class.
 * Which should find the root of a n-d linear problem in n+1 iterations.
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
#include "../../../Common/include/toolboxes/CQuasiNewtonInvLeastSquares.hpp"

struct Problem {
  static constexpr int N = 4;
  const passivedouble coeffs[N][N] = {
      {0.5, -0.7, 0.2, 3.0}, {1.0, -0.2, -0.6, 0.0}, {0.1, 0.2, 3.14, -1.0}, {-1.0, -0.4, 0.0, 1.6}};
  /*--- Row sum, sol should be {1.0}. ---*/
  const passivedouble rhs[N] = {3.0, 0.2, 2.44, 0.2};
  passivedouble sol[N] = {0.0};

  template <class T>
  void iterate(const T& x) {
    for (int i = 0; i < N; ++i) {
      sol[i] = x(i, 0) + rhs[i];
      for (int j = 0; j < N; ++j) sol[i] -= coeffs[i][j] * x(j, 0);
    }
  }
};

template <class P, class Q>
void iterate(P& p, Q& q) {
  p.iterate(q);
  for (int i = 0; i < P::N; ++i) q.FPresult(i, 0) = p.sol[i];
  q.compute();
}

TEST_CASE("QN-ILS", "[Toolboxes]") {
  Problem p;
  CQuasiNewtonInvLeastSquares<passivedouble> qnils(Problem::N + 1, Problem::N, 1);

  /*--- Solve ---*/
  for (int i = 0; i <= Problem::N; ++i) iterate(p, qnils);

  /*--- Check we solved in N+1 iterations. ---*/
  for (int i = 0; i < Problem::N; ++i) CHECK(qnils(i, 0) == Approx(1.0));

  /*--- Check we don't break a converged problem. ---*/
  iterate(p, qnils);

  for (int i = 0; i < Problem::N; ++i) CHECK(qnils(i, 0) == Approx(1.0));
}
