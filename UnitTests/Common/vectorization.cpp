/*!
 * \file vectorization.cpp
 * \brief Unit tests for the SIMD type and associated expression templates.
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
#include "../../Common/include/parallelization/vectorization.hpp"

using namespace std;

template <class T, class U>
struct arithmeticFun {
  static T f(T A, T B, T C, T D, U x, U y) { return pow((A + B - x * C) / y + pow(A * x - C / D, y), D); }
};

template <class T, class U>
struct logicFun {
  static T f(T A, T B, T C, T D, U x, U y) {
    // (B < A || B >= C) && ...
    return fmax(B < A, B >= fmin(C, -D)) * (abs(A) == abs(x)) * (abs(C) != abs(y));
  }
};

template <template <class, class> class Fun, class T, class U>
void computeAndCheck(T A, T B, T C, T D, U x, U y) {
  const auto result = Fun<T, U>::f(A, B, C, D, x, y);

  for (size_t k = 0; k < T::Size; ++k) {
    CHECK(result[k] == Fun<U, U>::f(A[k], B[k], C[k], D[k], x, y));
  }
}

TEST_CASE("SIMD INT", "[Vectorization]") {
  /*--- Integer types will use the expression templates. ---*/
  using Int = simd::Array<int>;

  Int A = 1, B = -3, C = 5, D = 2;
  int x = -1, y = 2;

  computeAndCheck<arithmeticFun>(A, B, C, D, x, y);
  computeAndCheck<logicFun>(A, B, C, D, x, y);

  Int t = sign(A) + sign(B + C);
  for (size_t k = 0; k < Int::Size; ++k) {
    CHECK(t[k] == 2);
  }
}

TEST_CASE("SIMD DOUBLE", "[Vectorization]") {
  /*--- Double use the explicitly vectorized template specializations. ---*/
  using Double = simd::Array<double>;

  Double A = 1, B = -3, C = 5, D = 2;
  double x = -1, y = 2;

  computeAndCheck<arithmeticFun>(A, B, C, D, x, y);
  computeAndCheck<logicFun>(A, B, C, D, x, y);

  Double t = sqrt(pow(B, 2) * C + D * y + A + x);
  for (size_t k = 0; k < Double::Size; ++k) {
    CHECK(t[k] == 7);
  }
}
