/*!
 * \file container_decorators_tests.cpp
 * \brief Unit tests for the container decorators.
 * \author K. R. Gulaj
 * \version 8.6.0 "Harrier"
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

#include <cstring>
#include <new>

#include "../../../Common/include/containers/container_decorators.hpp"

TEST_CASE("C3DContainerDecorator default state", "[Containers]") {
  /*--- The variable classes hold a C3DDoubleMatrix member that is default-initialized in memory
   * from operator new, which is not zeroed. The container is constructed here the same way, in
   * storage filled with a non-zero pattern, so that the test checks the default member initializers
   * and not zeroed memory ("new (buffer) C3DDoubleMatrix()" would value-initialize, which zeroes
   * every member first). ---*/
  alignas(C3DDoubleMatrix) unsigned char buffer[sizeof(C3DDoubleMatrix)];
  std::memset(buffer, 0xFF, sizeof(buffer));
  auto* mat = new (buffer) C3DDoubleMatrix;

  /*--- A container that was never resized is empty in every extent. ---*/
  CHECK(mat->size() == 0);
  CHECK(mat->length() == 0);
  CHECK(mat->rows() == 0);
  CHECK(mat->cols() == 0);

  mat->resize(3, 4, 2);
  CHECK(mat->size() == 24);
  CHECK(mat->length() == 3);
  CHECK(mat->rows() == 4);
  CHECK(mat->cols() == 2);

  mat->~C3DDoubleMatrix();
}
