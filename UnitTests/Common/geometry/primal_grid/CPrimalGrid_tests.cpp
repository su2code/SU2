/*!
 * \file CPrimalGrid_tests.cpp
 * \brief Unit tests for the primal grid classes
 * \author T. Albring
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
#include "../../../Common/include/geometry/primal_grid/CPrimalGrid.hpp"
#include "../../../Common/include/geometry/primal_grid/CHexahedron.hpp"

TEST_CASE("Center of gravity computation", "[Primal Grid]") {
  const int nDim = 3;

  su2double** coordinates = new su2double*[8];
  for (int i = 0; i < 8; i++) {
    coordinates[i] = new su2double[nDim];
  }
  coordinates[0][0] = 7.946516948817000e-01;
  coordinates[0][1] = 0.000000000000000e+00;
  coordinates[0][2] = -1.530741281331000e-03;
  coordinates[1][0] = 7.946538480939001e-01;
  coordinates[1][1] = 0.000000000000000e+00;
  coordinates[1][2] = -1.554823763559000e-03;
  coordinates[2][0] = 4.613089666440000e-01;
  coordinates[2][1] = 0.000000000000000e+00;
  coordinates[2][2] = -3.089699212314000e-03;
  coordinates[3][0] = 7.831152624057000e-01;
  coordinates[3][1] = 0.000000000000000e+00;
  coordinates[3][2] = -1.530741281331000e-03;
  coordinates[4][0] = 8.322953467911001e-01;
  coordinates[4][1] = 1.324931278110000e-01;
  coordinates[4][2] = -1.456502113362000e-03;
  coordinates[5][0] = 8.322984569865000e-01;
  coordinates[5][1] = 1.336230857244000e-01;
  coordinates[5][2] = -1.480923128397000e-03;
  coordinates[6][0] = 8.213630099601000e-01;
  coordinates[6][1] = 1.326360771765000e-01;
  coordinates[6][2] = -2.941176615903000e-03;
  coordinates[7][0] = 8.213597801418000e-01;
  coordinates[7][1] = 1.326371537826000e-01;
  coordinates[7][2] = -2.916814216089000e-03;

#define REQUIRE_CG(name, x, y, z)      \
  name.SetCoord_CG(nDim, coordinates); \
  REQUIRE(name.GetCG(0) == Approx(x)); \
  REQUIRE(name.GetCG(1) == Approx(y)); \
  REQUIRE(name.GetCG(2) == Approx(z));

  // It is sufficient to test the CG computation for Hexahedron.
  // Routine is the same for all elements (impl. in CPrimalGrid)
  CHexahedron hexa(0, 1, 2, 3, 4, 5, 6, 7);
  REQUIRE_CG(hexa, 0.7676307957, 0.0664236806, -0.0020626777);

  for (int i = 0; i < 8; i++) {
    delete[] coordinates[i];
  }
  delete[] coordinates;
}
