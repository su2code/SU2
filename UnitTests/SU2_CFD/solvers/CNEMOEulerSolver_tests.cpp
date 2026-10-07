/*!
 * \file CNEMOEulerSolver_tests.cpp
 * \brief Unit tests for the NEMO Euler solver.
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
#include "../../../SU2_CFD/include/solvers/CNEMOEulerSolver.hpp"

TEST_CASE("NEMO under-relaxation limits species and total-energy updates", "[NEMO][Solver]") {
  const unsigned short nSpecies = 2;
  const unsigned short nVar = 6;
  /*--- Two species, two momentum components, total energy, and vibrational energy. ---*/
  const su2double solution[nVar] = {0.75, 0.25, 1.0, 1.0, 10.0, 2.0};
  const su2double allowableRatio = 0.2;

  const struct {
    const char* name;
    su2double update[nVar];
    su2double expected;
  } cases[] = {
      {"zero update", {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 1.0},
      {"small energy update", {0.0, 0.0, 0.0, 0.0, 1.0, 0.0}, 1.0},
      {"positive energy excess", {0.0, 0.0, 0.0, 0.0, 8.0, 0.0}, 0.25},
      {"negative energy excess", {0.0, 0.0, 0.0, 0.0, -8.0, 0.0}, 0.25},
      {"opposing species updates", {0.4, -0.4, 0.0, 0.0, 0.0, 0.0}, 0.25},
      {"species sets the tighter limit", {0.8, -0.8, 0.0, 0.0, 8.0, 0.0}, 0.125},
      {"energy sets the tighter limit", {0.2, -0.2, 0.0, 0.0, 8.0, 0.0}, 0.25},
      {"momentum and vibrational energy are not limited", {0.0, 0.0, 100.0, -100.0, 0.0, 100.0}, 1.0},
      {"tiny factor cancels the update", {0.0, 0.0, 0.0, 0.0, 1e12, 0.0}, 0.0},
  };

  for (const auto& test : cases) {
    CAPTURE(test.name);
    const su2double factor =
        CNEMOEulerSolver::ComputeUnderRelaxationFactor(nSpecies, nVar, solution, test.update, allowableRatio);
    CHECK(factor == Approx(test.expected).margin(1e-12));
  }
}
