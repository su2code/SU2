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
#include "../../UnitQuadTestCase.hpp"
#include "../../../SU2_CFD/include/solvers/CNEMOEulerSolver.hpp"

/*!
 * \brief Build a NEMO Euler solver on the unit box and check which conserved
 * variables the implicit update limiter (MAX_UPDATE_FLOW) reacts to. The
 * total-energy check used to be nested inside the species block and could
 * never run, so the factor stayed at one for arbitrarily large energy updates.
 */
TEST_CASE("NEMO under-relaxation reacts to the total-energy update", "[NEMO][Solver]") {
  UnitQuadTestCase testCase;
  testCase.config_options =
      "SOLVER= NEMO_EULER\n"
      "FLUID_MODEL= SU2_NONEQ\n"
      "GAS_MODEL= AIR-5\n"
      "GAS_COMPOSITION= (0.77, 0.23, 0.0, 0.0, 0.0)\n"
      "MESH_FORMAT= BOX\n"
      "MESH_BOX_SIZE= 3,3,3\n"
      "MESH_BOX_LENGTH= 1,1,1\n"
      "MESH_BOX_OFFSET= 0,0,0\n"
      "INIT_OPTION= TD_CONDITIONS\n"
      "MACH_NUMBER= 5.0\n"
      "FREESTREAM_PRESSURE= 101325.0\n"
      "FREESTREAM_TEMPERATURE= 288.15\n"
      "FREESTREAM_TEMPERATURE_VE= 288.15\n"
      "MARKER_FAR= (x_minus, x_plus, y_minus, y_plus, z_minus, z_plus)\n"
      "TIME_DISCRE_FLOW= EULER_IMPLICIT\n"
      "MAX_UPDATE_FLOW= 0.2\n";
  testCase.InitConfig();
  testCase.InitGeometry();

  CConfig* config = testCase.config.get();
  CGeometry* geometry = testCase.geometry.get();

  /*--- Construct the class under test directly instead of going through
   * UnitQuadTestCase::InitSolver (CSolverFactory), so the test depends on
   * CNEMOEulerSolver alone. The constructor reports to cout; silence it the
   * way the harness helpers do. ---*/
  std::cout.rdbuf(nullptr);
  CNEMOEulerSolver solver(geometry, config, MESH_0);
  std::cout.rdbuf(testCase.orig_buf);

  const unsigned short nSpecies = config->GetnSpecies();
  const unsigned short nDim = geometry->GetnDim();
  const unsigned short nVar = solver.GetnVar();
  REQUIRE(nVar == nSpecies + nDim + 2);

  const unsigned long nPointDomain = geometry->GetnPointDomain();
  REQUIRE(nPointDomain > 1);
  const unsigned long iPoint = nPointDomain / 2;
  const unsigned long jPoint = (iPoint + 1) % nPointDomain;

  const su2double allowableRatio = config->GetMaxUpdateFractionFlow();
  REQUIRE(allowableRatio == Approx(0.2));

  /*--- Exceed the allowable update fraction by this factor for one variable at
   * one point. If the limiter sees the variable, the factor drops to 1/excess. ---*/
  const su2double excess = 4.0;
  CVariable* nodes = solver.GetNodes();

  auto factorWithUpdateOn = [&](unsigned short iVar, su2double reference) {
    solver.LinSysSol.SetValZero();
    solver.LinSysSol(iPoint, iVar) = excess * allowableRatio * reference;
    solver.ComputeUnderRelaxationFactor(config);
    return nodes->GetUnderRelaxation(iPoint);
  };

  const unsigned short rhoE = nVar - 2;
  REQUIRE(rhoE >= nSpecies);

  /*--- Control: the species check compares the summed species update with the
   * mixture density, so the reference for the first species is the density. ---*/
  su2double density = 0.0;
  for (unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies)
    density += fabs(nodes->GetSolution(iPoint, iSpecies));
  CHECK(factorWithUpdateOn(0, density + EPS) == Approx(1.0 / excess).margin(1e-12));
  CHECK(nodes->GetUnderRelaxation(jPoint) == 1.0);

  /*--- The total-energy update must be limited by the same rule. ---*/
  const su2double energy = fabs(nodes->GetSolution(iPoint, rhoE));
  CHECK(factorWithUpdateOn(rhoE, energy + EPS) == Approx(1.0 / excess).margin(1e-12));
  CHECK(nodes->GetUnderRelaxation(jPoint) == 1.0);

  /*--- An energy update within the allowable fraction leaves the factor at one. ---*/
  solver.LinSysSol.SetValZero();
  solver.LinSysSol(iPoint, rhoE) = 0.5 * allowableRatio * energy;
  solver.ComputeUnderRelaxationFactor(config);
  CHECK(nodes->GetUnderRelaxation(iPoint) == 1.0);
}
