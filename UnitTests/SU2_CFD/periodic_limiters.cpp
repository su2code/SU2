/*!
 * \file periodic_limiters.cpp
 * \brief Tests for periodic limiter storage and reconstruction helpers.
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

TEST_CASE("Configured zero angles select translation and tiny angles select rotation", "[PeriodicLimiter]") {
  const su2double noAngles[3] = {0.0, -0.0, 0.0}, negativeAngle[3] = {0.0, -0.5, 0.0};
  CHECK_FALSE(GeometryToolbox::HasRotation(noAngles));
  CHECK(GeometryToolbox::HasRotation(negativeAngle));
  const su2double deg2rad = PI_NUMBER / 180.0;
  for (int iDim = 0; iDim < 3; ++iDim) {
    su2double angles[3] = {0.0, 0.0, 0.0};
    angles[iDim] = su2double(0.0) * deg2rad;
    CHECK_FALSE(GeometryToolbox::HasRotation(angles));
    angles[iDim] *= -1.0;
    CHECK_FALSE(GeometryToolbox::HasRotation(angles));
    angles[iDim] = su2double(1.0e-18) * deg2rad;
    CHECK(GeometryToolbox::HasRotation(angles));
    angles[iDim] *= -1.0;
    CHECK(GeometryToolbox::HasRotation(angles));
  }
}

TEST_CASE("Rotation of component bounds encloses every corner", "[PeriodicLimiter]") {
  auto checkCorners = [](auto& rotation, auto& lower, auto& upper) {
    constexpr int nDim = sizeof(lower) / sizeof(lower[0]);
    su2double expectedMin[nDim], expectedMax[nDim];
    for (int iDim = 0; iDim < nDim; ++iDim) {
      expectedMin[iDim] = std::numeric_limits<passivedouble>::max();
      expectedMax[iDim] = -expectedMin[iDim];
    }
    const su2double origin[nDim] = {0.0};
    for (int corner = 0; corner < (1 << nDim); ++corner) {
      su2double point[nDim], rotated[nDim];
      for (int iDim = 0; iDim < nDim; ++iDim) point[iDim] = (corner & (1 << iDim)) ? upper[iDim] : lower[iDim];
      GeometryToolbox::Rotate(rotation, origin, point, rotated);
      for (int iDim = 0; iDim < nDim; ++iDim) {
        expectedMin[iDim] = std::min(expectedMin[iDim], rotated[iDim]);
        expectedMax[iDim] = std::max(expectedMax[iDim], rotated[iDim]);
      }
    }
    GeometryToolbox::RotateBox(rotation, lower, upper);
    for (int iDim = 0; iDim < nDim; ++iDim) {
      CHECK(lower[iDim] == Approx(expectedMin[iDim]));
      CHECK(upper[iDim] == Approx(expectedMax[iDim]));
    }
  };

  SECTION("2D oblique rotation") {
    su2double rotation[2][2], lower[2] = {-2.0, 1.0}, upper[2] = {3.0, 4.0};
    GeometryToolbox::RotationMatrix(su2double(PI_NUMBER / 4.0), rotation);
    checkCorners(rotation, lower, upper);
  }
  SECTION("3D rotation about all axes") {
    su2double rotation[3][3], lower[3] = {-2.0, 1.0, -3.0}, upper[3] = {3.0, 4.0, 2.0};
    GeometryToolbox::RotationMatrix(su2double(0.3), su2double(-0.7), su2double(1.2), rotation);
    checkCorners(rotation, lower, upper);
  }
}

TEST_CASE("Shared reconstruction increment has the MUSCL scaling", "[PeriodicLimiter]") {
  const su2double halfEdge[3] = {1.0, -1.0, 0.5};
  const su2double gradient[3] = {2.0, -1.0, 3.0};
  for (const auto kappa : {-1.0, 0.0, 0.5, 1.0}) {
    const auto increment = LimiterHelpers<>::reconstructionIncrement(3, halfEdge, gradient, su2double(7.0),
                                                                     su2double(15.0), su2double(kappa));
    CHECK(increment == Approx(4.5 - 0.5 * kappa));
    CHECK(LimiterHelpers<>::reconstructionIncrement(2, halfEdge, gradient, su2double(7.0), su2double(15.0),
                                                    su2double(kappa)) == Approx(3.0 + kappa));
    CHECK(LimiterHelpers<>::reconstructionIncrement(2, halfEdge, gradient, su2double(7.0), su2double(13.0),
                                                    su2double(kappa)) == Approx(3.0));
    const auto linearIncrement = LimiterHelpers<>::reconstructionIncrement(3, halfEdge, gradient, su2double(7.0),
                                                                           su2double(16.0), su2double(kappa));
    CHECK(linearIncrement == Approx(4.5));
  }
}

TEST_CASE("Periodic projection storage uses unique receive points", "[PeriodicLimiter]") {
  std::stringstream options;
  options << "SOLVER= EULER\n"
             "MARKER_PERIODIC= (per1, per2, 0,0,0, 0,0,45, 0,0,0)\n";
  CConfig config(options, SU2_COMPONENT::SU2_CFD, false);
  config.SetnMarker_All(1);
  config.SetMarker_All_KindBC(0, PERIODIC_BOUNDARY);
  config.SetMarker_All_TagBound(0, "per1");

  struct ProjectionSolver : CSolver {
    CVariable* GetBaseClassPointerToNodes() override { return nullptr; }
    ProjectionSolver() {
      nPoint = 1000000;
      nVar = 4;
      nPrimVarGrad = 4;
      rotate_periodic = true;
    }
  } solver;

  CGeometry geometry;
  geometry.nPeriodicRecv = 1;
  geometry.nPoint_PeriodicRecv = new int[2]{0, 4};
  geometry.Local_Point_PeriodicRecv = new unsigned long[4]{12, 27, 12, 99999};

  auto* storage = solver.GetPeriodicProjections(geometry, config);
  REQUIRE(storage != nullptr);
  CHECK(storage->rows() == 3);
  CHECK(storage->cols() == 8);
  CHECK(solver.GetPeriodicProjection(13) == nullptr);
  REQUIRE(solver.GetPeriodicProjection(12) != nullptr);
  REQUIRE(solver.GetPeriodicProjection(27) != nullptr);
  REQUIRE(solver.GetPeriodicProjection(99999) != nullptr);
  solver.GetPeriodicProjection(12)[0] = -2.0;
  CHECK(solver.GetPeriodicProjection(27) != solver.GetPeriodicProjection(12));
  CHECK(solver.GetPeriodicProjections(geometry, config) == storage);
  CHECK(solver.GetPeriodicProjection(12)[0] == -2.0);

  CGeometry noReceives;
  ProjectionSolver solverWithoutReceives;
  auto* emptyStorage = solverWithoutReceives.GetPeriodicProjections(noReceives, config);
  REQUIRE(emptyStorage != nullptr);
  CHECK(emptyStorage->empty());
  CHECK(solverWithoutReceives.GetPeriodicProjection(12) == nullptr);

  solver.SetRotatePeriodic(false);
  CHECK(solver.GetPeriodicProjections(geometry, config) == nullptr);
}
