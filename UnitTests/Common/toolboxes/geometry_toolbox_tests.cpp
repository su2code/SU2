/*!
 * \file geometry_toolbox_tests.cpp
 * \brief Unit tests for the geometry toolbox.
 * \author P. Gaur
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
#include "../../Common/include/toolboxes/geometry_toolbox.hpp"

TEST_CASE("PointInEllipsoid", "[Toolboxes]") {
  const double center[3] = {1.0, 2.0, 3.0};
  const double axes[3] = {2.0, 1.0, 0.5};
  const double pi = 3.14159265358979323846;

  /*--- 2D, unrotated: the x semi-axis is 2 and the y semi-axis is 1. ---*/
  const double inX[3] = {2.9, 2.0, 0.0};
  const double outY[3] = {1.0, 3.1, 0.0};
  CHECK(GeometryToolbox::PointInEllipsoid(2, center, center, axes, 0.0));
  CHECK(GeometryToolbox::PointInEllipsoid(2, inX, center, axes, 0.0));
  CHECK_FALSE(GeometryToolbox::PointInEllipsoid(2, outY, center, axes, 0.0));

  /*--- 2D, rotated by 90 degrees about z: the long axis now lies along y. ---*/
  CHECK_FALSE(GeometryToolbox::PointInEllipsoid(2, inX, center, axes, 0.5 * pi));
  const double inYrot[3] = {1.0, 3.9, 0.0};
  CHECK(GeometryToolbox::PointInEllipsoid(2, inYrot, center, axes, 0.5 * pi));

  /*--- 3D: the z coordinate is measured from the center and scaled by the z semi-axis. ---*/
  const double inZ[3] = {1.0, 2.0, 3.4};
  const double outZ[3] = {1.0, 2.0, 3.6};
  const double atOrigin[3] = {1.0, 2.0, 0.0};
  CHECK(GeometryToolbox::PointInEllipsoid(3, inZ, center, axes, 0.0));
  CHECK_FALSE(GeometryToolbox::PointInEllipsoid(3, outZ, center, axes, 0.0));
  CHECK_FALSE(GeometryToolbox::PointInEllipsoid(3, atOrigin, center, axes, 0.0));
}
