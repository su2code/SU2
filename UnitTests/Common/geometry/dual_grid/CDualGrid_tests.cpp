/*!
 * \file CDualGrid_tests.cpp
 * \brief Unit tests for the dual grid classes
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
#include <iomanip>
#include "../../../Common/include/geometry/dual_grid/CEdge.hpp"
#include "../../../Common/include/geometry/dual_grid/CVertex.hpp"

TEST_CASE("Volume Computation", "[Dual Grid]") {
  su2double Coord_FaceiPoint[3];
  su2double Coord_FaceElem_CG[3];
  su2double Coord_Elem_CG[3];
  su2double Coord_Edge_CG[3];
  su2double scaling = 100;

  Coord_FaceiPoint[0] = scaling * 0.664995;
  Coord_FaceiPoint[1] = scaling * 1.1462;
  Coord_FaceiPoint[2] = scaling * 0.00926223;
  Coord_FaceElem_CG[0] = scaling * 0.655997;
  Coord_FaceElem_CG[1] = scaling * 1.13054;
  Coord_FaceElem_CG[2] = scaling * 0.00945181;
  Coord_Elem_CG[0] = scaling * 0.653846;
  Coord_Elem_CG[1] = scaling * 1.12927;
  Coord_Elem_CG[2] = scaling * 0.00835789;
  Coord_Edge_CG[0] = scaling * 0.664943;
  Coord_Edge_CG[1] = scaling * 1.14623;
  Coord_Edge_CG[2] = scaling * 0.00935524;

  SECTION("2D Edge") {
    su2double volume = CEdge::GetVolume(Coord_FaceiPoint, Coord_Edge_CG, Coord_Elem_CG);
    REQUIRE(volume == Approx(0.00607415));
  }

  SECTION("3D Edge") {
    su2double volume = CEdge::GetVolume(Coord_FaceiPoint, Coord_Edge_CG, Coord_FaceElem_CG, Coord_Elem_CG);
    REQUIRE(volume == Approx(0.000546832));
  }

  CVertex vertex2d(0, 2);
  SECTION("2D Vertex") {
    vertex2d.SetNodes_Coord(Coord_Edge_CG, Coord_Elem_CG);
    REQUIRE(vertex2d.GetNormal()[0] == Approx(-1.696));
    REQUIRE(vertex2d.GetNormal()[1] == Approx(1.1097));
  }

  CVertex vertex3d(0, 3);
  SECTION("3D Vertex") {
    vertex3d.SetNodes_Coord(Coord_Edge_CG, Coord_FaceElem_CG, Coord_Elem_CG);
    REQUIRE(vertex3d.GetNormal()[0] == Approx(-0.0864312));
    REQUIRE(vertex3d.GetNormal()[1] == Approx(0.0499696));
    REQUIRE(vertex3d.GetNormal()[2] == Approx(0.111938));
  }
}
