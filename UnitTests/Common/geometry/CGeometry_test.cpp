/*!
 * \file CGeometry_tests.cpp
 * \brief Unit tests for CGeometry.
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
#include "../../UnitQuadTestCase.hpp"

std::unique_ptr<UnitQuadTestCase> TestCase;

TEST_CASE("Geometry constructor", "[Geometry]") {
  cout.rdbuf(nullptr);

  TestCase = std::unique_ptr<UnitQuadTestCase>(new UnitQuadTestCase());

  TestCase->InitConfig();

  auto aux_geometry = std::unique_ptr<CGeometry>(new CPhysicalGeometry(TestCase->config.get(), 0, 1));

  CHECK(aux_geometry->GetnPoint() == 125);
  CHECK(aux_geometry->GetnElem() == 64);
  CHECK(aux_geometry->GetnElemHexa() == 64);
  CHECK(aux_geometry->GetnEdge() == 0);
  CHECK(aux_geometry->GetnElem_Bound(0) == 16);
  CHECK(aux_geometry->GetnElem_Bound(5) == 16);

  TestCase->geometry = std::unique_ptr<CGeometry>(new CPhysicalGeometry(aux_geometry.get(), TestCase->config.get()));

  CHECK(TestCase->geometry->GetnPoint() == 125);
  CHECK(TestCase->geometry->GetnElem() == 64);
  CHECK(TestCase->geometry->GetnElemHexa() == 64);
  CHECK(TestCase->geometry->GetnEdge() == 0);
  CHECK(TestCase->geometry->GetnElem_Bound(0) == 16);
  CHECK(TestCase->geometry->GetnElem_Bound(5) == 16);

  cout.rdbuf(TestCase->orig_buf);
}

TEST_CASE("Set Send/Recv", "[Geometry]") {
  TestCase->geometry->SetSendReceive(TestCase->config.get());

  /*---- No check yet, since unit tests run in serial at the moment ---*/
}

TEST_CASE("Set Boundaries", "[Geometry]") {
  TestCase->geometry->SetBoundaries(TestCase->config.get());

  CHECK(TestCase->config->GetMarker_All_KindBC(0) == CUSTOM_BOUNDARY);
  CHECK(TestCase->config->GetMarker_All_KindBC(2) == HEAT_FLUX);
  CHECK(TestCase->config->GetSolid_Wall(2));
  CHECK(TestCase->config->GetSolid_Wall(3));
}

TEST_CASE("Set Point Connectivity", "[Geometry]") {
  TestCase->geometry->SetPoint_Connectivity();

  CHECK(TestCase->geometry->nodes->GetnElem(55) == 4);
  CHECK(TestCase->geometry->nodes->GetElem(30, 2) == 16);
  CHECK(TestCase->geometry->nodes->GetnNeighbor(3) == 4);
  CHECK(TestCase->geometry->nodes->GetPoint(99, 2) == 98);
}

TEST_CASE("Set elem connectivity", "[Geometry]") {
  TestCase->geometry->SetElement_Connectivity();

  CHECK(TestCase->geometry->elem[14]->GetnFaces() == 6);
  CHECK(TestCase->geometry->elem[14]->GetNeighbor_Elements(1) == 15);
}

TEST_CASE("Set bound volume", "[Geometry]") {
  TestCase->geometry->SetBoundVolume();

  CHECK(TestCase->geometry->bound[0][10]->GetDomainElement() == 40);
  CHECK(TestCase->geometry->bound[4][10]->GetDomainElement() == 10);
}

TEST_CASE("Set Edges", "[Geometry]") {
  TestCase->geometry->SetEdges();

  CHECK(TestCase->geometry->edges->GetnNodes() == 2);
  CHECK(TestCase->geometry->edges->GetNode(42, 0) == 15);
  CHECK(TestCase->geometry->edges->GetNode(87, 1) == 57);
}

TEST_CASE("Set vertex", "[Geometry]") {
  TestCase->geometry->SetVertex(TestCase->config.get());

  CHECK(TestCase->geometry->GetnVertex(0) == 25);
  CHECK(TestCase->geometry->vertex[0][20]->GetNode() == 100);
  CHECK(TestCase->geometry->nodes->GetVertex(100, 0) == 20);
  CHECK(TestCase->geometry->nodes->GetVertex(1, 0) == -1);
}

TEST_CASE("Set control volume", "[Geometry]") {
  TestCase->geometry->SetControlVolume(TestCase->config.get(), ALLOCATE);

  CHECK(TestCase->geometry->elem[42]->GetCG(0) == 0.625);
  CHECK(TestCase->geometry->elem[3]->GetCG(1) == 0.125);
  CHECK(TestCase->geometry->elem[25]->GetCG(2) == 0.375);

  CHECK(TestCase->geometry->nodes->GetVolume(42) == Approx(0.015625));

  CHECK(TestCase->geometry->edges->GetNormal(32)[0] == 0.03125);
  CHECK(TestCase->geometry->edges->GetNormal(5)[1] == 0.0);
  CHECK(TestCase->geometry->edges->GetNormal(10)[2] == 0.03125);

  CHECK(TestCase->config->GetDomainVolume() == Approx(1.0));
}

TEST_CASE("Set bound control volume", "[Geometry]") {
  TestCase->geometry->SetBoundControlVolume(TestCase->config.get(), ALLOCATE);

  CHECK(TestCase->geometry->bound[1][4]->GetCG(0) == 1.0);
  CHECK(TestCase->geometry->bound[3][2]->GetCG(1) == 1.0);
  CHECK(TestCase->geometry->bound[4][3]->GetCG(2) == 0.0);

  CHECK(TestCase->geometry->vertex[0][4]->GetNormal()[0] == -0.0625);
  CHECK(TestCase->geometry->vertex[3][2]->GetNormal()[1] == -0.0625);
  CHECK(TestCase->geometry->vertex[5][3]->GetNormal()[2] == 0.03125);
}
