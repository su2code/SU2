/*!
 * \file CGeometry_tests.cpp
 * \brief Unit tests for CGeometry.
 * \author T. Albring
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
#include "../../../Common/include/toolboxes/geometry_toolbox.hpp"
#include "../../UnitQuadTestCase.hpp"
#include "../../../Common/include/geometry/CMultiGridGeometry.hpp"

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

  CHECK(TestCase->geometry->edges->GetNormal(31)[0] == 0.03125);
  CHECK(TestCase->geometry->edges->GetNormal(5)[1] == 0.0);
  CHECK(TestCase->geometry->edges->GetNormal(11)[2] == 0.03125);

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

TEST_CASE("Periodic slip-wall normal", "[Periodic]") {
  const bool multigrid = GENERATE(false, true);
  UnitQuadTestCase field;
  const auto start = field.config_options.find("MARKER_HEATFLUX=");
  const auto end = field.config_options.find("VISCOSITY_MODEL=");
  field.config_options.replace(start, end - start, "MARKER_EULER= (y_minus,y_plus)\nMARKER_CUSTOM= (z_minus,z_plus)\n");
  field.AddOption("MARKER_PERIODIC= (x_minus,x_plus, 0,0,0, 30,0,0, 0,0,0)");
  if (multigrid) field.AddOption("MGLEVEL= 1");
  field.InitConfig();
  field.InitGeometry(true);
  /*--- Bend the BOX into an annular sector about x: y is radius, x is azimuth. ---*/
  for (auto i = 0ul; i < field.geometry->GetnPoint(); ++i) {
    const auto* x = field.geometry->nodes->GetCoord(i);
    const su2double angle = x[0] * PI_NUMBER / 6, radius = 1 + x[1], axial = x[2];
    const su2double coordinate[] = {axial, -radius * sin(angle), radius * cos(angle)};
    field.geometry->nodes->SetCoord(i, coordinate);
  }
  field.geometry->SetControlVolume(field.config.get(), UPDATE);
  field.geometry->SetBoundControlVolume(field.config.get(), UPDATE);
  field.geometry->MatchPeriodic(field.config.get(), 1);
  field.geometry->PreprocessPeriodicComms(field.geometry.get(), field.config.get());
  std::unique_ptr<CMultiGridGeometry> coarse;
  CGeometry* geometry = field.geometry.get();
  if (multigrid) {
    coarse.reset(new CMultiGridGeometry(geometry, field.config.get(), 1));
    coarse->SetPoint_Connectivity(geometry);
    coarse->SetEdges();
    coarse->SetVertex(geometry, field.config.get());
    coarse->SetControlVolume(geometry, ALLOCATE);
    coarse->SetBoundControlVolume(geometry, field.config.get(), ALLOCATE);
    coarse->SetCoord(geometry);
    coarse->SetMGLevel(1);
    coarse->MatchPeriodic(field.config.get(), 1);
    coarse->PreprocessPeriodicComms(coarse.get(), field.config.get());
    geometry = coarse.get();
  }
  su2double error = 0;
  unsigned long checked = 0;
  for (auto marker = 0u; marker < geometry->GetnMarker(); ++marker) {
    if (field.config->GetMarker_All_TagBound(marker) != "y_plus") continue;
    for (auto v = 0ul; v < geometry->GetnVertex(marker); ++v) {
      const auto i = geometry->vertex[marker][v]->GetNode();
      const auto* x = geometry->nodes->GetCoord(i);
      if (fabs(x[1]) > 1e-12 || fabs(x[2] - 2) > 1e-12 || !geometry->nodes->GetDomain(i)) continue;
      const auto it = geometry->symmetryNormals[marker].find(v);
      const auto* normal =
          it == geometry->symmetryNormals[marker].end() ? geometry->vertex[marker][v]->GetNormal() : it->second.data();
      error = std::max(error, fabs(normal[1]) / GeometryToolbox::Norm(3, normal));
      ++checked;
    }
  }
  unsigned long total = 0;
  SU2_MPI::Allreduce(&checked, &total, 1, MPI_UNSIGNED_LONG, MPI_SUM, SU2_MPI::GetComm());
  REQUIRE(total > 0);
  CHECK(error < 1e-12);
}

TEST_CASE("Periodic volume refresh after mesh deformation", "[Periodic]") {
  UnitQuadTestCase field;
  const auto start = field.config_options.find("MARKER_HEATFLUX=");
  const auto end = field.config_options.find("VISCOSITY_MODEL=");
  field.config_options.replace(start, end - start, "MARKER_CUSTOM= (y_minus,y_plus,z_plus,z_minus)\n");
  field.AddOption("MARKER_PERIODIC= (x_minus,x_plus, 0,0,0, 0,0,0, 1,0,0)");
  field.AddOption("DEFORM_MESH= YES");
  field.InitConfig();
  field.InitGeometry(true);
  field.geometry->MatchPeriodic(field.config.get(), 1);
  field.geometry->PreprocessPeriodicComms(field.geometry.get(), field.config.get());
  field.InitSolver();
  std::vector<std::array<su2double, 3>> original(field.geometry->GetnPoint());
  for (auto i = 0ul; i < original.size(); ++i)
    for (auto d = 0u; d < 3; ++d) original[i][d] = field.geometry->nodes->GetCoord(i, d);
  CGeometry* meshes[] = {field.geometry.get()};
  for (const auto amplitude : {0.05, -0.02, 0.0}) {
    for (auto i = 0ul; i < original.size(); ++i) {
      const auto& x = original[i];
      field.geometry->nodes->SetCoord(i, 0, x[0] + amplitude * sin(2 * PI_NUMBER * x[0]) * sin(PI_NUMBER * x[2]));
    }
    SU2_OMP_PARALLEL { CGeometry::UpdateGeometry(meshes, field.config.get()); }
    END_SU2_OMP_PARALLEL
    su2double expected = 0, stored = 0;
    for (auto i = 0ul; i < field.geometry->GetnPointDomain(); ++i) {
      const auto* x = field.geometry->nodes->GetCoord(i);
      if (fabs(x[1] - 0.5) > 1e-12 || fabs(x[2] - 0.5) > 1e-12) continue;
      if (fabs(x[0] - 1) < 1e-12) expected = field.geometry->nodes->GetVolume(i);
      if (fabs(x[0]) < 1e-12) stored = field.geometry->nodes->GetPeriodicVolume(i);
    }
    su2double totals[2] = {stored, expected}, global[2] = {};
    SU2_MPI::Allreduce(totals, global, 2, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());
    REQUIRE(global[1] > 0);
    CHECK(global[0] == Approx(global[1]).margin(1e-12));
    auto* solver = field.solver[FLOW_SOL];
    auto* nodes = solver->GetNodes();
    for (auto i = 0ul; i < field.geometry->GetnPoint(); ++i)
      for (auto v = 0u; v < solver->GetnPrimVarGrad(); ++v)
        nodes->SetPrimitive(i, v, 2 + field.geometry->nodes->GetCoord(i, 1));
    SU2_OMP_PARALLEL { solver->SetPrimitive_Gradient_GG(field.geometry.get(), field.config.get()); }
    END_SU2_OMP_PARALLEL
    su2double error = 0;
    for (auto i = 0ul; i < field.geometry->GetnPointDomain(); ++i) {
      const auto* x = field.geometry->nodes->GetCoord(i);
      if (fabs(x[1] - 0.5) > 1e-12 || fabs(x[2] - 0.5) > 1e-12) continue;
      if (fabs(x[0]) > 1e-12 && fabs(x[0] - 1) > 1e-12) continue;
      error = std::max(error, fabs(nodes->GetGradient_Primitive()(i, 0, 1) - 1));
    }
    CHECK(error < 1e-12);
  }
}
