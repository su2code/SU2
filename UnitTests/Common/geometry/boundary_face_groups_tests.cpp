/*!
 * \file boundary_face_groups_tests.cpp
 * \brief Unit tests of the face groups of boundary points (CGeometry::ComputeBoundaryFaceGroups).
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
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>
#include "../../../Common/include/CConfig.hpp"
#include "../../../Common/include/geometry/CPhysicalGeometry.hpp"
#include "../../../Common/include/toolboxes/geometry_toolbox.hpp"
#include "../../../SU2_CFD/include/solvers/CNSSolver.hpp"
#include "../../../SU2_CFD/include/variables/CNSVariable.hpp"

namespace {

struct TestGeometry {
  std::unique_ptr<CConfig> config;
  std::unique_ptr<CGeometry> geometry;
};

/*--- Finest grid as SU2_CFD preprocesses it (serial). ---*/
TestGeometry BuildGeometry(const std::string& options) {
  TestGeometry out;
  std::stringstream ss(options);
  out.config = std::unique_ptr<CConfig>(new CConfig(ss, SU2_COMPONENT::SU2_CFD, false));
  auto* config = out.config.get();
  {
    auto aux = std::unique_ptr<CGeometry>(new CPhysicalGeometry(config, 0, 1));
    aux->SetColorGrid_Parallel(config);
    out.geometry = std::unique_ptr<CGeometry>(new CPhysicalGeometry(aux.get(), config));
  }
  auto* geometry = out.geometry.get();
  geometry->SetSendReceive(config);
  geometry->SetBoundaries(config);
  geometry->SetPoint_Connectivity();
  geometry->SetElement_Connectivity();
  geometry->SetBoundVolume();
  geometry->Check_IntElem_Orientation(config);
  geometry->Check_BoundElem_Orientation(config);
  geometry->SetEdges();
  geometry->SetVertex(config);
  geometry->SetControlVolume(config, ALLOCATE);
  geometry->SetBoundControlVolume(config, ALLOCATE);
  geometry->FindNormal_Neighbor(config);
  return out;
}

using Key = std::array<long, 3>;

/*--- Integer key of a point (coordinates in units of h). ---*/
Key PointKey(const CGeometry& geometry, unsigned long iPoint, double h) {
  Key key = {0, 0, 0};
  for (unsigned short iDim = 0; iDim < geometry.GetnDim(); ++iDim)
    key[iDim] = std::lround(SU2_TYPE::GetValue(geometry.nodes->GetCoord(iPoint, iDim)) / h);
  return key;
}

/*--- Data of a boundary point over all markers of a family: per group (keyed by its unit normal) the total area (sum
 * of the shares of the point's vertices), the normal neighbour and its distance; the slip projector Q^T Q. ---*/
struct Group {
  double area = 0.0;
  Key neighbor = {0, 0, 0};
  double distance = 0.0;
};
struct PointGroups {
  std::map<Key, Group> groups;
  unsigned short slipRank = 0;
  std::array<double, 9> projector = {};
};

std::map<Key, PointGroups> CollectGroups(const CGeometry& geometry, const CConfig& config,
                                         const std::vector<std::string>& markers, double h) {
  std::map<Key, PointGroups> out;
  const auto nDim = geometry.GetnDim();
  for (unsigned short iMarker = 0; iMarker < geometry.GetnMarker(); ++iMarker) {
    const auto tag = config.GetMarker_All_TagBound(iMarker);
    if (std::find(markers.begin(), markers.end(), tag) == markers.end()) continue;
    for (unsigned long iVertex = 0; iVertex < geometry.GetnVertex(iMarker); ++iVertex) {
      const auto iPoint = geometry.vertex[iMarker][iVertex]->GetNode();
      auto& data = out[PointKey(geometry, iPoint, h)];
      for (size_t iGroup = 0; iGroup < geometry.GetnVertexGroups(iMarker, iVertex); ++iGroup) {
        su2double unitNormal[3] = {0.0, 0.0, 0.0}, area = 0.0, distance = 0.0;
        unsigned long neighbor = 0;
        geometry.GetVertexGroup(iMarker, iVertex, iGroup, unitNormal, area, neighbor, distance);
        Key normalKey = {0, 0, 0};
        for (unsigned short iDim = 0; iDim < nDim; ++iDim)
          normalKey[iDim] = std::lround(1e9 * SU2_TYPE::GetValue(unitNormal[iDim]));
        auto& group = data.groups[normalKey];
        group.area += SU2_TYPE::GetValue(area);
        group.neighbor = PointKey(geometry, neighbor, h);
        group.distance = SU2_TYPE::GetValue(distance);
      }
      su2double Q[3][3] = {};
      const auto rank = geometry.GetSlipBasis(iMarker, iVertex, Q);
      if (rank == 0) continue;
      data.slipRank = rank;
      for (unsigned short i = 0; i < 3; ++i)
        for (unsigned short j = 0; j < 3; ++j) {
          double p = 0.0;
          for (unsigned short r = 0; r < rank; ++r) p += SU2_TYPE::GetValue(Q[r][i] * Q[r][j]);
          data.projector[3 * i + j] = p;
        }
    }
  }
  return out;
}

/*--- Same groups, areas (1e-14 relative), neighbours and slip projectors at every point. A point without stored data
 * (rank 0) has one group, and then a basis of rank 1 is the projector of that group. ---*/
void CheckSameGroups(const std::map<Key, PointGroups>& a, const std::map<Key, PointGroups>& b) {
  REQUIRE(a.size() == b.size());
  for (const auto& item : a) {
    REQUIRE(b.count(item.first) == 1);
    const auto& x = item.second;
    const auto& y = b.at(item.first);
    REQUIRE(x.groups.size() == y.groups.size());
    for (const auto& group : x.groups) {
      REQUIRE(y.groups.count(group.first) == 1);
      CHECK(group.second.area == Approx(y.groups.at(group.first).area).epsilon(1e-14));
      CHECK(group.second.neighbor == y.groups.at(group.first).neighbor);
      CHECK(group.second.distance == Approx(y.groups.at(group.first).distance).epsilon(1e-14));
    }
    if (x.slipRank == y.slipRank) {
      for (int c = 0; c < 9; ++c) CHECK(x.projector[c] == Approx(y.projector[c]).margin(1e-14));
    } else {
      CHECK(x.groups.size() == 1);
      CHECK(std::max(x.slipRank, y.slipRank) == 1);
    }
  }
}

/*--- L-shaped 2D domain [0,2]^2 minus (1,2]x(1,2], triangles on a 0.5 grid. Wall edges: y = 0, x = 2, and the step
 * x = 1 (y >= 1), y = 1 (x >= 1); the others are far field. split = true puts the step on a second wall marker. The
 * wall corners (2,0), (2,1) and (1,1) (convex) and the far-field corner (0,2) have two face groups. ---*/
void WriteLMesh(const std::string& name, bool split) {
  const double h = 0.5;
  std::map<std::pair<int, int>, int> id;
  std::vector<std::array<double, 2>> pts;
  auto node = [&](int i, int j) {
    auto it = id.find({i, j});
    if (it != id.end()) return it->second;
    const int n = static_cast<int>(pts.size());
    id[{i, j}] = n;
    pts.push_back({i * h, j * h});
    return n;
  };
  std::vector<std::array<int, 3>> tri;
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j) {
      if (i >= 2 && j >= 2) continue;
      const int a = node(i, j), b = node(i + 1, j), c = node(i + 1, j + 1), d = node(i, j + 1);
      tri.push_back({a, b, c});
      tri.push_back({a, c, d});
    }
  std::map<std::pair<int, int>, int> count;
  for (const auto& t : tri)
    for (int e = 0; e < 3; ++e) {
      const int u = t[e], v = t[(e + 1) % 3];
      count[{std::min(u, v), std::max(u, v)}]++;
    }
  std::vector<std::array<int, 2>> wallA, wallB, far;
  for (const auto& t : tri)
    for (int e = 0; e < 3; ++e) {
      const int u = t[e], v = t[(e + 1) % 3];
      if (count[{std::min(u, v), std::max(u, v)}] != 1) continue;
      const double xm = 0.5 * (pts[u][0] + pts[v][0]), ym = 0.5 * (pts[u][1] + pts[v][1]);
      const bool outer = (std::fabs(ym) < 1e-12) || (std::fabs(xm - 2.0) < 1e-12);
      const bool step = (std::fabs(xm - 1.0) < 1e-12 && ym > 1.0) || (std::fabs(ym - 1.0) < 1e-12 && xm > 1.0);
      if (outer)
        wallA.push_back({u, v});
      else if (step)
        (split ? wallB : wallA).push_back({u, v});
      else
        far.push_back({u, v});
    }
  std::ofstream f(name);
  f << "NDIME= 2\nNELEM= " << tri.size() << "\n";
  for (size_t e = 0; e < tri.size(); ++e)
    f << "5 " << tri[e][0] << " " << tri[e][1] << " " << tri[e][2] << " " << e << "\n";
  f << "NPOIN= " << pts.size() << "\n";
  f.precision(17);
  for (size_t p = 0; p < pts.size(); ++p) f << pts[p][0] << " " << pts[p][1] << " " << p << "\n";
  f << "NMARK= " << (split ? 3 : 2) << "\n";
  auto marker = [&](const std::string& tag, const std::vector<std::array<int, 2>>& edges) {
    f << "MARKER_TAG= " << tag << "\nMARKER_ELEMS= " << edges.size() << "\n";
    for (const auto& e : edges) f << "3 " << e[0] << " " << e[1] << "\n";
  };
  marker("wall", wallA);
  if (split) marker("step", wallB);
  marker("far", far);
}

/*--- Unit cube [0,1]^3 with 2 x 2 x 2 hexahedra. The faces x = 0, y = 0 and z = 0 are the wall, split into 1 marker,
 * 2 markers (x = 0; y = 0 and z = 0), 3 markers (one per plane) or one marker per face (layouts 1 to 4); x = 1, y = 1,
 * z = 1 are far field. Returns the names of the wall markers. ---*/
std::vector<std::string> WriteCubeMesh(const std::string& name, int layout) {
  const int n = 2;
  auto id = [&](int i, int j, int k) { return i + (n + 1) * (j + (n + 1) * k); };
  std::ofstream f(name);
  f << "NDIME= 3\nNELEM= " << n * n * n << "\n";
  int e = 0;
  for (int k = 0; k < n; ++k)
    for (int j = 0; j < n; ++j)
      for (int i = 0; i < n; ++i)
        f << "12 " << id(i, j, k) << " " << id(i + 1, j, k) << " " << id(i + 1, j + 1, k) << " " << id(i, j + 1, k)
          << " " << id(i, j, k + 1) << " " << id(i + 1, j, k + 1) << " " << id(i + 1, j + 1, k + 1) << " "
          << id(i, j + 1, k + 1) << " " << e++ << "\n";
  f << "NPOIN= " << (n + 1) * (n + 1) * (n + 1) << "\n";
  for (int k = 0; k <= n; ++k)
    for (int j = 0; j <= n; ++j)
      for (int i = 0; i <= n; ++i)
        f << double(i) / n << " " << double(j) / n << " " << double(k) / n << " " << id(i, j, k) << "\n";

  /*--- Quads of the six planes: plane p = 0, 1, 2 (x, y, z) at index s = 0 or n. ---*/
  auto quads = [&](int p, int s) {
    std::vector<std::array<int, 4>> out;
    for (int b = 0; b < n; ++b)
      for (int a = 0; a < n; ++a) {
        auto node = [&](int u, int v) {
          if (p == 0) return id(s, u, v);
          if (p == 1) return id(u, s, v);
          return id(u, v, s);
        };
        out.push_back({node(a, b), node(a + 1, b), node(a + 1, b + 1), node(a, b + 1)});
      }
    return out;
  };
  std::vector<std::string> walls;
  std::map<std::string, std::vector<std::array<int, 4>>> elems;
  for (int p = 0; p < 3; ++p) {
    for (const auto& q : quads(p, 0)) {
      std::string tag = "wall";
      if (layout == 2) tag = (p == 0) ? "wallx" : "wallyz";
      if (layout == 3) tag += "xyz"[p];
      if (layout == 4) tag += std::to_string(walls.size());
      if (elems.count(tag) == 0) walls.push_back(tag);
      elems[tag].push_back(q);
    }
  }
  std::vector<std::pair<std::string, std::vector<std::array<int, 4>>>> markers;
  for (const auto& tag : walls) markers.push_back({tag, elems[tag]});
  markers.push_back({"far", {}});
  for (int p = 0; p < 3; ++p)
    for (const auto& q : quads(p, n)) markers.back().second.push_back(q);

  f << "NMARK= " << markers.size() << "\n";
  for (const auto& m : markers) {
    f << "MARKER_TAG= " << m.first << "\nMARKER_ELEMS= " << m.second.size() << "\n";
    for (const auto& q : m.second) f << "9 " << q[0] << " " << q[1] << " " << q[2] << " " << q[3] << "\n";
  }
  return walls;
}

std::string MarkerList(const std::vector<std::string>& markers, bool heatFlux) {
  std::string list;
  for (const auto& m : markers) list += (list.empty() ? "" : ", ") + m + (heatFlux ? ", 0.0" : "");
  return list;
}

}  // namespace

TEST_CASE("Boundary face groups: 2D corners of one marker and of two markers", "[BoundaryFaceGroups]") {
  if (SU2_MPI::GetSize() > 1) return;  // serial geometry test
  std::map<Key, PointGroups> wall[2], far[2];
  for (int split = 0; split < 2; ++split) {
    const std::string mesh = split ? "face_groups_L_split.su2" : "face_groups_L.su2";
    WriteLMesh(mesh, split);
    const std::string markers = split ? "wall, step" : "wall";
    const auto test = BuildGeometry("SOLVER= EULER\nMESH_FILENAME= " + mesh + "\nMESH_FORMAT= SU2\nMACH_NUMBER= 0.3\n" +
                                    "MARKER_EULER= ( " + markers + " )\nMARKER_FAR= ( far )\n");
    wall[split] = CollectGroups(*test.geometry, *test.config, {"wall", "step"}, 0.5);
    far[split] = CollectGroups(*test.geometry, *test.config, {"far"}, 0.5);
  }
  /*--- Slip bases exactly at the wall corners (2,0), (2,1) (convex) and (1,1), in units of 0.5: rank 2 (zero velocity),
   * two groups. The same groups, areas and bases with one and with two markers. ---*/
  for (int split = 0; split < 2; ++split) {
    for (const auto& item : wall[split]) {
      const bool corner = (item.first == Key{4, 0, 0}) || (item.first == Key{4, 2, 0}) || (item.first == Key{2, 2, 0});
      CHECK(item.second.slipRank == (corner ? 2 : 0));
      CHECK(item.second.groups.size() == (corner ? 2u : 1u));
    }
    CHECK(far[split].at({0, 4, 0}).groups.size() == 2);
  }
  CheckSameGroups(wall[0], wall[1]);
  CheckSameGroups(far[0], far[1]);
}

TEST_CASE("Boundary face groups: 3D cube corner with the wall on 1, 2, 3 and per-face markers",
          "[BoundaryFaceGroups]") {
  if (SU2_MPI::GetSize() > 1) return;  // serial geometry test
  for (const bool heatFlux : {false, true}) {
    std::map<Key, PointGroups> reference;
    for (int layout = 1; layout <= 4; ++layout) {
      const std::string mesh = "face_groups_cube_" + std::to_string(layout) + ".su2";
      const auto walls = WriteCubeMesh(mesh, layout);
      const std::string wallOption = heatFlux ? "SOLVER= NAVIER_STOKES\nREYNOLDS_NUMBER= 1000000\nMARKER_HEATFLUX= ( " +
                                                    MarkerList(walls, true) + " )\n"
                                              : "SOLVER= EULER\nMARKER_EULER= ( " + MarkerList(walls, false) + " )\n";
      const auto test = BuildGeometry(wallOption + "MESH_FILENAME= " + mesh +
                                      "\nMESH_FORMAT= SU2\nMACH_NUMBER= 0.3\nMARKER_FAR= ( far )\n");
      const auto groups = CollectGroups(*test.geometry, *test.config, walls, 0.5);
      if (layout > 1) {
        CheckSameGroups(reference, groups);
        continue;
      }
      reference = groups;

      /*--- One marker: 19 wall points; the corner (0,0,0) has three groups (rank 3), the other points of the three
       * edges two (rank 2), the rest one. Group areas: a quarter of a 0.5 x 0.5 face per incident face. Normals
       * point into the domain. Normal neighbours: off the walls along the group normal where such a neighbour exists
       * (face points), otherwise the largest-cosine neighbour; at the normal distance 0.5, so that the isothermal
       * Jacobian factor sum_q a_q / d_q of a point is twice its area. ---*/
      REQUIRE(groups.size() == 19);
      for (const auto& item : groups) {
        const auto& p = item.first;
        const int nZero = (p[0] == 0) + (p[1] == 0) + (p[2] == 0);
        CHECK(item.second.groups.size() == static_cast<size_t>(nZero));
        if (!heatFlux) CHECK(item.second.slipRank == ((nZero > 1) ? nZero : 0));
        for (const auto& group : item.second.groups) {
          int plane = 0;
          while (group.first[plane] == 0) ++plane;
          CHECK(group.first[plane] == 1000000000);
          int nFaces = 1;
          for (int d = 0; d < 3; ++d)
            if ((d != plane) && (p[d] == 1)) nFaces *= 2;
          CHECK(group.second.area == Approx(0.0625 * nFaces).epsilon(1e-14));
          Key expected = p;
          expected[plane] += 1;
          CHECK(group.second.neighbor == expected);
          CHECK(group.second.distance == Approx(0.5).epsilon(1e-14));
        }
      }
    }
  }
}

TEST_CASE("Boundary face groups: missing eligible neighbour uses vertex stencil", "[BoundaryFaceGroups]") {
  if (SU2_MPI::GetSize() > 1) return;
  const auto walls = WriteCubeMesh("face_groups_fallback.su2", 1);
  const auto test = BuildGeometry(
      "SOLVER= NAVIER_STOKES\nREYNOLDS_NUMBER= 1000000\nMARKER_HEATFLUX= ( wall, 0.0 )\n"
      "MESH_FILENAME= face_groups_fallback.su2\nMESH_FORMAT= SU2\n"
      "MACH_NUMBER= 0.3\nMARKER_FAR= ( far )\n");
  const auto& geometry = *test.geometry;
  bool found = false;
  for (unsigned short marker = 0; marker < geometry.GetnMarker(); ++marker) {
    if (test.config->GetMarker_All_TagBound(marker) != "wall") continue;
    for (unsigned long vertex = 0; vertex < geometry.GetnVertex(marker); ++vertex) {
      const auto point = geometry.vertex[marker][vertex]->GetNode();
      if (PointKey(geometry, point, 0.5) != Key{0, 0, 0}) continue;
      found = true;
      REQUIRE(geometry.HasVertexGroups(marker, vertex));
      CHECK_FALSE(geometry.HasEligibleGroupNeighbors(marker, vertex));
      CHECK(geometry.GetVertexNormalNeighbor(marker, vertex) == geometry.vertex[marker][vertex]->GetNormal_Neighbor());
    }
  }
  CHECK(found);
}

TEST_CASE("Boundary face groups: incompressible integrated flux uses vertex areas", "[BoundaryFaceGroups]") {
  if (SU2_MPI::GetSize() > 1) return;
  WriteCubeMesh("face_groups_incompressible.su2", 1);
  const auto test = BuildGeometry(
      "SOLVER= INC_NAVIER_STOKES\nVISCOSITY_MODEL= CONSTANT_VISCOSITY\nMU_CONSTANT= 1.8e-5\nMARKER_HEATFLUX= ( wall, "
      "0.0 )\n"
      "MESH_FILENAME= face_groups_incompressible.su2\nMESH_FORMAT= SU2\n"
      "MARKER_FAR= ( far )\n");
  auto& geometry = *test.geometry;
  REQUIRE_FALSE(CGeometry::UseBoundaryFaceGroups(test.config.get()));
  geometry.ComputeSurfaceAreaCfgFile(test.config.get());
  for (unsigned short marker = 0; marker < geometry.GetnMarker(); ++marker) {
    if (test.config->GetMarker_All_TagBound(marker) != "wall") continue;
    su2double area = 0.0;
    for (unsigned long vertex = 0; vertex < geometry.GetnVertex(marker); ++vertex) {
      const auto point = geometry.vertex[marker][vertex]->GetNode();
      if (!geometry.nodes->GetDomain(point)) continue;
      CHECK_FALSE(geometry.HasVertexGroups(marker, vertex));
      area += GeometryToolbox::Norm(geometry.GetnDim(), geometry.vertex[marker][vertex]->GetNormal());
    }
    CHECK(geometry.GetSurfaceArea(test.config.get(), marker) == Approx(SU2_TYPE::GetValue(area)));
  }
}

#ifdef HAVE_MPI
TEST_CASE("Boundary face groups: partitioned owned stars match serial", "[BoundaryFaceGroups]") {
  if (SU2_MPI::GetSize() != 2) return;
  const auto comm = SU2_MPI::GetComm();
  if (SU2_MPI::GetRank() == 0) WriteCubeMesh("face_groups_partition.su2", 1);
  SU2_MPI::Barrier(comm);
  const std::string options =
      "SOLVER= EULER\nMARKER_EULER= ( wall )\n"
      "MESH_FILENAME= face_groups_partition.su2\nMESH_FORMAT= SU2\n"
      "MACH_NUMBER= 0.3\nMARKER_FAR= ( far )\n";
  using FaceKey = std::vector<unsigned long>;
  using GroupKey = std::pair<std::vector<FaceKey>, unsigned long>;
  using Star = std::vector<GroupKey>;
  auto collect = [](const CGeometry& geometry, const CConfig& config) {
    std::map<unsigned long, Star> result;
    for (const auto& point : geometry.boundaryFaceGroupPoints) {
      if (config.GetMarker_All_TagBound(point.faces.front()[0]) != "wall") continue;
      Star groups(point.neighbor.size());
      for (size_t q = 0; q < groups.size(); ++q) groups[q].second = geometry.nodes->GetGlobalIndex(point.neighbor[q]);
      for (size_t f = 0; f < point.faces.size(); ++f) {
        const auto& face = point.faces[f];
        const auto* elem = geometry.bound[face[0]][face[1]];
        FaceKey key;
        for (unsigned short i = 0; i < elem->GetnNodes(); ++i)
          key.push_back(geometry.nodes->GetGlobalIndex(elem->GetNode(i)));
        std::sort(key.begin(), key.end());
        groups[point.faceGroup[f]].first.push_back(key);
      }
      for (auto& group : groups) std::sort(group.first.begin(), group.first.end());
      std::sort(groups.begin(), groups.end());
      result[geometry.nodes->GetGlobalIndex(point.point)] = groups;
    }
    for (unsigned short marker = 0; marker < geometry.GetnMarker(); ++marker) {
      if (config.GetMarker_All_TagBound(marker) != "wall") continue;
      for (unsigned long vertex = 0; vertex < geometry.GetnVertex(marker); ++vertex) {
        const auto point = geometry.vertex[marker][vertex]->GetNode();
        if (!geometry.nodes->GetDomain(point) || geometry.HasVertexGroups(marker, vertex)) continue;
        const auto id = geometry.nodes->GetGlobalIndex(point);
        result[id] = Star(1);
        result[id][0].second = geometry.nodes->GetGlobalIndex(geometry.GetVertexNormalNeighbor(marker, vertex));
      }
      for (unsigned long f = 0; f < geometry.nElem_Bound[marker]; ++f) {
        const auto* elem = geometry.bound[marker][f];
        FaceKey key;
        for (unsigned short i = 0; i < elem->GetnNodes(); ++i)
          key.push_back(geometry.nodes->GetGlobalIndex(elem->GetNode(i)));
        std::sort(key.begin(), key.end());
        for (unsigned short i = 0; i < elem->GetnNodes(); ++i) {
          const auto point = elem->GetNode(i);
          if (!geometry.nodes->GetDomain(point)) continue;
          const auto vertex = geometry.nodes->GetVertex(point, marker);
          if (vertex < 0 || geometry.HasVertexGroups(marker, vertex)) continue;
          result[geometry.nodes->GetGlobalIndex(point)][0].first.push_back(key);
        }
      }
    }
    for (auto& point : result)
      for (auto& group : point.second) std::sort(group.first.begin(), group.first.end());
    return result;
  };
  SU2_MPI::SetComm(MPI_COMM_SELF);
  const auto serial = BuildGeometry(options);
  const auto reference = collect(*serial.geometry, *serial.config);
  SU2_MPI::SetComm(comm);
  const auto parallel = BuildGeometry(options);
  const auto partitioned = collect(*parallel.geometry, *parallel.config);
  CHECK_FALSE(partitioned.empty());
  for (const auto& item : partitioned) {
    REQUIRE(reference.count(item.first) == 1);
    CHECK(item.second == reference.at(item.first));
  }
}
#endif

TEST_CASE("Boundary face groups: multigrid keeps legacy geometry", "[BoundaryFaceGroups]") {
  if (SU2_MPI::GetSize() > 1) return;
  WriteCubeMesh("face_groups_multigrid.su2", 1);
  const auto test = BuildGeometry(
      "SOLVER= EULER\nMARKER_EULER= ( wall )\nMGLEVEL= 1\n"
      "MESH_FILENAME= face_groups_multigrid.su2\nMESH_FORMAT= SU2\n"
      "MACH_NUMBER= 0.3\nMARKER_FAR= ( far )\n");
  CHECK_FALSE(CGeometry::UseBoundaryFaceGroups(test.config.get()));
  for (unsigned short marker = 0; marker < test.geometry->GetnMarker(); ++marker)
    for (unsigned long vertex = 0; vertex < test.geometry->GetnVertex(marker); ++vertex)
      CHECK_FALSE(test.geometry->HasVertexGroups(marker, vertex));
}

TEST_CASE("Boundary face groups: mixed thermal markers at one point", "[BoundaryFaceGroups]") {
  if (SU2_MPI::GetSize() > 1) return;
  WriteCubeMesh("face_groups_mixed_thermal.su2", 3);
  const auto test = BuildGeometry(
      "SOLVER= NAVIER_STOKES\nREYNOLDS_NUMBER= 1000000\nMARKER_HEATFLUX= ( wallx, 1.0 )\n"
      "MARKER_ISOTHERMAL= ( wally, 300.0 )\n"
      "MARKER_HEATTRANSFER= ( wallz, 10.0, 350.0 )\n"
      "MESH_FILENAME= face_groups_mixed_thermal.su2\nMESH_FORMAT= SU2\n"
      "MACH_NUMBER= 0.3\nMARKER_FAR= ( far )\n");
  const auto& geometry = *test.geometry;
  unsigned short markers = 0;
  double areaSum = 0.0;
  std::array<bool, 3> directions = {};
  for (unsigned short marker = 0; marker < geometry.GetnMarker(); ++marker) {
    const auto tag = test.config->GetMarker_All_TagBound(marker);
    if (tag != "wallx" && tag != "wally" && tag != "wallz") continue;
    for (unsigned long vertex = 0; vertex < geometry.GetnVertex(marker); ++vertex) {
      const auto point = geometry.vertex[marker][vertex]->GetNode();
      if (PointKey(geometry, point, 0.5) != Key{0, 0, 0}) continue;
      ++markers;
      REQUIRE(geometry.HasVertexGroups(marker, vertex));
      REQUIRE(geometry.GetnVertexGroups(marker, vertex) == 1);
      su2double normal[3] = {}, area = 0.0, distance = 0.0;
      unsigned long neighbor = 0;
      geometry.GetVertexGroup(marker, vertex, 0, normal, area, neighbor, distance);
      areaSum += SU2_TYPE::GetValue(area);
      for (unsigned short d = 0; d < 3; ++d)
        if (std::fabs(SU2_TYPE::GetValue(normal[d])) > 0.99) directions[d] = true;
    }
  }
  CHECK(markers == 3);
  CHECK(areaSum == Approx(3.0 / 16.0).margin(1e-14));
  for (bool present : directions) CHECK(present);
}

TEST_CASE("Boundary face groups: exactly cancelling cusp", "[BoundaryFaceGroups]") {
  if (SU2_MPI::GetSize() > 1) return;
  WriteLMesh("face_groups_exact_cusp.su2", false);
  const auto test = BuildGeometry(
      "SOLVER= EULER\nMARKER_EULER= ( wall )\n"
      "MESH_FILENAME= face_groups_exact_cusp.su2\nMESH_FORMAT= SU2\n"
      "MACH_NUMBER= 0.3\nMARKER_FAR= ( far )\n");
  auto& geometry = *test.geometry;
  unsigned long cusp = geometry.GetnPoint(), neighbor = geometry.GetnPoint();
  for (unsigned long point = 0; point < geometry.GetnPoint(); ++point) {
    if (PointKey(geometry, point, 0.5) == Key{4, 0, 0}) cusp = point;
    if (PointKey(geometry, point, 0.5) == Key{4, 1, 0}) neighbor = point;
  }
  REQUIRE(cusp < geometry.GetnPoint());
  REQUIRE(neighbor < geometry.GetnPoint());
  geometry.nodes->SetCoord(neighbor, 0, 1.5);
  geometry.nodes->SetCoord(neighbor, 1, 0.0);
  geometry.SetBoundControlVolume(test.config.get(), UPDATE);
  geometry.ComputeBoundaryFaceGroups(test.config.get());
  bool found = false;
  for (unsigned short marker = 0; marker < geometry.GetnMarker(); ++marker) {
    if (test.config->GetMarker_All_TagBound(marker) != "wall") continue;
    const auto vertex = geometry.nodes->GetVertex(cusp, marker);
    REQUIRE(vertex >= 0);
    const auto* normal = geometry.vertex[marker][vertex]->GetNormal();
    CHECK(normal[0] == EPS * EPS);
    CHECK(normal[1] == 0.0);
    CHECK(geometry.HasVertexGroups(marker, vertex));
    CHECK(geometry.GetnVertexGroups(marker, vertex) == 2);
    su2double basis[3][3] = {};
    CHECK(geometry.GetSlipBasis(marker, vertex, basis) == 1);
    found = true;
  }
  CHECK(found);
}

TEST_CASE("Boundary face groups: cancelling curved markers keep finite symmetry normals", "[BoundaryFaceGroups]") {
  if (SU2_MPI::GetSize() > 1) return;
  WriteLMesh("face_groups_split_cusp.su2", true);
  const auto test = BuildGeometry(
      "SOLVER= EULER\nMARKER_EULER= ( wall, step )\n"
      "MESH_FILENAME= face_groups_split_cusp.su2\nMESH_FORMAT= SU2\n"
      "MACH_NUMBER= 0.3\nMARKER_FAR= ( far )\n");
  auto& geometry = *test.geometry;
  unsigned long cusp = geometry.GetnPoint(), neighbor = geometry.GetnPoint();
  for (unsigned long point = 0; point < geometry.GetnPoint(); ++point) {
    if (PointKey(geometry, point, 0.5) == Key{4, 2, 0}) cusp = point;
    if (PointKey(geometry, point, 0.5) == Key{4, 1, 0}) neighbor = point;
  }
  REQUIRE(cusp < geometry.GetnPoint());
  REQUIRE(neighbor < geometry.GetnPoint());
  geometry.nodes->SetCoord(neighbor, 0, 1.5);
  geometry.nodes->SetCoord(neighbor, 1, 1.0);
  geometry.SetBoundControlVolume(test.config.get(), UPDATE);
  geometry.ComputeBoundaryFaceGroups(test.config.get());
  unsigned short found = 0;
  su2double normalSum[3] = {};
  for (unsigned short marker = 0; marker < geometry.GetnMarker(); ++marker) {
    const auto tag = test.config->GetMarker_All_TagBound(marker);
    if (tag != "wall" && tag != "step") continue;
    const auto vertex = geometry.nodes->GetVertex(cusp, marker);
    REQUIRE(vertex >= 0);
    REQUIRE(geometry.HasVertexGroups(marker, vertex));
    su2double basis[3][3] = {};
    CHECK(geometry.GetSlipBasis(marker, vertex, basis) == 1);
    for (unsigned short dim = 0; dim < geometry.GetnDim(); ++dim)
      normalSum[dim] += geometry.vertex[marker][vertex]->GetNormal(dim);
    const auto normal = geometry.symmetryNormals[marker].find(vertex);
    if (normal != geometry.symmetryNormals[marker].end())
      for (unsigned short dim = 0; dim < geometry.GetnDim(); ++dim)
        CHECK(std::isfinite(SU2_TYPE::GetValue(normal->second[dim])));
    ++found;
  }
  CHECK(found == 2);
  for (unsigned short dim = 0; dim < geometry.GetnDim(); ++dim) CHECK(normalSum[dim] == 0.0);
}

TEST_CASE("Grouped isothermal wall assembles its temperature linearization", "[BoundaryFaceGroups]") {
  if (SU2_MPI::GetSize() > 1) return;
  const std::string mesh = "face_groups_isothermal_assembly.su2";
  WriteLMesh(mesh, false);
  struct Values {
    double residual, jacobian, temperatureRef, density, gas, gamma;
  };
  auto evaluate = [&](double wallTemperature) {
    auto test = BuildGeometry(
        "SOLVER= NAVIER_STOKES\nREYNOLDS_NUMBER= 1000000\n"
        "MESH_FILENAME= " +
        mesh +
        "\nMESH_FORMAT= SU2\nMACH_NUMBER= 0.3\n"
        "TIME_DISCRE_FLOW= EULER_IMPLICIT\nMARKER_ISOTHERMAL= ( wall, " +
        std::to_string(wallTemperature) + " )\nMARKER_FAR= ( far )\n");
    auto* geometry = test.geometry.get();
    auto* config = test.config.get();
    geometry->SetGlobal_to_Local_Point();
    geometry->PreprocessP2PComms(geometry, config);
    CNSSolver solver(geometry, config, MESH_0);
    auto* nodes = static_cast<CNSVariable*>(solver.GetNodes());
    for (unsigned long point = 0; point < geometry->GetnPoint(); ++point) {
      nodes->SetPrimVar(point, 0.0, 0.0, solver.GetFluidModel());
      nodes->SetSecondaryVar(point, solver.GetFluidModel());
    }
    unsigned short wall = 0;
    while (wall < geometry->GetnMarker() && config->GetMarker_All_TagBound(wall) != "wall") ++wall;
    REQUIRE(wall < geometry->GetnMarker());
    unsigned long vertex = 0;
    while (vertex < geometry->GetnVertex(wall) &&
           !(PointKey(*geometry, geometry->vertex[wall][vertex]->GetNode(), 0.5) == Key{4, 2, 0}))
      ++vertex;
    REQUIRE(vertex < geometry->GetnVertex(wall));
    REQUIRE(geometry->HasVertexGroups(wall, vertex));
    REQUIRE(geometry->HasEligibleGroupNeighbors(wall, vertex));
    const auto point = geometry->vertex[wall][vertex]->GetNode();
    solver.LinSysRes.SetValZero();
    solver.Jacobian.SetValZero();
    CSolver* solvers[MAX_SOLS] = {};
    solvers[FLOW_SOL] = &solver;
    solver.BC_Isothermal_Wall(geometry, solvers, nullptr, nullptr, config, wall);
    const auto energy = geometry->GetnDim() + 1;
    const auto block = solver.Jacobian.GetBlockView(point, point);
    REQUIRE(block);
    return Values{SU2_TYPE::GetValue(solver.LinSysRes(point, energy)), SU2_TYPE::GetValue(block(energy, energy)),
                  SU2_TYPE::GetValue(config->GetTemperature_Ref()),    SU2_TYPE::GetValue(nodes->GetDensity(point)),
                  SU2_TYPE::GetValue(config->GetGas_ConstantND()),     SU2_TYPE::GetValue(config->GetGamma())};
  };
  const double wall = 300.0, step = 0.003;
  const auto center = evaluate(wall);
  REQUIRE(center.jacobian > 0.0);
  const auto plus = evaluate(wall + step);
  const auto minus = evaluate(wall - step);
  const double dResidualdT = (plus.residual - minus.residual) * center.temperatureRef / (2.0 * step);
  // The existing implicit wall linearization uses the negative temperature-residual derivative.
  CHECK(center.jacobian == Approx(-dResidualdT * (center.gamma - 1.0) / (center.gas * center.density)).epsilon(1e-7));
}
