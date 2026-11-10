/*!
 * \file sharp_edge_normals_tests.cpp
 * \brief Unit tests of the grouping of the dual area pieces of a boundary point at sharp edges.
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
#include <random>
#include <vector>
#include "../../../Common/include/geometry/CGeometry.hpp"
#include "../../../Common/include/toolboxes/sharp_edge_normals.hpp"
#include "../../../SU2_CFD/include/gradients/correctGradientsSymmetry.hpp"

namespace {

/*--- Fan of six triangles around the origin: ring node k at angle 60 k deg; faces of the upper half (y > 0) are rotated
 * about the x axis by `fold` (0: flat). Returns the corner vectors and the edge neighbours of the origin. ---*/
void Fan(double fold, double (*N)[3], unsigned long (*edge)[2], unsigned short* nEdge) {
  const double pi = 4.0 * std::atan(1.0);
  double x[7][3] = {{0.0}};
  for (int k = 0; k < 6; ++k) {
    const double c = std::cos(pi * k / 3.0), s = std::sin(pi * k / 3.0);
    x[k + 1][0] = c;
    x[k + 1][1] = (s > 1e-12) ? s * std::cos(fold) : s;
    x[k + 1][2] = (s > 1e-12) ? s * std::sin(fold) : 0.0;
  }
  for (int f = 0; f < 6; ++f) {
    const int a = f + 1, b = (f + 1) % 6 + 1;
    const double* u = x[a];
    const double* v = x[b];
    N[f][0] = (u[1] * v[2] - u[2] * v[1]) / 6.0;
    N[f][1] = (u[2] * v[0] - u[0] * v[2]) / 6.0;
    N[f][2] = (u[0] * v[1] - u[1] * v[0]) / 6.0;
    edge[f][0] = a;
    edge[f][1] = b;
    nEdge[f] = 2;
  }
}

}  // namespace

TEST_CASE("Sharp-edge normals: face groups", "[SharpEdgeNormals]") {
  namespace se = SharpEdgeNormals;
  const double pi = 4.0 * std::atan(1.0);
  double N[6][3], g[32][3];
  unsigned long edge[6][2];
  unsigned short nEdge[6];

  /*--- Flat fan and a moderate fold: one group. ---*/
  Fan(0.0, N, edge, nEdge);
  CHECK(se::FaceGroups(3, 6, N, edge, nEdge, g) == 1);
  Fan(40.0 * pi / 180.0, N, edge, nEdge);
  CHECK(se::FaceGroups(3, 6, N, edge, nEdge, g) == 1);

  /*--- Thin wedge (164 deg fold, a sharp trailing edge): two groups, each parallel to its faces; the vertex normal
   * (sum) leaves a slip plane that crosses the faces almost normally, the two group normals restrict the velocity to
   * the edge. ---*/
  Fan(164.0 * pi / 180.0, N, edge, nEdge);
  unsigned short faceGroup[6];
  REQUIRE(se::FaceGroups(3, 6, N, edge, nEdge, g, faceGroup) == 2);
  for (int f = 0; f < 6; ++f) {
    const auto& q = g[faceGroup[f]];
    const double c = (q[0] * N[f][0] + q[1] * N[f][1] + q[2] * N[f][2]) /
                     (std::sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2]) *
                      std::sqrt(N[f][0] * N[f][0] + N[f][1] * N[f][1] + N[f][2] * N[f][2]));
    CHECK(c == Approx(1.0).margin(1e-12));
  }
  double unit[2][3], Q[3][3] = {{0.0}};
  for (int q = 0; q < 2; ++q) {
    const double n = std::sqrt(g[q][0] * g[q][0] + g[q][1] * g[q][1] + g[q][2] * g[q][2]);
    for (int d = 0; d < 3; ++d) unit[q][d] = g[q][d] / n;
  }
  REQUIRE(se::Basis(3, 2, unit, Q) == 2);
  double m[5] = {0.0, 0.3, -0.7, 0.4, 0.0};
  se::Project(3, 2, Q, 1, m);
  CHECK(std::fabs(m[2]) < 1e-14);
  CHECK(std::fabs(m[3]) < 1e-14);
  CHECK(m[1] == Approx(0.3).margin(1e-14));

  /*--- A chain of positive adjacent products with a negative pair (0, 2): face 2 stays apart, and every group normal
   * lies in the dual cone of its faces. ---*/
  {
    const double r3 = std::sqrt(3.0);
    double M[3][3] = {{0.0, 0.0, 1.0}, {0.5 * r3, 0.0, 0.5}, {r3, 3.0, -2.0}};
    unsigned long e3[3][2] = {{1, 2}, {2, 3}, {3, 4}};
    unsigned short n3[3] = {2, 2, 2};
    unsigned short fg[3];
    REQUIRE(se::FaceGroups(3, 3, M, e3, n3, g, fg) == 2);
    for (int f = 0; f < 3; ++f) {
      const auto& q = g[fg[f]];
      CHECK(q[0] * M[f][0] + q[1] * M[f][1] + q[2] * M[f][2] >=
            M[f][0] * M[f][0] + M[f][1] * M[f][1] + M[f][2] * M[f][2] - 1e-14);
    }
  }

  /*--- 90 deg fold: split; faces sharing only a vertex that is not an edge neighbour: never joined. ---*/
  {
    double M[2][3] = {{0.0, 0.0, 1.0}, {0.0, 1.0, 0.0}};
    unsigned long e2[2][2] = {{1, 2}, {2, 3}};
    unsigned short n2[2] = {2, 2};
    CHECK(se::FaceGroups(3, 2, M, e2, n2, g) == 2);
    M[1][1] = 0.6;
    M[1][2] = 0.8;
    CHECK(se::FaceGroups(3, 2, M, e2, n2, g) == 1);
    e2[1][0] = 5;
    e2[1][1] = 6;
    CHECK(se::FaceGroups(3, 2, M, e2, n2, g) == 2);
  }

  /*--- 2D: the two segments of a point; smooth (one group), thin wedge (two). ---*/
  {
    double M[2][3] = {{0.0, -1.0, 0.0}, {std::sin(0.1), -std::cos(0.1), 0.0}};
    unsigned long e2[2][2] = {{0, 0}, {0, 0}};
    unsigned short n2[2] = {0, 0};
    CHECK(se::FaceGroups(2, 2, M, e2, n2, g) == 1);
    const double w = 164.0 * pi / 180.0;
    M[1][0] = std::sin(w);
    M[1][1] = -std::cos(w);
    CHECK(se::FaceGroups(2, 2, M, e2, n2, g) == 2);
  }
}

namespace {

using Vec3 = std::array<double, 3>;

/*--- Groups of the faces (corner vectors N, edge neighbours, optional link classes) for every order of the faces: all
 * n! orders up to n = 7, otherwise 2000 random orders (fixed seed). Checks that the partition (groups numbered in the
 * order of their first face in the original order) is the expected one in every order, that every group normal lies in
 * the dual cone of its faces (n_q . N_f >= |N_f|^2) and that the basis of the group normals has the expected rank.
 * ---*/
void CheckAllOrders(unsigned short nDim, const std::vector<Vec3>& N,
                    const std::vector<std::array<unsigned long, 2>>& edge, const std::vector<unsigned short>& linkClass,
                    const std::vector<unsigned short>& expected, unsigned short expectedRank) {
  const auto n = static_cast<unsigned short>(N.size());
  std::vector<unsigned short> perm(n);
  for (unsigned short f = 0; f < n; ++f) perm[f] = f;
  std::mt19937 gen(2026);
  bool samePartition = true, dualCone = true, sameRank = true;
  for (int trial = 0;; ++trial) {
    su2double Np[32][3] = {{0.0}}, groupN[32][3], unit[32][3], Q[3][3];
    unsigned long ep[32][2] = {{0}};
    unsigned short nep[32], lc[32], fg[32];
    for (unsigned short a = 0; a < n; ++a) {
      for (unsigned short d = 0; d < 3; ++d) Np[a][d] = N[perm[a]][d];
      ep[a][0] = edge[perm[a]][0];
      ep[a][1] = edge[perm[a]][1];
      nep[a] = (nDim == 3) ? 2 : 0;
      lc[a] = linkClass.empty() ? 0 : linkClass[perm[a]];
    }
    const auto nGroups =
        SharpEdgeNormals::FaceGroups(nDim, n, Np, ep, nep, groupN, fg, linkClass.empty() ? nullptr : lc);
    std::vector<unsigned short> byFace(n), canon(n), slot(32, 99);
    for (unsigned short a = 0; a < n; ++a) byFace[perm[a]] = fg[a];
    unsigned short next = 0;
    for (unsigned short f = 0; f < n; ++f) {
      if (slot[byFace[f]] == 99) slot[byFace[f]] = next++;
      canon[f] = slot[byFace[f]];
    }
    samePartition = samePartition && (canon == expected);
    for (unsigned short a = 0; a < n; ++a) {
      double nf = 0.0, ff = 0.0;
      for (unsigned short d = 0; d < nDim; ++d) {
        nf += groupN[fg[a]][d] * Np[a][d];
        ff += Np[a][d] * Np[a][d];
      }
      dualCone = dualCone && (nf >= ff * (1.0 - 1e-12));
    }
    for (unsigned short q = 0; q < nGroups; ++q) {
      double len = 0.0;
      for (unsigned short d = 0; d < nDim; ++d) len += groupN[q][d] * groupN[q][d];
      for (unsigned short d = 0; d < 3; ++d) unit[q][d] = groupN[q][d] / std::sqrt(len);
    }
    sameRank = sameRank && (SharpEdgeNormals::Basis(nDim, nGroups, unit, Q) == expectedRank);
    if (n <= 7) {
      if (!std::next_permutation(perm.begin(), perm.end())) break;
    } else {
      if (trial == 2000) break;
      std::shuffle(perm.begin(), perm.end(), gen);
    }
  }
  CHECK(samePartition);
  CHECK(dualCone);
  CHECK(sameRank);
}

/*--- Edge neighbours of a closed fan of n faces around a point: face f has the other nodes f + 1 and f + 2 (mod n), so
 * face f shares an edge with faces f - 1 and f + 1. ---*/
std::vector<std::array<unsigned long, 2>> Ring(unsigned long n) {
  std::vector<std::array<unsigned long, 2>> e(n);
  for (unsigned long f = 0; f < n; ++f) e[f] = {f + 1, (f + 1) % n + 1};
  return e;
}

/*--- Open chain of n faces: face f shares an edge with faces f - 1 and f + 1 only. ---*/
std::vector<std::array<unsigned long, 2>> Chain(unsigned long n) {
  std::vector<std::array<unsigned long, 2>> e(n);
  for (unsigned long f = 0; f < n; ++f) e[f] = {f + 1, f + 2};
  return e;
}

/*--- Unit normal rotated by t degrees from z about the x axis. ---*/
Vec3 RotX(double t) {
  const double r = t * 4.0 * std::atan(1.0) / 180.0;
  return {0.0, std::sin(r), std::cos(r)};
}

Vec3 Scaled(const Vec3& v, double s) { return {s * v[0], s * v[1], s * v[2]}; }

}  // namespace

TEST_CASE("Sharp-edge normals: order-independent groups", "[SharpEdgeNormals]") {
  const double pi = 4.0 * std::atan(1.0), deg = pi / 180.0;
  const std::vector<unsigned short> none;

  /*--- Box corner: three planes at 90 deg, unequal faces (3 + 3 + 2, more than 7 faces: random orders). Three groups,
   * rank 3 (zero velocity at a convex corner). ---*/
  {
    const std::vector<Vec3> N = {{1.0, 0, 0}, {0.4, 0, 0}, {0.8, 0, 0}, {0, 0.7, 0},
                                 {0, 2.0, 0}, {0, 0.1, 0}, {0, 0, 0.3}, {0, 0, 1.1}};
    CheckAllOrders(3, N, Ring(8), none, {0, 0, 0, 1, 1, 1, 2, 2}, 3);
  }

  /*--- Apex of a pyramid with 4 faces, normals (+-c, 0, s) and (0, +-c, s), s = sin(b): adjacent faces always
   * compatible (cosine s^2), opposite faces at 180 - 2b deg. b = 30 deg (opposite faces 120 deg apart): the ring is not
   * pairwise compatible and all its links have exactly the same cosine, so all are removed: four groups, rank 3. b = 60
   * deg (opposite faces 60 deg apart): one group. ---*/
  for (const double b : {30.0, 60.0}) {
    const double c = std::cos(b * deg), s = std::sin(b * deg);
    const std::vector<Vec3> N = {{c, 0, s}, {0, c, s}, {-c, 0, s}, {0, -c, s}};
    if (b < 45.0)
      CheckAllOrders(3, N, Ring(4), none, {0, 1, 2, 3}, 3);
    else
      CheckAllOrders(3, N, Ring(4), none, {0, 0, 0, 0}, 1);
  }

  /*--- Apex of an irregular pyramid with 5 faces (azimuths 0, 64, 132, 204, 280 deg: gaps of 64, 68, 72, 76, 80 deg),
   * unequal face areas. b = 30 deg: faces two apart are incompatible; the links break in the order of decreasing gap
   * until every component is compatible, faces 0 and 1 (64 deg gap) stay together: four groups, rank 3. b = 60 deg: all
   * faces pairwise compatible, one group. ---*/
  for (const double b : {30.0, 60.0}) {
    const double az[5] = {0.0, 64.0, 132.0, 204.0, 280.0};
    std::vector<Vec3> N;
    for (int k = 0; k < 5; ++k)
      N.push_back(Scaled(
          {std::cos(b * deg) * std::cos(az[k] * deg), std::cos(b * deg) * std::sin(az[k] * deg), std::sin(b * deg)},
          1.0 + 0.3 * k));
    if (b < 45.0)
      CheckAllOrders(3, N, Ring(5), none, {0, 0, 1, 2, 3}, 3);
    else
      CheckAllOrders(3, N, Ring(5), none, {0, 0, 0, 0, 0}, 1);
  }

  /*--- Polyhedron corner with folds below and above 90 deg between adjacent faces: top T (normal z), chamfer C (40 deg
   * from T), side S (normal y, 50 deg from C), front F (100 deg from S), bevel B (36 deg from F, 55 deg from T). The
   * component is not compatible (T . S = 0); its weakest link B-T (55 deg) breaks first, then C-S (50 deg) in {T, C,
   * S}: groups {T, C}, {S}, {F, B}, rank 3. ---*/
  {
    const Vec3 F = {std::cos(10 * deg), -std::sin(10 * deg), 0.0}, B = {std::sin(55 * deg), 0.0, std::cos(55 * deg)};
    const std::vector<Vec3> N = {RotX(0.0), Scaled(RotX(40.0), 0.5), Scaled(RotX(90.0), 2.0), Scaled(F, 0.7),
                                 Scaled(B, 1.3)};
    CheckAllOrders(3, N, Ring(5), none, {0, 0, 1, 2, 2}, 3);
  }

  /*--- A face compatible with two mutually incompatible groups (corner of a wing tip): upper (0, 1) and lower (4, 5)
   * faces of a 10 deg wedge, a cap (2, 3) at 90 - e deg to both halves. Cap at 85 deg to both (exactly equal cosines):
   * its own group, rank 3; tilted by 1 deg towards the upper half: it joins that half, rank 2; cap at 95 deg: three
   * groups. ---*/
  for (const double capDeg : {85.0, 95.0}) {
    const double w = 5.0 * deg, e = (90.0 - capDeg) * deg;
    const Vec3 U = {std::sin(w) * std::cos(e), std::sin(e), std::cos(w) * std::cos(e)};
    std::vector<Vec3> N = {
        U, Scaled(U, 0.6), {0, 0.8, 0}, {0, 0.8, 0}, Scaled({U[0], U[1], -U[2]}, 0.6), {U[0], U[1], -U[2]}};
    CheckAllOrders(3, N, Ring(6), none, {0, 0, 1, 1, 2, 2}, 3);
    if (capDeg < 90.0) {
      N[2] = N[3] = {0, 0.8 * std::cos(deg), 0.8 * std::sin(deg)};
      CheckAllOrders(3, N, Ring(6), none, {0, 0, 0, 0, 1, 1}, 2);
    }
  }

  /*--- A face X compatible with three mutually incompatible groups A, B, C (normals 85 deg from X at azimuths 60, 180,
   * -60 deg; X shares an edge with A and C). Exactly equal cosines: X forms its own group (4 groups, rank 3); X tilted
   * by 1 deg towards A: it joins A (3 groups, rank 3). ---*/
  {
    const double s = std::sin(85 * deg), c = std::cos(85 * deg), r3 = std::sqrt(3.0);
    const Vec3 A = {0.5 * s, 0.5 * r3 * s, c}, B = {-s, 0.0, c}, C = {0.5 * s, -0.5 * r3 * s, c};
    std::vector<Vec3> N = {{0, 0, 1}, A, Scaled(A, 0.6), B, Scaled(B, 1.2), Scaled(C, 0.8), C};
    CheckAllOrders(3, N, Ring(7), none, {0, 1, 1, 2, 2, 3, 3}, 3);
    N[0] = {0.5 * std::sin(deg), 0.5 * r3 * std::sin(deg), std::cos(deg)};
    CheckAllOrders(3, N, Ring(7), none, {0, 0, 0, 1, 1, 2, 2}, 3);
  }

  /*--- Longer component (chain A' - A - X - B, normals at -60, -20, 33, 93 deg about one axis): X is compatible with A
   * (cosine 0.60) and B (0.50), which are incompatible, but X does not join A: the weakest link X-B breaks first, and
   * in {A', A, X} (A' . X < 0) the link A-X is then the weakest. Groups {A', A}, {X}, {B}. ---*/
  {
    const std::vector<Vec3> N = {RotX(-60.0), RotX(-20.0), RotX(33.0), RotX(93.0)};
    CheckAllOrders(3, N, Chain(4), none, {0, 0, 1, 2}, 2);
  }

  /*--- Link classes: a flat fan whose faces 0-2 and 3-5 are in different classes forms two groups (parallel normals:
   * rank 1); without classes one group. ---*/
  {
    const std::vector<Vec3> N = {{0, 0, 1}, {0, 0, 2}, {0, 0, 1}, {0, 0, 1}, {0, 0, 0.5}, {0, 0, 1}};
    CheckAllOrders(3, N, Ring(6), {0, 0, 0, 7, 7, 7}, {0, 0, 0, 1, 1, 1}, 1);
    CheckAllOrders(3, N, Ring(6), none, {0, 0, 0, 0, 0, 0}, 1);
  }

  /*--- Compatibility threshold: cosine 0.0005 (below COMPATIBLE_COS) is a sharp edge, cosine 0.002 is not. ---*/
  for (const double t : {0.0005, 0.002}) {
    const std::vector<Vec3> N = {{0.0, 0.0, 1.0}, {0.0, std::sqrt(1.0 - t * t), t}};
    CheckAllOrders(3, N, Chain(2), none,
                   t < 0.001 ? std::vector<unsigned short>{0, 1} : std::vector<unsigned short>{0, 0},
                   t < 0.001 ? 2 : 1);
  }

  /*--- ONERA M6 wing (hybrid mesh), two points of the closed wing tip: point 46 (trailing edge on the tip ridge z = 0)
   * and point 47 (tip ridge). Lower faces (z component < 0) and upper faces fold by 108 deg across the trailing edge
   * and by 88 deg across the ridge. Groups: lower and upper halves in every order (rank 2: velocity along the sharp
   * edge). A greedy merge in face order gives 2 or 3 groups at point 46 depending on the order. ---*/
  {
    const std::vector<Vec3> N = {{-3.4527637101296927e-06, 0.00013549416073123523, -0.00015240393609880425},
                                 {-3.425224979111915e-06, 0.00013441347941039864, 0.00023411606631673615},
                                 {-5.6299153361818569e-06, 0.0001216527217636005, -0.00011034177632880386},
                                 {-5.8212976089669138e-06, 0.0001257881612846936, 0.00012987386080468954},
                                 {1.2126984884427711e-05, 0.00012108851495268514, -0.00012929906974200344},
                                 {5.2382358223712074e-05, 0.00012513052349801823, 0.00017125460926028838}};
    const std::vector<std::array<unsigned long, 2>> edge = {{45, 1441}, {3123, 45},   {1407, 47},
                                                            {47, 3121}, {1441, 1407}, {3121, 3123}};
    CheckAllOrders(3, N, edge, none, {0, 1, 0, 1, 0, 1}, 2);
  }
  {
    const std::vector<Vec3> N = {{-5.6299153361818536e-06, 0.00012165272176360047, -0.00011034177632880187},
                                 {-5.8212976089685045e-06, 0.00012578816128469322, 0.00012987386080468927},
                                 {-5.8941065175494447e-06, 8.69062243930457e-05, -6.1310200144633052e-05},
                                 {-5.53290568833489e-06, 8.158046378434997e-05, 4.8596030875096337e-05},
                                 {7.18835195206529e-06, 9.514498481312483e-05, -7.9450342581273144e-05},
                                 {5.3368426594275388e-05, 5.4546638530220041e-05, 0.00010323437735996845},
                                 {-4.3552878196415469e-05, 9.095605988952612e-05, 0.0001297251427599337}};
    const std::vector<std::array<unsigned long, 2>> edge = {{46, 1407},   {3121, 46},   {1391, 48},  {48, 3091},
                                                            {1407, 1391}, {3091, 3093}, {3093, 3121}};
    CheckAllOrders(3, N, edge, none, {0, 1, 0, 1, 0, 1, 1}, 2);
  }

  /*--- 2D, the two segments of a point: smooth (10 deg turn, one group), 90 deg corner and 164 deg wedge (two groups,
   * rank 2: zero velocity at the corner). ---*/
  for (const double turn : {10.0, 90.0, 164.0}) {
    const std::vector<Vec3> N = {{0.0, -1.0, 0.0}, Scaled({std::sin(turn * deg), -std::cos(turn * deg), 0.0}, 0.3)};
    const std::vector<std::array<unsigned long, 2>> edge = {{0, 0}, {0, 0}};
    if (turn < 45.0)
      CheckAllOrders(2, N, edge, none, {0, 0}, 1);
    else
      CheckAllOrders(2, N, edge, none, {0, 1}, 2);
  }
}

TEST_CASE("Sharp-edge normals: gradient correction over the basis", "[SharpEdgeNormals]") {
  std::mt19937 gen(515);
  std::uniform_real_distribution<double> uni(-1.0, 1.0);
  for (int trial = 0; trial < 20; ++trial) {
    double n1[3] = {uni(gen), uni(gen), uni(gen)}, n2[3] = {uni(gen), uni(gen), uni(gen)};
    double list[2][3], Q[3][3] = {{0.0}};
    for (int d = 0; d < 3; ++d) {
      list[0][d] = n1[d];
      list[1][d] = n2[d];
    }
    REQUIRE(SharpEdgeNormals::Basis(3, 2, list, Q) == 2);
    double A[5][3], B[5][3];
    for (int v = 0; v < 5; ++v)
      for (int d = 0; d < 3; ++d) A[v][d] = B[v][d] = uni(gen);
    detail::correctGradient<3>(0, 5, 1, Q[0], A);
    detail::correctGradient<3>(0, 5, 1, Q[1], A);
    detail::correctGradient<3>(0, 5, 1, Q[1], B);
    detail::correctGradient<3>(0, 5, 1, Q[0], B);
    for (int v = 0; v < 5; ++v)
      for (int d = 0; d < 3; ++d) CHECK(A[v][d] == Approx(B[v][d]).margin(1e-13));
    /*--- Idempotent. ---*/
    detail::correctGradient<3>(0, 5, 1, Q[0], B);
    detail::correctGradient<3>(0, 5, 1, Q[1], B);
    for (int v = 0; v < 5; ++v)
      for (int d = 0; d < 3; ++d) CHECK(A[v][d] == Approx(B[v][d]).margin(1e-13));
    /*--- Both reflections (Blazek eq. 8.40 for q1 and for q2): scalars have no derivative along either, and the
     * velocity gradient has no entry coupling a crease direction with another direction. ---*/
    double t[3] = {Q[0][1] * Q[1][2] - Q[0][2] * Q[1][1], Q[0][2] * Q[1][0] - Q[0][0] * Q[1][2],
                   Q[0][0] * Q[1][1] - Q[0][1] * Q[1][0]};
    const double* dirs[3] = {Q[0], Q[1], t};
    for (const int v : {0, 4})
      for (int r = 0; r < 2; ++r) CHECK(std::fabs(A[v][0] * Q[r][0] + A[v][1] * Q[r][1] + A[v][2] * Q[r][2]) < 1e-13);
    for (int x = 0; x < 3; ++x)
      for (int y = 0; y < 3; ++y) {
        if (x == y) continue;
        double e = 0.0;
        for (int i = 0; i < 3; ++i)
          for (int j = 0; j < 3; ++j) e += dirs[x][i] * A[1 + i][j] * dirs[y][j];
        CHECK(std::fabs(e) < 1e-13);
      }
  }
}
