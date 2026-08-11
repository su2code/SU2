/*!
 * \file CTrapezoidalMap.hpp
 * \brief Memory-efficient trapezoidal map for 2D lookup table queries,
 *        based on the LUT implementation of P. Gomes (https://github.com/pcarruscag/LUT).
 * \author T. Kiymaz
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

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <numeric>
#include <tuple>
#include <utility>
#include <vector>

#include "../basic_types/datatype_structure.hpp"

namespace su2_lut {

using IntT = int32_t;
using RealT = su2double;

/*--- Simple row-major matrix with a fixed number of columns. ---*/
template <typename T, size_t N>
struct Matrix {
  std::vector<T> data;

  void resize(size_t rows, size_t) { data.resize(rows * N); }
  size_t rows() const { return data.size() / N; }

  const T& operator()(size_t i, size_t j) const { return data[i * N + j]; }
  T& operator()(size_t i, size_t j) { return data[i * N + j]; }
};

using Matrix2i = Matrix<IntT, 2>;
using Matrix3i = Matrix<IntT, 3>;
using VectorInt = std::vector<IntT>;
using VectorReal = std::vector<RealT>;

/*--- The map is defined by the limits of the bands in the x direction and a CSR of
 * the edge IDs in each band, sorted by the edge y position at the band midpoint. ---*/
struct TrapezoidalMap {
  VectorInt offsets, edge_id;
  VectorReal x_bands, edge_y;
};

/*!
 * \brief Orders points by ascending x coordinates and updates triangle indices.
 */
inline void ReorderPoints(Matrix3i& triangles, VectorReal& x, VectorReal& y) {
  const IntT n_pts = static_cast<IntT>(x.size());

  std::vector<IntT> perm(n_pts);
  std::iota(perm.begin(), perm.end(), 0);
  std::sort(perm.begin(), perm.end(),
            [&x, &y](const auto i, const auto j) { return x[i] != x[j] ? x[i] < x[j] : y[i] < y[j]; });

  auto reorder = [n_pts, &perm](const auto& v) {
    VectorReal tmp(n_pts);
    for (IntT i = 0; i < n_pts; ++i) {
      tmp[i] = v[perm[i]];
    }
    return tmp;
  };
  x = reorder(x);
  y = reorder(y);

  std::vector<IntT> inv_perm(n_pts);
  for (IntT i = 0; i < n_pts; ++i) {
    inv_perm[perm[i]] = i;
  }
  for (IntT i = 0; i < static_cast<IntT>(triangles.rows()); ++i) {
    for (IntT j = 0; j < 3; ++j) {
      triangles(i, j) = inv_perm[triangles(i, j)];
    }
  }
}

/*!
 * \brief Extracts unique edges from triangles. Edges are defined by two point IDs and
 * up to two adjacent triangles (boundary edges have the second triangle ID < 0).
 */
inline void ExtractEdges(const Matrix3i& triangles, Matrix2i& edge_pts, Matrix2i& edge_faces) {
  std::vector<std::array<IntT, 3>> edges;
  edges.resize(3 * triangles.rows());

  for (IntT i_tri = 0; i_tri < static_cast<IntT>(triangles.rows()); ++i_tri) {
    for (IntT i = 0; i < 3; ++i) {
      const IntT j = (i + 1) % 3;
      const IntT i_pt = std::min(triangles(i_tri, i), triangles(i_tri, j));
      const IntT j_pt = std::max(triangles(i_tri, i), triangles(i_tri, j));
      edges[3 * i_tri + i] = {i_pt, j_pt, i_tri};
    }
  }

  /*--- Sort to identify duplicates. ---*/
  std::sort(edges.begin(), edges.end(),
            [](const auto& a, const auto& b) { return a[0] != b[0] ? (a[0] < b[0]) : (a[1] < b[1]); });

  auto is_equal = [](const auto& a, const auto& b) { return a[0] == b[0] && a[1] == b[1]; };

  IntT n_edges = 1;
  for (IntT i = 1; i < static_cast<IntT>(edges.size()); ++i) {
    n_edges += static_cast<IntT>(!is_equal(edges[i], edges[i - 1]));
  }

  edge_pts.resize(n_edges, 2);
  edge_faces.resize(n_edges, 2);
  IntT pos = 0;

  auto new_edge = [&](const auto& edge) {
    edge_pts(pos, 0) = edge[0];
    edge_pts(pos, 1) = edge[1];
    edge_faces(pos, 0) = edge[2];
    edge_faces(pos, 1) = -1;
    ++pos;
  };

  new_edge(edges[0]);
  for (IntT i = 1; i < static_cast<IntT>(edges.size()); ++i) {
    if (is_equal(edges[i], edges[i - 1])) {
      edge_faces(pos - 1, 1) = edges[i][2];
    } else {
      new_edge(edges[i]);
    }
  }
}

/*!
 * \brief Detects the x bands of the map. One band per unique x coordinate is used unless
 * that would exceed max_bands, in which case equal-width bands are used to limit memory.
 * \return Tuple of (n_bands, x_bands).
 */
inline auto DetectBands(const VectorReal& x, IntT max_bands = 0) {
  if (max_bands <= 0) {
    max_bands = std::min(IntT{5000}, static_cast<IntT>(4.0 * std::sqrt(static_cast<double>(x.size()))));
  }

  IntT n_unique = 1;
  for (IntT i = 1; i < static_cast<IntT>(x.size()); ++i) {
    if (x[i] != x[i - 1]) n_unique++;
  }

  if (n_unique <= max_bands) {
    const IntT n_bands = n_unique - 1;
    VectorReal x_bands(n_unique);
    IntT pos = 0;
    x_bands[pos] = x[0];
    for (IntT i = 1; i < static_cast<IntT>(x.size()); ++i) {
      if (x[i] != x_bands[pos]) {
        x_bands[++pos] = x[i];
      }
    }
    return std::make_tuple(n_bands, std::move(x_bands));
  }

  const RealT x_min = x.front();
  const RealT x_max = x.back();
  const RealT band_width = (x_max - x_min) / max_bands;

  VectorReal x_bands(max_bands + 1);
  for (IntT i = 0; i <= max_bands; ++i) {
    x_bands[i] = x_min + i * band_width;
  }
  x_bands[max_bands] = x_max;

  return std::make_tuple(max_bands, std::move(x_bands));
}

/*!
 * \brief Builds the trapezoidal map for a set of edges (points must be ordered by x).
 */
inline void BuildTrapezoidalMap(const Matrix2i& edge_pts, const VectorReal& x, const VectorReal& y,
                                TrapezoidalMap& map) {
  auto& x_bands = map.x_bands;
  auto& offsets = map.offsets;
  auto& edge_id = map.edge_id;
  auto& edge_y = map.edge_y;

  auto clear_map = [&]() {
    x_bands.clear();
    offsets.clear();
    edge_id.clear();
    edge_y.clear();
  };

  const auto [n_bands, bands] = DetectBands(x);
  x_bands = std::move(bands);

  if (n_bands <= 0) {
    clear_map();
    return;
  }

  auto find_band = [&x_bands, n_bands = n_bands](RealT x_val) -> IntT {
    auto it = std::lower_bound(x_bands.begin(), x_bands.end(), x_val);
    const IntT idx = static_cast<IntT>(it - x_bands.begin());
    return std::min(std::max(IntT{0}, idx - 1), n_bands - 1);
  };

  /*--- Count edges per band. Each edge is stored in every band between the bands of its
   * two endpoints (inclusive), a superset of the bands it overlaps. ---*/
  auto& counts = offsets;
  counts.clear();
  counts.resize(n_bands + 1, 0);

  for (IntT i = 0; i < static_cast<IntT>(edge_pts.rows()); ++i) {
    const IntT band_0 = find_band(x[edge_pts(i, 0)]);
    const IntT band_1 = find_band(x[edge_pts(i, 1)]);

    for (IntT j = std::min(band_0, band_1); j <= std::max(band_0, band_1); ++j) {
      ++counts[j + 1];
    }
  }

  /*--- Convert counts to offsets (CSR format). ---*/
  for (IntT i = 2; i < static_cast<IntT>(offsets.size()); ++i) {
    offsets[i] += offsets[i - 1];
  }

  /*--- Give up (build failure) rather than allocating an excessive amount of memory. ---*/
  const size_t memory_mb = static_cast<size_t>(offsets.back()) * (sizeof(IntT) + sizeof(RealT)) / (1024 * 1024);
  if (memory_mb > 2048) {
    clear_map();
    return;
  }

  edge_id.resize(offsets.back());
  edge_y.resize(offsets.back());
  auto pos = offsets;

  for (IntT i_edge = 0; i_edge < static_cast<IntT>(edge_pts.rows()); ++i_edge) {
    const IntT pt_0 = edge_pts(i_edge, 0);
    const IntT pt_1 = edge_pts(i_edge, 1);
    const RealT x_0 = x[pt_0], y_0 = y[pt_0];
    const RealT x_1 = x[pt_1], y_1 = y[pt_1];

    const IntT band_0 = find_band(x_0);
    const IntT band_1 = find_band(x_1);

    const RealT dx = x_1 - x_0;
    const bool vertical = std::abs(SU2_TYPE::GetValue(dx)) < 1e-30;
    const RealT dy_dx = vertical ? RealT{0} : (y_1 - y_0) / dx;

    for (IntT j = std::min(band_0, band_1); j <= std::max(band_0, band_1); ++j) {
      edge_id[pos[j]] = i_edge;
      const RealT x_mid = (x_bands[j] + x_bands[j + 1]) / 2;
      edge_y[pos[j]] = vertical ? RealT((y_0 + y_1) / 2) : RealT(y_0 + dy_dx * (x_mid - x_0));
      ++pos[j];
    }
  }

  /*--- Sort the edges in each band by y coordinate. ---*/
  std::vector<std::pair<IntT, RealT>> tmp;
  for (IntT i = 0; i < n_bands; ++i) {
    const IntT begin = offsets[i];
    const IntT end = offsets[i + 1];
    if (begin >= end) continue;

    tmp.resize(end - begin);
    for (auto k = begin; k < end; ++k) {
      tmp[k - begin] = {edge_id[k], edge_y[k]};
    }
    std::sort(tmp.begin(), tmp.end(), [](const auto& a, const auto& b) { return a.second < b.second; });
    for (auto k = begin; k < end; ++k) {
      edge_id[k] = tmp[k - begin].first;
      edge_y[k] = tmp[k - begin].second;
    }
  }
}

/*!
 * \brief Returns the IDs of the edges directly below and above a query point
 * (either ID can be -1 if the point is at a boundary).
 */
inline auto QueryTrapezoidalMap(const TrapezoidalMap& map, const Matrix2i& edge_pts, const VectorReal& x_coords,
                                const VectorReal& y_coords, const RealT& x, const RealT& y) {
  if (map.x_bands.size() < 2 || map.offsets.empty()) {
    return std::make_pair(IntT{-1}, IntT{-1});
  }

  const auto& x_bands = map.x_bands;
  const IntT n_bands = static_cast<IntT>(x_bands.size()) - 1;
  auto it = std::lower_bound(x_bands.begin(), x_bands.end(), x);
  const IntT d = static_cast<IntT>(it - x_bands.begin());
  const IntT band_idx = std::min(std::max(IntT{0}, d - 1), n_bands - 1);

  RealT best_y_below = -1e300;
  RealT best_y_above = 1e300;
  IntT edge_below = -1;
  IntT edge_above = -1;

  const IntT begin = map.offsets[band_idx];
  const IntT end = map.offsets[band_idx + 1];

  for (IntT k = begin; k < end; ++k) {
    const IntT e_id = map.edge_id[k];

    const IntT p0 = edge_pts(e_id, 0);
    const IntT p1 = edge_pts(e_id, 1);
    const RealT x0 = x_coords[p0], y0 = y_coords[p0];
    const RealT x1 = x_coords[p1], y1 = y_coords[p1];

    if (x < std::min(x0, x1) - 1e-10 || x > std::max(x0, x1) + 1e-10) {
      continue;
    }

    /*--- y position of the edge at the query x. ---*/
    RealT edge_y_at_x;
    const RealT dx = x1 - x0;
    if (std::abs(SU2_TYPE::GetValue(dx)) < 1e-30) {
      edge_y_at_x = (y0 + y1) / 2.0;
    } else {
      edge_y_at_x = y0 + (x - x0) / dx * (y1 - y0);
    }

    if (edge_y_at_x <= y + 1e-10 && edge_y_at_x > best_y_below) {
      best_y_below = edge_y_at_x;
      edge_below = e_id;
    }
    if (edge_y_at_x >= y - 1e-10 && edge_y_at_x < best_y_above) {
      best_y_above = edge_y_at_x;
      edge_above = e_id;
    }
  }

  return std::make_pair(edge_below, edge_above);
}

/*!
 * \brief Returns the IDs of the triangles adjacent to two query edges (up to 3 triangles).
 */
inline auto AdjacentTriangles(const IntT edge_0, const IntT edge_1, const Matrix2i& edge_faces) {
  std::array<IntT, 3> tris = {-1, -1, -1};
  IntT pos = 0;

  auto insert = [&tris, &pos](const IntT t) {
    if (t < 0) return;
    for (IntT i = 0; i < pos; ++i) {
      if (t == tris[i]) return;
    }
    tris[pos++] = t;
  };

  auto get_tris = [&edge_faces](const IntT e) {
    if (e < 0) return std::array{IntT{-1}, IntT{-1}};
    return std::array{edge_faces(e, 0), edge_faces(e, 1)};
  };

  for (const auto e : {edge_0, edge_1}) {
    for (const auto t : get_tris(e)) {
      insert(t);
    }
  }
  return tris;
}

/*!
 * \brief Computes the barycentric coordinates of point (x_q, y_q) in a triangle.
 */
inline auto TriangleCoords(const IntT i_tri, const Matrix3i& triangles, const VectorReal& x, const VectorReal& y,
                           const RealT x_q, const RealT y_q) {
  const IntT p0 = triangles(i_tri, 0);
  const IntT p1 = triangles(i_tri, 1);
  const IntT p2 = triangles(i_tri, 2);

  const RealT x0 = x[p0], y0 = y[p0];
  const RealT x1 = x[p1], y1 = y[p1];
  const RealT x2 = x[p2], y2 = y[p2];

  const RealT dx1 = x1 - x0, dy1 = y1 - y0;
  const RealT dx2 = x2 - x0, dy2 = y2 - y0;

  auto cross = [](const RealT ux, const RealT uy, const RealT vx, const RealT vy) { return ux * vy - uy * vx; };

  const RealT det = cross(dx1, dy1, dx2, dy2);
  if (std::abs(SU2_TYPE::GetValue(det)) < 1e-30) {
    return std::array<RealT, 3>{RealT{0}, RealT{0}, RealT{0}};
  }

  const RealT inv_det = 1.0 / det;
  const RealT a = (cross(x_q, y_q, dx2, dy2) - cross(x0, y0, dx2, dy2)) * inv_det;
  const RealT b = (cross(x0, y0, dx1, dy1) - cross(x_q, y_q, dx1, dy1)) * inv_det;

  return std::array<RealT, 3>{1 - a - b, a, b};
}

/*!
 * \brief Checks if a point is inside a triangle based on its barycentric coordinates.
 */
inline bool InTriangle(const std::array<RealT, 3>& coords, const RealT tol = 0.0) {
  return coords[0] >= -tol && coords[1] >= -tol && coords[2] >= -tol;
}

/*!
 * \brief Finds the triangle containing a point using the trapezoidal map.
 */
inline IntT FindTriangle(const TrapezoidalMap& map, const Matrix3i& triangles, const Matrix2i& edge_pts,
                         const Matrix2i& edge_faces, const VectorReal& x, const VectorReal& y, const RealT x_q,
                         const RealT y_q, std::array<RealT, 3>& bary_out) {
  const auto [e_below, e_above] = QueryTrapezoidalMap(map, edge_pts, x, y, x_q, y_q);
  const auto candidates = AdjacentTriangles(e_below, e_above, edge_faces);

  const RealT tol = 1e-12;
  for (const auto t : candidates) {
    if (t < 0) continue;

    const auto coords = TriangleCoords(t, triangles, x, y, x_q, y_q);
    if (InTriangle(coords, tol)) {
      bary_out = coords;
      return t;
    }
  }

  bary_out = {0.0, 0.0, 0.0};
  return -1;
}

}  // namespace su2_lut

/*!
 * \class CTrapezoidalMap
 * \ingroup LookUpInterp
 * \brief Trapezoidal map for finding the triangle containing a query point in a 2D triangulation.
 */
class CTrapezoidalMap {
 private:
  su2_lut::Matrix3i triangles;
  su2_lut::Matrix2i edge_pts, edge_faces;
  su2_lut::VectorReal x_coords, y_coords;
  su2_lut::TrapezoidalMap map;

  unsigned long n_points = 0;
  unsigned long n_triangles = 0;

 public:
  CTrapezoidalMap() = default;

  /*!
   * \brief Build the trapezoidal map from a triangulation.
   * \return True on success.
   */
  bool Build(unsigned long num_points, unsigned long num_triangles, const su2double* x, const su2double* y,
             const unsigned long* connectivity) {
    n_points = num_points;
    n_triangles = num_triangles;

    if (num_points == 0 || num_triangles == 0) return false;

    x_coords.assign(x, x + num_points);
    y_coords.assign(y, y + num_points);

    triangles.resize(num_triangles, 3);
    for (size_t i = 0; i < 3 * num_triangles; ++i) {
      triangles.data[i] = static_cast<su2_lut::IntT>(connectivity[i]);
    }

    su2_lut::ReorderPoints(triangles, x_coords, y_coords);
    su2_lut::ExtractEdges(triangles, edge_pts, edge_faces);
    su2_lut::BuildTrapezoidalMap(edge_pts, x_coords, y_coords, map);

    return !map.x_bands.empty() && !map.offsets.empty();
  }

  /*!
   * \brief Find the triangle containing a query point.
   * \return True if the point is inside the triangulation.
   */
  bool FindTriangle(su2double val_x, su2double val_y, unsigned long& triangle_id,
                    std::array<su2double, 3>& bary_coords) const {
    if (n_triangles == 0 || n_points == 0 || map.x_bands.empty()) {
      bary_coords = {0.0, 0.0, 0.0};
      return false;
    }

    std::array<su2_lut::RealT, 3> bary;
    const su2_lut::IntT tri_id =
        su2_lut::FindTriangle(map, triangles, edge_pts, edge_faces, x_coords, y_coords, val_x, val_y, bary);

    if (tri_id < 0) return false;

    triangle_id = static_cast<unsigned long>(tri_id);
    bary_coords = {bary[0], bary[1], bary[2]};
    return true;
  }

  /*!
   * \brief Get the memory footprint of the map in MB.
   */
  double GetMemoryFootprint() const {
    const size_t bytes =
        (map.edge_id.size() + map.offsets.size() + edge_pts.data.size() + edge_faces.data.size() +
         triangles.data.size()) *
            sizeof(su2_lut::IntT) +
        (map.edge_y.size() + map.x_bands.size() + x_coords.size() + y_coords.size()) * sizeof(su2_lut::RealT);
    return double(bytes) / (1024.0 * 1024.0);
  }
};
