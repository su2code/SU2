/*!
 * \file CTrapezoidalMapFAST.hpp
 * \brief Fast trapezoidal map implementation for 2D LUT queries in SU2.
 *        Based on Pedro Gomes' (pcarruscag) LUT implementation:
 *        https://github.com/pcarruscag/LUT
 */

#ifndef CTRAPEZOIDALMAP_FAST_HPP
#define CTRAPEZOIDALMAP_FAST_HPP

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <numeric>
#include <tuple>
#include <utility>
#include <vector>

namespace su2_lut {

using IntT = int32_t;
using RealT = su2double;

/*--- Simple Matrix container (row-major) ---*/
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

/*--- Trapezoidal Map structure ---*/
struct TrapezoidalMap {
    VectorInt offsets, edge_id;
    VectorReal x_bands, edge_y;
};

/*!
 * \brief Orders points by ascending x coordinates and updates triangle indices.
 */
inline void ReorderPoints(Matrix3i& triangles, VectorReal& x, VectorReal& y) {
    const IntT n_pts = static_cast<IntT>(x.size());

    // Create sort permutation based on (x, y) coordinates
    std::vector<IntT> perm(n_pts);
    std::iota(perm.begin(), perm.end(), 0);
    std::sort(perm.begin(), perm.end(), [&x, &y](const auto i, const auto j) {
        return x[i] != x[j] ? x[i] < x[j] : y[i] < y[j];
    });

    // Reorder coordinate arrays
    auto reorder = [n_pts, &perm](const auto& v) {
        VectorReal tmp(n_pts);
        for (IntT i = 0; i < n_pts; ++i) {
            tmp[i] = v[perm[i]];
        }
        return tmp;
    };
    x = reorder(x);
    y = reorder(y);

    // Build inverse permutation to update triangle indices
    std::vector<IntT> inv_perm(n_pts);
    for (IntT i = 0; i < n_pts; ++i) {
        inv_perm[perm[i]] = i;
    }

    // Update triangle connectivity
    for (IntT i = 0; i < static_cast<IntT>(triangles.rows()); ++i) {
        for (IntT j = 0; j < 3; ++j) {
            triangles(i, j) = inv_perm[triangles(i, j)];
        }
    }
}

/*!
 * \brief Check if two edges are equal (same point IDs).
 */
inline bool EdgesEqual(const std::array<IntT, 3>& a, const std::array<IntT, 3>& b) {
    return a[0] == b[0] && a[1] == b[1];
}

/*!
 * \brief Extract unique edges from triangles.
 */
inline void ExtractEdges(const Matrix3i& triangles, Matrix2i& edge_pts,
                         Matrix2i& edge_faces) {
    // Extract edges from triangles
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

    // Sort to identify duplicates
    std::sort(edges.begin(), edges.end(), [](const auto& a, const auto& b) {
        return a[0] != b[0] ? (a[0] < b[0]) : (a[1] < b[1]);
    });

    // Count unique edges
    IntT n_edges = 1;
    for (IntT i = 1; i < static_cast<IntT>(edges.size()); ++i) {
        n_edges += static_cast<IntT>(!EdgesEqual(edges[i], edges[i - 1]));
    }

    // Map edge points and edge faces
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
        if (EdgesEqual(edges[i], edges[i - 1])) {
            edge_faces(pos - 1, 1) = edges[i][2];
        } else {
            new_edge(edges[i]);
        }
    }
}

/*!
 * \brief Detects bands based on unique x-coordinates.
 * Uses coarse banding to limit memory for large tables.
 * Returns tuple of (n_bands, pt_to_band, x_bands).
 */
inline auto DetectBands(const VectorReal& x, IntT max_bands = 0) {
    // Safety check: empty input
    if (x.empty()) {
        return std::make_tuple(IntT{0}, std::vector<IntT>{}, VectorReal{});
    }
    
    // First pass: count unique x-coordinates
    IntT n_unique = 1;
    for (IntT i = 1; i < static_cast<IntT>(x.size()); ++i) {
        if (x[i] != x[i - 1]) n_unique++;
    }
    
    // Determine if we need to use coarser bands
    if (max_bands <= 0) {
        max_bands = std::min(IntT{5000}, static_cast<IntT>(4.0 * std::sqrt(static_cast<double>(x.size()))));
    }
    
    const bool use_coarse_bands = (n_unique > max_bands);
    
    if (!use_coarse_bands) {
        // Original algorithm: use all unique x-coordinates as band boundaries
        std::vector<IntT> pt_to_band(x.size());
        pt_to_band[0] = 0;
        IntT n_bands = 0;
        
        for (IntT i = 1; i < static_cast<IntT>(x.size()); ++i) {
            n_bands += static_cast<IntT>(x[i] != x[i - 1]);
            pt_to_band[i] = n_bands;
        }

        VectorReal x_bands(n_bands + 1);
        IntT pos = 0;
        x_bands[pos] = x[0];
        for (IntT i = 1; i < static_cast<IntT>(x.size()); ++i) {
            if (x[i] != x_bands[pos]) {
                x_bands[++pos] = x[i];
            }
        }
        
        return std::make_tuple(n_bands, std::move(pt_to_band), std::move(x_bands));
    } else {
        // Coarse banding: divide x-range into max_bands equal intervals
        const RealT x_min = x.front();
        const RealT x_max = x.back();
        const RealT x_range = x_max - x_min;
        
        // Handle degenerate case where all x are the same
        if (x_range < 1e-30) {
            std::vector<IntT> pt_to_band(x.size(), 0);
            VectorReal x_bands = {x_min, x_max};
            return std::make_tuple(IntT{1}, std::move(pt_to_band), std::move(x_bands));
        }
        
        const IntT n_bands = max_bands;
        const RealT band_width = x_range / max_bands;
        
        // Build band boundaries
        VectorReal x_bands(max_bands + 1);
        for (IntT i = 0; i <= max_bands; ++i) {
            x_bands[i] = x_min + i * band_width;
        }
        x_bands[max_bands] = x_max;
        
        // Map each point to its band
        std::vector<IntT> pt_to_band(x.size());
        for (IntT i = 0; i < static_cast<IntT>(x.size()); ++i) {
            IntT band = static_cast<IntT>((x[i] - x_min) / band_width);
            band = std::max(IntT{0}, std::min(band, max_bands - 1));
            pt_to_band[i] = band;
        }
        
        return std::make_tuple(n_bands, std::move(pt_to_band), std::move(x_bands));
    }
}

/*!
 * \brief Builds the trapezoidal map for efficient 2D queries.
 */
inline void BuildTrapezoidalMap(const Matrix2i& edge_pts, const VectorReal& x,
                                const VectorReal& y, TrapezoidalMap& map) {
    auto& x_bands = map.x_bands;
    auto& offsets = map.offsets;
    auto& edge_id = map.edge_id;
    auto& edge_y = map.edge_y;

    const auto [n_bands, pt_to_band, bands] = DetectBands(x);
    x_bands = std::move(bands);
    
    // Safety check
    if (n_bands <= 0 || x_bands.size() < 2) {
        x_bands.clear();
        offsets.clear();
        edge_id.clear();
        edge_y.clear();
        return;
    }

    // Helper function to find band index for an x-coordinate
    auto find_band = [&x_bands, n_bands](RealT x_val) -> IntT {
        auto it = std::lower_bound(x_bands.begin(), x_bands.end(), x_val);
        IntT idx = static_cast<IntT>(it - x_bands.begin());
        idx = std::max(IntT{0}, idx - 1);
        idx = std::min(idx, n_bands - 1);
        return idx;
    };

    // Count edges per band
    auto& counts = offsets;
    counts.clear();
    counts.resize(n_bands + 1, 0);
    
    for (IntT i = 0; i < static_cast<IntT>(edge_pts.rows()); ++i) {
        const IntT pt_0 = edge_pts(i, 0);
        const IntT pt_1 = edge_pts(i, 1);
        const RealT x_0 = x[pt_0];
        const RealT x_1 = x[pt_1];
        
        const IntT band_0 = find_band(x_0);
        const IntT band_1 = find_band(x_1);
        
        const IntT min_band = std::min(band_0, band_1);
        const IntT max_band = std::max(band_0, band_1);
        
        for (IntT j = min_band; j <= max_band && j < n_bands; ++j) {
            ++counts[j + 1];
        }
    }

    // Convert counts to offsets (CSR format)
    for (IntT i = 2; i < static_cast<IntT>(offsets.size()); ++i) {
        offsets[i] += offsets[i - 1];
    }

    // Check memory requirement
    const size_t total_entries = static_cast<size_t>(offsets.back());
    const size_t memory_mb = total_entries * (sizeof(IntT) + sizeof(RealT)) / (1024 * 1024);
    
    if (memory_mb > 2048) {
        x_bands.clear();
        offsets.clear();
        edge_id.clear();
        edge_y.clear();
        return;
    }

    // Allocate storage
    edge_id.resize(offsets.back());
    edge_y.resize(offsets.back());
    auto pos = offsets;

    // Store edges in each band
    for (IntT i_edge = 0; i_edge < static_cast<IntT>(edge_pts.rows()); ++i_edge) {
        const IntT pt_0 = edge_pts(i_edge, 0);
        const IntT pt_1 = edge_pts(i_edge, 1);
        const RealT x_0 = x[pt_0], y_0 = y[pt_0];
        const RealT x_1 = x[pt_1], y_1 = y[pt_1];
        
        const IntT band_0 = find_band(x_0);
        const IntT band_1 = find_band(x_1);
        
        if (band_0 == band_1) continue;
        
        const IntT min_band = std::min(band_0, band_1);
        const IntT max_band = std::max(band_0, band_1);
        const RealT dy_dx = (y_1 - y_0) / (x_1 - x_0);
        
        for (IntT j = min_band; j < max_band && j < n_bands; ++j) {
            if (pos[j] < static_cast<IntT>(edge_id.size())) {
                edge_id[pos[j]] = i_edge;
                const RealT x_mid = (x_bands[j] + x_bands[j + 1]) / 2;
                edge_y[pos[j]] = y_0 + dy_dx * (x_mid - x_0);
                ++pos[j];
            }
        }
    }

    // Sort edges in each band by y coordinate
    std::vector<std::pair<IntT, RealT>> tmp;
    for (IntT i = 0; i < n_bands; ++i) {
        const IntT begin = offsets[i];
        const IntT end = offsets[i + 1];
        if (begin >= end) continue;
        
        tmp.resize(end - begin);
        for (auto k = begin; k < end; ++k) {
            tmp[k - begin] = {edge_id[k], edge_y[k]};
        }
        std::sort(tmp.begin(), tmp.end(), [](const auto& a, const auto& b) {
            return a.second < b.second;
        });
        for (auto k = begin; k < end; ++k) {
            edge_id[k] = tmp[k - begin].first;
            edge_y[k] = tmp[k - begin].second;
        }
    }
}

/*!
 * \brief Query the trapezoidal map for edges bounding a point.
 * Returns pair of (edge_below, edge_above). Either can be -1 if at boundary.
 */
inline auto QueryTrapezoidalMap(const TrapezoidalMap& map, const Matrix2i& edge_pts,
                                const VectorReal& x_coords, const VectorReal& y_coords,
                                const RealT& x, const RealT& y) {
    // Safety check
    if (map.x_bands.size() < 2 || map.offsets.empty()) {
        return std::make_pair(IntT{-1}, IntT{-1});
    }
    
    const IntT n_bands = static_cast<IntT>(map.x_bands.size()) - 1;
    
    // Find the band containing x
    const auto& x_bands = map.x_bands;
    auto it = std::lower_bound(x_bands.begin(), x_bands.end(), x);
    IntT d = static_cast<IntT>(it - x_bands.begin());
    IntT center_band = std::min(std::max(IntT{0}, d - 1), n_bands - 1);
    
    // Search in current band AND adjacent bands
    RealT best_y_below = -1e300;
    RealT best_y_above = 1e300;
    IntT edge_below = -1;
    IntT edge_above = -1;
    
    for (IntT band_offset = -1; band_offset <= 1; ++band_offset) {
        IntT band_idx = center_band + band_offset;
        if (band_idx < 0 || band_idx >= n_bands) continue;
        if (band_idx + 1 >= static_cast<IntT>(map.offsets.size())) continue;
        
        const IntT begin = map.offsets[band_idx];
        const IntT end = map.offsets[band_idx + 1];
        
        if (begin < 0 || end > static_cast<IntT>(map.edge_id.size()) || begin >= end) {
            continue;
        }
        
        for (IntT k = begin; k < end; ++k) {
            const IntT e_id = map.edge_id[k];
            
            const IntT p0 = edge_pts(e_id, 0);
            const IntT p1 = edge_pts(e_id, 1);
            const RealT x0 = x_coords[p0], y0 = y_coords[p0];
            const RealT x1 = x_coords[p1], y1 = y_coords[p1];
            
            // Check if query x is within edge's x-extent
            if (x < std::min(x0, x1) - 1e-10 || x > std::max(x0, x1) + 1e-10) {
                continue;
            }
            
            // Compute y-position of edge at query x
            RealT edge_y_at_x;
            const RealT dx = x1 - x0;
            if (std::abs(dx) < 1e-30) {
                edge_y_at_x = (y0 + y1) / 2.0;
            } else {
                const RealT t = (x - x0) / dx;
                edge_y_at_x = y0 + t * (y1 - y0);
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
    }
    
    return std::make_pair(edge_below, edge_above);
}

/*!
 * \brief Get triangles adjacent to two query edges (up to 3 triangles).
 */
inline auto AdjacentTriangles(const IntT edge_0, const IntT edge_1,
                              const Matrix2i& edge_faces) {
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
 * \brief Compute barycentric coordinates of point (x, y) in triangle.
 */
inline auto TriangleCoords(const IntT i_tri, const Matrix3i& triangles,
                           const VectorReal& x, const VectorReal& y,
                           const RealT x_q, const RealT y_q) {
    const IntT p0 = triangles(i_tri, 0);
    const IntT p1 = triangles(i_tri, 1);
    const IntT p2 = triangles(i_tri, 2);
    
    const RealT x0 = x[p0], y0 = y[p0];
    const RealT x1 = x[p1], y1 = y[p1];
    const RealT x2 = x[p2], y2 = y[p2];
    
    const RealT dx1 = x1 - x0, dy1 = y1 - y0;
    const RealT dx2 = x2 - x0, dy2 = y2 - y0;
    
    auto cross = [](const RealT ux, const RealT uy, const RealT vx, const RealT vy) {
        return ux * vy - uy * vx;
    };
    
    const RealT det = cross(dx1, dy1, dx2, dy2);
    if (std::abs(det) < 1e-30) {
        return std::array{RealT{0}, RealT{0}, RealT{0}};
    }
    
    const RealT inv_det = 1.0 / det;
    const RealT a = (cross(x_q, y_q, dx2, dy2) - cross(x0, y0, dx2, dy2)) * inv_det;
    const RealT b = (cross(x0, y0, dx1, dy1) - cross(x_q, y_q, dx1, dy1)) * inv_det;
    
    return std::array{1 - a - b, a, b};
}

/*!
 * \brief Check if point is inside triangle.
 */
inline bool InTriangle(const std::array<RealT, 3>& coords, const RealT tol = 0.0) {
    return coords[0] >= -tol && coords[1] >= -tol && coords[2] >= -tol;
}

/*!
 * \brief Find triangle containing point using the trapezoidal map.
 */
inline IntT FindTriangle(const TrapezoidalMap& map, const Matrix3i& triangles,
                         const Matrix2i& edge_pts, const Matrix2i& edge_faces, 
                         const VectorReal& x, const VectorReal& y,
                         const RealT x_q, const RealT y_q, std::array<RealT, 3>& bary_out) {
    
    const auto [e_below, e_above] = QueryTrapezoidalMap(map, edge_pts, x, y, x_q, y_q);
    const auto candidates = AdjacentTriangles(e_below, e_above, edge_faces);
    
    constexpr RealT tol = 1e-12;
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
 * \class CTrapezoidalMapFAST
 * \brief Wrapper class for SU2 integration.
 */
class CTrapezoidalMapFAST {
private:
    su2_lut::Matrix3i triangles;
    su2_lut::Matrix2i edge_pts, edge_faces;
    su2_lut::VectorReal x_coords, y_coords;
    su2_lut::TrapezoidalMap map;
    
    unsigned long n_points;
    unsigned long n_triangles;
    bool build_successful;

public:
    CTrapezoidalMapFAST() : n_points(0), n_triangles(0), build_successful(false) {}

    /*!
     * \brief Build the trapezoidal map from mesh data.
     */
    bool Build(unsigned long num_points, unsigned long num_triangles,
               const su2double* x, const su2double* y,
               const unsigned long* connectivity) {
        
        n_points = num_points;
        n_triangles = num_triangles;
        build_successful = false;

        if (num_points == 0 || num_triangles == 0) {
            x_coords.clear();
            y_coords.clear();
            triangles.resize(0, 3);
            edge_pts.resize(0, 2);
            edge_faces.resize(0, 2);
            map.x_bands.clear();
            map.offsets.clear();
            map.edge_id.clear();
            map.edge_y.clear();
            return false;
        }

        // Copy coordinates
        x_coords.assign(x, x + num_points);
        y_coords.assign(y, y + num_points);

        // Copy triangle connectivity
        triangles.resize(num_triangles, 3);
        for (size_t i = 0; i < num_triangles; ++i) {
            unsigned long v0 = connectivity[3*i + 0];
            unsigned long v1 = connectivity[3*i + 1];
            unsigned long v2 = connectivity[3*i + 2];
            
            if (v0 >= num_points || v1 >= num_points || v2 >= num_points) {
                triangles(i, 0) = 0;
                triangles(i, 1) = 0;
                triangles(i, 2) = 0;
            } else {
                triangles(i, 0) = static_cast<su2_lut::IntT>(v0);
                triangles(i, 1) = static_cast<su2_lut::IntT>(v1);
                triangles(i, 2) = static_cast<su2_lut::IntT>(v2);
            }
        }

        // Build the map
        su2_lut::ReorderPoints(triangles, x_coords, y_coords);
        su2_lut::ExtractEdges(triangles, edge_pts, edge_faces);
        su2_lut::BuildTrapezoidalMap(edge_pts, x_coords, y_coords, map);
        
        if (map.x_bands.empty() || map.offsets.empty()) {
            build_successful = false;
            return false;
        }
        
        build_successful = true;
        return true;
    }
    
    bool IsBuilt() const { return build_successful; }

    /*!
     * \brief Find the triangle containing a query point.
     */
    bool FindTriangle(su2double val_x, su2double val_y,
                      unsigned long& triangle_id,
                      std::array<su2double, 3>& bary_coords) const {
        
        if (n_triangles == 0 || n_points == 0 || map.x_bands.empty()) {
            bary_coords = {0.0, 0.0, 0.0};
            return false;
        }
        
        std::array<su2_lut::RealT, 3> bary;
        su2_lut::IntT tri_id = su2_lut::FindTriangle(
            map, triangles, edge_pts, edge_faces, x_coords, y_coords, val_x, val_y, bary);
        
        if (tri_id >= 0) {
            triangle_id = static_cast<unsigned long>(tri_id);
            bary_coords = {bary[0], bary[1], bary[2]};
            return true;
        }
        return false;
    }

    unsigned long GetNTriangles() const { return n_triangles; }
    unsigned long GetNPoints() const { return n_points; }
    unsigned long GetNEdges() const { return edge_pts.rows(); }
    unsigned long GetNBands() const { return map.x_bands.size() > 0 ? map.x_bands.size() - 1 : 0; }
};

#endif  // CTRAPEZOIDALMAP_FAST_HPP
