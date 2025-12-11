/*!
 * \file CTrapezoidalMapFAST.hpp
 * \brief Fast trapezoidal map implementation for 2D LUT queries in SU2.
 *        Based on Pedro Gomes' (pcarruscag) LUT implementation:
 *        https://github.com/pcarruscag/LUT
 *
 *
 */

#ifndef CTRAPEZOIDALMAP_FAST_HPP
#define CTRAPEZOIDALMAP_FAST_HPP

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <numeric>
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

/*--- Comparison functors ---*/
struct PointComparator {
    const VectorReal& x;
    const VectorReal& y;
    
    PointComparator(const VectorReal& x_ref, const VectorReal& y_ref) 
        : x(x_ref), y(y_ref) {}
    
    bool operator()(const IntT i, const IntT j) const {
        return x[i] != x[j] ? x[i] < x[j] : y[i] < y[j];
    }
};

struct EdgeComparator {
    bool operator()(const std::array<IntT, 3>& a, const std::array<IntT, 3>& b) const {
        return a[0] != b[0] ? (a[0] < b[0]) : (a[1] < b[1]);
    }
};

struct PairSecondComparator {
    bool operator()(const std::pair<IntT, RealT>& a, const std::pair<IntT, RealT>& b) const {
        return a.second < b.second;
    }
};

/*!
 * \brief Orders points by ascending x coordinates and updates triangle indices.
 */
inline void ReorderPoints(Matrix3i& triangles, VectorReal& x, VectorReal& y) {
    const IntT n_pts = static_cast<IntT>(x.size());

    // Create sort permutation based on (x, y) coordinates
    std::vector<IntT> perm(n_pts);
    for (IntT i = 0; i < n_pts; ++i) {
        perm[i] = i;
    }
    std::sort(perm.begin(), perm.end(), PointComparator(x, y));

    // Reorder coordinate arrays
    VectorReal x_tmp(n_pts), y_tmp(n_pts);
    for (IntT i = 0; i < n_pts; ++i) {
        x_tmp[i] = x[perm[i]];
        y_tmp[i] = y[perm[i]];
    }
    x = x_tmp;
    y = y_tmp;

    // Build inverse permutation to update triangle indices
    std::vector<IntT> inv_perm(n_pts);
    for (IntT i = 0; i < n_pts; ++i) {
        inv_perm[perm[i]] = i;
    }

    // Update all triangle vertex indices to use new point ordering
    for (IntT i = 0; i < static_cast<IntT>(triangles.rows()); ++i) {
        for (IntT j = 0; j < 3; ++j) {
            triangles(i, j) = inv_perm[triangles(i, j)];
        }
    }
}

/*!
 * \brief Helper to check if two edges are equal (same endpoints).
 */
inline bool EdgesEqual(const std::array<IntT, 3>& a, const std::array<IntT, 3>& b) {
    return a[0] == b[0] && a[1] == b[1];
}

/*!
 * \brief Extracts unique edges from triangles with face adjacency information.
 * Each edge stores which triangles (up to 2) share it.
 * Boundary edges have edge_faces(i,1) = -1.
 */
inline void ExtractEdges(const Matrix3i& triangles, Matrix2i& edge_pts,
                         Matrix2i& edge_faces) {
    // Extract all edges (3 per triangle), storing (min_pt, max_pt, triangle_id)
    std::vector<std::array<IntT, 3> > edges;
    edges.resize(3 * triangles.rows());
    
    for (IntT i_tri = 0; i_tri < static_cast<IntT>(triangles.rows()); ++i_tri) {
        for (IntT i = 0; i < 3; ++i) {
            const IntT j = (i + 1) % 3;
            const IntT i_pt = std::min(triangles(i_tri, i), triangles(i_tri, j));
            const IntT j_pt = std::max(triangles(i_tri, i), triangles(i_tri, j));
            std::array<IntT, 3> edge_data = {{i_pt, j_pt, i_tri}};
            edges[3 * i_tri + i] = edge_data;
        }
    }

    // Sort to identify duplicate edges (shared between triangles)
    std::sort(edges.begin(), edges.end(), EdgeComparator());

    // Count unique edges
    IntT n_edges = 1;
    for (IntT i = 1; i < static_cast<IntT>(edges.size()); ++i) {
        n_edges += static_cast<IntT>(!EdgesEqual(edges[i], edges[i - 1]));
    }

    // Build edge_pts and edge_faces arrays
    edge_pts.resize(n_edges, 2);
    edge_faces.resize(n_edges, 2);
    
    IntT pos = 0;
    
    // First edge
    edge_pts(pos, 0) = edges[0][0];
    edge_pts(pos, 1) = edges[0][1];
    edge_faces(pos, 0) = edges[0][2];
    edge_faces(pos, 1) = -1;
    ++pos;
    
    for (IntT i = 1; i < static_cast<IntT>(edges.size()); ++i) {
        if (EdgesEqual(edges[i], edges[i - 1])) {
            // Duplicate edge - record second adjacent triangle
            edge_faces(pos - 1, 1) = edges[i][2];
        } else {
            // New edge
            edge_pts(pos, 0) = edges[i][0];
            edge_pts(pos, 1) = edges[i][1];
            edge_faces(pos, 0) = edges[i][2];
            edge_faces(pos, 1) = -1;
            ++pos;
        }
    }
}

/*!
 * \brief Band detection result structure.
 */
struct BandResult {
    IntT n_bands;
    std::vector<IntT> pt_to_band;
    VectorReal x_bands;
};

/*!
 * \brief Detects bands based on unique x-coordinates.
 * Uses coarse banding to limit memory for large tables.
 */
inline BandResult DetectBands(const VectorReal& x, IntT max_bands = 0) {
    BandResult result;
    
    // Safety check: empty input
    if (x.empty()) {
        result.n_bands = 0;
        return result;
    }
    
    // First pass: count unique x-coordinates
    IntT n_unique = 1;
    for (IntT i = 1; i < static_cast<IntT>(x.size()); ++i) {
        if (x[i] != x[i - 1]) n_unique++;
    }
    
    // Determine if we need to use coarser bands
    // Default max_bands: limit to sqrt(n_points) * 4, capped at 5000
    if (max_bands <= 0) {
        max_bands = std::min(IntT(5000), static_cast<IntT>(4.0 * std::sqrt(static_cast<double>(x.size()))));
    }
    
    const bool use_coarse_bands = (n_unique > max_bands);
    
    if (!use_coarse_bands) {
        // Original algorithm: use all unique x-coordinates as band boundaries
        result.pt_to_band.resize(x.size());
        result.pt_to_band[0] = 0;
        result.n_bands = 0;
        
        for (IntT i = 1; i < static_cast<IntT>(x.size()); ++i) {
            result.n_bands += static_cast<IntT>(x[i] != x[i - 1]);
            result.pt_to_band[i] = result.n_bands;
        }

        result.x_bands.resize(result.n_bands + 1);
        IntT pos = 0;
        result.x_bands[pos] = x[0];
        for (IntT i = 1; i < static_cast<IntT>(x.size()); ++i) {
            if (x[i] != result.x_bands[pos]) {
                result.x_bands[++pos] = x[i];
            }
        }
    } else {
        // Coarse banding: divide x-range into max_bands equal intervals
        const RealT x_min = x.front();  // x is sorted
        const RealT x_max = x.back();
        const RealT x_range = x_max - x_min;
        
        // Handle degenerate case where all x are the same
        if (x_range < 1e-30) {
            result.n_bands = 1;
            result.pt_to_band.resize(x.size(), 0);
            result.x_bands.resize(2);
            result.x_bands[0] = x_min;
            result.x_bands[1] = x_max;
            return result;
        }
        
        result.n_bands = max_bands;
        const RealT band_width = x_range / max_bands;
        
        // Build band boundaries
        result.x_bands.resize(max_bands + 1);
        for (IntT i = 0; i <= max_bands; ++i) {
            result.x_bands[i] = x_min + i * band_width;
        }
        // Ensure last band includes x_max exactly
        result.x_bands[max_bands] = x_max;
        
        // Map each point to its band
        result.pt_to_band.resize(x.size());
        for (IntT i = 0; i < static_cast<IntT>(x.size()); ++i) {
            // Calculate band index: floor((x - x_min) / band_width)
            IntT band = static_cast<IntT>((x[i] - x_min) / band_width);
            // Clamp to valid range [0, n_bands-1]
            band = std::max(IntT(0), std::min(band, max_bands - 1));
            result.pt_to_band[i] = band;
        }
    }

    return result;
}

/*!
 * \brief Builds the trapezoidal map for efficient 2D queries.
 * For each band, stores edges crossing it, sorted by y-coordinate at band midpoint.
 */
inline void BuildTrapezoidalMap(const Matrix2i& edge_pts, const VectorReal& x,
                                const VectorReal& y, TrapezoidalMap& map) {
    VectorReal& x_bands = map.x_bands;
    VectorInt& offsets = map.offsets;
    VectorInt& edge_id = map.edge_id;
    VectorReal& edge_y = map.edge_y;

    BandResult band_result = DetectBands(x);
    const IntT n_bands = band_result.n_bands;
    x_bands = band_result.x_bands;
    
    // Safety check
    if (n_bands <= 0 || x_bands.size() < 2) {
        x_bands.clear();
        offsets.clear();
        edge_id.clear();
        edge_y.clear();
        return;
    }

    // Helper function to find band index for an x-coordinate
    // Uses binary search on x_bands
    auto find_band = [&x_bands, n_bands](RealT x_val) -> IntT {
        // Binary search for the band containing x_val
        VectorReal::const_iterator it = std::lower_bound(x_bands.begin(), x_bands.end(), x_val);
        IntT idx = static_cast<IntT>(it - x_bands.begin());
        // lower_bound returns iterator to first element >= x_val
        // We want the band to the left, so subtract 1 (but not below 0)
        idx = std::max(IntT(0), idx - 1);
        // Also cap at n_bands - 1
        idx = std::min(idx, n_bands - 1);
        return idx;
    };

    // Count edges per band - use actual x-coordinates to determine band range
    VectorInt& counts = offsets;
    counts.clear();
    counts.resize(n_bands + 1, 0);
    
    for (IntT i = 0; i < static_cast<IntT>(edge_pts.rows()); ++i) {
        const IntT pt_0 = edge_pts(i, 0);
        const IntT pt_1 = edge_pts(i, 1);
        const RealT x_0 = x[pt_0];
        const RealT x_1 = x[pt_1];
        
        // Find bands for each endpoint using actual x-coordinates
        const IntT band_0 = find_band(x_0);
        const IntT band_1 = find_band(x_1);
        
        // Edge spans from min_band to max_band
        const IntT min_band = std::min(band_0, band_1);
        const IntT max_band = std::max(band_0, band_1);
        
        // Count edge in each band it crosses
        for (IntT j = min_band; j <= max_band && j < n_bands; ++j) {
            ++counts[j + 1];
        }
    }

    // Convert counts to offsets (CSR format)
    for (IntT i = 2; i < static_cast<IntT>(offsets.size()); ++i) {
        offsets[i] += offsets[i - 1];
    }

    // Check total memory requirement before allocation
    const size_t total_entries = static_cast<size_t>(offsets.back());
    const size_t memory_bytes = total_entries * (sizeof(IntT) + sizeof(RealT));
    const size_t memory_mb = memory_bytes / (1024 * 1024);
    
    // Memory limit: 2GB for edge storage (adjustable)
    const size_t MAX_MEMORY_MB = 2048;
    if (memory_mb > MAX_MEMORY_MB) {
        // Too much memory - clear and return empty map
        x_bands.clear();
        offsets.clear();
        edge_id.clear();
        edge_y.clear();
        return;
    }

    // Fill edge data for each band
    edge_id.resize(total_entries);
    edge_y.resize(total_entries);
    VectorInt pos = offsets;  // Working copy of offsets

    for (IntT i_edge = 0; i_edge < static_cast<IntT>(edge_pts.rows()); ++i_edge) {
        const IntT pt_0 = edge_pts(i_edge, 0);
        const IntT pt_1 = edge_pts(i_edge, 1);
        const RealT x_0 = x[pt_0];
        const RealT x_1 = x[pt_1];
        const RealT y_0 = y[pt_0];
        const RealT y_1 = y[pt_1];
        
        // Find bands for each endpoint
        const IntT band_0 = find_band(x_0);
        const IntT band_1 = find_band(x_1);
        const IntT min_band = std::min(band_0, band_1);
        const IntT max_band = std::max(band_0, band_1);
        
        // Calculate slope (handle vertical edges)
        const RealT dx = x_1 - x_0;
        const bool is_vertical = (std::abs(dx) < 1e-30);
        const RealT dy_dx = is_vertical ? 0.0 : (y_1 - y_0) / dx;

        // Add edge to each band it crosses
        for (IntT j = min_band; j <= max_band && j < n_bands; ++j) {
            edge_id[pos[j]] = i_edge;
            const RealT x_mid = (x_bands[j] + x_bands[j + 1]) / 2.0;
            // Calculate y at band midpoint
            if (is_vertical) {
                edge_y[pos[j]] = (y_0 + y_1) / 2.0;  // Use average y for vertical edges
            } else {
                edge_y[pos[j]] = y_0 + dy_dx * (x_mid - x_0);
            }
            ++pos[j];
        }
    }

    // Sort edges within each band by y-coordinate
    std::vector<std::pair<IntT, RealT> > tmp;
    for (IntT i = 0; i < n_bands; ++i) {
        const IntT begin = offsets[i];
        const IntT end = offsets[i + 1];
        tmp.resize(end - begin);
        
        for (IntT k = begin; k < end; ++k) {
            tmp[k - begin] = std::make_pair(edge_id[k], edge_y[k]);
        }
        std::sort(tmp.begin(), tmp.end(), PairSecondComparator());
        
        for (IntT k = begin; k < end; ++k) {
            edge_id[k] = tmp[k - begin].first;
            edge_y[k] = tmp[k - begin].second;
        }
    }
}

/*!
 * \brief Query result structure.
 */
struct QueryResult {
    IntT edge_below;
    IntT edge_above;
};

/*!
 * \brief Query the trapezoidal map for edges bounding a point.
 * Returns edges below and above the point. Either can be -1 if point is at boundary.
 * With coarse banding, recomputes edge y-positions at the actual query x.
 * Also searches adjacent bands to avoid missing candidates.
 */
inline QueryResult QueryTrapezoidalMap(const TrapezoidalMap& map, const Matrix2i& edge_pts,
                                       const VectorReal& x_coords, const VectorReal& y_coords,
                                       const RealT& x, const RealT& y) {
    QueryResult result;
    result.edge_below = -1;
    result.edge_above = -1;
    
    // Safety check: empty map
    if (map.x_bands.size() < 2 || map.offsets.empty()) {
        return result;
    }
    
    const IntT n_bands = static_cast<IntT>(map.x_bands.size()) - 1;
    
    // Find the band containing x
    const VectorReal& x_bands = map.x_bands;
    VectorReal::const_iterator it = std::lower_bound(x_bands.begin(), x_bands.end(), x);
    IntT d = static_cast<IntT>(it - x_bands.begin());
    IntT center_band = std::min(std::max(IntT(0), d - 1), n_bands - 1);
    
    // Search in current band AND adjacent bands to catch edges near boundaries
    RealT best_y_below = -1e300;
    RealT best_y_above = 1e300;
    
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
            
            // Get edge endpoints
            const IntT p0 = edge_pts(e_id, 0);
            const IntT p1 = edge_pts(e_id, 1);
            const RealT x0 = x_coords[p0], y0 = y_coords[p0];
            const RealT x1 = x_coords[p1], y1 = y_coords[p1];
            
            // Check if query x is within the edge's x-extent
            const RealT x_min_edge = std::min(x0, x1);
            const RealT x_max_edge = std::max(x0, x1);
            if (x < x_min_edge - 1e-10 || x > x_max_edge + 1e-10) {
                continue;  // Query x is outside this edge's range
            }
            
            // Compute y-position of edge at query x
            RealT edge_y_at_x;
            const RealT dx = x1 - x0;
            if (std::abs(dx) < 1e-30) {
                // Vertical edge - use average y
                edge_y_at_x = (y0 + y1) / 2.0;
            } else {
                // Linear interpolation
                const RealT t = (x - x0) / dx;
                edge_y_at_x = y0 + t * (y1 - y0);
            }
            
            // Check if this edge is below or above the query point
            if (edge_y_at_x <= y + 1e-10) {
                if (edge_y_at_x > best_y_below) {
                    best_y_below = edge_y_at_x;
                    result.edge_below = e_id;
                }
            }
            if (edge_y_at_x >= y - 1e-10) {
                if (edge_y_at_x < best_y_above) {
                    best_y_above = edge_y_at_x;
                    result.edge_above = e_id;
                }
            }
        }
    }
    
    return result;
}

/*!
 * \brief Get triangles adjacent to two query edges (up to 3 triangles).
 */
inline std::array<IntT, 3> AdjacentTriangles(const IntT edge_0, const IntT edge_1,
                                              const Matrix2i& edge_faces) {
    std::array<IntT, 3> tris = {{-1, -1, -1}};
    IntT pos = 0;

    // Helper to insert unique triangle
    if (edge_0 >= 0) {
        IntT t0 = edge_faces(edge_0, 0);
        IntT t1 = edge_faces(edge_0, 1);
        if (t0 >= 0) {
            tris[pos++] = t0;
        }
        if (t1 >= 0) {
            bool found = false;
            for (IntT k = 0; k < pos; ++k) {
                if (t1 == tris[k]) { found = true; break; }
            }
            if (!found) tris[pos++] = t1;
        }
    }
    
    if (edge_1 >= 0) {
        IntT t0 = edge_faces(edge_1, 0);
        IntT t1 = edge_faces(edge_1, 1);
        if (t0 >= 0) {
            bool found = false;
            for (IntT k = 0; k < pos; ++k) {
                if (t0 == tris[k]) { found = true; break; }
            }
            if (!found && pos < 3) tris[pos++] = t0;
        }
        if (t1 >= 0) {
            bool found = false;
            for (IntT k = 0; k < pos; ++k) {
                if (t1 == tris[k]) { found = true; break; }
            }
            if (!found && pos < 3) tris[pos++] = t1;
        }
    }

    return tris;
}

/*!
 * \brief Compute barycentric coordinates and check if point is inside triangle.
 */
inline std::array<RealT, 3> TriangleCoords(const IntT i_tri, const Matrix3i& triangles,
                                           const VectorReal& x, const VectorReal& y,
                                           const RealT x_q, const RealT y_q) {
    const IntT i0 = triangles(i_tri, 0);
    const IntT i1 = triangles(i_tri, 1);
    const IntT i2 = triangles(i_tri, 2);

    const RealT x0 = x[i0], y0 = y[i0];
    const RealT x1 = x[i1], y1 = y[i1];
    const RealT x2 = x[i2], y2 = y[i2];

    const RealT det = (y1 - y2) * (x0 - x2) + (x2 - x1) * (y0 - y2);
    const RealT l0 = ((y1 - y2) * (x_q - x2) + (x2 - x1) * (y_q - y2)) / det;
    const RealT l1 = ((y2 - y0) * (x_q - x2) + (x0 - x2) * (y_q - y2)) / det;
    const RealT l2 = 1.0 - l0 - l1;

    std::array<RealT, 3> coords = {{l0, l1, l2}};
    return coords;
}

/*!
 * \brief Check if point is inside triangle (barycentric coords all >= 0).
 */
inline bool InTriangle(const std::array<RealT, 3>& coords, const RealT tol = 0.0) {
    return coords[0] >= -tol && coords[1] >= -tol && coords[2] >= -tol;
}

/*!
 * \brief Find triangle containing point (x_q, y_q) using the trapezoidal map.
 * Returns triangle index and barycentric coordinates.
 * Returns -1 if point is outside the mesh.
 */
inline IntT FindTriangle(const TrapezoidalMap& map, const Matrix3i& triangles,
                         const Matrix2i& edge_pts, const Matrix2i& edge_faces, 
                         const VectorReal& x, const VectorReal& y,
                         const RealT x_q, const RealT y_q, std::array<RealT, 3>& bary_out) {
    
    QueryResult qr = QueryTrapezoidalMap(map, edge_pts, x, y, x_q, y_q);
    std::array<IntT, 3> candidates = AdjacentTriangles(qr.edge_below, qr.edge_above, edge_faces);
    
    const RealT tol = 1e-12;
    for (IntT i = 0; i < 3; ++i) {
        if (candidates[i] < 0) continue;
        
        std::array<RealT, 3> coords = TriangleCoords(candidates[i], triangles, x, y, x_q, y_q);
        if (InTriangle(coords, tol)) {
            bary_out = coords;
            return candidates[i];
        }
    }
    
    // Not found - return -1
    bary_out[0] = bary_out[1] = bary_out[2] = 0.0;
    return -1;
}

}  // namespace su2_lut


/*!
 * \class CTrapezoidalMapFAST
 * \brief Wrapper class for SU2 integration of Pedro's trapezoidal map algorithm.
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
     * \brief Build interface matching CLookUpTable.cpp's expected signature.
     * \param[in] num_points - Number of mesh points
     * \param[in] num_triangles - Number of triangles
     * \param[in] x - X coordinates of all points
     * \param[in] y - Y coordinates of all points
     * \param[in] connectivity - Triangle connectivity (flat array: tri0_p0, tri0_p1, tri0_p2, ...)
     * \return true if map was built successfully, false if memory limit exceeded
     */
    bool Build(unsigned long num_points, unsigned long num_triangles,
               const su2double* x, const su2double* y,
               const unsigned long* connectivity) {
        
        n_points = num_points;
        n_triangles = num_triangles;
        build_successful = false;

        // Safety check: handle empty input
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

        // Copy triangle connectivity with bounds checking
        triangles.resize(num_triangles, 3);
        for (size_t i = 0; i < num_triangles; ++i) {
            unsigned long v0 = connectivity[3*i + 0];
            unsigned long v1 = connectivity[3*i + 1];
            unsigned long v2 = connectivity[3*i + 2];
            
            // Bounds check
            if (v0 >= num_points || v1 >= num_points || v2 >= num_points) {
                // Invalid triangle - skip by setting to valid indices
                triangles(i, 0) = 0;
                triangles(i, 1) = 0;
                triangles(i, 2) = 0;
            } else {
                triangles(i, 0) = static_cast<su2_lut::IntT>(v0);
                triangles(i, 1) = static_cast<su2_lut::IntT>(v1);
                triangles(i, 2) = static_cast<su2_lut::IntT>(v2);
            }
        }

        // Build the map using Pedro's algorithm
        su2_lut::ReorderPoints(triangles, x_coords, y_coords);
        su2_lut::ExtractEdges(triangles, edge_pts, edge_faces);
        su2_lut::BuildTrapezoidalMap(edge_pts, x_coords, y_coords, map);
        
        // Check if map was successfully built (not empty due to memory limit)
        if (map.x_bands.empty() || map.offsets.empty()) {
            build_successful = false;
            return false;
        }
        
        build_successful = true;
        return true;
    }
    
    /*!
     * \brief Check if the map was successfully built.
     */
    bool IsBuilt() const { return build_successful; }

    /*!
     * \brief Find the triangle containing a query point.
     * \param[in] val_x - X coordinate of query point
     * \param[in] val_y - Y coordinate of query point
     * \param[out] triangle_id - ID of the containing triangle
     * \param[out] bary_coords - Barycentric coordinates of query point in triangle
     * \return true if point is inside the mesh, false otherwise
     */
    bool FindTriangle(su2double val_x, su2double val_y,
                      unsigned long& triangle_id,
                      std::array<su2double, 3>& bary_coords) const {
        
        // Safety check: empty map
        if (n_triangles == 0 || n_points == 0 || map.x_bands.empty()) {
            bary_coords[0] = bary_coords[1] = bary_coords[2] = 0.0;
            return false;
        }
        
        std::array<su2_lut::RealT, 3> bary;
        su2_lut::IntT tri_id = su2_lut::FindTriangle(
            map, triangles, edge_pts, edge_faces, x_coords, y_coords, val_x, val_y, bary);
        
        if (tri_id >= 0) {
            triangle_id = static_cast<unsigned long>(tri_id);
            bary_coords[0] = bary[0];
            bary_coords[1] = bary[1];
            bary_coords[2] = bary[2];
            return true;
        }
        return false;
    }

    /*--- Accessors ---*/
    unsigned long GetNTriangles() const { return n_triangles; }
    unsigned long GetNPoints() const { return n_points; }
    unsigned long GetNEdges() const { return edge_pts.rows(); }
    unsigned long GetNBands() const { return map.x_bands.size() > 0 ? map.x_bands.size() - 1 : 0; }
};

#endif  // CTRAPEZOIDALMAP_FAST_HPP
