/*!
 * \file sharp_edge_normals.hpp
 * \brief Wall normals at sharp edges inside boundary markers (option WALL_SHARP_EDGE_NORMALS): grouping of the dual
 *        area pieces of a boundary point and orthonormal bases of the resulting normals.
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

#pragma once

#include <algorithm>
#include <cmath>
#include <utility>
#include <vector>
#include "../basic_types/datatype_structure.hpp"

/*--- The vertex normal of a boundary point is the sum of the dual area pieces ("corner vectors" N_f) of its incident
 * boundary faces. Where the wall folds sharply at the point (e.g. a sharp trailing edge, whose upper and lower faces have
 * nearly opposite normals) these vectors almost cancel and their sum is a small residue whose direction is set by the
 * local face sizes, not by the geometry: a slip condition along it lets the flow cross the faces. The faces of a point
 * are therefore grouped: two groups are merged across a shared boundary edge only if every pair of faces of the merged
 * group has N_f . N_g > 0. Inside a group the summed normal n satisfies n . N_f >= |N_f|^2 > 0 for every face of the
 * group, whatever the face areas; across a fold of 90 degrees or more this can fail. A point with one group keeps its
 * vertex normal; a point with two or more groups is a sharp-edge point and uses the group normals. ---*/

namespace SharpEdgeNormals {

/*! \brief Largest number of incident faces handled at a point (more: the vertex normal is kept). */
constexpr unsigned short MAX_FACES = 32;

/*! \brief SU2's parallel test of ComputeModifiedSymmetryNormals, 1 - |cos| < 0.001, as a squared sine: 2 t - t^2. */
constexpr passivedouble PARALLEL_SIN2 = 2.0 * 0.001 - 0.001 * 0.001;

/*!
 * \brief Groups of the nFaces incident faces of a point.
 * \param[in] nDim - Number of dimensions.
 * \param[in] nFaces - Number of faces (at most MAX_FACES).
 * \param[in] N - Corner area vectors of the faces at the point, consistently oriented.
 * \param[in] edge - For each face the nodes joined to the point by an edge of the face (2 in 3D; unused in 2D, where
 *            the two segments of a point are adjacent).
 * \param[in] nEdge - Number of entries of edge per face.
 * \param[out] groupN - Summed corner vector of every group, groups in the order of their first face.
 * \param[out] faceGroup - Group of every face (optional).
 * \return Number of groups.
 */
inline unsigned short FaceGroups(unsigned short nDim, unsigned short nFaces, const su2double (*N)[3],
                                 const unsigned long (*edge)[2], const unsigned short* nEdge, su2double (*groupN)[3],
                                 unsigned short* faceGroup = nullptr) {
  unsigned short group[MAX_FACES];
  for (unsigned short f = 0; f < nFaces; ++f) group[f] = f;
  auto dot = [&](unsigned short f, unsigned short g) {
    su2double d = 0.0;
    for (unsigned short k = 0; k < nDim; ++k) d += N[f][k] * N[g][k];
    return d;
  };
  for (unsigned short f = 0; f < nFaces; ++f) {
    for (unsigned short g = f + 1; g < nFaces; ++g) {
      if (group[f] == group[g]) continue;
      bool adjacent = (nDim == 2);
      for (unsigned short a = 0; a < nEdge[f] && !adjacent; ++a)
        for (unsigned short b = 0; b < nEdge[g]; ++b) adjacent = adjacent || (edge[f][a] == edge[g][b]);
      if (!adjacent) continue;
      bool positive = true;
      for (unsigned short p = 0; p < nFaces && positive; ++p) {
        if (group[p] != group[f]) continue;
        for (unsigned short q = 0; q < nFaces && positive; ++q)
          if (group[q] == group[g]) positive = (dot(p, q) > 0.0);
      }
      if (!positive) continue;
      const auto keep = std::min(group[f], group[g]), drop = std::max(group[f], group[g]);
      for (unsigned short p = 0; p < nFaces; ++p)
        if (group[p] == drop) group[p] = keep;
    }
  }
  unsigned short nGroups = 0, slotOf[MAX_FACES];
  for (unsigned short f = 0; f < nFaces; ++f) {
    if (group[f] == f) {
      slotOf[f] = nGroups;
      for (unsigned short d = 0; d < 3; ++d) groupN[nGroups][d] = 0.0;
      ++nGroups;
    }
    for (unsigned short d = 0; d < nDim; ++d) groupN[slotOf[group[f]]][d] += N[f][d];
    if (faceGroup) faceGroup[f] = slotOf[group[f]];
  }
  return nGroups;
}

/*!
 * \brief Orthonormal basis Q (rows) of the span of nNormals normals: Gram-Schmidt in the given order, a normal dropped when
 *        its orthogonalized squared norm (relative to its own) is below PARALLEL_SIN2.
 * \return Rank (number of rows of Q, at most nDim).
 */
inline unsigned short Basis(unsigned short nDim, unsigned short nNormals, const su2double (*normals)[3],
                            su2double (*Q)[3]) {
  unsigned short rank = 0;
  for (unsigned short a = 0; a < nNormals && rank < nDim; ++a) {
    su2double v[3] = {0.0, 0.0, 0.0}, norm2 = 0.0;
    for (unsigned short d = 0; d < nDim; ++d) {
      v[d] = normals[a][d];
      norm2 += v[d] * v[d];
    }
    if (!(norm2 > 0.0)) continue;
    for (unsigned short r = 0; r < rank; ++r) {
      su2double proj = 0.0;
      for (unsigned short d = 0; d < nDim; ++d) proj += Q[r][d] * v[d];
      for (unsigned short d = 0; d < nDim; ++d) v[d] -= proj * Q[r][d];
    }
    su2double rest2 = 0.0;
    for (unsigned short d = 0; d < nDim; ++d) rest2 += v[d] * v[d];
    if (rest2 < PARALLEL_SIN2 * norm2) continue;
    const su2double inv = 1.0 / sqrt(rest2);
    for (unsigned short d = 0; d < nDim; ++d) Q[rank][d] = v[d] * inv;
    ++rank;
  }
  return rank;
}

/*!
 * \brief vec[iVel + d] <- ((I - Q^T Q) vec)[iVel + d], the components of vec[iVel..iVel+nDim-1] along the orthonormal rows
 *        of Q removed.
 */
template <class VecType>
inline void Project(unsigned short nDim, unsigned short rank, const su2double (*Q)[3], unsigned short iVel,
                    VecType& vec) {
  for (unsigned short r = 0; r < rank; ++r) {
    su2double proj = 0.0;
    for (unsigned short d = 0; d < nDim; ++d) proj += Q[r][d] * vec[iVel + d];
    for (unsigned short d = 0; d < nDim; ++d) vec[iVel + d] -= proj * Q[r][d];
  }
}

/*! \brief Per face group of a no-slip wall vertex: wetted area A_q, normal neighbour and its normal distance d_q. */
using Group = std::pair<su2double, std::pair<unsigned long, su2double>>;

/*!
 * \brief Isothermal wall flux of a sharp-edge vertex: flux = sum_q k (Twall - T_q) / d_q A_q (the conducted heat, same sign
 *        convention as k dT/dn Area of the isothermal wall), areaOverDist = sum_q A_q / d_q (for the Jacobian).
 * \param[in] T - Callable returning the temperature of a point.
 * \return false if a group has no eligible neighbour (d_q <= 0): the caller keeps its single-normal stencil.
 */
template <class TemperatureFn>
inline bool IsothermalFlux(const std::vector<Group>& groups, su2double k, su2double Twall, const TemperatureFn& T,
                           su2double& flux, su2double& areaOverDist) {
  flux = 0.0;
  areaOverDist = 0.0;
  for (const auto& g : groups)
    if (!(g.second.second > 0.0)) return false;
  for (const auto& g : groups) {
    flux += k * (Twall - T(g.second.first)) / g.second.second * g.first;
    areaOverDist += g.first / g.second.second;
  }
  return true;
}

}  // namespace SharpEdgeNormals
