/*!
 * \file sharp_edge_normals.hpp
 * \brief Grouping of the dual area pieces of a boundary point at sharp edges, and orthonormal bases of the resulting
 *        group normals.
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
#include "../basic_types/datatype_structure.hpp"

/*--- The vertex normal of a boundary point is the sum of the dual area pieces ("corner vectors" N_f) of its incident
 * boundary faces. Where the wall folds sharply at the point (e.g. a sharp trailing edge, whose upper and lower faces
 * have nearly opposite normals) these vectors almost cancel and their sum is a small residue whose direction is set by
 * the local face sizes, not by the geometry. The faces of a point are therefore grouped.
 *
 * Compatibility: N_f . N_g > c0 |N_f| |N_g|, c0 = COMPATIBLE_COS. For a group G of nonzero, pairwise compatible corner
 * vectors the summed normal n_G = sum_{g in G} N_g satisfies, for every face f of G and any number of faces,
 *   n_G . N_f = |N_f|^2 + sum_{g != f} N_f . N_g >= |N_f|^2 > 0        (inside the dual cone of the faces), and
 *   |n_G|^2 = sum_g |N_g|^2 + sum_{f != g} N_f . N_g >= sum_g |N_g|^2   (no cancellation),
 * whatever the face areas. For a pair with N_f . N_g < 0 both fail for some area ratio (scaling N_g by
 * t > |N_f|^2 / (-N_f . N_g) gives n . N_f < 0), so 90 degrees is the bound; c0 > 0 is a round-off guard above it, so
 * that a nominal right angle is a sharp edge whatever the round-off.
 *
 * Rule: links are the pairs of faces that are adjacent (they share an edge through the point; in 2D the two segments of
 * a point), permitted (same link class) and compatible, weighted by their cosine. The groups are the maximal pairwise
 * compatible clusters of the single-linkage hierarchy of the links: the connected components of the links, where every
 * component that is not pairwise compatible loses all its weakest links (smallest cosine, exact ties together) and is
 * split into its components, recursively. Each pass removes at least one link of every incompatible component and a
 * single face is compatible, so this terminates, with any number of groups (1 to nFaces). The result depends only on
 * the set of faces, their adjacency, link classes and cosines (each computed by the same symmetric expression for every
 * pair), not on the order of the faces; only the numbering of the groups and the summation order of their normals
 * follow that order. Where every link component is pairwise compatible the groups are the components (the result of a
 * greedy merge with the same compatibility test). A face compatible with two mutually incompatible neighbours of an
 * otherwise isolated chain joins the one with the larger cosine, and forms its own group if the cosines are equal; in a
 * longer component weaker links elsewhere are removed first.
 *
 * Implied feature angle: pairwise compatible faces (less than 90 degrees apart) are averaged into one group, exactly as
 * the vertex normal averages all faces of a point. 90 degrees is not a tuned threshold but the largest fold angle for
 * which the averaged normal can never cancel and stays inside the dual cone for any face areas. ---*/

namespace SharpEdgeNormals {

/*! \brief Largest number of incident faces handled at a point (more: the caller keeps the vertex normal). */
constexpr unsigned short MAX_FACES = 32;

/*! \brief Compatibility threshold c0 of the cosine of two corner vectors (the parallel tolerance of
 *         CGeometry::ComputeModifiedSymmetryNormals). */
constexpr passivedouble COMPATIBLE_COS = 0.001;

/*! \brief SU2's parallel test of ComputeModifiedSymmetryNormals, 1 - |cos| < 0.001, as a squared sine: 2 t - t^2. */
constexpr passivedouble PARALLEL_SIN2 = 2.0 * 0.001 - 0.001 * 0.001;

/*!
 * \brief Groups of the nFaces incident faces of a point (see the rule above).
 * \param[in] nDim - Number of dimensions.
 * \param[in] nFaces - Number of faces (at most MAX_FACES).
 * \param[in] N - Nonzero corner area vectors of the faces at the point, consistently oriented.
 * \param[in] edge - For each face the nodes joined to the point by an edge of the face (2 in 3D; unused in 2D, where
 *            the two segments of a point are adjacent).
 * \param[in] nEdge - Number of entries of edge per face.
 * \param[out] groupN - Summed corner vector of every group, groups numbered in the order of their first face.
 * \param[out] faceGroup - Group of every face (optional).
 * \param[in] linkClass - Link class of every face (optional): two faces can only be linked if their classes are equal.
 * \return Number of groups.
 */
inline unsigned short FaceGroups(unsigned short nDim, unsigned short nFaces, const su2double (*N)[3],
                                 const unsigned long (*edge)[2], const unsigned short* nEdge, su2double (*groupN)[3],
                                 unsigned short* faceGroup = nullptr, const unsigned short* linkClass = nullptr) {
  constexpr unsigned short MAX_LINKS = MAX_FACES * (MAX_FACES - 1) / 2;

  /*--- Products and lengths; the test is symmetric in f and g (c0 multiplies the product of the lengths). ---*/
  su2double dotFG[MAX_FACES][MAX_FACES], len[MAX_FACES];
  for (unsigned short f = 0; f < nFaces; ++f) {
    for (unsigned short g = 0; g < nFaces; ++g) {
      su2double d = 0.0;
      for (unsigned short k = 0; k < nDim; ++k) d += N[f][k] * N[g][k];
      dotFG[f][g] = d;
    }
    len[f] = sqrt(dotFG[f][f]);
  }
  auto compatible = [&](unsigned short f, unsigned short g) {
    return dotFG[f][g] > COMPATIBLE_COS * (len[f] * len[g]);
  };

  /*--- Links: adjacent, permitted and compatible pairs, with their cosine. ---*/
  unsigned short linkA[MAX_LINKS], linkB[MAX_LINKS], nLinks = 0;
  su2double linkCos[MAX_LINKS];
  bool active[MAX_LINKS];
  for (unsigned short f = 0; f < nFaces; ++f) {
    for (unsigned short g = f + 1; g < nFaces; ++g) {
      if (linkClass && (linkClass[f] != linkClass[g])) continue;
      bool adjacent = (nDim == 2);
      for (unsigned short a = 0; a < nEdge[f] && !adjacent; ++a)
        for (unsigned short b = 0; b < nEdge[g]; ++b) adjacent = adjacent || (edge[f][a] == edge[g][b]);
      if (!adjacent || !compatible(f, g)) continue;
      linkA[nLinks] = f;
      linkB[nLinks] = g;
      linkCos[nLinks] = dotFG[f][g] / (len[f] * len[g]);
      active[nLinks] = true;
      ++nLinks;
    }
  }

  /*--- Components of the active links, each labelled by its smallest face. Every component that is not pairwise
   * compatible loses all of its weakest links; repeat until all components are compatible. ---*/
  unsigned short group[MAX_FACES];
  while (true) {
    for (unsigned short f = 0; f < nFaces; ++f) group[f] = f;
    for (bool changed = true; changed;) {
      changed = false;
      for (unsigned short l = 0; l < nLinks; ++l) {
        if (!active[l] || (group[linkA[l]] == group[linkB[l]])) continue;
        const auto keep = std::min(group[linkA[l]], group[linkB[l]]);
        const auto drop = std::max(group[linkA[l]], group[linkB[l]]);
        for (unsigned short p = 0; p < nFaces; ++p)
          if (group[p] == drop) group[p] = keep;
        changed = true;
      }
    }
    bool split = false;
    for (unsigned short r = 0; r < nFaces; ++r) {
      if (group[r] != r) continue;
      bool clique = true;
      for (unsigned short p = r; p < nFaces && clique; ++p) {
        if (group[p] != r) continue;
        for (unsigned short q = p + 1; q < nFaces && clique; ++q)
          if (group[q] == r) clique = compatible(p, q);
      }
      if (clique) continue;
      su2double weakest = 2.0;
      for (unsigned short l = 0; l < nLinks; ++l)
        if (active[l] && (group[linkA[l]] == r)) weakest = std::min(weakest, linkCos[l]);
      for (unsigned short l = 0; l < nLinks; ++l)
        if (active[l] && (group[linkA[l]] == r) && (linkCos[l] <= weakest)) active[l] = false;
      split = true;
    }
    if (!split) break;
  }

  /*--- Number the groups in the order of their first face and sum their corner vectors. ---*/
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
 * \brief Orthonormal basis Q (rows) of the span of nNormals normals: Gram-Schmidt in the given order, a normal dropped
 * when its orthogonalized squared norm (relative to its own) is below PARALLEL_SIN2. \return Rank (number of rows of Q,
 * at most nDim).
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
 * \brief vec[iVel + d] <- ((I - Q^T Q) vec)[iVel + d], the components of vec[iVel..iVel+nDim-1] along the orthonormal
 * rows of Q removed.
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

}  // namespace SharpEdgeNormals
