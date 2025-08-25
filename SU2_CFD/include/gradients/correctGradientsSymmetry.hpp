/*!
 * \file correctGradientsSymmetry.hpp
 * \brief Implements the symmetry boundary conditions for the gradient computations.
 * \author N. Beishuizen
 * \version 8.2.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2025, SU2 Contributors (cf. AUTHORS.md)
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

#include <vector>
#include <algorithm>

#include "../../../Common/include/parallelization/omp_structure.hpp"
#include "../../../Common/include/toolboxes/geometry_toolbox.hpp"

namespace detail {
/*!
 * \brief Correct the gradient on a symmetry plane.
 * \ingroup FvmAlgos
 * \param[in] varBegin - Start of the variables.
 * \param[in] varEnd - End of the variables.
 * \param[in] idxVel - Variable index where velocity gradients start (-1 if all variables are scalars).
 * \param[in] n - Normal direction.
 * \param[in, out] gradients - The gradients to be modified.
 */
template <size_t nDim, class Matrix, class Scalar>
inline void correctGradient(const size_t varBegin, const size_t varEnd, const int idxVel, const Scalar* n,
                            Matrix&& gradients) {

  /*--- First we correct the part that involves velocities. ---*/

  if (idxVel >= static_cast<int>(varBegin) && idxVel + nDim <= varEnd) {
    /*--- Normal gradient of velocity components and gradient of normal velocity. ---*/
    su2double normalGrad[nDim] = {}, gradNormalVel[nDim] = {};

    for (size_t iDim = 0; iDim < nDim; iDim++) {
      normalGrad[iDim] = GeometryToolbox::DotProduct(nDim, gradients[idxVel + iDim], n);

      for (size_t jDim = 0; jDim < nDim; jDim++) {
        gradNormalVel[jDim] += n[iDim] * gradients[idxVel + iDim][jDim];
      }
    }

    /*--- Normal gradient of the normal velocity. ---*/
    const su2double normalGradNormalVel = GeometryToolbox::DotProduct(nDim, n, gradNormalVel);

    /*--- Remove the tangential projection (I - n.n^T) of the normal gradients.
     * And the normal projection (n.n^T) of the tangential gradients.
     * dV = dV - (I - n.n^T) dV n.n^T - n.n^T dV (I - n.n^T) ---*/

    for (size_t iDim = 0; iDim < nDim; iDim++) {
      for (size_t jDim = 0; jDim < nDim; jDim++) {
        gradients[idxVel + iDim][jDim] -= normalGrad[iDim] * n[jDim] + n[iDim] * gradNormalVel[jDim];
        gradients[idxVel + iDim][jDim] += 2 * n[iDim] * normalGradNormalVel * n[jDim];
      }
    }
  }

  /*--- Remove the normal component for all scalars (excluding velocities). ---*/

  for (auto iVar = varBegin; iVar < varEnd; iVar++) {
    if (idxVel != -1 && static_cast<int>(iVar) >= idxVel && iVar < idxVel + nDim) continue;

    const su2double normalGrad = GeometryToolbox::DotProduct(nDim, n, gradients[iVar]);
    for (size_t iDim = 0; iDim < nDim; iDim++) {
      gradients[iVar][iDim] -= normalGrad * n[iDim];
    }
  }
}
}

/*!
 * \brief Correct gradients on symmetry and Euler (slip) markers to respect the conditions:
 *        1. n.grad(phi) = 0
 *        2. n.grad(v.t) = 0
 *        3. t.grad(v.n) = 0
 * \note See Blazek eq. 8.40.
 * \ingroup FvmAlgos
 * \param[in] geometry - Geometric grid properties.
 * \param[in] config - Configuration of the problem, used to identify types of boundaries.
 * \param[in] varBegin - Index of first variable for which to compute the gradient.
 * \param[in] varEnd - Index of last variable for which to compute the gradient.
 * \param[in] idxVel - Index to velocity, -1 if no velocity is present in the solver.
 * \param[in,out] gradient - Generic object implementing operator (iPoint, iVar, iDim).
 */
template <size_t nDim, class GradientType>
void correctGradientsSymmetry(CGeometry& geometry, const CConfig& config, const size_t varBegin,
                              const size_t varEnd, const int idxVel, GradientType& gradient) {

  /*--- Check how many symmetry planes there are. ---*/
  std::vector<unsigned short> symMarkers;
  for (auto iMarker = 0u; iMarker < geometry.GetnMarker(); ++iMarker) {
    if (config.GetMarker_All_KindBC(iMarker) == SYMMETRY_PLANE ||
        config.GetMarker_All_KindBC(iMarker) == EULER_WALL) {
      symMarkers.push_back(iMarker);
    }
  }

  for (const auto iMarker : symMarkers) {
    SU2_OMP_FOR_STAT(32)
    for (size_t iVertex = 0; iVertex < geometry.GetnVertex(iMarker); ++iVertex) {

      const auto iPoint = geometry.vertex[iMarker][iVertex]->GetNode();

      /*--- Get the normal of the current symmetry. This may be the original normal of the vertex
       * or a modified normal if there are intersecting symmetries. ---*/

      su2double unitNormal[nDim] = {};
      const auto it = geometry.symmetryNormals[iMarker].find(iVertex);

      if (it != geometry.symmetryNormals[iMarker].end()) {
        for (auto iDim = 0u; iDim < nDim; iDim++) unitNormal[iDim] = it->second[iDim];
      } else {
        geometry.vertex[iMarker][iVertex]->GetNormal(unitNormal);
        const su2double area = GeometryToolbox::Norm(nDim, unitNormal);
        for (auto iDim = 0u; iDim < nDim; iDim++) unitNormal[iDim] /= area;
      }

      detail::correctGradient<nDim>(varBegin, varEnd, idxVel, unitNormal, gradient[iPoint]);

    } // loop over vertices
    END_SU2_OMP_FOR
  } // markers

}

