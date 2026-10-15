/*!
 * \file correctGradientsWall.hpp
 * \brief Imposes the known tangential derivatives on gradients at no-slip walls.
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

#include "WallGradientVars.hpp"
#include "../../../Common/include/parallelization/omp_structure.hpp"
#include "../../../Common/include/toolboxes/geometry_toolbox.hpp"

namespace detail {

/*!
 * \brief Whether the condition on a boundary marker is a no-slip wall.
 */
inline bool isNoSlipWall(const CConfig& config, unsigned short iMarker) {
  switch (config.GetMarker_All_KindBC(iMarker)) {
    case ISOTHERMAL:
    case HEAT_FLUX:
    case HEAT_TRANSFER:
    case CHT_WALL_INTERFACE:
      return true;
    default:
      return false;
  }
}

/*!
 * \brief Whether the vertex also lies on another boundary, where the wall normal is ambiguous.
 */
inline bool isSharedWithOtherBoundary(const CGeometry& geometry, const CConfig& config, unsigned long iPoint,
                                      unsigned short iMarker) {
  for (auto jMarker = 0u; jMarker < geometry.GetnMarker(); ++jMarker) {
    if (jMarker == iMarker || geometry.nodes->GetVertex(iPoint, jMarker) < 0) continue;
    const auto kind = config.GetMarker_All_KindBC(jMarker);
    if (kind != INTERNAL_BOUNDARY && kind != NEARFIELD_BOUNDARY && kind != PERIODIC_BOUNDARY) return true;
  }
  return false;
}

/*!
 * \brief Remove the tangential part of the gradient of a scalar that is uniform along the wall.
 */
template <size_t nDim, class Vector>
inline void removeTangentialGradient(const su2double* n, Vector&& grad) {
  const su2double normalGrad = GeometryToolbox::DotProduct(nDim, n, grad);
  for (size_t iDim = 0; iDim < nDim; iDim++) grad[iDim] = normalGrad * n[iDim];
}

}  // namespace detail

/*!
 * \brief Impose the known derivatives at no-slip walls: zero tangential derivative of the velocity
 *        components, the temperature (isothermal walls) and wall-vanishing variables (smooth walls).
 * \note With divFreeVel the normal derivative of the normal velocity is also removed (div(u) = 0).
 *       Vertices shared with another boundary and moving walls are left untouched.
 * \ingroup FvmAlgos
 * \param[in] geometry - Geometric grid properties.
 * \param[in] config - Configuration of the problem, used to identify types of boundaries.
 * \param[in] varBegin - Index of first variable for which the gradient was computed.
 * \param[in] varEnd - Index of last variable for which the gradient was computed.
 * \param[in] vars - Variables with known wall behavior.
 * \param[in,out] gradient - Generic object implementing operator (iPoint, iVar, iDim).
 */
template <size_t nDim, class GradientType>
void correctGradientsWall(CGeometry& geometry, const CConfig& config, const size_t varBegin, const size_t varEnd,
                          const WallGradientVars& vars, GradientType& gradient) {
  if (!vars.Active()) return;

  const auto inRange = [&](int i) { return i >= static_cast<int>(varBegin) && i < static_cast<int>(varEnd); };
  const bool velocity = vars.idxVel >= 0 && !config.GetDynamic_Grid() &&
                        inRange(vars.idxVel) && inRange(vars.idxVel + static_cast<int>(nDim) - 1);
  const bool temperature = vars.idxTemp >= 0 && inRange(vars.idxTemp);
  const bool zeroVars = vars.zeroEnd > vars.zeroBegin && inRange(vars.zeroBegin) && inRange(vars.zeroEnd - 1);

  for (auto iMarker = 0u; iMarker < geometry.GetnMarker(); ++iMarker) {
    if (!detail::isNoSlipWall(config, iMarker)) continue;

    const bool isothermal = config.GetMarker_All_KindBC(iMarker) == ISOTHERMAL && !config.GetMarker_All_PyCustom(iMarker);
    const bool smooth = !config.GetWall_Functions() &&
                        config.GetWallRoughnessProperties(config.GetMarker_All_TagBound(iMarker)).first == WALL_TYPE::SMOOTH;

    SU2_OMP_FOR_STAT(32)
    for (size_t iVertex = 0; iVertex < geometry.GetnVertex(iMarker); ++iVertex) {
      const auto iPoint = geometry.vertex[iMarker][iVertex]->GetNode();

      if (detail::isSharedWithOtherBoundary(geometry, config, iPoint, iMarker)) continue;

      su2double n[nDim] = {};
      geometry.vertex[iMarker][iVertex]->GetNormal(n);
      const su2double area = GeometryToolbox::Norm(nDim, n);
      if (area < EPS) continue;
      for (size_t iDim = 0; iDim < nDim; iDim++) n[iDim] /= area;

      auto&& grad = gradient[iPoint];

      if (velocity) {
        /*--- The velocity gradient reduces to dU/dn (x) n, with no normal-velocity component if div(u) = 0. ---*/
        su2double dUdn[nDim] = {};
        for (size_t iDim = 0; iDim < nDim; iDim++) dUdn[iDim] = GeometryToolbox::DotProduct(nDim, n, grad[vars.idxVel + iDim]);

        if (vars.divFreeVel) {
          const su2double dUndn = GeometryToolbox::DotProduct(nDim, n, dUdn);
          for (size_t iDim = 0; iDim < nDim; iDim++) dUdn[iDim] -= dUndn * n[iDim];
        }
        for (size_t iDim = 0; iDim < nDim; iDim++)
          for (size_t jDim = 0; jDim < nDim; jDim++) grad[vars.idxVel + iDim][jDim] = dUdn[iDim] * n[jDim];
      }

      if (temperature && isothermal) detail::removeTangentialGradient<nDim>(n, grad[vars.idxTemp]);

      if (zeroVars && smooth)
        for (int iVar = vars.zeroBegin; iVar < vars.zeroEnd; ++iVar) detail::removeTangentialGradient<nDim>(n, grad[iVar]);
    }
    END_SU2_OMP_FOR
  }
}
