/*!
 * \file computeGradientsSymmetry.hpp
 * \brief Implements the symmetry boundary conditions for the gradient computations.
 * \author N. Beishuizen
 * \version 8.0.1 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2024, SU2 Contributors (cf. AUTHORS.md)
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
 * \param[in] nDim - Number of dimensions, 2 or 3.
 * \param[in] varBegin - Start of the variables.
 * \param[in] varEnd - End of the variables.
 * \param[in] n - Normal direction.
 * \param[in] idxVel - Variable index where velocity gradients start (-1 if all variables are scalars).
 * \param[in, out] gradients - The gradients to be modified.
 */
template <class Matrix, class Scalar>
inline void CorrectGradient(const int nDim_, const int varBegin, const int varEnd, const Scalar* n,
                            const int idxVel, Matrix&& gradients) {
  static constexpr size_t MAXNDIM = 3;
  const int nDim = nDim_ == 2 ? 2 : 3;

  /*--- First we correct the part that involves velocities. ---*/

  if (idxVel != -1) {
    /*--- Normal gradient of velocity components and gradient of normal velocity. ---*/
    su2double normalGrad[MAXNDIM] = {}, gradNormalVel[MAXNDIM] = {};

    for (auto iDim = 0; iDim < nDim; iDim++) {
      normalGrad[iDim] = GeometryToolbox::DotProduct(nDim, gradients[idxVel + iDim], n);

      for (auto jDim = 0; jDim < nDim; jDim++) {
        gradNormalVel[jDim] += n[iDim] * gradients[idxVel + iDim][jDim];
      }
    }

    /*--- Normal gradient of the normal velocity. ---*/
    const su2double normalGradNormalVel = GeometryToolbox::DotProduct(nDim, n, gradNormalVel);

    /*--- Remove the tangential projection (I - n.n^T) of the normal gradients.
     * And the normal projection (n.n^T) of the tangential gradients.
     * dV = dV - n.n^T dV (I - n.n^T) - (I - n.n^T) dV n.n^T ---*/

    for (auto iDim = 0; iDim < nDim; iDim++) {
      for (auto jDim = 0; jDim < nDim; jDim++) {
        gradients[idxVel + iDim][jDim] -= normalGrad[iDim] * n[jDim] + n[iDim] * gradNormalVel[jDim];
        gradients[idxVel + iDim][jDim] += 2 * n[iDim] * normalGradNormalVel * n[jDim];
      }
    }
  }

  /*--- Remove the normal component for all scalars (excluding velocities). ---*/

  for (auto iVar = varBegin; iVar < varEnd; iVar++) {
    if (idxVel != -1 && iVar >= idxVel && iVar < idxVel + nDim) continue;

    const su2double normalGrad = GeometryToolbox::DotProduct(nDim, n, gradients[iVar]);
    for (auto iDim = 0; iDim < nDim; iDim++) {
      gradients[iVar][iDim] -= normalGrad * n[iDim];
    }
  }
}

}


template <class FieldType, class GradientType>
void computeGradientsSymmetry(unsigned short nDim, CSolver* solver, MPI_QUANTITIES kindMpiComm, PERIODIC_QUANTITIES kindPeriodicComm,
                                CGeometry& geometry, const CConfig& config, const FieldType& field, size_t varBegin,
                                size_t varEnd, GradientType& gradient, short idxVel) {

  static constexpr size_t MAXNDIM = 3;

  /* For symmetry planes (and Euler walls), we need to impose the conditions (Blazek eq. 8.40):
   * 1. n.grad(phi) = 0
   * 2. n.grad(v.t) = 0
   * 3. t.grad(v.n) = 0
   */

  /*--- Check how many symmetry planes there are ---*/
  unsigned short nSym = 0;
  for (size_t iMarker = 0; iMarker < geometry.GetnMarker(); ++iMarker) {
    if ((config.GetMarker_All_KindBC(iMarker) == SYMMETRY_PLANE) ||
       (config.GetMarker_All_KindBC(iMarker) == EULER_WALL)) {
      nSym++;
    }
  }

  /*--- No symmetry or Euler walls are present. ---*/
  if (nSym == 0) return;

  for (size_t iMarker = 0; iMarker < geometry.GetnMarker(); ++iMarker) {

    if ((config.GetMarker_All_KindBC(iMarker) == SYMMETRY_PLANE) || (config.GetMarker_All_KindBC(iMarker) == EULER_WALL)) {
      SU2_OMP_FOR_STAT(32)
      for (size_t iVertex = 0; iVertex < geometry.GetnVertex(iMarker); ++iVertex) {

        size_t iPoint = geometry.vertex[iMarker][iVertex]->GetNode();

        /*--- Normal vector for this vertex (negate for outward convention). ---*/
        const su2double* VertexNormal = geometry.vertex[iMarker][iVertex]->GetNormal();

        const auto NormArea = GeometryToolbox::Norm(nDim, VertexNormal);

        su2double UnitNormal[MAXNDIM] = {0.0};
        for (size_t iDim = 0; iDim < nDim; iDim++)
          UnitNormal[iDim] = VertexNormal[iDim] / NormArea;

        /*--- When we have more than 1 symmetry or Euler wall, check if there are shared nodes.
              Then correct the normal at those node ---*/
        if (nSym > 1) {
          if (geometry.symmetryNormals.size() > 0) {
            auto it =  std::find_if(
              geometry.symmetryNormals[iMarker].begin(),
              geometry.symmetryNormals[iMarker].end(),
              findSymNormalIndex(iPoint));
            if (it != geometry.symmetryNormals[iMarker].end()) {
              for (auto iDim = 0u; iDim < nDim; iDim++) UnitNormal[iDim] = it->normal[iDim];
            }
          }
        }

        detail::CorrectGradient(nDim, varBegin, varEnd, UnitNormal, idxVel, gradient[iPoint]);

      } // loop over vertices
      END_SU2_OMP_FOR
    } // symmetry
  } // markers

}

