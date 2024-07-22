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
 * \brief Correct the gradient on a symmetry by using a tensor mapping.
 * \ingroup FvmAlgos
 * \param[in] nDim - number of dimensions, 2 or 3.
 * \param[in] varBegin - start of the variables.
 * \param[in] varEnd - end of the variables.
 * \param[in] isFlowSolver - are we using the flow solver.
 * \param[in] TensorMAp - the tensor map to map to rotated base.
 * \param[out] Gradients_iPoint - the gradient for the point.
 */
template <typename Int, class Matrix, class Scalar>
inline void CorrectGradient(Int nDim, size_t& varBegin, size_t& varEnd, const Scalar* n,
                            Matrix& Gradients_iPoint, short idxVel) {
  static constexpr size_t MAXNDIM = 3;

  su2activematrix Gradients_Velocity(nDim, nDim);
  su2activematrix Gradients_Velocity_Reflected(nDim, nDim);
  su2double gradPhi[MAXNDIM] = {0.0};

  /*--- n.n^T ---*/
  su2activematrix nn(nDim, nDim);
  su2activematrix Vnn(nDim, nDim);
  su2activematrix nnV(nDim, nDim);



  /*--- n.n^T ---*/
  for (auto jDim = 0u; jDim < nDim; jDim++) {
    for (auto iDim = 0u; iDim < nDim; iDim++) {
      nn[iDim][jDim] = n[iDim]*n[jDim];
    }
  }

  /*--- First we correct the part that involves velocities ---*/
  if (idxVel != -1) {
    /*--- Get gradients of primitives of boundary cell ---*/
    for (auto iVar = 0u; iVar < nDim; iVar++) {
      for (auto iDim = 0u; iDim < nDim; iDim++) {
        Gradients_Velocity[iVar][iDim] = Gradients_iPoint[idxVel + iVar][iDim];
        Vnn[iVar][iDim] = 0.0;
        nnV[iVar][iDim] = 0.0;
      }
    }

    for (auto iDim = 0u; iDim < nDim; iDim++) {
      for (auto jDim = 0u; jDim < nDim; jDim++) {
        for (auto kDim = 0u; kDim < nDim; kDim++) {
          Vnn[iDim][jDim] += Gradients_Velocity[iDim][kDim]*nn[kDim][jDim];
          nnV[iDim][jDim] += nn[iDim][kDim]*Gradients_Velocity[kDim][jDim];
        }
      }
    }

    for (auto iDim = 0u; iDim < nDim; iDim++) {
      for (auto jDim = 0u; jDim < nDim; jDim++) {
        Gradients_iPoint[iDim + idxVel][jDim] -= (Vnn[iDim][jDim] + nnV[iDim][jDim]);
        for (auto kDim = 0u; kDim < nDim; kDim++) {
          Gradients_iPoint[iDim + idxVel][jDim] += 2.0*nnV[iDim][kDim]*nn[kDim][jDim];
        }
      }
    }
  }

  /*-------------------------------------------------------------------*/
  /*--- Reflect the gradients for all scalars (we exclude velocity). --*/
  for (auto iVar = varBegin; iVar < varEnd; iVar++) {
    if ((idxVel == -1) || ((idxVel != -1) && (iVar == 0 || iVar > nDim))) {

       for (auto iDim = 0u; iDim < nDim; iDim++) {
         gradPhi[iDim] = Gradients_iPoint[iVar][iDim];
       }

      // I - n'.n.gradphi
      for (auto iDim = 0u; iDim < nDim; iDim++) {
        for (auto jDim = 0u; jDim < nDim; jDim++) {
          Gradients_iPoint[iVar][iDim] -= nn[iDim][jDim] * gradPhi[jDim];
        }
      }

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

        su2activematrix Gradients_iPoint(varEnd-varBegin,nDim);

        /*--- Fill the local gradient tensor ---*/
        for (auto iVar = varBegin; iVar < varEnd; iVar++) {
          for (auto iDim = 0u; iDim < nDim; iDim++) {
            Gradients_iPoint[iVar][iDim] = gradient(iPoint, iVar, iDim);
          }
        }

        detail::CorrectGradient(nDim, varBegin,varEnd, UnitNormal, Gradients_iPoint, idxVel);

        /*--- Write the corrected gradient tensor ---*/
        for (auto iVar = varBegin; iVar < varEnd; iVar++) {
          for (auto iDim = 0u; iDim < nDim; iDim++) {
            gradient(iPoint, iVar, iDim) = Gradients_iPoint[iVar][iDim];
          }
        }

      } // loop over vertices
      END_SU2_OMP_FOR
    } // symmetry
  } // markers

}

