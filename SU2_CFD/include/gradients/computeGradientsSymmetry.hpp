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
 * \brief Reflect a gradient using a tensor mapping. Used for symmetry reflection.
 * \ingroup FvmAlgos
 * \param[in] nDim - number of dimensions, 2 or 3.
 * \param[in] varBegin - start of the variables.
 * \param[in] varEnd - end of the variables.
 * \param[in] isFlowSolver - are we using the flow solver.
 * \param[in] TensorMAp - the tensor map to map to rotated base.
 * \param[out] Gradients_iPoint - the gradient for the point.
 */


template <typename Int, class Matrix>
inline void ReflectGradient(Int nDim, size_t& varBegin, size_t& varEnd, Matrix& TensorMap,
                            Matrix& Gradients_iPoint, int idx_vel) {
  static constexpr size_t MAXNDIM = 3;

  su2activematrix Gradients_Velocity(nDim, nDim);
  su2activematrix Gradients_Velocity_Reflected(nDim, nDim);
  su2double gradPhi[MAXNDIM] = {0.0};
  su2double gradPhiReflected[MAXNDIM] = {0.0};

  if (idx_vel != -1) {
    /*--- Get gradients of primitives of boundary cell ---*/
    for (auto iVar = 0u; iVar < nDim; iVar++) {
      for (auto iDim = 0u; iDim < nDim; iDim++) {
        // todo: 1 ->idx.velocity
        Gradients_Velocity[iVar][iDim] = Gradients_iPoint[idx_vel + iVar][iDim];
        Gradients_Velocity_Reflected[iVar][iDim] = 0.0;
      }
    }

    /*--- Q' = L^T*Q*T ---*/
    for (auto iDim = 0u; iDim < nDim; iDim++) {
      for (auto jDim = 0u; jDim < nDim; jDim++) {
        for (auto kDim = 0u; kDim < nDim; kDim++) {
          for (auto mDim = 0u; mDim < nDim; mDim++) {
            Gradients_Velocity_Reflected[iDim][jDim] +=
                TensorMap[iDim][mDim] * TensorMap[jDim][kDim] * Gradients_Velocity[mDim][kDim];
          }
        }
      }
    }

    /*--- we have aligned such that U is the direction of the normal
     *    in 2D: dU/dy = dV/dx = 0
     *    in 3D: dU/dy = dV/dx = 0
     *           dU/dz = dW/dx = 0 ---*/
    for (auto iDim = 1u; iDim < nDim; iDim++) {
      Gradients_Velocity_Reflected[0][iDim] = 0.0;
      Gradients_Velocity_Reflected[iDim][0] = 0.0;
    }

    for (auto iDim = 0u; iDim < nDim; iDim++) {
      for (auto jDim = 0u; jDim < nDim; jDim++) {
        Gradients_Velocity[iDim][jDim] = 0.0;
      }
    }

    /*--- now transform back the corrected velocity gradients by taking the inverse again
     * T = (L^-1)*T' ---*/
    for (auto iDim = 0u; iDim < nDim; iDim++) {
      for (auto jDim = 0u; jDim < nDim; jDim++) {
        for (auto kDim = 0u; kDim < nDim; kDim++) {
          for (auto mDim = 0u; mDim < nDim; mDim++) {
            Gradients_Velocity[iDim][jDim] +=
                TensorMap[mDim][iDim] * TensorMap[kDim][jDim] * Gradients_Velocity_Reflected[mDim][kDim];
          }
        }
      }
    }

    for (auto iDim = 0u; iDim < nDim; iDim++) {
      for (auto jDim = 0u; jDim < nDim; jDim++) {
        // todo: 1->idx.velocity
        Gradients_iPoint[iDim + idx_vel][jDim] = Gradients_Velocity[iDim][jDim];
      }
    }
  }

  /*--- Reflect the gradients for all scalars (we exclude velocity). --*/
  for (auto iVar = varBegin; iVar < varEnd; iVar++) {
    if ((idx_vel == -1) || ((idx_vel != -1) && (iVar == 0 || iVar > nDim))) {
      /*--- project to symmetry aligned base ---*/
      for (auto iDim = 0u; iDim < nDim; iDim++) {
        gradPhi[iDim] = Gradients_iPoint[iVar][iDim];
        gradPhiReflected[iDim] = 0.0;
      }

      for (auto jDim = 0u; jDim < nDim; jDim++) {
        for (auto iDim = 0u; iDim < nDim; iDim++) {
          /*--- map transpose T' * grad(phi) ---*/
          gradPhiReflected[jDim] += TensorMap[jDim][iDim] * Gradients_iPoint[iVar][iDim];
        }
      }

      for (auto iDim = 0u; iDim < nDim; iDim++) gradPhi[iDim] = 0.0;

      /*--- gradient in direction normal to symmetry is cancelled ---*/
      gradPhiReflected[0] = 0.0;

      /*--- Now transform back ---*/
      for (auto jDim = 0u; jDim < nDim; jDim++) {
        for (auto iDim = 0u; iDim < nDim; iDim++) {
          gradPhi[jDim] += TensorMap[iDim][jDim] * gradPhiReflected[iDim];
        }
      }

      for (auto iDim = 0u; iDim < nDim; iDim++) Gradients_iPoint[iVar][iDim] = gradPhi[iDim];
    }
  }
}

/*! \brief Construct a 2D or 3D base given a normal vector.
           Constructs 1 (2D) or 2 (3D) additional vectors orthogonal to the normal to form a base. */
template <class Matrix, class Scalar, typename Int>
inline void BaseFromNormal(Int nDim, const Scalar* UnitNormal, Matrix& TensorMap) {
  /*--- Preprocessing: Compute unit tangential, the direction is arbitrary as long as
        t*n=0 && |t|_2 = 1 ---*/
  Scalar Tangential[3] = {0.0};
  Scalar Orthogonal[3] = {0.0};
  switch (nDim) {
    case 2: {
      Tangential[0] = -UnitNormal[1];
      Tangential[1] = UnitNormal[0];
      for (auto iDim = 0u; iDim < nDim; iDim++) {
        TensorMap[0][iDim] = UnitNormal[iDim];
        TensorMap[1][iDim] = Tangential[iDim];
      }
      break;
    }
    case 3: {
      /*--- n = ai + bj + ck, if |b| > |c| ---*/
      if (abs(UnitNormal[1]) > abs(UnitNormal[2])) {
        /*--- t = bi + (c-a)j - bk  ---*/
        Tangential[0] = UnitNormal[1];
        Tangential[1] = UnitNormal[2] - UnitNormal[0];
        Tangential[2] = -UnitNormal[1];
      } else {
        /*--- t = ci - cj + (b-a)k  ---*/
        Tangential[0] = UnitNormal[2];
        Tangential[1] = -UnitNormal[2];
        Tangential[2] = UnitNormal[1] - UnitNormal[0];
      }
      /*--- Make it a unit vector. ---*/
      Scalar TangentialNorm = GeometryToolbox::Norm(3, Tangential);
      Tangential[0] = Tangential[0] / TangentialNorm;
      Tangential[1] = Tangential[1] / TangentialNorm;
      Tangential[2] = Tangential[2] / TangentialNorm;

      /*--- Compute 3rd direction of the base using cross product ---*/
      GeometryToolbox::CrossProduct(UnitNormal, Tangential, Orthogonal);

      // now we construct the tensor mapping T, note that its inverse is the transpose of T
      for (auto iDim = 0u; iDim < nDim; iDim++) {
        TensorMap[0][iDim] = UnitNormal[iDim];
        TensorMap[1][iDim] = Tangential[iDim];
        TensorMap[2][iDim] = Orthogonal[iDim];
      }
      break;
    }
  }  // switch
}


}


template <class FieldType, class GradientType>
void computeGradientsSymmetry(unsigned short nDim, CSolver* solver, ENUM_MPI_QUANTITIES kindMpiComm, PERIODIC_QUANTITIES kindPeriodicComm,
                                CGeometry& geometry, const CConfig& config, const FieldType& field, size_t varBegin,
                                size_t varEnd, GradientType& gradient, int idx_vel) {

  static constexpr size_t MAXNDIM = 3;
  static constexpr size_t MAXNSYMS = 100;

 /* For symmetry planes, we need to impose the conditions (Blazek eq. 8.40):
   * 1. n.grad(phi) = 0
   * 2. n.grad(v.t) = 0
   * 3. t.grad(v.n) = 0
   */
  /*--- Check how many symmetry planes there are ---*/
  unsigned short Syms[MAXNSYMS] = {0};
  unsigned short nSym = 0;
  for (size_t iMarker = 0; iMarker < geometry.GetnMarker(); ++iMarker) {
    if ((config.GetMarker_All_KindBC(iMarker) == SYMMETRY_PLANE) ||
       (config.GetMarker_All_KindBC(iMarker) == EULER_WALL)) {
    Syms[nSym] = iMarker;
    nSym++;
    }
  }

  if (nSym==0) return;

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

        /*--- Normal of the primary symmetry plane ---*/
        su2double NormalPrim[MAXNDIM] = {0.0}, UnitNormalPrim[MAXNDIM] = {0.0};

        /*--- At this point we can find out if the node is shared with another symmetry.
         * Step 1: do we have other symmetries? ---*/
        if (nSym>1) {
          /*--- Step 2: are we on a shared node? ---*/
          for (auto jMarker=0;jMarker<nSym;jMarker++) {
            /*--- we do not need the current symmetry ---*/
            if (iMarker!= Syms[jMarker]) {
              /*--- Loop over all points on the other symmetry and check if current iPoint is on the symmetry ---*/
              for (auto jVertex = 0ul; jVertex < geometry.nVertex[Syms[jMarker]]; jVertex++) {
                const auto jPoint = geometry.vertex[Syms[jMarker]][jVertex]->GetNode();
                if (iPoint==jPoint) {
                  /*--- Does the other symmetry have a lower ID? Then that is the primary symmetry ---*/
                  if (Syms[jMarker]<iMarker) {
                    /*--- So whe have to get the normal of that other marker ---*/
                    geometry.vertex[Syms[jMarker]][jVertex]->GetNormal(NormalPrim);
                    su2double AreaPrim = GeometryToolbox::Norm(nDim, NormalPrim);
                    for(unsigned short iDim = 0; iDim < nDim; iDim++) {
                      UnitNormalPrim[iDim] = NormalPrim[iDim] / AreaPrim;
                    }

                    /*--- Correct the current normal as n2_new = n2 - (n2.n1)n1 ---*/
                    su2double ProjNorm = 0.0;
                    for (auto iDim = 0u; iDim < nDim; iDim++) ProjNorm += UnitNormal[iDim] * UnitNormalPrim[iDim];
                    /*--- We check if the normal of the 2 planes coincide.
                     * We only update the normal if the normals of the symmetry planes are different. ---*/
                    if (fabs(1.0-ProjNorm)>EPS) {
                      for (auto iDim = 0u; iDim < nDim; iDim++) UnitNormal[iDim] -= ProjNorm * UnitNormalPrim[iDim];
                      /* Make normalized vector ---*/
                      su2double newarea=GeometryToolbox::Norm(nDim, UnitNormal);
                      for (auto iDim = 0u; iDim < nDim; iDim++) UnitNormal[iDim] = UnitNormal[iDim]/newarea;
                    }
                  }
                }
              }
            }
          }
        }

        /*--- Tensor mapping from global Cartesian to local Cartesian aligned with symmetry plane ---*/
        su2activematrix TensorMap(nDim,nDim);

        /*--- Compute a new base for TensorMap aligned with the unit normal ---*/
        detail::BaseFromNormal(nDim,UnitNormal,TensorMap);

        su2activematrix Gradients_iPoint(varEnd-varBegin,nDim);

        for (auto iVar = varBegin; iVar < varEnd; iVar++) {
          for (auto iDim = 0u; iDim < nDim; iDim++) {
            Gradients_iPoint[iVar][iDim] = gradient(iPoint, iVar, iDim);
          }
        }

        detail::ReflectGradient(nDim, varBegin,varEnd, TensorMap, Gradients_iPoint, idx_vel);

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

