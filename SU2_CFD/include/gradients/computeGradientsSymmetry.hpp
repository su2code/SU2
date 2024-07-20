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
template <typename Int, class Matrix>
inline void CorrectGradient(Int nDim, size_t& varBegin, size_t& varEnd, Matrix& TensorMap,
                            Matrix& Gradients_iPoint, short idxVel) {
  static constexpr size_t MAXNDIM = 3;

  su2activematrix Gradients_Velocity(nDim, nDim);
  su2activematrix Gradients_Velocity_Reflected(nDim, nDim);
  su2double gradPhi[MAXNDIM] = {0.0};

  /*--- normal vector---*/
  su2double n[MAXNDIM] = {0.0};
  /*--- n.n^T ---*/
  su2activematrix nn(nDim, nDim);
  su2activematrix Vnn(nDim, nDim);
  su2activematrix nnV(nDim, nDim);
  //su2activematrix nnVnn(nDim, nDim);

  // get the normal vector of the symmetry/Euler
  for (auto iDim = 0u; iDim < nDim; iDim++) {
    n[iDim] = TensorMap[0][iDim];
  }
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
        //Gradients_Velocity_Reflected[iVar][iDim] = 0.0;
        nn[iVar][iDim] = 0.0;
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
        for (auto kDim = 0u; kDim < nDim; kDim++) {
          Gradients_Velocity[iDim][jDim] += 2.0*nnV[iDim][kDim]*nn[kDim][jDim];
        }
      }
   }

   for (auto iDim = 0u; iDim < nDim; iDim++) {
      for (auto jDim = 0u; jDim < nDim; jDim++) {
        Gradients_Velocity[iDim][jDim] -= Vnn[iDim][jDim] - nnV[iDim][jDim];
      }
   }



    for (auto iDim = 0u; iDim < nDim; iDim++) {
      for (auto jDim = 0u; jDim < nDim; jDim++) {
        Gradients_iPoint[iDim + idxVel][jDim] = Gradients_Velocity[iDim][jDim];
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
          Gradients_iPoint[iVar][iDim] -= TensorMap[0][iDim] * TensorMap[0][jDim] * gradPhi[jDim];
        }
      }

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

      /*--- now we construct the tensor mapping T, note that its inverse is the transpose of T ---*/
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
//cout << "ipoint="<<iPoint << endl;
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

        /*--- Tensor mapping from global Cartesian to local Cartesian aligned with symmetry plane ---*/
        su2activematrix TensorMap(nDim,nDim);

        /*--- Compute a new base for TensorMap aligned with the unit normal ---*/
        detail::BaseFromNormal(nDim,UnitNormal,TensorMap);

        su2activematrix Gradients_iPoint(varEnd-varBegin,nDim);

        /*--- Fill the local gradient tensor ---*/
        for (auto iVar = varBegin; iVar < varEnd; iVar++) {
          for (auto iDim = 0u; iDim < nDim; iDim++) {
            Gradients_iPoint[iVar][iDim] = gradient(iPoint, iVar, iDim);
          }
        }

        detail::CorrectGradient(nDim, varBegin,varEnd, TensorMap, Gradients_iPoint, idxVel);

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

