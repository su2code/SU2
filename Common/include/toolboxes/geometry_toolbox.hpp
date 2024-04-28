/*!
 * \file geometry_toolbox.hpp
 * \brief Collection of common lightweight geometry-oriented methods.
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

#include <cmath>

// for su2activematrix
#include "../containers/C2DContainer.hpp"

namespace GeometryToolbox {
/// \addtogroup GeometryToolbox
/// @{

/*! \return ||a-b||^2 */
template <class T, class U, typename Int>
inline T SquaredDistance(Int nDim, const T* a, const U* b) {
  T d(0);
  for (Int i = 0; i < nDim; i++) d += pow(a[i] - b[i], 2);
  return d;
}

/*! \return ||a-b|| */
template <class T, class U, typename Int>
inline T Distance(Int nDim, const T* a, const U* b) {
  return sqrt(SquaredDistance(nDim, a, b));
}

/*! \brief d = a-b */
template <class T, typename Int>
inline void Distance(Int nDim, const T* a, const T* b, T* d) {
  for (Int i = 0; i < nDim; i++) d[i] = a[i] - b[i];
}

/*! \brief Reflect a at b: c = 2*b - a
 */
template <class T, typename Int>
inline void PointPointReflect(Int nDim, const T* a, const T* b, T* d) {
  for (Int i = 0; i < nDim; i++) d[i] = 2 * b[i] - a[i];
}

/*! \return a.b */
template <class T, typename Int>
inline T DotProduct(Int nDim, const T* a, const T* b) {
  T d(0);
  for (Int i = 0; i < nDim; ++i) d += a[i] * b[i];
  return d;
}

/*! \return ||a||^2 */
template <class T, typename Int>
inline T SquaredNorm(Int nDim, const T* a) {
  return DotProduct(nDim, a, a);
}

/*! \return ||a|| */
template <class T, typename Int>
inline T Norm(Int nDim, const T* a) {
  return sqrt(SquaredNorm(nDim, a));
}

/*! \brief c = a x b */
template <class T>
inline void CrossProduct(const T* a, const T* b, T* c) {
  c[0] = a[1] * b[2] - a[2] * b[1];
  c[1] = a[2] * b[0] - a[0] * b[2];
  c[2] = a[0] * b[1] - a[1] * b[0];
}

/*!
 * \brief Compute the coordinate (c) where the line defined by coordinate l0 and
 *        direction d intersects the plane defined by point p0 and normal n.
 * \return The intersection distance.
 */
template <class T, int nDim>
inline T LinePlaneIntersection(const T* l0, const T* d, const T* p0, const T* n, T* c) {
  T dist[nDim] = {0.0};
  Distance(nDim, p0, l0, dist);
  T alpha = DotProduct(nDim, dist, n) / DotProduct(nDim, d, n);
  for (int iDim = 0; iDim < nDim; ++iDim) c[iDim] = l0[iDim] + alpha * d[iDim];
  return fabs(alpha) * Norm(nDim, d);
}

/*!
 * \brief Compute the coordinate (c) where point p1 intersects the plane defined
 *        by point p0 and normal n if projected perpendicular to it.
 * \return The normal distance.
 */
template <class T, int nDim>
inline T PointPlaneProjection(const T* p1, const T* p0, const T* n, T* c) {
  return LinePlaneIntersection<T, nDim>(p1, n, p0, n, c);
}

/*! \brief Set U as the normal to a 2D line defined by coords[iPoint][iDim]. */
template <class T, class U>
inline void LineNormal(const T& coords, U* normal) {
  normal[0] = coords[0][1] - coords[1][1];
  normal[1] = coords[1][0] - coords[0][0];
}

/*! \brief Normal vector of a triangle, cross product of two sides. */
template <class T, class U>
inline void TriangleNormal(const T& coords, U* normal) {
  U a[3], b[3];

  for (int iDim = 0; iDim < 3; iDim++) {
    a[iDim] = coords[1][iDim] - coords[0][iDim];
    b[iDim] = coords[2][iDim] - coords[0][iDim];
  }

  CrossProduct(a, b, normal);
  normal[0] *= 0.5;
  normal[1] *= 0.5;
  normal[2] *= 0.5;
}

/*! \brief Normal vector of a quadrilateral, cross product of the two diagonals. */
template <class T, class U>
inline void QuadrilateralNormal(const T& coords, U* normal) {
  U a[3], b[3];

  for (int iDim = 0; iDim < 3; iDim++) {
    a[iDim] = coords[2][iDim] - coords[0][iDim];
    b[iDim] = coords[3][iDim] - coords[1][iDim];
  }

  CrossProduct(a, b, normal);
  normal[0] *= 0.5;
  normal[1] *= 0.5;
  normal[2] *= 0.5;
}

/*! \brief Signed distance from a point to a plane defined by 3 coordinates. */
template <class T, class U>
inline U PointToPlaneDistance(const T& coords, const U* point) {
  U normal[3], distance[3];
  TriangleNormal(coords, normal);
  Distance(3, point, coords[0], distance);
  return DotProduct(3, distance, normal) / Norm(3, normal);
}

/*!
 * \brief Compute a 3D rotation matrix.
 * \note The implicit ordering is rotation about the x, y, and then z axis.
 */
template <class Scalar, class Matrix>
inline void RotationMatrix(Scalar theta, Scalar phi, Scalar psi, Matrix& mat) {
  Scalar cosTheta = cos(theta);
  Scalar cosPhi = cos(phi);
  Scalar cosPsi = cos(psi);
  Scalar sinTheta = sin(theta);
  Scalar sinPhi = sin(phi);
  Scalar sinPsi = sin(psi);

  mat[0][0] = cosPhi * cosPsi;
  mat[1][0] = cosPhi * sinPsi;
  mat[2][0] = -sinPhi;

  mat[0][1] = sinTheta * sinPhi * cosPsi - cosTheta * sinPsi;
  mat[1][1] = sinTheta * sinPhi * sinPsi + cosTheta * cosPsi;
  mat[2][1] = sinTheta * cosPhi;

  mat[0][2] = cosTheta * sinPhi * cosPsi + sinTheta * sinPsi;
  mat[1][2] = cosTheta * sinPhi * sinPsi - sinTheta * cosPsi;
  mat[2][2] = cosTheta * cosPhi;
}

/*! \brief Compute a 2D rotation matrix. */
template <class Scalar, class Matrix>
inline void RotationMatrix(Scalar psi, Matrix& mat) {
  Scalar cosPsi = cos(psi);
  Scalar sinPsi = sin(psi);

  mat[0][0] = cosPsi;
  mat[0][1] = -sinPsi;
  mat[1][0] = sinPsi;
  mat[1][1] = cosPsi;
}

/*! \brief Apply a rotation matrix (R) about origin (O) to a point at
 *         distance (d) from it to obtain new coordinate (c). */
template <class Scalar, int nDim>
inline void Rotate(const Scalar R[][nDim], const Scalar* O, const Scalar* d, Scalar* c) {
  for (int iDim = 0; iDim < nDim; ++iDim) {
    c[iDim] = O[iDim];
    for (int k = 0; k < nDim; ++k) c[iDim] += R[iDim][k] * d[k];
  }
}

/*! \brief Tangent projection  */
template <class Mat, class Scalar, class Int>
inline void TangentProjection(Int nDim, const Mat& tensor, const Scalar* vector, Scalar* proj) {
  for (Int iDim = 0; iDim < nDim; iDim++) proj[iDim] = DotProduct(nDim, tensor[iDim], vector);

  auto normalProj = DotProduct(nDim, proj, vector);

  for (Int iDim = 0; iDim < nDim; iDim++) proj[iDim] -= normalProj * vector[iDim];
}

/*! \brief Reflect a gradient using a tensor mapping. Used for symmetry reflection. */
template <typename Int, class Matrix>
inline void ReflectGradient(Int nDim, size_t& varBegin, size_t& varEnd, bool isFlowSolver, Matrix& TensorMap,
                            Matrix& Gradients_iPoint) {
  su2activematrix Gradients_Velocity(nDim, nDim);
  su2activematrix Gradients_Velocity_Reflected(nDim, nDim);
  static constexpr size_t MAXNDIM = 3;
  su2double gradPhi[MAXNDIM] = {0.0};
  su2double gradPhiReflected[MAXNDIM] = {0.0};

  if (isFlowSolver == true) {
    /*--- Get gradients of primitives of boundary cell ---*/
    for (auto iVar = 0u; iVar < nDim; iVar++) {
      for (auto iDim = 0u; iDim < nDim; iDim++) {
        // Gradients_Velocity[iVar][iDim] = gradient(iPoint, 1 + iVar, iDim);
        Gradients_Velocity[iVar][iDim] = Gradients_iPoint[1 + iVar][iDim];
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
        // gradient(iPoint,iDim+1,jDim) = Gradients_Velocity[iDim][jDim];
        Gradients_iPoint[iDim + 1][jDim] = Gradients_Velocity[iDim][jDim];
      }
    }
  }

  /*--- Reflect the gradients for all scalars (we exclude velocity). --*/
  for (auto iVar = varBegin; iVar < varEnd; iVar++) {
    if ((isFlowSolver == false) || ((isFlowSolver == true) && (iVar == 0 || iVar > nDim))) {
      /*--- project to symmetry aligned base ---*/
      for (auto iDim = 0u; iDim < nDim; iDim++) {
        // gradPhi[iDim] = gradient(iPoint, iVar, iDim);
        gradPhi[iDim] = Gradients_iPoint[iVar][iDim];
        gradPhiReflected[iDim] = 0.0;
      }

      for (auto jDim = 0u; jDim < nDim; jDim++) {
        for (auto iDim = 0u; iDim < nDim; iDim++) {
          /*--- map transpose T' * grad(phi) ---*/
          // gradPhiReflected[jDim] += TensorMap[jDim][iDim]*gradient(iPoint,iVar,iDim);
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

      for (auto iDim = 0u; iDim < nDim; iDim++)
        // gradient(iPoint,iVar,iDim) = gradPhi[iDim];
        Gradients_iPoint[iVar][iDim] = gradPhi[iDim];
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
      Scalar TangentialNorm = Norm(3, Tangential);
      Tangential[0] = Tangential[0] / TangentialNorm;
      Tangential[1] = Tangential[1] / TangentialNorm;
      Tangential[2] = Tangential[2] / TangentialNorm;

      /*--- Compute 3rd direction of the base using cross product ---*/
      CrossProduct(UnitNormal, Tangential, Orthogonal);

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

/// @}
}  // namespace GeometryToolbox
