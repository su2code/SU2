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
// Added by max

/*!
 * \brief Check if the edge intersects the plane
 * \param[in] nDim - Number of dimensions of the particular problem.
 * \param[in] p0 - Point defining the plane.
 * \param[in] n - Plane normal vector.
 * \param[in] p1,p2 - Points defining the edge
 * \return Returns <code>TRUE</code> if the edge intersect the plane, <code>FALSE</code> it doesn't.
 */
template <class T, class Int>
inline bool IntersectEdge(Int nDim, const T* p0, const T* n, const T* p1, const T* p2) {
  T d = -DotProduct(nDim, n, p0);
  auto a = DotProduct(nDim, n, p1) + d;
  auto b = DotProduct(nDim, n, p2) + d;
  if (a * b <= 0 && (a != 0 || b != 0))
    return true;
  else
    return false;
}

/*!
 * \brief Check if a point is inside a polygon
 * \param[in] nDim - Number of dimensions of the particular problem.
 * \param[in] vector - Vector to normalize.
 */
template <class T, class Int>
inline void NormalizeVector(Int nDim, T* vector) {
  T norm = Norm(nDim, vector);
  for (Int i = 0; i < nDim; i++) vector[i] /= norm;
}

/*!
 * \brief Check if a point is inside a polygon
 * \param[in] nDim - Number of dimensions of the particular problem.
 * \param[in] pVert - Point defining the polygon.
 * \param[in] p0 - Point to test
 * \param[in] nVert - Number of points that form the polygon.
 * \return Returns TRUE if the edge intersect the plane, FALSE it doesn't.
 */
template <class Mat, class T, class Int>
inline bool PointInConvexPolygon(Int nDim, const Mat& pVert, const T* p0, int nVert) {
  unsigned short i = 0, j = 1;
  unsigned short idxPoint1 = 0, idxPoint2 = 1;
  unsigned short nIntersections = 0;
  unsigned short proj_index = 0;
  T max_proj = 0;

  /* Check which of the x,y,z planes is the most parallel to the polygon */

  if (nDim == 3) {
    T plane_norm[3];
    TriangleNormal(pVert, plane_norm);
    NormalizeVector(3, plane_norm);

    for (unsigned short iDim = 0; iDim < nDim; iDim++) {
      if (abs(plane_norm[iDim]) > max_proj) {
        proj_index = iDim;
        max_proj = abs(plane_norm[iDim]);
      }
    }

    switch (proj_index) {
      case 0:
        i = 1;
        j = 2;
        break;
      case 1:
        i = 2;
        j = 0;
        break;
      default:
        i = 0;
        j = 1;
        break;
    }
  }

  /* Loop across the polygon edges and check if the ray intersect a vertex */

  for (unsigned short iVert = 1; iVert < nVert + 1; iVert++) {
    idxPoint2 = iVert % nVert;

    if (((p0[j] < pVert[idxPoint1][j]) != (p0[j] < pVert[idxPoint2][j])) &&
        (p0[i] < pVert[idxPoint1][i] + ((p0[j] - pVert[idxPoint1][j]) / (pVert[idxPoint2][j] - pVert[idxPoint1][j])) *
                                           (pVert[idxPoint2][i] - pVert[idxPoint1][i]))) {
      nIntersections++;
    }
    idxPoint1 = idxPoint2;
  }

  return nIntersections % 2 == 1;
}

// end added by max
/// @}
}  // namespace GeometryToolbox
