/*!
 * \file geometry_toolbox.hpp
 * \brief Collection of common lightweight geometry-oriented methods.
 * \version 7.0.3 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2020, SU2 Contributors (cf. AUTHORS.md)
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

namespace GeometryToolbox {

/*! \return ||a-b||^2 */
template<class T, typename Int>
inline T SquaredDistance(Int nDim, const T* a, const T* b) {
  T d(0);
  for(Int i = 0; i < nDim; i++) d += pow(a[i]-b[i], 2);
  return d;
}

/*! \return ||a-b|| */
template<class T, typename Int>
inline T Distance(Int nDim, const T* a, const T* b) {
  return sqrt(SquaredDistance(nDim, a, b));
}

/*! \brief d = a-b */
template<class T, typename Int>
inline void Distance(Int nDim, const T* a, const T* b, T* d) {
  for(Int i = 0; i < nDim; i++) d[i] = a[i] - b[i];
}

/*! \return a.b */
template<class T, typename Int>
inline T DotProduct(Int nDim, const T* a, const T* b) {
  T d(0);
  for (Int i = 0; i < nDim; ++i) d += a[i]*b[i];
  return d;
}

/*! \return ||a||^2 */
template<class T, typename Int>
inline T SquaredNorm(Int nDim, const T* a) {
  return DotProduct(nDim, a, a);
}

/*! \return ||a|| */
template<class T, typename Int>
inline T Norm(Int nDim, const T* a) {
  return sqrt(SquaredNorm(nDim, a));
}

/*! \brief c = a x b */
template<class T>
inline void CrossProduct(const T* a, const T* b, T* c) {
  c[0] = a[1]*b[2] - a[2]*b[1];
  c[1] = a[2]*b[0] - a[0]*b[2];
  c[2] = a[0]*b[1] - a[1]*b[0];
}

/*! \brief Set U as the normal to a 2D line defined by coords[iPoint][iDim]. */
template<class T, class U>
inline void LineNormal(const T& coords, U* normal) {
  normal[0] = coords[0][1] - coords[1][1];
  normal[1] = coords[1][0] - coords[0][0];
}

/*! \brief Normal vector of a triangle, cross product of two sides. */
template<class T, class U>
inline void TriangleNormal(const T& coords, U* normal) {

  U a[3], b[3];

  for (int iDim = 0; iDim < 3; iDim++) {
    a[iDim] = coords[1][iDim] - coords[0][iDim];
    b[iDim] = coords[2][iDim] - coords[0][iDim];
  }

  CrossProduct(a, b, normal);
  normal[0] *= 0.5; normal[1] *= 0.5; normal[2] *= 0.5;
}

/*! \brief Normal vector of a quadrilateral, cross product of the two diagonals. */
template<class T, class U>
inline void QuadrilateralNormal(const T& coords, U* normal) {

  U a[3], b[3];

  for (int iDim = 0; iDim < 3; iDim++) {
    a[iDim] = coords[2][iDim] - coords[0][iDim];
    b[iDim] = coords[3][iDim] - coords[1][iDim];
  }

  CrossProduct(a, b, normal);
  normal[0] *= 0.5; normal[1] *= 0.5; normal[2] *= 0.5;
}

}
