/*!
 * \file CEdge.cpp
 * \brief Implementation of the edge class.
 * \author F. Palacios, T. Economon
 * \version 8.3.0 "Harrier"
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

#include "../../../include/geometry/dual_grid/CEdge.hpp"
#include "../../../include/toolboxes/geometry_toolbox.hpp"
#include "../../../include/parallelization/omp_structure.hpp"

using namespace GeometryToolbox;

CEdge::CEdge(unsigned long nEdge_, unsigned long nDim)
    : nEdge(nEdge_), nEdgeSIMD(nextMultiple(nEdge_, simd::preferredLen<su2double>())) {
  /*--- Allocate with padding. ---*/
  Nodes.resize(nEdgeSIMD, 2) = 0;
  Normal.resize(nEdgeSIMD, nDim) = su2double(0.0);

  a.resize(nEdgeSIMD, nDim) = su2double(0.0);
  b.resize(nEdgeSIMD, nDim) = su2double(0.0);
  c.resize(nEdgeSIMD, nDim) = su2double(0.0);
}

void CEdge::SetZeroValues() { Normal = su2double(0.0); a = su2double(0.0);  b = su2double(0.0);  c = su2double(0.0);}

su2double CEdge::GetVolume(const su2double* coord_Edge_CG, const su2double* coord_FaceElem_CG,
                           const su2double* coord_Elem_CG, const su2double* coord_Point) {
  constexpr unsigned long nDim = 3;

  su2double vec_a[nDim] = {0.0}, vec_b[nDim] = {0.0}, vec_c[nDim] = {0.0}, vec_d[nDim] = {0.0};

  Distance(nDim, coord_Edge_CG, coord_Point, vec_a);
  Distance(nDim, coord_FaceElem_CG, coord_Point, vec_b);
  Distance(nDim, coord_Elem_CG, coord_Point, vec_c);

  CrossProduct(vec_a, vec_b, vec_d);

  return fabs(DotProduct(nDim, vec_c, vec_d)) / 6.0;
}

su2double CEdge::GetVolume(const su2double* coord_Edge_CG, const su2double* coord_Elem_CG,
                           const su2double* coord_Point) {
  constexpr unsigned long nDim = 2;

  su2double vec_a[nDim] = {0.0}, vec_b[nDim] = {0.0};

  Distance(nDim, coord_Elem_CG, coord_Point, vec_a);
  Distance(nDim, coord_Edge_CG, coord_Point, vec_b);

  return 0.5 * fabs(vec_a[0] * vec_b[1] - vec_a[1] * vec_b[0]);
}

void CEdge::SetNodes_Coord(unsigned long iEdge, const su2double* coord_Edge_CG, const su2double* coord_FaceElem_CG,
                           const su2double* coord_Elem_CG) {
  constexpr unsigned long nDim = 3;

  su2double vec_a[nDim] = {0.0}, vec_b[nDim] = {0.0}, Dim_Normal[nDim];

  Distance(nDim, coord_Elem_CG, coord_Edge_CG, vec_a);
  Distance(nDim, coord_FaceElem_CG, coord_Edge_CG, vec_b);

  CrossProduct(vec_a, vec_b, Dim_Normal);

  for (auto iDim = 0ul; iDim < nDim; ++iDim) Normal(iEdge, iDim) += 0.5 * Dim_Normal[iDim];
}

void CEdge::SetNodes_CoordCorrection(unsigned long iEdge, const su2double* coord_Edge_CG, const su2double* coord_FaceElem_CG,
                           const su2double* coord_Elem_CG, const su2double* quadraturePoint) {
  constexpr unsigned long nDim = 3;

  su2double vec_a[nDim] = {0.0}, vec_b[nDim] = {0.0}, NT[nDim] = {0.0};

  Distance(nDim, coord_Elem_CG, coord_Edge_CG, vec_a);
  Distance(nDim, coord_FaceElem_CG, coord_Edge_CG, vec_b);

  CrossProduct(vec_a, vec_b, NT);

  for (auto iDim = 0ul; iDim < nDim; ++iDim) NT[iDim] *= 0.5;

  su2double CG[nDim];
  for (int iDim = 0; iDim < nDim; ++iDim)
    CG[iDim] = (coord_Edge_CG[iDim] + coord_FaceElem_CG[iDim] + coord_Elem_CG[iDim]) / 3.0;

  for (int iDim = 0; iDim < nDim; ++iDim) {
    a(iEdge,iDim) += NT[iDim] * (CG[0] - quadraturePoint[0]);
    b(iEdge,iDim) += NT[iDim] * (CG[1] - quadraturePoint[1]);
    c(iEdge,iDim) += NT[iDim] * (CG[2] - quadraturePoint[2]);
  }

}

void CEdge::SetNodes_Coord(unsigned long iEdge, const su2double* coord_Edge_CG, const su2double* coord_Elem_CG) {
  Normal(iEdge, 0) += coord_Elem_CG[1] - coord_Edge_CG[1];
  Normal(iEdge, 1) -= coord_Elem_CG[0] - coord_Edge_CG[0];
}
void CEdge::SetNodes_CoordCorrection(unsigned long iEdge, const su2double* coord_Edge_CG,
                                     const su2double* coord_Elem_CG, const su2double* quadraturePoint) {

  const su2double CG[2] = {0.5*(coord_Elem_CG[0] + coord_Edge_CG[0]), 0.5*(coord_Elem_CG[1] + coord_Edge_CG[1])};
  const su2double NT[2] = {coord_Elem_CG[1] - coord_Edge_CG[1], coord_Edge_CG[0] - coord_Elem_CG[0]};
//  su2double S[2][2] = {0.0};
//
//  for (int iDim = 0; iDim < 2; ++iDim) {
//    for (int jDim = 0; jDim < 2; ++jDim) {
//      S[iDim][jDim] = NT[iDim] * (CG[jDim] - quadraturePoint[jDim]);
//    }
//  }

  for (int iDim = 0; iDim < 2; ++iDim) {
    a(iEdge,iDim) += NT[iDim] * (CG[0] - quadraturePoint[0]);
    b(iEdge,iDim) += NT[iDim] * (CG[1] - quadraturePoint[1]);
//    c(iEdge,jDim) += NT[2] * (CG[jDim] - quadraturePoint[jDim]);
  }

}
