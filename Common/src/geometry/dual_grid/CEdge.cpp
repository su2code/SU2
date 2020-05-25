/*!
 * \file CEdge.cpp
 * \brief Implementation of the edge class.
 * \author F. Palacios, T. Economon
 * \version 7.0.4 "Blackbird"
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

#include "../../../include/geometry/dual_grid/CEdge.hpp"
#include "../../../include/toolboxes/geometry_toolbox.hpp"

using namespace GeometryToolbox;


CEdge::CEdge(unsigned long nEdge, unsigned long nDim) :
  Nodes(nEdge,2), Normal(nEdge,nDim), Coord_CG(nEdge,nDim) {
  Normal = su2double(0.0);
  Coord_CG = su2double(0.0);
}

void CEdge::SetZeroValues(void) {
  Normal = su2double(0.0);
}

su2double CEdge::GetVolume(const su2double *coord_Edge_CG,
                           const su2double *coord_FaceElem_CG,
                           const su2double *coord_Elem_CG,
                           const su2double *coord_Point) {

  constexpr unsigned long nDim = 3;

  su2double vec_a[nDim] = {0.0}, vec_b[nDim] = {0.0}, vec_c[nDim] = {0.0}, vec_d[nDim] = {0.0};

  AD::StartPreacc();
  AD::SetPreaccIn(coord_Edge_CG, nDim);
  AD::SetPreaccIn(coord_Elem_CG, nDim);
  AD::SetPreaccIn(coord_FaceElem_CG, nDim);
  AD::SetPreaccIn(coord_Point, nDim);

  Distance(nDim, coord_Edge_CG,     coord_Point, vec_a);
  Distance(nDim, coord_FaceElem_CG, coord_Point, vec_b);
  Distance(nDim, coord_Elem_CG,     coord_Point, vec_c);

  CrossProduct(vec_a, vec_b, vec_d);

  su2double Local_Volume = fabs(DotProduct(nDim, vec_c, vec_d)) / 6.0;

  AD::SetPreaccOut(Local_Volume);
  AD::EndPreacc();

  return Local_Volume;
}

su2double CEdge::GetVolume(const su2double *coord_Edge_CG,
                           const su2double *coord_Elem_CG,
                           const su2double *coord_Point) {

  constexpr unsigned long nDim = 2;

  su2double vec_a[nDim] = {0.0}, vec_b[nDim] = {0.0};

  AD::StartPreacc();
  AD::SetPreaccIn(coord_Edge_CG, nDim);
  AD::SetPreaccIn(coord_Elem_CG, nDim);
  AD::SetPreaccIn(coord_Point, nDim);

  Distance(nDim, coord_Elem_CG, coord_Point, vec_a);
  Distance(nDim, coord_Edge_CG, coord_Point, vec_b);

  su2double Local_Volume = 0.5 * fabs(vec_a[0]*vec_b[1] - vec_a[1]*vec_b[0]);

  AD::SetPreaccOut(Local_Volume);
  AD::EndPreacc();

  return Local_Volume;
}

void CEdge::SetNodes_Coord(unsigned long iEdge,
                           const su2double *coord_Edge_CG,
                           const su2double *coord_FaceElem_CG,
                           const su2double *coord_Elem_CG) {

  constexpr unsigned long nDim = 3;

  su2double vec_a[nDim] = {0.0}, vec_b[nDim] = {0.0}, Dim_Normal[nDim];

  AD::StartPreacc();
  AD::SetPreaccIn(coord_Edge_CG, nDim);
  AD::SetPreaccIn(coord_Elem_CG, nDim);
  AD::SetPreaccIn(coord_FaceElem_CG, nDim);
  AD::SetPreaccIn(Normal[iEdge], nDim);

  Distance(nDim, coord_Elem_CG, coord_Edge_CG, vec_a);
  Distance(nDim, coord_FaceElem_CG, coord_Edge_CG, vec_b);

  CrossProduct(vec_a, vec_b, Dim_Normal);

  for (auto iDim = 0ul; iDim < nDim; ++iDim)
    Normal(iEdge,iDim) += 0.5 * Dim_Normal[iDim];

  AD::SetPreaccOut(Normal[iEdge], nDim);
  AD::EndPreacc();
}

void CEdge::SetNodes_Coord(unsigned long iEdge,
                           const su2double *coord_Edge_CG,
                           const su2double *coord_Elem_CG) {

  constexpr unsigned long nDim = 2;

  AD::StartPreacc();
  AD::SetPreaccIn(coord_Elem_CG, nDim);
  AD::SetPreaccIn(coord_Edge_CG, nDim);
  AD::SetPreaccIn(Normal[iEdge], nDim);

  Normal(iEdge,0) += coord_Elem_CG[1] - coord_Edge_CG[1];
  Normal(iEdge,1) -= coord_Elem_CG[0] - coord_Edge_CG[0];

  AD::SetPreaccOut(Normal[iEdge], nDim);
  AD::EndPreacc();

}
