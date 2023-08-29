/*!
 * \file CVertex.cpp
 * \brief Main classes for defining the vertices of the dual grid
 * \author F. Palacios, T. Economon
 * \version 8.0.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2023, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../../include/geometry/dual_grid/CVertex.hpp"
#include "../../../include/toolboxes/geometry_toolbox.hpp"

using namespace GeometryToolbox;

CVertex::CVertex(unsigned long val_point, unsigned short val_nDim) : CDualGrid(val_nDim) { Nodes[0] = val_point; }

void CVertex::SetNodes_Coord(const su2double* coord_Edge_CG, const su2double* coord_FaceElem_CG,
                             const su2double* coord_Elem_CG) {
  constexpr unsigned long nDim = 3;
  su2double vec_a[nDim] = {0.0}, vec_b[nDim] = {0.0}, Dim_Normal[nDim];

  Distance(nDim, coord_Elem_CG, coord_Edge_CG, vec_a);
  Distance(nDim, coord_FaceElem_CG, coord_Edge_CG, vec_b);

  CrossProduct(vec_a, vec_b, Dim_Normal);

  for (auto iDim = 0ul; iDim < nDim; ++iDim) Normal[iDim] += 0.5 * Dim_Normal[iDim];
}

void CVertex::SetNodes_Coord(const su2double* val_coord_Edge_CG, const su2double* val_coord_Elem_CG) {
  Normal[0] += val_coord_Elem_CG[1] - val_coord_Edge_CG[1];
  Normal[1] -= val_coord_Elem_CG[0] - val_coord_Edge_CG[0];
}
