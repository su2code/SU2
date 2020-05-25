/*!
 * \file CEdge.cpp
 * \brief Main classes for defining the edges of the dual grid
 * \author F. Palacios, T. Economon
 * \version 7.0.1 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2019, SU2 Contributors (cf. AUTHORS.md)
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

CEdge::CEdge(unsigned long val_iPoint, unsigned long val_jPoint, unsigned short val_nDim) : CDualGrid(val_nDim) {

  unsigned short iDim;

  /*--- Pointers initialization ---*/
  Coord_CG = NULL;
  Normal   = NULL;
  Nodes    = NULL;

  /*--- Allocate center of gravity coordinates, nodes, and face normal ---*/
  Coord_CG = new su2double [nDim];
  Normal   = new su2double [nDim];
  Nodes    = new unsigned long[2];

  /*--- Initializate the structure ---*/
  for (iDim = 0; iDim < nDim; iDim++) {
    Coord_CG[iDim] = 0.0;
    Normal[iDim]   = 0.0;
  }

  Nodes[0] = val_iPoint;
  Nodes[1] = val_jPoint;

}

CEdge::~CEdge() {

  if (Coord_CG != NULL) delete[] Coord_CG;
  if (Normal   != NULL) delete[] Normal;
  if (Nodes    != NULL) delete[] Nodes;

}

void CEdge::SetCoord_CG(su2double **val_coord) {

  unsigned short iDim, iNode;

  for (iDim = 0; iDim < nDim; iDim++) {
    Coord_CG[iDim] = 0.0;
    for (iNode = 0; iNode < 2;  iNode++)
      Coord_CG[iDim] += val_coord[iNode][iDim] / 2.0;
  }

}

su2double CEdge::GetVolume(su2double *val_coord_Edge_CG, su2double *val_coord_FaceElem_CG, su2double *val_coord_Elem_CG, su2double *val_coord_Point) const {

  unsigned short iDim;
  su2double vec_a[3] = {0.0,0.0,0.0}, vec_b[3] = {0.0,0.0,0.0}, vec_c[3] = {0.0,0.0,0.0}, vec_d[3] = {0.0,0.0,0.0}, Local_Volume;

  AD::StartPreacc();
  AD::SetPreaccIn(val_coord_Edge_CG, nDim);
  AD::SetPreaccIn(val_coord_Elem_CG, nDim);
  AD::SetPreaccIn(val_coord_FaceElem_CG, nDim);
  AD::SetPreaccIn(val_coord_Point, nDim);

  for (iDim = 0; iDim < nDim; iDim++) {
    vec_a[iDim] = val_coord_Edge_CG[iDim]     - val_coord_Point[iDim];
    vec_b[iDim] = val_coord_FaceElem_CG[iDim] - val_coord_Point[iDim];
    vec_c[iDim] = val_coord_Elem_CG[iDim]     - val_coord_Point[iDim];
  }

  vec_d[0] =   vec_a[1] * vec_b[2] - vec_a[2] * vec_b[1];
  vec_d[1] = -(vec_a[0] * vec_b[2] - vec_a[2] * vec_b[0]);
  vec_d[2] =   vec_a[0] * vec_b[1] - vec_a[1] * vec_b[0];

  Local_Volume = fabs( vec_c[0] * vec_d[0] + vec_c[1] * vec_d[1] + vec_c[2] * vec_d[2] ) / 6.0;

  AD::SetPreaccOut(Local_Volume);
  AD::EndPreacc();

  return Local_Volume;

}

su2double CEdge::GetVolume(su2double *val_coord_Edge_CG, su2double *val_coord_Elem_CG, su2double *val_coord_Point) const {

  unsigned short iDim;
  su2double vec_a[2] = {0.0,0.0}, vec_b[2] = {0.0,0.0}, Local_Volume;

  AD::StartPreacc();
  AD::SetPreaccIn(val_coord_Edge_CG, nDim);
  AD::SetPreaccIn(val_coord_Elem_CG, nDim);
  AD::SetPreaccIn(val_coord_Point, nDim);

  for (iDim = 0; iDim < nDim; iDim++) {
    vec_a[iDim] = val_coord_Elem_CG[iDim] - val_coord_Point[iDim];
    vec_b[iDim] = val_coord_Edge_CG[iDim] - val_coord_Point[iDim];
  }

  Local_Volume = 0.5 * fabs( vec_a[0] * vec_b[1] - vec_a[1] * vec_b[0] );

  AD::SetPreaccOut(Local_Volume);
  AD::EndPreacc();

  return Local_Volume;

}

void CEdge::SetNodes_Coord(su2double *val_coord_Edge_CG, su2double *val_coord_FaceElem_CG, su2double *val_coord_Elem_CG) {

  unsigned short iDim;
  su2double vec_a[3] = {0.0,0.0,0.0}, vec_b[3] = {0.0,0.0,0.0}, Dim_Normal[3];

  AD::StartPreacc();
  AD::SetPreaccIn(val_coord_Edge_CG, nDim);
  AD::SetPreaccIn(val_coord_Elem_CG, nDim);
  AD::SetPreaccIn(val_coord_FaceElem_CG, nDim);
  AD::SetPreaccIn(Normal, nDim);

  for (iDim = 0; iDim < nDim; iDim++) {
    vec_a[iDim] = val_coord_Elem_CG[iDim]-val_coord_Edge_CG[iDim];
    vec_b[iDim] = val_coord_FaceElem_CG[iDim]-val_coord_Edge_CG[iDim];
  }

  Dim_Normal[0] =  0.5 * ( vec_a[1] * vec_b[2] - vec_a[2] * vec_b[1] );
  Dim_Normal[1] = -0.5 * ( vec_a[0] * vec_b[2] - vec_a[2] * vec_b[0] );
  Dim_Normal[2] =  0.5 * ( vec_a[0] * vec_b[1] - vec_a[1] * vec_b[0] );

  Normal[0] += Dim_Normal[0];
  Normal[1] += Dim_Normal[1];
  Normal[2] += Dim_Normal[2];

  AD::SetPreaccOut(Normal, nDim);
  AD::EndPreacc();
}

void CEdge::SetNodes_Coord(su2double *val_coord_Edge_CG, su2double *val_coord_Elem_CG) {

  su2double Dim_Normal[2];

  AD::StartPreacc();
  AD::SetPreaccIn(val_coord_Elem_CG, nDim);
  AD::SetPreaccIn(val_coord_Edge_CG, nDim);
  AD::SetPreaccIn(Normal, nDim);

  Dim_Normal[0] =   val_coord_Elem_CG[1] - val_coord_Edge_CG[1];
  Dim_Normal[1] = -(val_coord_Elem_CG[0] - val_coord_Edge_CG[0]);

  Normal[0] += Dim_Normal[0];
  Normal[1] += Dim_Normal[1];

  AD::SetPreaccOut(Normal, nDim);
  AD::EndPreacc();

}
