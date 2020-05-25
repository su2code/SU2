/*!
 * \file CVertex.cpp
 * \brief Main classes for defining the vertices of the dual grid
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

#include "../../../include/geometry/dual_grid/CVertex.hpp"

CVertex::CVertex(unsigned long val_point, unsigned short val_nDim) : CDualGrid(val_nDim) {

  unsigned short iDim;

  /*--- Set periodic points to zero ---*/

  PeriodicPoint[0] = -1; PeriodicPoint[1] = -1; PeriodicPoint[2] = -1;
  PeriodicPoint[3] = -1; PeriodicPoint[4] = -1;

  /*--- Identify the points at the perimeter of the actuatrod disk ---*/

  ActDisk_Perimeter = false;

  /*--- Initializate the structure ---*/

  Nodes[0] = val_point;

  for (iDim = 0; iDim < 3; iDim ++) Normal[iDim] = 0.0;

  /*--- Set to zero the variation of the coordinates ---*/

  for (iDim = 0; iDim < 3; iDim ++) VarCoord[iDim] = 0.0;

  /*--- Set to nullptr variation of the rotation  ---*/

  VarRot = nullptr;

  /*--- Set to nullptr donor arrays for interpolation ---*/

  Donor_Points  = nullptr;
  Donor_Proc    = nullptr;
  Donor_Coeff   = nullptr;
  nDonor_Points = 1;

}

CVertex::~CVertex() {

  delete[] VarRot;
  delete[] Donor_Coeff;
  delete[] Donor_Proc;
  delete[] Donor_Points;

}

void CVertex::SetNodes_Coord(su2double *val_coord_Edge_CG, su2double *val_coord_FaceElem_CG, su2double *val_coord_Elem_CG) {

  su2double vec_a[3] = {0.0,0.0,0.0}, vec_b[3] = {0.0,0.0,0.0};
  unsigned short iDim;

  assert(nDim == 3);

  AD::StartPreacc();
  AD::SetPreaccIn(val_coord_Edge_CG, nDim);
  AD::SetPreaccIn(val_coord_Elem_CG, nDim);
  AD::SetPreaccIn(val_coord_FaceElem_CG, nDim);
  AD::SetPreaccIn(Normal, nDim);

  for (iDim = 0; iDim < nDim; iDim++) {
    vec_a[iDim] = val_coord_Elem_CG[iDim]-val_coord_Edge_CG[iDim];
    vec_b[iDim] = val_coord_FaceElem_CG[iDim]-val_coord_Edge_CG[iDim];
  }

  Normal[0] += 0.5 * ( vec_a[1] * vec_b[2] - vec_a[2] * vec_b[1]);
  Normal[1] -= 0.5 * ( vec_a[0] * vec_b[2] - vec_a[2] * vec_b[0]);
  Normal[2] += 0.5 * ( vec_a[0] * vec_b[1] - vec_a[1] * vec_b[0]);

  AD::SetPreaccOut(Normal, nDim);
  AD::EndPreacc();

}

void CVertex::SetNodes_Coord(su2double *val_coord_Edge_CG, su2double *val_coord_Elem_CG) {

  AD::StartPreacc();
  AD::SetPreaccIn(val_coord_Elem_CG, nDim);
  AD::SetPreaccIn(val_coord_Edge_CG, nDim);
  AD::SetPreaccIn(Normal, nDim);

  Normal[0] += val_coord_Elem_CG[1]-val_coord_Edge_CG[1];
  Normal[1] -= (val_coord_Elem_CG[0]-val_coord_Edge_CG[0]);

  AD::SetPreaccOut(Normal, nDim);
  AD::EndPreacc();

}

void CVertex::Allocate_DonorInfo(unsigned short nDonor) {

  nDonor_Points = nDonor;

  delete [] Donor_Points;
  delete [] Donor_Proc;
  delete [] Donor_Coeff;

  Donor_Points = new unsigned long [nDonor_Points];
  Donor_Proc   = new unsigned long [nDonor_Points];
  Donor_Coeff  = new su2double [nDonor_Points];

}
