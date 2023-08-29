/*!
 * \file CTRIA3.cpp
 * \brief Definition of the 3-node triangular element with three Gauss points.
 * \author T. Dick
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

#include "../../../include/geometry/elements/CElement.hpp"

CTRIA3::CTRIA3() : CElementWithKnownSizes<NGAUSS, NNODE, NDIM>() {
  /*--- Gauss coordinates and weights ---*/

  GaussCoord[0][0] = 0.66666666666666666667;
  GaussCoord[0][1] = 0.16666666666666666667;
  GaussWeight(0) = 0.33333333333333333333;
  GaussCoord[1][0] = 0.16666666666666666667;
  GaussCoord[1][1] = 0.66666666666666666667;
  GaussWeight(1) = 0.33333333333333333333;
  GaussCoord[2][0] = 0.16666666666666666667;
  GaussCoord[2][1] = 0.16666666666666666667;
  GaussWeight(2) = 0.33333333333333333333;

  /*--- Store the values of the shape functions and their derivatives ---*/

  unsigned short iGauss;
  su2double Xi, Eta, val_Ni;

  for (iGauss = 0; iGauss < nGaussPoints; iGauss++) {
    Xi = GaussCoord[iGauss][0];
    Eta = GaussCoord[iGauss][1];

    val_Ni = 1 - Xi - Eta;
    GaussPoint[iGauss].SetNi(val_Ni, 0);
    val_Ni = Xi;
    GaussPoint[iGauss].SetNi(val_Ni, 1);
    val_Ni = Eta;
    GaussPoint[iGauss].SetNi(val_Ni, 2);

    /*--- dN/d xi, dN/d eta ---*/

    dNiXj[iGauss][0][0] = -1.0;
    dNiXj[iGauss][0][1] = -1.0;
    dNiXj[iGauss][1][0] = 1.0;
    dNiXj[iGauss][1][1] = 0.0;
    dNiXj[iGauss][2][0] = 0.0;
    dNiXj[iGauss][2][1] = 1.0;
  }

  /*--- Shape functions evaluated at the nodes for extrapolation of the stresses at the Gaussian Points ---*/
  /*--- The stress is constant over a TRIA1 element ---*/

  NodalExtrap[0][0] = 1.0;
  NodalExtrap[1][0] = 1.0;
  NodalExtrap[2][0] = 1.0;
}

su2double CTRIA3::ComputeArea(const FrameType mode) const {
  unsigned short iDim;
  su2double a[3] = {0.0, 0.0, 0.0}, b[3] = {0.0, 0.0, 0.0};
  su2double Area = 0.0;

  /*--- Select the appropriate source for the nodal coordinates depending on the frame requested
        for the gradient computation, REFERENCE (undeformed) or CURRENT (deformed) ---*/
  const su2activematrix& Coord = (mode == REFERENCE) ? RefCoord : CurrentCoord;

  for (iDim = 0; iDim < NDIM; iDim++) {
    a[iDim] = Coord[0][iDim] - Coord[2][iDim];
    b[iDim] = Coord[1][iDim] - Coord[2][iDim];
  }

  Area = 0.5 * fabs(a[0] * b[1] - a[1] * b[0]);

  return Area;
}
