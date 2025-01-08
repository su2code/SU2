/*!
 * \file CTETRA4.cpp
 * \brief Definition of 4-node tetrahedral element with 4 Gauss point.
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
#include "../../../include/toolboxes/geometry_toolbox.hpp"

CTETRA4::CTETRA4() : CElementWithKnownSizes<NGAUSS, NNODE, NDIM>() {
  /*--- Gauss coordinates and weights ---*/

  su2double r = ((5.0 - sqrt(5.0)) / 20);
  su2double s = ((5.0 + 3 * sqrt(5.0)) / 20);
  GaussCoord[0][0] = r;
  GaussCoord[0][1] = r;
  GaussCoord[0][2] = r;
  GaussWeight(0) = 1.0 / 24.0;
  GaussCoord[0][0] = s;
  GaussCoord[0][1] = r;
  GaussCoord[0][2] = r;
  GaussWeight(1) = 1.0 / 24.0;
  GaussCoord[0][0] = r;
  GaussCoord[0][1] = s;
  GaussCoord[0][2] = r;
  GaussWeight(2) = 1.0 / 24.0;
  GaussCoord[0][0] = r;
  GaussCoord[0][1] = r;
  GaussCoord[0][2] = s;
  GaussWeight(3) = 1.0 / 24.0;

  /*--- Store the values of the shape functions and their derivatives ---*/

  unsigned short iGauss;
  su2double Xi, Eta, Zeta, val_Ni;

  for (iGauss = 0; iGauss < nGaussPoints; iGauss++) {
    Xi = GaussCoord[iGauss][0];
    Eta = GaussCoord[iGauss][1];
    Zeta = GaussCoord[iGauss][2];

    val_Ni = Xi;
    GaussPoint[iGauss].SetNi(val_Ni, 0);
    val_Ni = Eta;
    GaussPoint[iGauss].SetNi(val_Ni, 1);
    val_Ni = 1.0 - Xi - Eta - Zeta;
    GaussPoint[iGauss].SetNi(val_Ni, 2);
    val_Ni = Zeta;
    GaussPoint[iGauss].SetNi(val_Ni, 3);

    /*--- dN/d xi, dN/d eta, dN/d zeta ---*/

    dNiXj[iGauss][0][0] = 1.0;
    dNiXj[iGauss][0][1] = 0.0;
    dNiXj[iGauss][0][2] = 0.0;
    dNiXj[iGauss][1][0] = 0.0;
    dNiXj[iGauss][1][1] = 1.0;
    dNiXj[iGauss][1][2] = 0.0;
    dNiXj[iGauss][2][0] = -1.0;
    dNiXj[iGauss][2][1] = -1.0;
    dNiXj[iGauss][2][2] = -1.0;
    dNiXj[iGauss][3][0] = 0.0;
    dNiXj[iGauss][3][1] = 0.0;
    dNiXj[iGauss][3][2] = 1.0;
  }

  /*--- Shape functions evaluated at the nodes for extrapolation of the stresses at the Gaussian Points ---*/
  /*--- TODO: Use correct extrapolation, currently assumes stress is constant similar to TETRA1 element ---*/

  NodalExtrap[0][0] = 1.0;
  NodalExtrap[1][0] = 1.0;
  NodalExtrap[2][0] = 1.0;
  NodalExtrap[3][0] = 1.0;
}

su2double CTETRA4::ComputeVolume(const FrameType mode) const {
  unsigned short iDim;
  su2double r1[3] = {0.0, 0.0, 0.0}, r2[3] = {0.0, 0.0, 0.0}, r3[3] = {0.0, 0.0, 0.0},
            CrossProduct[3] = {0.0, 0.0, 0.0};
  su2double Volume = 0.0;

  /*--- Select the appropriate source for the nodal coordinates depending on the frame requested
        for the gradient computation, REFERENCE (undeformed) or CURRENT (deformed)---*/
  const su2activematrix& Coord = (mode == REFERENCE) ? RefCoord : CurrentCoord;

  for (iDim = 0; iDim < NDIM; iDim++) {
    r1[iDim] = Coord[1][iDim] - Coord[0][iDim];
    r2[iDim] = Coord[2][iDim] - Coord[0][iDim];
    r3[iDim] = Coord[3][iDim] - Coord[0][iDim];
  }

  GeometryToolbox::CrossProduct(r1, r2, CrossProduct);
  Volume = fabs(GeometryToolbox::DotProduct(3, CrossProduct, r3)) / 6.0;

  return Volume;
}
