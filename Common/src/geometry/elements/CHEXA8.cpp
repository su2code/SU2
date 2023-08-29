/*!
 * \file CHEXA8.cpp
 * \brief Definition of 8-node hexa element with 8 Gauss points.
 * \author R. Sanchez
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

CHEXA8::CHEXA8() : CElementWithKnownSizes<NGAUSS, NNODE, NDIM>() {
  /*--- Gauss coordinates and weights ---*/

  su2double oneOnSqrt3 = 0.577350269189626;

  GaussCoord[0][0] = -oneOnSqrt3;
  GaussCoord[0][1] = -oneOnSqrt3;
  GaussCoord[0][2] = -oneOnSqrt3;
  GaussWeight(0) = 1.0;
  GaussCoord[1][0] = oneOnSqrt3;
  GaussCoord[1][1] = -oneOnSqrt3;
  GaussCoord[1][2] = -oneOnSqrt3;
  GaussWeight(1) = 1.0;
  GaussCoord[2][0] = oneOnSqrt3;
  GaussCoord[2][1] = oneOnSqrt3;
  GaussCoord[2][2] = -oneOnSqrt3;
  GaussWeight(2) = 1.0;
  GaussCoord[3][0] = -oneOnSqrt3;
  GaussCoord[3][1] = oneOnSqrt3;
  GaussCoord[3][2] = -oneOnSqrt3;
  GaussWeight(3) = 1.0;
  GaussCoord[4][0] = -oneOnSqrt3;
  GaussCoord[4][1] = -oneOnSqrt3;
  GaussCoord[4][2] = oneOnSqrt3;
  GaussWeight(4) = 1.0;
  GaussCoord[5][0] = oneOnSqrt3;
  GaussCoord[5][1] = -oneOnSqrt3;
  GaussCoord[5][2] = oneOnSqrt3;
  GaussWeight(5) = 1.0;
  GaussCoord[6][0] = oneOnSqrt3;
  GaussCoord[6][1] = oneOnSqrt3;
  GaussCoord[6][2] = oneOnSqrt3;
  GaussWeight(6) = 1.0;
  GaussCoord[7][0] = -oneOnSqrt3;
  GaussCoord[7][1] = oneOnSqrt3;
  GaussCoord[7][2] = oneOnSqrt3;
  GaussWeight(7) = 1.0;

  /*--- Store the values of the shape functions and their derivatives ---*/

  unsigned short iNode, iGauss;
  su2double Xi, Eta, Zeta, val_Ni;

  for (iGauss = 0; iGauss < NGAUSS; iGauss++) {
    Xi = GaussCoord[iGauss][0];
    Eta = GaussCoord[iGauss][1];
    Zeta = GaussCoord[iGauss][2];

    val_Ni = 0.125 * (1.0 - Xi) * (1.0 - Eta) * (1.0 - Zeta);
    GaussPoint[iGauss].SetNi(val_Ni, 0);
    val_Ni = 0.125 * (1.0 + Xi) * (1.0 - Eta) * (1.0 - Zeta);
    GaussPoint[iGauss].SetNi(val_Ni, 1);
    val_Ni = 0.125 * (1.0 + Xi) * (1.0 + Eta) * (1.0 - Zeta);
    GaussPoint[iGauss].SetNi(val_Ni, 2);
    val_Ni = 0.125 * (1.0 - Xi) * (1.0 + Eta) * (1.0 - Zeta);
    GaussPoint[iGauss].SetNi(val_Ni, 3);
    val_Ni = 0.125 * (1.0 - Xi) * (1.0 - Eta) * (1.0 + Zeta);
    GaussPoint[iGauss].SetNi(val_Ni, 4);
    val_Ni = 0.125 * (1.0 + Xi) * (1.0 - Eta) * (1.0 + Zeta);
    GaussPoint[iGauss].SetNi(val_Ni, 5);
    val_Ni = 0.125 * (1.0 + Xi) * (1.0 + Eta) * (1.0 + Zeta);
    GaussPoint[iGauss].SetNi(val_Ni, 6);
    val_Ni = 0.125 * (1.0 - Xi) * (1.0 + Eta) * (1.0 + Zeta);
    GaussPoint[iGauss].SetNi(val_Ni, 7);

    /*--- dN/d xi ---*/

    dNiXj[iGauss][0][0] = -0.125 * (1.0 - Eta) * (1.0 - Zeta);
    dNiXj[iGauss][1][0] = 0.125 * (1.0 - Eta) * (1.0 - Zeta);
    dNiXj[iGauss][2][0] = 0.125 * (1.0 + Eta) * (1.0 - Zeta);
    dNiXj[iGauss][3][0] = -0.125 * (1.0 + Eta) * (1.0 - Zeta);
    dNiXj[iGauss][4][0] = -0.125 * (1.0 - Eta) * (1.0 + Zeta);
    dNiXj[iGauss][5][0] = 0.125 * (1.0 - Eta) * (1.0 + Zeta);
    dNiXj[iGauss][6][0] = 0.125 * (1.0 + Eta) * (1.0 + Zeta);
    dNiXj[iGauss][7][0] = -0.125 * (1.0 + Eta) * (1.0 + Zeta);

    /*--- dN/d eta ---*/

    dNiXj[iGauss][0][1] = -0.125 * (1.0 - Xi) * (1.0 - Zeta);
    dNiXj[iGauss][1][1] = -0.125 * (1.0 + Xi) * (1.0 - Zeta);
    dNiXj[iGauss][2][1] = 0.125 * (1.0 + Xi) * (1.0 - Zeta);
    dNiXj[iGauss][3][1] = 0.125 * (1.0 - Xi) * (1.0 - Zeta);
    dNiXj[iGauss][4][1] = -0.125 * (1.0 - Xi) * (1.0 + Zeta);
    dNiXj[iGauss][5][1] = -0.125 * (1.0 + Xi) * (1.0 + Zeta);
    dNiXj[iGauss][6][1] = 0.125 * (1.0 + Xi) * (1.0 + Zeta);
    dNiXj[iGauss][7][1] = 0.125 * (1.0 - Xi) * (1.0 + Zeta);

    /*--- dN/d zeta ---*/

    dNiXj[iGauss][0][2] = -0.125 * (1.0 - Xi) * (1.0 - Eta);
    dNiXj[iGauss][1][2] = -0.125 * (1.0 + Xi) * (1.0 - Eta);
    dNiXj[iGauss][2][2] = -0.125 * (1.0 + Xi) * (1.0 + Eta);
    dNiXj[iGauss][3][2] = -0.125 * (1.0 - Xi) * (1.0 + Eta);
    dNiXj[iGauss][4][2] = 0.125 * (1.0 - Xi) * (1.0 - Eta);
    dNiXj[iGauss][5][2] = 0.125 * (1.0 + Xi) * (1.0 - Eta);
    dNiXj[iGauss][6][2] = 0.125 * (1.0 + Xi) * (1.0 + Eta);
    dNiXj[iGauss][7][2] = 0.125 * (1.0 - Xi) * (1.0 + Eta);
  }

  /*--- Store the extrapolation functions ---*/

  su2double ExtrapCoord[8][3], sqrt3 = 1.732050807568877;

  ExtrapCoord[0][0] = -sqrt3;
  ExtrapCoord[0][1] = -sqrt3;
  ExtrapCoord[0][2] = -sqrt3;
  ExtrapCoord[1][0] = sqrt3;
  ExtrapCoord[1][1] = -sqrt3;
  ExtrapCoord[1][2] = -sqrt3;
  ExtrapCoord[2][0] = sqrt3;
  ExtrapCoord[2][1] = sqrt3;
  ExtrapCoord[2][2] = -sqrt3;
  ExtrapCoord[3][0] = -sqrt3;
  ExtrapCoord[3][1] = sqrt3;
  ExtrapCoord[3][2] = -sqrt3;
  ExtrapCoord[4][0] = -sqrt3;
  ExtrapCoord[4][1] = -sqrt3;
  ExtrapCoord[4][2] = sqrt3;
  ExtrapCoord[5][0] = sqrt3;
  ExtrapCoord[5][1] = -sqrt3;
  ExtrapCoord[5][2] = sqrt3;
  ExtrapCoord[6][0] = sqrt3;
  ExtrapCoord[6][1] = sqrt3;
  ExtrapCoord[6][2] = sqrt3;
  ExtrapCoord[7][0] = -sqrt3;
  ExtrapCoord[7][1] = sqrt3;
  ExtrapCoord[7][2] = sqrt3;

  for (iNode = 0; iNode < NNODE; iNode++) {
    Xi = ExtrapCoord[iNode][0];
    Eta = ExtrapCoord[iNode][1];
    Zeta = ExtrapCoord[iNode][2];

    NodalExtrap[iNode][0] = 0.125 * (1.0 - Xi) * (1.0 - Eta) * (1.0 - Zeta);
    NodalExtrap[iNode][1] = 0.125 * (1.0 + Xi) * (1.0 - Eta) * (1.0 - Zeta);
    NodalExtrap[iNode][2] = 0.125 * (1.0 + Xi) * (1.0 + Eta) * (1.0 - Zeta);
    NodalExtrap[iNode][3] = 0.125 * (1.0 - Xi) * (1.0 + Eta) * (1.0 - Zeta);
    NodalExtrap[iNode][4] = 0.125 * (1.0 - Xi) * (1.0 - Eta) * (1.0 + Zeta);
    NodalExtrap[iNode][5] = 0.125 * (1.0 + Xi) * (1.0 - Eta) * (1.0 + Zeta);
    NodalExtrap[iNode][6] = 0.125 * (1.0 + Xi) * (1.0 + Eta) * (1.0 + Zeta);
    NodalExtrap[iNode][7] = 0.125 * (1.0 - Xi) * (1.0 + Eta) * (1.0 + Zeta);
  }
}

su2double CHEXA8::ComputeVolume(const FrameType mode) const {
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
    r3[iDim] = Coord[5][iDim] - Coord[0][iDim];
  }

  CrossProduct[0] = (r1[1] * r2[2] - r1[2] * r2[1]) * r3[0];
  CrossProduct[1] = (r1[2] * r2[0] - r1[0] * r2[2]) * r3[1];
  CrossProduct[2] = (r1[0] * r2[1] - r1[1] * r2[0]) * r3[2];

  Volume = fabs(CrossProduct[0] + CrossProduct[1] + CrossProduct[2]) / 6.0;

  for (iDim = 0; iDim < NDIM; iDim++) {
    r1[iDim] = Coord[2][iDim] - Coord[0][iDim];
    r2[iDim] = Coord[7][iDim] - Coord[0][iDim];
    r3[iDim] = Coord[5][iDim] - Coord[0][iDim];
  }

  CrossProduct[0] = (r1[1] * r2[2] - r1[2] * r2[1]) * r3[0];
  CrossProduct[1] = (r1[2] * r2[0] - r1[0] * r2[2]) * r3[1];
  CrossProduct[2] = (r1[0] * r2[1] - r1[1] * r2[0]) * r3[2];

  Volume += fabs(CrossProduct[0] + CrossProduct[1] + CrossProduct[2]) / 6.0;

  for (iDim = 0; iDim < NDIM; iDim++) {
    r1[iDim] = Coord[2][iDim] - Coord[0][iDim];
    r2[iDim] = Coord[3][iDim] - Coord[0][iDim];
    r3[iDim] = Coord[7][iDim] - Coord[0][iDim];
  }

  CrossProduct[0] = (r1[1] * r2[2] - r1[2] * r2[1]) * r3[0];
  CrossProduct[1] = (r1[2] * r2[0] - r1[0] * r2[2]) * r3[1];
  CrossProduct[2] = (r1[0] * r2[1] - r1[1] * r2[0]) * r3[2];

  Volume += fabs(CrossProduct[0] + CrossProduct[1] + CrossProduct[2]) / 6.0;

  for (iDim = 0; iDim < NDIM; iDim++) {
    r1[iDim] = Coord[5][iDim] - Coord[0][iDim];
    r2[iDim] = Coord[7][iDim] - Coord[0][iDim];
    r3[iDim] = Coord[4][iDim] - Coord[0][iDim];
  }

  CrossProduct[0] = (r1[1] * r2[2] - r1[2] * r2[1]) * r3[0];
  CrossProduct[1] = (r1[2] * r2[0] - r1[0] * r2[2]) * r3[1];
  CrossProduct[2] = (r1[0] * r2[1] - r1[1] * r2[0]) * r3[2];

  Volume += fabs(CrossProduct[0] + CrossProduct[1] + CrossProduct[2]) / 6.0;

  for (iDim = 0; iDim < NDIM; iDim++) {
    r1[iDim] = Coord[7][iDim] - Coord[2][iDim];
    r2[iDim] = Coord[5][iDim] - Coord[2][iDim];
    r3[iDim] = Coord[6][iDim] - Coord[2][iDim];
  }

  CrossProduct[0] = (r1[1] * r2[2] - r1[2] * r2[1]) * r3[0];
  CrossProduct[1] = (r1[2] * r2[0] - r1[0] * r2[2]) * r3[1];
  CrossProduct[2] = (r1[0] * r2[1] - r1[1] * r2[0]) * r3[2];

  Volume += fabs(CrossProduct[0] + CrossProduct[1] + CrossProduct[2]) / 6.0;

  return Volume;
}
