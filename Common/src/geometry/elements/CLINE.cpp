/*!
 * \file CLINE.cpp
 * \brief Definition of the 2-node line element with two Gauss points.
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

CLINE::CLINE() : CElementWithKnownSizes<NGAUSS, NNODE, NDIM>() {
  su2double Xi, val_Ni;

  /*--- Gauss coordinates and weights ---*/

  su2double oneOnTwoSqrt3 = 0.288675134594813;
  GaussCoord[0][0] = 0.5 - oneOnTwoSqrt3;
  GaussWeight(0) = 0.5;
  GaussCoord[1][0] = 0.5 + oneOnTwoSqrt3;
  GaussWeight(1) = 0.5;

  /*--- Store the values of the shape functions and their derivatives ---*/

  unsigned short iGauss;
  for (iGauss = 0; iGauss < nGaussPoints; iGauss++) {
    Xi = GaussCoord[iGauss][0];

    val_Ni = 1.0 - Xi;
    GaussPoint[iGauss].SetNi(val_Ni, 0);
    val_Ni = Xi;
    GaussPoint[iGauss].SetNi(val_Ni, 1);

    /*--- dN/d xi ---*/

    dNiXj[iGauss][0][0] = -1.0;
    dNiXj[iGauss][1][0] = 1.0;
  }
}

su2double CLINE::ComputeLength(const FrameType mode) const {
  /*--- Select the appropriate source for the nodal coordinates depending on the frame requested
        for the gradient computation, REFERENCE (undeformed) or CURRENT (deformed)---*/
  const su2activematrix& Coord = (mode == REFERENCE) ? RefCoord : CurrentCoord;
  return fabs(Coord[1][0] - Coord[0][0]);
}
