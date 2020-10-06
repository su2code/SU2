/*!
 * \file CFEMStandardQuadGrid.cpp
 * \brief Functions for the class CFEMStandardQuadGrid.
 * \author E. van der Weide
 * \version 7.0.6 "Blackbird"
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

#include "../../include/fem/CFEMStandardQuadGrid.hpp"

/*----------------------------------------------------------------------------------*/
/*             Public member functions of CFEMStandardQuadGrid.                     */
/*----------------------------------------------------------------------------------*/

CFEMStandardQuadGrid::CFEMStandardQuadGrid(const unsigned short val_nPoly,
                                           const unsigned short val_orderExact,
                                           const bool           val_surfElement)
  : CFEMStandardQuad(val_nPoly, val_orderExact) {

  /*--- Determine the number of space dimensions, which is 3 if this standard
        element is a surface element and 2 when it is a value element. ---*/
  nDim = val_surfElement ? 3 : 2;

  /*--- Compute the values of the 1D Lagrangian basis functions in the integration
        points for both the equidistant and LGL point distribution. ---*/
  LagBasisIntPointsLine(rLineDOFsEqui, rLineInt, lagBasisLineIntEqui);
  LagBasisIntPointsLine(rLineDOFsLGL,  rLineInt, lagBasisLineIntLGL);

  /*--- Compute the values of the derivatives of the 1D Lagrangian basis functions in
        the integration points for both the equidistant and LGL point distribution. ---*/
  DerLagBasisIntPointsLine(rLineDOFsEqui, rLineInt, derLagBasisLineIntEqui);
  DerLagBasisIntPointsLine(rLineDOFsLGL,  rLineInt, derLagBasisLineIntLGL);
}

void CFEMStandardQuadGrid::DerivativesCoorIntPoints(const bool                         LGLDistribution,
                                                    ColMajorMatrix<su2double>          &matCoor,
                                                    vector<ColMajorMatrix<su2double> > &matDerCoor) {

  /*--- Check for which point distribution the derivatives must be computed. ---*/
  if( LGLDistribution ) {

    /*--- LGL distribution. Call the function TensorProductIntegrationPoints 2 times to compute the
          derivatives of the Cartesian coordinates w.r.t. the two parametric coordinates. The first
          argument in the function call is nDim, which corresponds to the number of Cartesian
          coordinates (3 for a surface element and 2 for a volume element). ---*/
    TensorProductIntegrationPoints(nDim, derLagBasisLineIntLGL, lagBasisLineIntLGL,
                                   matCoor, matDerCoor[0], nullptr);
    TensorProductIntegrationPoints(nDim, lagBasisLineIntLGL, derLagBasisLineIntLGL,
                                   matCoor, matDerCoor[1], nullptr);
  }
  else {

    /*--- Equidistant distribution. Call the function TensorProductIntegrationPoints 3 times to compute the
          derivatives of the Cartesian coordinates w.r.t. the three parametric coordinates. The first
          argument in the function call is nDim, which corresponds to the number of Cartesian
          coordinates (3 for a surface element and 2 for a volume element). ---*/
    TensorProductIntegrationPoints(nDim, derLagBasisLineIntEqui, lagBasisLineIntEqui,
                                   matCoor, matDerCoor[0], nullptr);
    TensorProductIntegrationPoints(nDim, lagBasisLineIntEqui, derLagBasisLineIntEqui,
                                   matCoor, matDerCoor[1], nullptr);
  }
}
