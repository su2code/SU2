/*!
 * \file CFEMStandardLineGrid.cpp
 * \brief Functions for the class CFEMStandardLineGrid.
 * \author E. van der Weide
 * \version 7.0.7 "Blackbird"
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

#include "../../include/fem/CFEMStandardLineGrid.hpp"

/*----------------------------------------------------------------------------------*/
/*            Public member functions of CFEMStandardLineGrid.                      */
/*----------------------------------------------------------------------------------*/

CFEMStandardLineGrid::CFEMStandardLineGrid(const unsigned short val_nPoly,
                                           const unsigned short val_orderExact)
  : CFEMStandardLine(val_nPoly, val_orderExact) {

  /*--- Compute the values of the 1D Lagrangian basis functions in the integration
        points for both the equidistant and LGL point distribution. ---*/
  LagBasisIntPointsLine(rLineDOFsEqui, rLineInt, lagBasisLineIntEqui);
  LagBasisIntPointsLine(rLineDOFsLGL,  rLineInt, lagBasisLineIntLGL);

  /*--- Compute the values of the derivatives of the 1D Lagrangian basis functions in
        the integration points for both the equidistant and LGL point distribution. ---*/
  DerLagBasisIntPointsLine(rLineDOFsEqui, rLineInt, derLagBasisLineIntEqui);
  DerLagBasisIntPointsLine(rLineDOFsLGL,  rLineInt, derLagBasisLineIntLGL);

  /*--- Set up the jitted gemm call, if supported. For this particular standard
        element the derivatives of the coordinates are computed, which is 2,
        because the line element is only present as a surface element in a
        2D simulation. ---*/
  SetUpJittedGEMM(nIntegrationPad, 2, nDOFs);
}

void CFEMStandardLineGrid::DerivativesCoorIntPoints(const bool                         LGLDistribution,
                                                    ColMajorMatrix<su2double>          &matCoor,
                                                    vector<ColMajorMatrix<su2double> > &matDerCoor) {

  /*--- Check for which point distribution the derivatives must be computed. ---*/
  if( LGLDistribution ) {

    /*--- LGL distribution. Call the function OwnGemm to compute the derivatives
          of the Cartesian coordinates w.r.t. the parametric coordinate. ---*/
    OwnGemm(nIntegrationPad, 2, nDOFs, derLagBasisLineIntLGL, matCoor, matDerCoor[0], nullptr);
  }
  else {

    /*--- Equidistant distribution. Call the function OwnGemm to compute the derivatives
          of the Cartesian coordinates w.r.t. the parametric coordinate. ---*/
    OwnGemm(nIntegrationPad, 2, nDOFs, derLagBasisLineIntEqui, matCoor, matDerCoor[0], nullptr);
  }
}
