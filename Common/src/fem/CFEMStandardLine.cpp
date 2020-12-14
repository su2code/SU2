/*!
 * \file CFEMStandardLine.cpp
 * \brief Functions for the class CFEMStandardLine.
 * \author E. van der Weide
 * \version 7.0.8 "Blackbird"
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

#include "../../include/fem/CFEMStandardLine.hpp"
#include "../../include/fem/fem_gauss_jacobi_quadrature.hpp"

/*----------------------------------------------------------------------------------*/
/*              Public member functions of CFEMStandardLine.                        */
/*----------------------------------------------------------------------------------*/

CFEMStandardLine::CFEMStandardLine(const unsigned short val_nPoly,
                                   const unsigned short val_orderExact) {

  /*--- Store the command line arguments. ---*/
  VTK_Type   = LINE;
  nPoly      = val_nPoly;
  orderExact = val_orderExact;

  /*--- Determine the total number of DOFs and its padded version. ---*/
  nDOFs    = (nPoly+1);
  nDOFsPad = ((nDOFs+baseVectorLen-1)/baseVectorLen)*baseVectorLen;

  /*--- Determine the 1D parametric locations of the grid DOFs. ---*/
  Location1DGridDOFsEquidistant(rLineDOFsEqui);
  Location1DGridDOFsLGL(rLineDOFsLGL);

  /*--- Determine the total number of integration points
        and its padded version. ---*/
  nIntegration    = orderExact/2 + 1;
  nIntegrationPad = ((nIntegration+baseVectorLen-1)/baseVectorLen)*baseVectorLen;

  /*--- Determine the location and the weights of the 1D integration points. ---*/
  rLineInt.resize(nIntegration);
  wLineInt.resize(nIntegration);

  CGaussJacobiQuadrature GaussJacobi;
  GaussJacobi.GetQuadraturePoints(0.0, 0.0, -1.0, 1.0, rLineInt, wLineInt);

  /*--- Allocate the memory for the padded number of integration points and
        initialize the weights to zero. This is done such that the padded
        values are initialized appropriately. ---*/
  wIntegration.resize(nIntegrationPad);
  wIntegration.setConstant(0.0);

  /*--- Copy the values from wLineInt to wIntegration. ---*/
  for(unsigned short i=0; i<nIntegration; ++i)
    wIntegration(i) = wLineInt[i];
}
