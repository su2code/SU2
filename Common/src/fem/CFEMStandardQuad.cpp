/*!
 * \file CFEMStandardQuad.cpp
 * \brief Functions for the class CFEMStandardQuad.
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

#include "../../include/fem/CFEMStandardQuad.hpp"
#include "../../include/fem/fem_gauss_jacobi_quadrature.hpp"

/*----------------------------------------------------------------------------------*/
/*          Public member functions of CFEMStandardVolumeQuadGrid.                  */
/*----------------------------------------------------------------------------------*/

CFEMStandardQuad::CFEMStandardQuad(const unsigned short val_nPoly,
                                   const unsigned short val_orderExact) {

  /*--- Store the command line arguments. ---*/
  VTK_Type   = QUADRILATERAL;
  nPoly      = val_nPoly;
  orderExact = val_orderExact;

  /*--- Determine the number of DOFs in 1D and the total number of DOFs.
        Also determine the padded value of the latter. ---*/
  nDOFs1D  = nPoly + 1;
  nDOFs    = nDOFs1D*nDOFs1D;
  nDOFsPad = ((nDOFs+baseVectorLen-1)/baseVectorLen)*baseVectorLen;

  /*--- Determine the 1D parametric locations of the grid DOFs. 1D is enough,
        because a tensor product is used to obtain the 2D coordinates. ---*/
  Location1DGridDOFsEquidistant(rLineDOFsEqui);
  Location1DGridDOFsLGL(rLineDOFsLGL);

  /*--- Determine the number of integration points in 1D as well as the total number
        of integration points. Also determine the padded value of the latter. ---*/
  nInt1D          = orderExact/2 + 1;
  nIntegration    = nInt1D*nInt1D;
  nIntegrationPad = ((nIntegration+baseVectorLen-1)/baseVectorLen)*baseVectorLen;

  /*--- Determine the location and the weights of the 1D integration points.
        The 2D integration points are obtained via a tensor product. ---*/
  rLineInt.resize(nInt1D);
  wLineInt.resize(nInt1D);

  CGaussJacobiQuadrature GaussJacobi;
  GaussJacobi.GetQuadraturePoints(0.0, 0.0, -1.0, 1.0, rLineInt, wLineInt);

  /*--- Allocate the memory for the padded number of integration points and
        initialize the weights to zero. This is done such that the padded
        values are initialized appropriately. ---*/
  wIntegration.resize(nIntegrationPad);
  wIntegration.setConstant(0.0);

  /*--- Double loop to set the integration weights in the 2D integration points. ---*/
  unsigned short ii = 0;
  for(unsigned short j=0; j<nInt1D; ++j)
    for(unsigned short i=0; i<nInt1D; ++i, ++ii)
      wIntegration(ii) = wLineInt[i]*wLineInt[j];
}
