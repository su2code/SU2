/*!
 * \file CFEMStandardPrism.cpp
 * \brief Functions for the class CFEMStandardPrism.
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

#include "../../include/fem/CFEMStandardPrism.hpp"
#include "../../include/fem/fem_gauss_jacobi_quadrature.hpp"

/*----------------------------------------------------------------------------------*/
/*                 Public member functions of CFEMStandardPrism.                    */
/*----------------------------------------------------------------------------------*/

CFEMStandardPrism::CFEMStandardPrism(const unsigned short val_nPoly,
                                     const unsigned short val_orderExact) {

  /*--- Store the command line arguments. ---*/
  VTK_Type   = PRISM;
  nPoly      = val_nPoly;
  orderExact = val_orderExact;

  /*--- Determine the number of DOFs in 1D, the number of DOFs of the base triangle
        and the total number of DOFs.  Also determine the padded value of the 
        total number of DOFs. ---*/
  nDOFs1D       = nPoly + 1;
  nDOFsTriangle = nDOFs1D*(nDOFs1D+1)/2;
  nDOFs         = nDOFs1D*nDOFsTriangle;
  nDOFsPad      = ((nDOFs+baseVectorLen-1)/baseVectorLen)*baseVectorLen;

  /*--- Determine the 1D parametric locations of the grid DOFs. These are needed
        as the 3D grid DOFs are obtained by taking the tensor product of the
        1D grid DOFs and the grid DOFs of the base triangle. ---*/
  Location1DGridDOFsEquidistant(rLineDOFsEqui);
  Location1DGridDOFsLGL(rLineDOFsLGL);

  /*--- Determine the parametric locations of the grid DOFs of the triangle. ---*/
  LocationTriangleGridDOFsEquidistant(rTriangleDOFsEqui, sTriangleDOFsEqui);
  LocationTriangleGridDOFsLGL(rTriangleDOFsLGL, sTriangleDOFsLGL);

  /*--- Determine the 1D integration points of a line, which corresponds to the
        direction normal to the base triangle of the prism. ---*/
  nInt1D = orderExact/2 + 1;
  rLineInt.resize(nInt1D);
  wLineInt.resize(nInt1D);

  CGaussJacobiQuadrature GaussJacobi;
  GaussJacobi.GetQuadraturePoints(0.0, 0.0, -1.0, 1.0, rLineInt, wLineInt);

  /*--- Determine the parametric location and weights of the
        integration rule of the base triangle. ---*/
  IntegrationPointsTriangle(rTriangleInt, sTriangleInt, wTriangleInt);
  nIntTriangle = rTriangleInt.size();

  /*--- The 3D quadrature rule is a tensor product of the 1D Gauss-Legendre
        quadrature rule and the integration rule of the triangle. Determine
        the total number of integration points and its padded value. ---*/
  nIntegration    = nInt1D*nIntTriangle;
  nIntegrationPad = ((nIntegration+baseVectorLen-1)/baseVectorLen)*baseVectorLen;

  /*--- Allocate the memory for the padded number of integration points and
        initialize the weights to zero. This is done such that the padded
        values are initialized appropriately. ---*/
  wIntegration.resize(nIntegrationPad);
  wIntegration.setConstant(0.0);

  /*--- Determine the integration weights of the prism. ---*/
  unsigned short ii = 0;
  for(unsigned short j=0; j<nInt1D; ++j)
    for(unsigned short i=0; i<nIntTriangle; ++i, ++ii)
        wIntegration(ii) = wTriangleInt[i]*wLineInt[j];
}

void CFEMStandardPrism::LocationAllIntegrationPoints(vector<passivedouble> &rInt,
                                                     vector<passivedouble> &sInt,
                                                     vector<passivedouble> &tInt) {

  /*--- Allocate the memory for rInt, sInt and tInt. ---*/
  rInt.resize(nIntegration);
  sInt.resize(nIntegration);
  tInt.resize(nIntegration);

  /*--- Determine the location of all the integration points. ---*/
  unsigned short ii = 0;
  for(unsigned short k=0; k<rLineInt.size(); ++k) {
    for(unsigned short j=0; j<rTriangleInt.size(); ++j, ++ii) {
      rInt[ii] = rTriangleInt[j];
      sInt[ii] = sTriangleInt[j];
      tInt[ii] = rLineInt[k];
    }
  }
}
