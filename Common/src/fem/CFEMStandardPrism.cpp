/*!
 * \file CFEMStandardPrism.cpp
 * \brief Functions for the class CFEMStandardPrism.
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

#include "../../include/fem/CFEMStandardPrism.hpp"

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
  nDOFsPad      = ((nDOFs+vecLen-1)/vecLen)*vecLen;

  /*--- Determine the parametric location and weights of the
        integration rule of the base triangle of the prism. ---*/
  vector<passivedouble> rTriangle, sTriangle, wTriangle;
  IntegrationPointsTriangle(rTriangle, sTriangle, wTriangle);

  nIntTriangle = rTriangle.size();

  /*--- The 3D quadrature rule is a tensor product of the 1D Gauss-Legendre
        quadrature rule and the integration rule of the triangle. Determine
        the number of integration points in 1D and the total number.
        Also determine the padded value of the latter. ---*/
  nInt1D          = orderExact/2 + 1;
  nIntegration    = nInt1D*nIntTriangle;
  nIntegrationPad = ((nIntegration+vecLen-1)/vecLen)*vecLen;
}
