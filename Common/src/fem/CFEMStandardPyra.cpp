/*!
 * \file CFEMStandardPyra.cpp
 * \brief Functions for the class CFEMStandardPyra.
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

#include "../../include/fem/CFEMStandardPyra.hpp"

/*----------------------------------------------------------------------------------*/
/*              Public member functions of CFEMStandardPyra.                        */
/*----------------------------------------------------------------------------------*/

CFEMStandardPyra::CFEMStandardPyra(const unsigned short val_nPoly,
                                   const unsigned short val_orderExact) {

  /*--- Store the command line arguments. ---*/
  VTK_Type   = PYRAMID;
  nPoly      = val_nPoly;
  orderExact = val_orderExact;

  /*--- Determine the number of DOFs in 1D of the base quad and the total number of DOFs.
        Also determine the padded value of the latter. ---*/
  nDOFs1D  = nPoly + 1;
  nDOFs    = nDOFs1D*(nDOFs1D+1)*(2*nDOFs1D+1)/6;
  nDOFsPad = ((nDOFs+vecLen-1)/vecLen)*vecLen;

  /*--- The 3D quadrature rule for a pyramid is obtained by transforming the
        standard pyramid into a standard hexahedron by means of the Duffy
        transformation. On the created hexahedron a tensor product rule is used,
        but in the t-direction a Gauss-Jacobi rule with alpha = 2, beta = 0 is
        used to account for the transformation from the pyramid to the
        hexahedron. The determinant of the Jacobian of this transformation
        is 0.25*(1-t)^2, which explains the alpha = 2 in the Gauss-Jacobi.
        Determine the number of integration points in 1D and the total
        number. Also determine the padded value of the latter. ---*/
  nInt1D          = orderExact/2 + 1;
  nIntegration    = nInt1D*nInt1D*nInt1D;
  nIntegrationPad = ((nIntegration+vecLen-1)/vecLen)*vecLen;
}
