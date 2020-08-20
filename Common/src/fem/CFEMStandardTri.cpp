/*!
 * \file CFEMStandardTri.cpp
 * \brief Functions for the class CFEMStandardTri.
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

#include "../../include/fem/CFEMStandardTri.hpp"

/*----------------------------------------------------------------------------------*/
/*               Public member functions of CFEMStandardTri.                        */
/*----------------------------------------------------------------------------------*/

CFEMStandardTri::CFEMStandardTri(const unsigned short val_nPoly,
                                 const unsigned short val_orderExact) {

  /*--- Store the command line arguments. ---*/
  VTK_Type   = TRIANGLE;
  nPoly      = val_nPoly;
  orderExact = val_orderExact;

  /*--- Determine the total number of DOFs and its padded version. ---*/
  nDOFs    = (nPoly+1)*(nPoly+2)/2;
  nDOFsPad = ((nDOFs+vecLen-1)/vecLen)*vecLen;

  /*--- Determine the parametric location and weights of the
        integration rule of the base triangle of the prism. ---*/
  vector<su2double> rTriangle, sTriangle, wTriangle;
  IntegrationPointsTriangle(rTriangle, sTriangle, wTriangle);

  /*--- Determine the total number of integration points
        and its padded version. ---*/
  nIntegration    = rTriangle.size();
  nIntegrationPad = ((nIntegration+vecLen-1)/vecLen)*vecLen;
}
