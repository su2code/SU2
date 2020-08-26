/*!
 * \file CFEMStandardTet.cpp
 * \brief Functions for the class CFEMStandardTet.
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

#include "../../include/fem/CFEMStandardTet.hpp"

/*----------------------------------------------------------------------------------*/
/*               Public member functions of CFEMStandardTet.                        */
/*----------------------------------------------------------------------------------*/

CFEMStandardTet::CFEMStandardTet(const unsigned short val_nPoly,
                                 const unsigned short val_orderExact) {

  /*--- Store the command line arguments. ---*/
  VTK_Type   = TETRAHEDRON;
  nPoly      = val_nPoly;
  orderExact = val_orderExact;

  /*--- Determine the total number of DOFs and its padded version. ---*/
  nDOFs    = (nPoly+1)*(nPoly+2)*(nPoly+3)/6;
  nDOFsPad = ((nDOFs+vecLen-1)/vecLen)*vecLen;

  /*--- Determine the parametric locations of the grid DOFs of the tetrahedron. ---*/
  LocationTetGridDOFsEquidistant();
  LocationTetGridDOFsLGL();

  /*--- Determine the parametric location and weights of the
        integration rule of the tetrahedron. ---*/
  IntegrationPointsTetrahedron(rTetInt, sTetInt, tTetInt, wTetInt);

  /*--- Determine the total number of integration points
        and its padded version. ---*/
  nIntegration    = rTetInt.size();
  nIntegrationPad = ((nIntegration+vecLen-1)/vecLen)*vecLen;

  /*--- Allocate the memory for the padded number of integration points and
        initialize the weights to zero. This is done such that the padded
        values are initialized appropriately. ---*/
  wIntegration.resize(nIntegrationPad);
  wIntegration.setConstant(0.0);

  /*--- Copy the values from wTetInt to wIntegration. ---*/
  for(unsigned short i=0; i<nIntegration; ++i)
    wIntegration(i) = wTetInt[i];
}

void CFEMStandardTet::LocationTetGridDOFsEquidistant() {

  /*--- For a tetrahedron it is not possible to apply
        a tensor product and therefore all DOFs are
        simply stored. Allocate the memory. ---*/
  rTetDOFsEqui.resize(nDOFs);
  sTetDOFsEqui.resize(nDOFs);
  tTetDOFsEqui.resize(nDOFs);

  /*--- Determine the equidistant spacing along an edge. ---*/
  const passivedouble dh = 2.0/nPoly;

  /*--- Triple loop to compute the location of the grid DOFs. ---*/
  unsigned short ii = 0;
  for(unsigned short k=0; k<=nPoly; ++k) {
    const passivedouble t = -1.0 + k*dh;
    const unsigned short uppBoundJ = nPoly - k;
    for(unsigned short j=0; j<=uppBoundJ; ++j) {
      const passivedouble s = -1.0 + j*dh;
      const unsigned short uppBoundI = nPoly - k - j;
      for(unsigned short i=0; i<=uppBoundI; ++i, ++ii) {
        rTetDOFsEqui[ii] = -1.0 + i*dh;
        sTetDOFsEqui[ii] =  s;
        tTetDOFsEqui[ii] =  t;
      }
    }
  }
}

void CFEMStandardTet::LocationTetGridDOFsLGL() {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}
