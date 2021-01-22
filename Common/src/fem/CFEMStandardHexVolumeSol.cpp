/*!
 * \file CFEMStandardHexVolumeSol.cpp
 * \brief Functions for the class CFEMStandardHexVolumeSol.
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

#include "../../include/fem/CFEMStandardHexVolumeSol.hpp"

/*----------------------------------------------------------------------------------*/
/*             Public member functions of CFEMStandardHexVolumeSol.                 */
/*----------------------------------------------------------------------------------*/

CFEMStandardHexVolumeSol::CFEMStandardHexVolumeSol(const unsigned short val_nPoly,
                                                   const unsigned short val_orderExact,
                                                   const unsigned short val_locGridDOFs,
                                                   const unsigned short val_nVar)
  : CFEMStandardHexBase(val_nPoly,val_orderExact) {

  /*--- Compute the 1D parametric coordinates of the solution DOFs. Only
        different from grid DOFs when a different polynomial degree is
        used for the grid and solution. ---*/
  if(val_locGridDOFs == LGL) Location1DGridDOFsLGL(nPoly, rLineSolDOFs);
  else                       Location1DGridDOFsEquidistant(nPoly, rLineSolDOFs);

  /*--- Compute the 1D Legendre basis functions and its first
        and second derivatives in the integration points. ---*/
  const unsigned short nInt1DPad = PaddedValue(nInt1D);
  legBasisLineInt.resize(nInt1DPad, nDOFs1D);     legBasisLineInt.setConstant(0.0);
  derLegBasisLineInt.resize(nInt1DPad, nDOFs1D);  derLegBasisLineInt.setConstant(0.0);
  hesLegBasisLineInt.resize(nInt1DPad, nDOFs1D);  hesLegBasisLineInt.setConstant(0.0);

  Vandermonde1D(nPoly, rLineInt, legBasisLineInt);
  GradVandermonde1D(nPoly, rLineInt, derLegBasisLineInt);
  HesVandermonde1D(nPoly, rLineInt, hesLegBasisLineInt);

  /*--- Compute the 1D Legendre basis functions and its first
        derivatives in the solution DOFs. ---*/
  const unsigned short nDOFs1DPad = PaddedValue(nDOFs1D);
  legBasisLineSolDOFs.resize(nDOFs1DPad, nDOFs1D);    legBasisLineSolDOFs.setConstant(0.0);
  derLegBasisLineSolDOFs.resize(nDOFs1DPad, nDOFs1D); derLegBasisLineSolDOFs.setConstant(0.0);

  Vandermonde1D(nPoly, rLineSolDOFs, legBasisLineSolDOFs);
  GradVandermonde1D(nPoly, rLineSolDOFs, derLegBasisLineSolDOFs);

  /*--- Determine the local subconnectivity of this standard element when split
        in several linear elements. Used for a.o. plotting and searching. ---*/
  CFEMStandardHexBase::SubConnLinearElements();

  /*--- Set the function pointers for the tensor product multiplications to
        determine the data in the volume integration points and volume
        nodal solution DOFs. ---*/
  SetFunctionPointerVolumeDataHex(nDOFs1D, nInt1D, TensorProductDataVolIntPoints);
  SetFunctionPointerVolumeDataHex(nDOFs1D, nDOFs1D, TensorProductDataVolSolDOFs);
}

void CFEMStandardHexVolumeSol::BasisFunctionsInPoints(const vector<vector<passivedouble> > &parCoor,
                                                      ColMajorMatrix<passivedouble>        &matBasis) {

  /*--- Determine the 1D basis functions for all parametric coordinates in
        the three directions. ---*/
  const unsigned short nCoor = parCoor[0].size();

  ColMajorMatrix<passivedouble> legR(nCoor,nDOFs1D), legS(nCoor,nDOFs1D), legT(nCoor,nDOFs1D);

  Vandermonde1D(nPoly, parCoor[0], legR);
  Vandermonde1D(nPoly, parCoor[1], legS);
  Vandermonde1D(nPoly, parCoor[2], legT);

  /*--- Allocate the memory for matBasis and set the values. ---*/
  matBasis.resize(nCoor, nDOFs);

  unsigned short ii = 0;
  for(unsigned short k=0; k<nDOFs1D; ++k)
    for(unsigned short j=0; j<nDOFs1D; ++j)
      for(unsigned short i=0; i<nDOFs1D; ++i, ++ii)
        for(unsigned short l=0; l<nCoor; ++l)
          matBasis(l,ii) = legR(l,i)*legS(l,j)*legT(l,k);
}
