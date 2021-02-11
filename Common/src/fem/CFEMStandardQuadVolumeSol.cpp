/*!
 * \file CFEMStandardQuadVolumeSol.cpp
 * \brief Functions for the class CFEMStandardQuadVolumeSol.
 * \author E. van der Weide
 * \version 7.1.0 "Blackbird"
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

#include "../../include/fem/CFEMStandardQuadVolumeSol.hpp"
#include "../../include/toolboxes/CSquareMatrixCM.hpp"

/*----------------------------------------------------------------------------------*/
/*             Public member functions of CFEMStandardQuadVolumeSol.                */
/*----------------------------------------------------------------------------------*/

CFEMStandardQuadVolumeSol::CFEMStandardQuadVolumeSol(const unsigned short val_nPoly,
                                                     const unsigned short val_orderExact,
                                                     const unsigned short val_locGridDOFs,
                                                     const unsigned short val_nVar)
  : CFEMStandardQuadBase(val_nPoly, val_orderExact) {

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

  /*--- Compute the 1D Legendre basis functions in the solution DOFs.
        Store it in a square matrix as also the inverse is needed. ---*/
  CSquareMatrixCM Vtmp(nDOFs1D);
  Vandermonde1D(nPoly, rLineSolDOFs, Vtmp.GetMat());

  /*--- Store the contents of Vtmp in legBasisLineSolDOFs. Note that
        the first dimension of legBasisLineSolDOFs is padded. ---*/
  const unsigned short nDOFs1DPad = PaddedValue(nDOFs1D);
  legBasisLineSolDOFs.resize(nDOFs1DPad, nDOFs1D);
  legBasisLineSolDOFs.setConstant(0.0);

  for(unsigned short j=0; j<nDOFs1D; ++j)
    for(unsigned short i=0; i<nDOFs1D; ++i)
      legBasisLineSolDOFs(i,j) = Vtmp(i,j);

  /*--- Compute the inverse of Vtmp and store the contents in legBasisLineSolDOFsInv.
        Note that the first dimension of legBasisLineSolDOFsInv is padded. ---*/
  Vtmp.Invert();

  legBasisLineSolDOFsInv.resize(nDOFs1DPad, nDOFs1D);
  legBasisLineSolDOFsInv.setConstant(0.0);

  for(unsigned short j=0; j<nDOFs1D; ++j)
    for(unsigned short i=0; i<nDOFs1D; ++i)
      legBasisLineSolDOFsInv(i,j) = Vtmp(i,j);

  /*--- Compute the first derivatives of the 1D Legendre basis
        functions in the solution DOFs. ---*/
  derLegBasisLineSolDOFs.resize(nDOFs1DPad, nDOFs1D);
  derLegBasisLineSolDOFs.setConstant(0.0);

  GradVandermonde1D(nPoly, rLineSolDOFs, derLegBasisLineSolDOFs);

  /*--- Determine the local subconnectivity of this standard element when split
        in several linear elements. Used for a.o. plotting and searching. ---*/
  CFEMStandardQuadBase::SubConnLinearElements();

  /*--- Set the function pointers for the tensor product multiplications to
        determine the data in the volume integration points and volume
        nodal solution DOFs. ---*/
  SetFunctionPointerVolumeDataQuad(nDOFs1D, nInt1D, TensorProductDataVolIntPoints);
  SetFunctionPointerVolumeDataQuad(nDOFs1D, nDOFs1D, TensorProductDataVolSolDOFs);
}

void CFEMStandardQuadVolumeSol::BasisFunctionsInPoints(const vector<vector<passivedouble> > &parCoor,
                                                       ColMajorMatrix<passivedouble>        &matBasis) {

  /*--- Determine the 1D basis functions for all parametric coordinates in
        the two directions. ---*/
  const unsigned short nCoor = parCoor[0].size();

  ColMajorMatrix<passivedouble> legR(nCoor,nDOFs1D), legS(nCoor,nDOFs1D);

  Vandermonde1D(nPoly, parCoor[0], legR);
  Vandermonde1D(nPoly, parCoor[1], legS);

  /*--- Allocate the memory for matBasis and set the values. ---*/
  matBasis.resize(nCoor, nDOFs);

  unsigned short ii = 0;
  for(unsigned short j=0; j<nDOFs1D; ++j)
    for(unsigned short i=0; i<nDOFs1D; ++i, ++ii)
      for(unsigned short l=0; l<nCoor; ++l)
        matBasis(l,ii) = legR(l,i)*legS(l,j);
}

void CFEMStandardQuadVolumeSol::ModalToNodal(ColMajorMatrix<su2double> &solDOFs) {

  /*--- Copy solDOFs into tmp and carry out the tensor product for
        the conversion to the nodal formulation. ---*/
  const ColMajorMatrix<su2double> tmp = solDOFs;

  TensorProductVolumeDataQuad(TensorProductDataVolSolDOFs, solDOFs.cols(), nDOFs1D,
                             nDOFs1D, legBasisLineSolDOFs, legBasisLineSolDOFs,
                             tmp, solDOFs, nullptr);
}

void CFEMStandardQuadVolumeSol::NodalToModal(ColMajorMatrix<su2double> &solDOFs) {

  /*--- Copy solDOFs into tmp and carry out the tensor product for
        the conversion to the modal formulation. ---*/
  const ColMajorMatrix<su2double> tmp = solDOFs;

  TensorProductVolumeDataQuad(TensorProductDataVolSolDOFs, solDOFs.cols(), nDOFs1D,
                              nDOFs1D, legBasisLineSolDOFsInv, legBasisLineSolDOFsInv,
                              tmp, solDOFs, nullptr);
}

passivedouble CFEMStandardQuadVolumeSol::ValBasis0(void) {

  /*--- Easier storage of the 1D Legendre basis function
        and return the value of the 2D basis function.  ---*/
  const passivedouble leg1D = legBasisLineSolDOFs(0,0);
  return leg1D*leg1D;
}
