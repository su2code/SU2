/*!
 * \file CFEMStandardQuadVolumeSol.cpp
 * \brief Functions for the class CFEMStandardQuadVolumeSol.
 * \author E. van der Weide
 * \version 7.1.1 "Blackbird"
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

  /*--- Create the transpose of legBasisLineInt and derLegBasisLineInt. ---*/
  const unsigned short nDOFs1DPad = PaddedValue(nDOFs1D);
  legBasisLineIntTranspose.resize(nDOFs1DPad, nInt1D);    legBasisLineIntTranspose.setConstant(0.0);
  derLegBasisLineIntTranspose.resize(nDOFs1DPad, nInt1D); derLegBasisLineIntTranspose.setConstant(0.0);

  for(unsigned short j=0; j<nInt1D; ++j) {
    for(unsigned short i=0; i<nDOFs1D; ++i) {
      legBasisLineIntTranspose(i,j)    = legBasisLineInt(j,i);
      derLegBasisLineIntTranspose(i,j) = derLegBasisLineInt(j,i);
    }
  }

  Vandermonde1D(nPoly, rLineInt, legBasisLineInt);
  GradVandermonde1D(nPoly, rLineInt, derLegBasisLineInt);
  HesVandermonde1D(nPoly, rLineInt, hesLegBasisLineInt);

  /*--- Compute the 1D Legendre basis functions in the solution DOFs.
        Store it in a square matrix as also the inverse is needed. ---*/
  CSquareMatrixCM Vtmp(nDOFs1D);
  Vandermonde1D(nPoly, rLineSolDOFs, Vtmp.GetMat());

  /*--- Store the contents of Vtmp in legBasisLineSolDOFs. Note that
        the first dimension of legBasisLineSolDOFs is padded. ---*/
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

  /*--- Determine the correction factors for the inviscid and viscous
        spectral radii for the high order element. These factors depend
        on the polynomial degree and the element type. ---*/
  switch( nPoly ) {
    case 0: factInviscidRad =  2.0; factViscousRad =     6.0; break;
    case 1: factInviscidRad =  6.0; factViscousRad =    36.0; break;
    case 2: factInviscidRad = 12.0; factViscousRad =   150.0; break;
    case 3: factInviscidRad = 20.0; factViscousRad =   420.0; break;
    case 4: factInviscidRad = 28.0; factViscousRad =   980.0; break;
    case 5: factInviscidRad = 38.0; factViscousRad =  1975.0; break;
    case 6: factInviscidRad = 50.0; factViscousRad =  3575.0; break;
    case 7: factInviscidRad = 64.0; factViscousRad =  7000.0; break;
    case 8: factInviscidRad = 80.0; factViscousRad = 14000.0; break;
    case 9: factInviscidRad = 98.0; factViscousRad = 28000.0; break;
    default:
      SU2_MPI::Error(string("Polynomial order not foreseen"), CURRENT_FUNCTION);
  }
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
                             tmp, solDOFs, true, nullptr);
}

void CFEMStandardQuadVolumeSol::NodalToModal(ColMajorMatrix<su2double> &solDOFs) {

  /*--- Copy solDOFs into tmp and carry out the tensor product for
        the conversion to the modal formulation. ---*/
  const ColMajorMatrix<su2double> tmp = solDOFs;

  TensorProductVolumeDataQuad(TensorProductDataVolSolDOFs, solDOFs.cols(), nDOFs1D,
                              nDOFs1D, legBasisLineSolDOFsInv, legBasisLineSolDOFsInv,
                              tmp, solDOFs, true, nullptr);
}

void CFEMStandardQuadVolumeSol::GradSolIntPoints(ColMajorMatrix<su2double>          &matSolDOF,
                                                 vector<ColMajorMatrix<su2double> > &matGradSolInt) {

  /*--- Call the function TensorProductVolumeDataHex 2 times to compute the derivatives
        of the solution coordinates w.r.t. the two parametric coordinates. ---*/
  TensorProductVolumeDataQuad(TensorProductDataVolIntPoints, matSolDOF.cols(), nDOFs1D, nInt1D,
                              derLegBasisLineInt, legBasisLineInt, matSolDOF, matGradSolInt[0],
                              true, nullptr);
  TensorProductVolumeDataQuad(TensorProductDataVolIntPoints, matSolDOF.cols(), nDOFs1D, nInt1D,
                              legBasisLineInt, derLegBasisLineInt, matSolDOF, matGradSolInt[1],
                              true, nullptr);

  /*--- Fill the padded data to avoid problems. ---*/
  for(unsigned short j=0; j<matSolDOF.cols(); ++j) {
    for(unsigned short i=nIntegration; i<nIntegrationPad; ++i) {
      matGradSolInt[0](i,j) = matGradSolInt[0](0,j);
      matGradSolInt[1](i,j) = matGradSolInt[1](0,j);
    }
  }
}

void CFEMStandardQuadVolumeSol::SolIntPoints(ColMajorMatrix<su2double> &matSolDOF,
                                             ColMajorMatrix<su2double> &matSolInt) {

  /*--- Call TensorProductVolumeDataQuad with the appropriate arguments
        to carry out the actual job. ---*/
  TensorProductVolumeDataQuad(TensorProductDataVolIntPoints, matSolDOF.cols(), nDOFs1D, nInt1D,
                              legBasisLineInt, legBasisLineInt, matSolDOF, matSolInt, true, nullptr);

  /*--- Fill the padded data to avoid problems. ---*/
  for(unsigned short j=0; j<matSolInt.cols(); ++j)
    for(unsigned short i=nIntegration; i<nIntegrationPad; ++i)
      matSolInt(i,j) = matSolInt(0,j);
}

void CFEMStandardQuadVolumeSol::SolIntPointsDOFsPadded(ColMajorMatrix<su2double> &matSolDOF,
                                                       ColMajorMatrix<su2double> &matSolInt) {
  SolIntPoints(matSolDOF, matSolInt);
}

void CFEMStandardQuadVolumeSol::ResidualBasisFunctions(ColMajorMatrix<su2double> &scalarDataInt,
                                                       ColMajorMatrix<su2double> &resDOFs) {

  /*--- Call TensorProductVolumeDataQuad with the appropriate arguments
        to carry out the actual job. ---*/
  TensorProductVolumeDataQuad(TensorProductResVolDOFs, resDOFs.cols(), nInt1D, nDOFs1D,
                              legBasisLineIntTranspose, legBasisLineIntTranspose,
                              scalarDataInt, resDOFs, false, nullptr);
}

void CFEMStandardQuadVolumeSol::ResidualGradientBasisFunctions(vector<ColMajorMatrix<su2double> > &vectorDataInt,
                                                               ColMajorMatrix<su2double>          &resDOFs) {

  /*--- Call TensorProductVolumeDataHex 2 times with the appropriate arguments
        to carry out the actual job. ---*/
  TensorProductVolumeDataQuad(TensorProductResVolDOFs, resDOFs.cols(), nInt1D, nDOFs1D,
                              derLegBasisLineIntTranspose, legBasisLineIntTranspose,
                              vectorDataInt[0], resDOFs, false, nullptr);
  TensorProductVolumeDataQuad(TensorProductResVolDOFs, resDOFs.cols(), nInt1D, nDOFs1D,
                              legBasisLineIntTranspose, derLegBasisLineIntTranspose,
                              vectorDataInt[1], resDOFs, false, nullptr);
}

passivedouble CFEMStandardQuadVolumeSol::ValBasis0(void) {

  /*--- Easier storage of the 1D Legendre basis function
        and return the value of the 2D basis function.  ---*/
  const passivedouble leg1D = legBasisLineSolDOFs(0,0);
  return leg1D*leg1D;
}
