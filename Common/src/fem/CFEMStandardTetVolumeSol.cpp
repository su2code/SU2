/*!
 * \file CFEMStandardTetVolumeSol.cpp
 * \brief Functions for the class CFEMStandardTetVolumeSol.
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

#include "../../include/fem/CFEMStandardTetVolumeSol.hpp"
#include "../../include/toolboxes/CSquareMatrixCM.hpp"

/*----------------------------------------------------------------------------------*/
/*             Public member functions of CFEMStandardTetVolumeSol.                 */
/*----------------------------------------------------------------------------------*/

CFEMStandardTetVolumeSol::CFEMStandardTetVolumeSol(const unsigned short val_nPoly,
                                                   const unsigned short val_orderExact,
                                                   const unsigned short val_locGridDOFs,
                                                   const unsigned short val_nVar)
  : CFEMStandardTetBase(val_nPoly, val_orderExact) {

  /*--- Determine the location of the nodal solution DOFs. ---*/
  if(val_locGridDOFs == LGL)
    LocationTetGridDOFsLGL(nPoly, rTetSolDOFs, sTetSolDOFs, tTetSolDOFs);
  else
    LocationTetGridDOFsEquidistant(nPoly, rTetSolDOFs, sTetSolDOFs, tTetSolDOFs);

  /*--- Allocate the memory for the Legendre basis functions and its
        1st and 2nd derivatives in the integration points. ---*/
  legBasisInt.resize(nIntegrationPad, nDOFs); legBasisInt.setConstant(0.0);

  derLegBasisInt.resize(3);
  derLegBasisInt[0].resize(nIntegrationPad, nDOFs); derLegBasisInt[0].setConstant(0.0);
  derLegBasisInt[1].resize(nIntegrationPad, nDOFs); derLegBasisInt[1].setConstant(0.0);
  derLegBasisInt[2].resize(nIntegrationPad, nDOFs); derLegBasisInt[2].setConstant(0.0);

  hesLegBasisInt.resize(6);
  hesLegBasisInt[0].resize(nIntegrationPad, nDOFs); hesLegBasisInt[0].setConstant(0.0);
  hesLegBasisInt[1].resize(nIntegrationPad, nDOFs); hesLegBasisInt[1].setConstant(0.0);
  hesLegBasisInt[2].resize(nIntegrationPad, nDOFs); hesLegBasisInt[2].setConstant(0.0);
  hesLegBasisInt[3].resize(nIntegrationPad, nDOFs); hesLegBasisInt[3].setConstant(0.0);
  hesLegBasisInt[4].resize(nIntegrationPad, nDOFs); hesLegBasisInt[4].setConstant(0.0);
  hesLegBasisInt[5].resize(nIntegrationPad, nDOFs); hesLegBasisInt[5].setConstant(0.0);

  /*--- Compute the Legendre basis functions and its first
        and second derivatives in the integration points. ---*/
  VandermondeTetrahedron(nPoly, rTetInt, sTetInt, tTetInt, legBasisInt);
  GradVandermondeTetrahedron(nPoly, rTetInt, sTetInt, tTetInt, derLegBasisInt[0],
                             derLegBasisInt[1], derLegBasisInt[2]);
  HesVandermondeTetrahedron(nPoly, rTetInt, sTetInt, tTetInt, hesLegBasisInt[0],
                            hesLegBasisInt[1], hesLegBasisInt[2], hesLegBasisInt[3],
                            hesLegBasisInt[4], hesLegBasisInt[5]);

  /*--- Fill the padded entries of the matrices just computed. ---*/
  for(unsigned short j=0; j<nDOFs; ++j) {
    for(unsigned short i=nIntegration; i<nIntegrationPad; ++i) {
      legBasisInt(i,j)       = legBasisInt(0,j);

      derLegBasisInt[0](i,j) = derLegBasisInt[0](0,j);
      derLegBasisInt[1](i,j) = derLegBasisInt[1](0,j);
      derLegBasisInt[2](i,j) = derLegBasisInt[2](0,j);

      hesLegBasisInt[0](i,j) = hesLegBasisInt[0](0,j);
      hesLegBasisInt[1](i,j) = hesLegBasisInt[1](0,j);
      hesLegBasisInt[2](i,j) = hesLegBasisInt[2](0,j);
      hesLegBasisInt[3](i,j) = hesLegBasisInt[3](0,j);
      hesLegBasisInt[4](i,j) = hesLegBasisInt[4](0,j);
      hesLegBasisInt[5](i,j) = hesLegBasisInt[5](0,j);
    }
  }

  /*--- Create the transpose of legBasisInt and derLegBasisInt. ---*/
  legBasisIntTranspose.resize(nDOFsPad, nIntegration); legBasisIntTranspose.setConstant(0.0);

  derLegBasisIntTranspose.resize(3);
  derLegBasisIntTranspose[0].resize(nDOFsPad, nIntegration); derLegBasisIntTranspose[0].setConstant(0.0);
  derLegBasisIntTranspose[1].resize(nDOFsPad, nIntegration); derLegBasisIntTranspose[1].setConstant(0.0);
  derLegBasisIntTranspose[2].resize(nDOFsPad, nIntegration); derLegBasisIntTranspose[2].setConstant(0.0);

  for(unsigned short j=0; j<nIntegration; ++j) {
    for(unsigned short i=0; i<nDOFs; ++i) {
      legBasisIntTranspose(i,j)       = legBasisInt(j,i);
      derLegBasisIntTranspose[0](i,j) = derLegBasisInt[0](j,i);
      derLegBasisIntTranspose[1](i,j) = derLegBasisInt[1](j,i);
      derLegBasisIntTranspose[2](i,j) = derLegBasisInt[2](j,i);
    }
  }

  /*--- Compute the Legendre basis functions in the solution DOFs.
        Store it in a square matrix as also the inverse is needed. ---*/
  CSquareMatrixCM Vtmp(nDOFs);
  VandermondeTetrahedron(nPoly, rTetSolDOFs, sTetSolDOFs, tTetSolDOFs, Vtmp.GetMat());

  /*--- Store the contents of Vtmp in legBasisLineSolDOFs. ---*/
  legBasisSolDOFs.resize(nDOFs, nDOFs);

  for(unsigned short j=0; j<nDOFs; ++j)
    for(unsigned short i=0; i<nDOFs; ++i)
      legBasisSolDOFs(i,j) = Vtmp(i,j);

  /*--- Compute the inverse of Vtmp and store the contents in legBasisLineSolDOFsInv. ---*/
  Vtmp.Invert();
  legBasisSolDOFsInv.resize(nDOFs, nDOFs);

  for(unsigned short j=0; j<nDOFs; ++j)
    for(unsigned short i=0; i<nDOFs; ++i)
      legBasisSolDOFsInv(i,j) = Vtmp(i,j);

  /*--- Compute the first derivatives of the basis
        functions in the solution DOFs. ---*/
  derLegBasisSolDOFs.resize(3);
  derLegBasisSolDOFs[0].resize(nDOFs, nDOFs);
  derLegBasisSolDOFs[1].resize(nDOFs, nDOFs);
  derLegBasisSolDOFs[2].resize(nDOFs, nDOFs);

  GradVandermondeTetrahedron(nPoly, rTetSolDOFs, sTetSolDOFs, tTetSolDOFs,
                             derLegBasisSolDOFs[0], derLegBasisSolDOFs[1],
                             derLegBasisSolDOFs[2]);

  /*--- Determine the local subconnectivity of this standard element when split
        in several linear elements. Used for a.o. plotting and searching. ---*/
  CFEMStandardTetBase::SubConnLinearElements();

  /*--- Set up the jitted gemm calls, if supported. ---*/
  SetUpJittedGEMM(nIntegrationPad, val_nVar, nDOFs, nIntegrationPad,
                  nDOFs, nIntegrationPad, true, jitterDOFs2Int, gemmDOFs2Int);
  SetUpJittedGEMM(nIntegrationPad, val_nVar, nDOFs, nIntegrationPad,
                  nDOFsPad, nIntegrationPad, true, jitterDOFsPad2Int, gemmDOFsPad2Int);
  SetUpJittedGEMM(nDOFs, val_nVar, nDOFs, nDOFs, nDOFs,
                  nDOFs, true, jitterDOFs2SolDOFs, gemmDOFs2SolDOFs);
  SetUpJittedGEMM(nDOFsPad, val_nVar, nIntegration, nDOFsPad, nIntegrationPad,
                  nDOFsPad, false, jitterInt2DOFs, gemmInt2DOFs);
}

CFEMStandardTetVolumeSol::~CFEMStandardTetVolumeSol() {

#if defined(PRIMAL_SOLVER) && defined(HAVE_MKL)
  if( jitterDOFs2Int ) {
    mkl_jit_destroy(jitterDOFs2Int);
    jitterDOFs2Int = nullptr;
  }

  if( jitterDOFsPad2Int ) {
    mkl_jit_destroy(jitterDOFsPad2Int);
    jitterDOFsPad2Int = nullptr;
  }

  if( jitterDOFs2SolDOFs ) {
    mkl_jit_destroy(jitterDOFs2SolDOFs);
    jitterDOFs2SolDOFs = nullptr;
  }

  if( jitterInt2DOFs ) {
    mkl_jit_destroy(jitterInt2DOFs);
    jitterInt2DOFs = nullptr;
  }
#endif
}

passivedouble CFEMStandardTetVolumeSol::GetFactorInviscidSpectralRadius(const unsigned short nPolyElem) const {

  passivedouble factInviscidRad = 0.0;
  switch( nPolyElem ) {
    case 0: factInviscidRad =  3.0; break;
    case 1: factInviscidRad =  6.0; break;
    case 2: factInviscidRad = 12.0; break;
    case 3: factInviscidRad = 20.0; break;
    case 4: factInviscidRad = 28.0; break;
    case 5: factInviscidRad = 38.0; break;
    case 6: factInviscidRad = 50.0; break;
    case 7: factInviscidRad = 64.0; break;
    case 8: factInviscidRad = 80.0; break;
    case 9: factInviscidRad = 98.0; break;
    default:
      SU2_MPI::Error(string("Polynomial order not foreseen"), CURRENT_FUNCTION);
  }
  return factInviscidRad;
}

passivedouble CFEMStandardTetVolumeSol::GetFactorViscousSpectralRadius(const unsigned short nPolyElem) const {
  passivedouble factViscousRad = 0.0;
  switch( nPolyElem ) {
    case 0: factViscousRad =     6.0; break;
    case 1: factViscousRad =    36.0; break;
    case 2: factViscousRad =   150.0; break;
    case 3: factViscousRad =   420.0; break;
    case 4: factViscousRad =   980.0; break;
    case 5: factViscousRad =  1975.0; break;
    case 6: factViscousRad =  3575.0; break;
    case 7: factViscousRad =  7000.0; break;
    case 8: factViscousRad = 14000.0; break;
    case 9: factViscousRad = 28000.0; break;
    default:
      SU2_MPI::Error(string("Polynomial order not foreseen"), CURRENT_FUNCTION);
  }
  return factViscousRad;
}

void CFEMStandardTetVolumeSol::BasisFunctionsInPoints(const vector<vector<passivedouble> > &parCoor,
                                                      ColMajorMatrix<passivedouble>        &matBasis) {

  /*--- Allocate the memory for matBasis. ---*/
  const unsigned short nCoor = parCoor[0].size();
  matBasis.resize(nCoor, nDOFs);

  /*--- Determine the values of the basis functions in the given
        parametric coordinates. ---*/
  VandermondeTetrahedron(nPoly, parCoor[0], parCoor[1], parCoor[2], matBasis);
}

void CFEMStandardTetVolumeSol::ModalToNodal(ColMajorMatrix<su2double> &solDOFs) {

  /*--- Copy solDOFs into tmp and carry out the GEMM call for
        the conversion to the nodal formulation. ---*/
  ColMajorMatrix<su2double> tmp = solDOFs;

  OwnGemm(gemmDOFs2SolDOFs, jitterDOFs2SolDOFs, nDOFs, tmp.cols(), nDOFs,
          nDOFs, nDOFs, nDOFs, true, legBasisSolDOFs, tmp, solDOFs, nullptr);
}

void CFEMStandardTetVolumeSol::NodalToModal(ColMajorMatrix<su2double> &solDOFs) {

  /*--- Copy solDOFs into tmp and carry out the GEMM call for
        the conversion to the modal formulation. ---*/
  ColMajorMatrix<su2double> tmp = solDOFs;

  OwnGemm(gemmDOFs2SolDOFs, jitterDOFs2SolDOFs, nDOFs, tmp.cols(), nDOFs,
          nDOFs, nDOFs, nDOFs, true, legBasisSolDOFsInv, tmp, solDOFs, nullptr);
}

void CFEMStandardTetVolumeSol::GradSolIntPoints(ColMajorMatrix<su2double>          &matSolDOF,
                                                vector<ColMajorMatrix<su2double> > &matGradSolInt) {

  /*--- Call OwnGemm with the appropriate arguments to compute the data. ---*/
  for(unsigned short nn=0; nn<3; ++nn)
    OwnGemm(gemmDOFs2Int, jitterDOFs2Int, nIntegrationPad, matSolDOF.cols(), nDOFs,
            nIntegrationPad, nDOFs, nIntegrationPad, true,
            derLegBasisInt[nn], matSolDOF, matGradSolInt[nn], nullptr);
}

void CFEMStandardTetVolumeSol::SolIntPoints(ColMajorMatrix<su2double> &matSolDOF,
                                            ColMajorMatrix<su2double> &matSolInt) {

  /*--- Call OwnGemm with the appropriate arguments to carry out the actual job. ---*/
  OwnGemm(gemmDOFs2Int, jitterDOFs2Int, nIntegrationPad, matSolDOF.cols(), nDOFs,
          nIntegrationPad, nDOFs, nIntegrationPad, true,
          legBasisInt, matSolDOF, matSolInt, nullptr);
}

void CFEMStandardTetVolumeSol::SolIntPointsDOFsPadded(ColMajorMatrix<su2double> &matSolDOF,
                                                      ColMajorMatrix<su2double> &matSolInt) {

  /*--- Call OwnGemm with the appropriate arguments to carry out the actual job. ---*/
  OwnGemm(gemmDOFsPad2Int, jitterDOFsPad2Int, nIntegrationPad, matSolDOF.cols(),
          nDOFs, nIntegrationPad, nDOFsPad, nIntegrationPad, true,
          legBasisInt, matSolDOF, matSolInt, nullptr);
}

void CFEMStandardTetVolumeSol::ResidualBasisFunctions(ColMajorMatrix<su2double> &scalarDataInt,
                                                      ColMajorMatrix<su2double> &resDOFs) {

  /*--- Call OwnGemm with the appropriate arguments to carry out the actual job. ---*/
  OwnGemm(gemmInt2DOFs, jitterInt2DOFs, nDOFsPad, resDOFs.cols(), nIntegration,
          nDOFsPad, nIntegrationPad, nDOFsPad, false,
          legBasisIntTranspose, scalarDataInt, resDOFs, nullptr);
}

void CFEMStandardTetVolumeSol::ResidualGradientBasisFunctions(vector<ColMajorMatrix<su2double> > &vectorDataInt,
                                                              ColMajorMatrix<su2double>          &resDOFs) {

  /*--- Call OwnGemm 3 times with the appropriate arguments to carry out the actual job. ---*/
  OwnGemm(gemmInt2DOFs, jitterInt2DOFs, nDOFsPad, resDOFs.cols(), nIntegration,
          nDOFsPad, nIntegrationPad, nDOFsPad, false,
          derLegBasisIntTranspose[0], vectorDataInt[0], resDOFs, nullptr);
  OwnGemm(gemmInt2DOFs, jitterInt2DOFs, nDOFsPad, resDOFs.cols(), nIntegration,
          nDOFsPad, nIntegrationPad, nDOFsPad, false,
          derLegBasisIntTranspose[1], vectorDataInt[1], resDOFs, nullptr);
  OwnGemm(gemmInt2DOFs, jitterInt2DOFs, nDOFsPad, resDOFs.cols(), nIntegration,
          nDOFsPad, nIntegrationPad, nDOFsPad, false,
          derLegBasisIntTranspose[2], vectorDataInt[2], resDOFs, nullptr);
}
