/*!
 * \file CFEMStandardTriVolumeSol.cpp
 * \brief Functions for the class CFEMStandardTriVolumeSol.
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

#include "../../include/fem/CFEMStandardTriVolumeSol.hpp"
#include "../../include/toolboxes/CSquareMatrixCM.hpp"

/*----------------------------------------------------------------------------------*/
/*             Public member functions of CFEMStandardTriVolumeSol.                 */
/*----------------------------------------------------------------------------------*/

CFEMStandardTriVolumeSol::CFEMStandardTriVolumeSol(const unsigned short val_nPoly,
                                                   const unsigned short val_orderExact,
                                                   const unsigned short val_locGridDOFs,
                                                   const unsigned short val_nVar)
  : CFEMStandardTriBase(val_nPoly, val_orderExact) {

  /*--- Determine the location of the nodal solution DOFs. ---*/
  if(val_locGridDOFs == LGL)
    LocationTriangleGridDOFsLGL(nPoly, rTriangleSolDOFs, sTriangleSolDOFs);
  else
    LocationTriangleGridDOFsEquidistant(nPoly, rTriangleSolDOFs, sTriangleSolDOFs);

  /*--- Allocate the memory for the Legendre basis functions and its
        1st and 2nd derivatives in the integration points. ---*/
  legBasisInt.resize(nIntegrationPad, nDOFs);

  derLegBasisInt.resize(2);
  derLegBasisInt[0].resize(nIntegrationPad, nDOFs);
  derLegBasisInt[1].resize(nIntegrationPad, nDOFs);

  hesLegBasisInt.resize(3);
  hesLegBasisInt[0].resize(nIntegrationPad, nDOFs);
  hesLegBasisInt[1].resize(nIntegrationPad, nDOFs);
  hesLegBasisInt[2].resize(nIntegrationPad, nDOFs);

  /*--- Compute the Legendre basis functions and its first
        and second derivatives in the integration points. ---*/
  VandermondeTriangle(nPoly, rTriangleInt, sTriangleInt, legBasisInt);
  GradVandermondeTriangle(nPoly, rTriangleInt, sTriangleInt, derLegBasisInt[0], derLegBasisInt[1]);
  HesVandermondeTriangle(nPoly, rTriangleInt, sTriangleInt, hesLegBasisInt[0],
                         hesLegBasisInt[1], hesLegBasisInt[2]);

  /*--- Fill the padded entries of the matrices just computed. ---*/
  for(unsigned short j=0; j<nDOFs; ++j) {
    for(unsigned short i=nIntegration; i<nIntegrationPad; ++i) {
      legBasisInt(i,j)       = legBasisInt(0,j);

      derLegBasisInt[0](i,j) = derLegBasisInt[0](0,j);
      derLegBasisInt[1](i,j) = derLegBasisInt[1](0,j);

      hesLegBasisInt[0](i,j) = hesLegBasisInt[0](0,j);
      hesLegBasisInt[1](i,j) = hesLegBasisInt[1](0,j);
      hesLegBasisInt[2](i,j) = hesLegBasisInt[2](0,j);
    }
  }

  /*--- Create the transpose of legBasisInt and derLegBasisInt. ---*/
  legBasisIntTranspose.resize(nDOFsPad, nIntegration); legBasisIntTranspose.setConstant(0.0);

  derLegBasisIntTranspose.resize(2);
  derLegBasisIntTranspose[0].resize(nDOFsPad, nIntegration); derLegBasisIntTranspose[0].setConstant(0.0);
  derLegBasisIntTranspose[1].resize(nDOFsPad, nIntegration); derLegBasisIntTranspose[1].setConstant(0.0);

  for(unsigned short j=0; j<nIntegration; ++j) {
    for(unsigned short i=0; i<nDOFs; ++i) {
      legBasisIntTranspose(i,j)       = legBasisInt(j,i);
      derLegBasisIntTranspose[0](i,j) = derLegBasisInt[0](j,i);
      derLegBasisIntTranspose[1](i,j) = derLegBasisInt[1](j,i);
    }
  }

  /*--- Compute the Legendre basis functions in the solution DOFs.
        Store it in a square matrix as also the inverse is needed. ---*/
  CSquareMatrixCM Vtmp(nDOFs);
  VandermondeTriangle(nPoly, rTriangleSolDOFs, sTriangleSolDOFs, Vtmp.GetMat());

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
  derLegBasisSolDOFs.resize(2);
  derLegBasisSolDOFs[0].resize(nDOFs, nDOFs);
  derLegBasisSolDOFs[1].resize(nDOFs, nDOFs);

  GradVandermondeTriangle(nPoly, rTriangleSolDOFs, sTriangleSolDOFs,
                          derLegBasisSolDOFs[0], derLegBasisSolDOFs[1]);

  /*--- Determine the local subconnectivity of this standard element when split
        in several linear elements. Used for a.o. plotting and searcing. ---*/
  CFEMStandardTriBase::SubConnLinearElements();

  /*--- Set up the jitted gemm calls, if supported. ---*/
  SetUpJittedGEMM(nIntegrationPad, val_nVar, nDOFs, nIntegrationPad,
                  nDOFs, nIntegrationPad, true, jitterDOFs2Int, gemmDOFs2Int);
  SetUpJittedGEMM(nDOFs, val_nVar, nDOFs, nDOFs, nDOFs,
                  nDOFs, true, jitterDOFs2SolDOFs, gemmDOFs2SolDOFs);
  SetUpJittedGEMM(nDOFsPad, val_nVar, nIntegration, nDOFsPad, nIntegrationPad,
                  nDOFsPad, false, jitterInt2DOFs, gemmInt2DOFs);

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

CFEMStandardTriVolumeSol::~CFEMStandardTriVolumeSol() {

#if defined(PRIMAL_SOLVER) && defined(HAVE_MKL)
if( jitterDOFs2Int ) {
    mkl_jit_destroy(jitterDOFs2Int);
    jitterDOFs2Int = nullptr;
  }

  if( jitterDOFs2SolDOFs ) {
    mkl_jit_destroy(jitterDOFs2SolDOFs);
    jitterDOFs2SolDOFs = nullptr;
  }  
#endif
}

void CFEMStandardTriVolumeSol::BasisFunctionsInPoints(const vector<vector<passivedouble> > &parCoor,
                                                      ColMajorMatrix<passivedouble>        &matBasis) {

  /*--- Allocate the memory for matBasis. ---*/
  const unsigned short nCoor = parCoor[0].size();
  matBasis.resize(nCoor, nDOFs);

  /*--- Determine the values of the basis functions in the given
        parametric coordinates. ---*/
  VandermondeTriangle(nPoly, parCoor[0], parCoor[1], matBasis);
}

void CFEMStandardTriVolumeSol::ModalToNodal(ColMajorMatrix<su2double> &solDOFs) {

  /*--- Copy solDOFs into tmp and carry out the GEMM call for
        the conversion to the nodal formulation. ---*/
  ColMajorMatrix<su2double> tmp = solDOFs;

  OwnGemm(gemmDOFs2SolDOFs, jitterDOFs2SolDOFs, nDOFs, tmp.cols(), nDOFs,
          nDOFs, nDOFs, nDOFs, true, legBasisSolDOFs, tmp, solDOFs, nullptr);
}

void CFEMStandardTriVolumeSol::NodalToModal(ColMajorMatrix<su2double> &solDOFs) {

  /*--- Copy solDOFs into tmp and carry out the GEMM call for
        the conversion to the modal formulation. ---*/
  ColMajorMatrix<su2double> tmp = solDOFs;

  OwnGemm(gemmDOFs2SolDOFs, jitterDOFs2SolDOFs, nDOFs, tmp.cols(), nDOFs,
          nDOFs, nDOFs, nDOFs, true, legBasisSolDOFsInv, tmp, solDOFs, nullptr);
}

void CFEMStandardTriVolumeSol::GradSolIntPoints(ColMajorMatrix<su2double>          &matSolDOF,
                                                vector<ColMajorMatrix<su2double> > &matGradSolInt) {

  /*--- Call OwnGemm with the appropriate arguments to compute the data. ---*/
  for(unsigned short nn=0; nn<2; ++nn)
    OwnGemm(gemmDOFs2Int, jitterDOFs2Int, nIntegrationPad, matSolDOF.cols(), nDOFs,
            nIntegrationPad, nDOFs, nIntegrationPad, true,
            derLegBasisInt[nn], matSolDOF, matGradSolInt[nn], nullptr);
}

void CFEMStandardTriVolumeSol::SolIntPoints(ColMajorMatrix<su2double> &matSolDOF,
                                            ColMajorMatrix<su2double> &matSolInt) {

  /*--- Call OwnGemm with the appropriate arguments to carry out the actual job. ---*/
  OwnGemm(gemmDOFs2Int, jitterDOFs2Int, nIntegrationPad, matSolDOF.cols(), nDOFs,
          nIntegrationPad, nDOFs, nIntegrationPad, true,
          legBasisInt, matSolDOF, matSolInt, nullptr);
}

void CFEMStandardTriVolumeSol::ResidualBasisFunctions(ColMajorMatrix<su2double> &scalarDataInt,
                                                      ColMajorMatrix<su2double> &resDOFs) {

  /*--- Call OwnGemm with the appropriate arguments to carry out the actual job. ---*/
  OwnGemm(gemmInt2DOFs, jitterInt2DOFs, nDOFsPad, resDOFs.cols(), nIntegration,
          nDOFsPad, nIntegrationPad, nDOFsPad, false,
          legBasisIntTranspose, scalarDataInt, resDOFs, nullptr);
}

void CFEMStandardTriVolumeSol::ResidualGradientBasisFunctions(vector<ColMajorMatrix<su2double> > &vectorDataInt,
                                                              ColMajorMatrix<su2double>          &resDOFs) {

  /*--- Call OwnGemm 2 times with the appropriate arguments to carry out the actual job. ---*/
  OwnGemm(gemmInt2DOFs, jitterInt2DOFs, nDOFsPad, resDOFs.cols(), nIntegration,
          nDOFsPad, nIntegrationPad, nDOFsPad, false,
          derLegBasisIntTranspose[0], vectorDataInt[0], resDOFs, nullptr);
  OwnGemm(gemmInt2DOFs, jitterInt2DOFs, nDOFsPad, resDOFs.cols(), nIntegration,
          nDOFsPad, nIntegrationPad, nDOFsPad, false,
          derLegBasisIntTranspose[1], vectorDataInt[1], resDOFs, nullptr);
}
