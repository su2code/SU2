/*!
 * \file CFEMStandardTriAdjacentPrismSol.cpp
 * \brief Functions for the class CFEMStandardTriAdjacentPrismSol.
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

#include "../../include/fem/CFEMStandardTriAdjacentPrismSol.hpp"

/*----------------------------------------------------------------------------------*/
/*            Public member functions of CFEMStandardTriAdjacentPrismSol.           */
/*----------------------------------------------------------------------------------*/

CFEMStandardTriAdjacentPrismSol::CFEMStandardTriAdjacentPrismSol(const unsigned short val_nPoly,
                                                                 const unsigned short val_orderExact,
                                                                 const unsigned short val_faceID_Elem,
                                                                 const unsigned short val_orientation,
                                                                 CGemmBase           *val_gemm_1,
                                                                 CGemmBase           *val_gemm_2)
  : CFEMStandardPrismBase(),
    CFEMStandardTriBase(val_nPoly, val_orderExact) {

  /*--- Store the faceID of the element and the orientation. ---*/
  faceID_Elem = val_faceID_Elem;
  orientation = val_orientation;

  /*--- Convert the pointers for the gemm functionalities. ---*/
  gemmDOFs2Int = dynamic_cast<CGemmStandard *> (val_gemm_1);
  gemmInt2DOFs = dynamic_cast<CGemmStandard *> (val_gemm_2);
  if(!gemmDOFs2Int || !gemmInt2DOFs)
    SU2_MPI::Error(string("Dynamic cast failure. This should not happen"), CURRENT_FUNCTION);

  /*--- Convert the 2D parametric coordinates of the integration points of the
        triangular face to the 3D parametric coordinates of the adjacent prism. ---*/
  vector<passivedouble> rTrianglePrism, sTrianglePrism, rLinePrism;
  ConvertCoor2DTriFaceTo3DPrism(rTriangleInt, sTriangleInt, val_faceID_Elem, val_orientation,
                                rTrianglePrism, sTrianglePrism, rLinePrism);

  /*--- Create the vector to store all parametric t-coordinates of the integration
        points of the triangular face. ---*/
  vector<passivedouble> tTrianglePrism(rTrianglePrism.size(), rLinePrism[0]);

  /*--- Allocate the memory for the Legendre basis functions and its
        1st derivatives in the integration points. ---*/
  nDOFs = (nPoly+1)*(nPoly+1)*(nPoly+2)/2;
  legBasisInt.resize(nIntegrationPad, nDOFs); legBasisInt.setConstant(0.0);

  derLegBasisInt.resize(3);
  derLegBasisInt[0].resize(nIntegrationPad, nDOFs); derLegBasisInt[0].setConstant(0.0);
  derLegBasisInt[1].resize(nIntegrationPad, nDOFs); derLegBasisInt[1].setConstant(0.0);
  derLegBasisInt[2].resize(nIntegrationPad, nDOFs); derLegBasisInt[2].setConstant(0.0);

  /*--- Compute the Legendre basis functions and its first
        derivatives in the integration points. ---*/
  VandermondePrism(nPoly, rTrianglePrism, sTrianglePrism, tTrianglePrism, legBasisInt);
  GradVandermondePrism(nPoly, rTrianglePrism, sTrianglePrism, tTrianglePrism,
                       derLegBasisInt[0], derLegBasisInt[1], derLegBasisInt[2]);

  /*--- Make sure that the padded values of legBasisInt are initialized properly
        to avoid problems. The gradients are set to zero, which is fine. ---*/
  for(unsigned short i=nIntegration; i<nIntegrationPad; ++i)
    for(unsigned short j=0; j<nDOFs; ++j)
      legBasisInt(i,j) = legBasisInt(0,j);

  /*--- Create the transpose of legBasisInt and derLegBasisInt. ---*/
  nDOFsPad = PaddedValue(nDOFs);
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

  /*--- The following two functions are mainly for writing surface output.
        Determine first the local connectivity of all faces ---*/
  CFEMStandardPrismBase::LocalGridConnFaces();

  /*--- Then determine local subconnectivity at the corresponding face ID.
        This will also instantiate the variable subConn1ForPlotting. ---*/
  CFEMStandardTriBase::SubConnLinearElementsFace(faceID_Elem);
}

void CFEMStandardTriAdjacentPrismSol::GradSolIntPoints(ColMajorMatrix<su2double>          &matSolDOF,
                                                       vector<ColMajorMatrix<su2double> > &matGradSolInt) {

  /*--- Call the general functionality of gemmDOFs2Int with the appropriate
        arguments to compute the gradient of the solution in the integration
        points of the face. ---*/
  gemmDOFs2Int->gemm(derLegBasisInt[0], matSolDOF.cols(), matSolDOF, matGradSolInt[0], nullptr);
  gemmDOFs2Int->gemm(derLegBasisInt[1], matSolDOF.cols(), matSolDOF, matGradSolInt[1], nullptr);
  gemmDOFs2Int->gemm(derLegBasisInt[2], matSolDOF.cols(), matSolDOF, matGradSolInt[2], nullptr);
}

void CFEMStandardTriAdjacentPrismSol::SolIntPoints(ColMajorMatrix<su2double> &matSolDOF,
                                                   ColMajorMatrix<su2double> &matSolInt) {

  /*--- Call the general functionality of gemmDOFs2Int with the appropriate
        arguments to compute the solution in the integration points
        of the face. ---*/
  gemmDOFs2Int->gemm(legBasisInt, matSolDOF.cols(), matSolDOF, matSolInt, nullptr);
}

void CFEMStandardTriAdjacentPrismSol::ResidualBasisFunctions(ColMajorMatrix<su2double> &scalarDataInt,
                                                             ColMajorMatrix<su2double> &resDOFs) {

  /*--- Call the generic functionality of gemmInt2DOFs with the appropriate
        arguments to compute the residuals in the DOFs of the adjacent element. ---*/
  gemmInt2DOFs->gemm(legBasisIntTranspose, scalarDataInt.cols(), scalarDataInt, resDOFs, nullptr);
}

void CFEMStandardTriAdjacentPrismSol::ResidualGradientBasisFunctions(vector<ColMajorMatrix<su2double> > &vectorDataInt,
                                                                     ColMajorMatrix<su2double>          &resDOFs) {

  /*--- Call the generic functionality of gemmInt2DOFs 3 times with the appropriate
        arguments to compute the residuals in the DOFs of the adjacent element. ---*/
  gemmInt2DOFs->gemm(derLegBasisIntTranspose[0], vectorDataInt[0].cols(), vectorDataInt[0], resDOFs, nullptr);
  gemmInt2DOFs->gemm(derLegBasisIntTranspose[1], vectorDataInt[1].cols(), vectorDataInt[1], resDOFs, nullptr);
  gemmInt2DOFs->gemm(derLegBasisIntTranspose[2], vectorDataInt[2].cols(), vectorDataInt[2], resDOFs, nullptr);
}
