/*!
 * \file CFEMStandardLineAdjacentTriSol.cpp
 * \brief Functions for the class CFEMStandardLineAdjacentTriSol.
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

#include "../../include/fem/CFEMStandardLineAdjacentTriSol.hpp"

/*----------------------------------------------------------------------------------*/
/*            Public member functions of CFEMStandardLineAdjacentTriSol.            */
/*----------------------------------------------------------------------------------*/

CFEMStandardLineAdjacentTriSol::CFEMStandardLineAdjacentTriSol(const unsigned short val_nPoly,
                                                               const unsigned short val_orderExact,
                                                               const unsigned short val_faceID_Elem,
                                                               const unsigned short val_orientation,
                                                               CGemmBase           *val_gemm_1,
                                                               CGemmBase           *val_gemm_2)
  : CFEMStandardTriBase(),
    CFEMStandardLineBase(val_nPoly, val_orderExact) {

  /*--- Store the faceID of the element and the orientation. ---*/
  faceID_Elem = val_faceID_Elem;
  orientation = val_orientation;

  /*--- Convert the pointers for the gemm functionalities. ---*/
  gemmDOFs2Int = dynamic_cast<CGemmStandard *> (val_gemm_1);
  gemmInt2DOFs = dynamic_cast<CGemmStandard *> (val_gemm_2);
  if(!gemmDOFs2Int || !gemmInt2DOFs)
    SU2_MPI::Error(string("Dynamic cast failure. This should not happen"), CURRENT_FUNCTION);

  /*--- Convert the 1D parametric coordinates of the integration points of the
        face to the 2D parametric coordinates of the adjacent triangle. ---*/
  ConvertCoor1DFaceTo2DTriangle(rLineInt, val_faceID_Elem, val_orientation,
                                rTriangleInt, sTriangleInt);

  /*--- Allocate the memory for the Legendre basis functions and its
        1st derivatives in the integration points. ---*/
  nDOFs = (nPoly+1)*(nPoly+2)/2;
  legBasisInt.resize(nIntegrationPad, nDOFs); legBasisInt.setConstant(0.0);

  derLegBasisInt.resize(2);
  derLegBasisInt[0].resize(nIntegrationPad, nDOFs); derLegBasisInt[0].setConstant(0.0);
  derLegBasisInt[1].resize(nIntegrationPad, nDOFs); derLegBasisInt[1].setConstant(0.0);

  /*--- Compute the Legendre basis functions and its first
        derivatives in the integration points. ---*/
  VandermondeTriangle(nPoly, rTriangleInt, sTriangleInt, legBasisInt);
  GradVandermondeTriangle(nPoly, rTriangleInt, sTriangleInt, derLegBasisInt[0], derLegBasisInt[1]);

  /*--- Make sure that the padded values of legBasisInt are initialized properly
        to avoid problems. The gradients are set to zero, which is fine. ---*/
  for(unsigned short i=nIntegration; i<nIntegrationPad; ++i)
    for(unsigned short j=0; j<nDOFs; ++j)
      legBasisInt(i,j) = legBasisInt(0,j);
}

void CFEMStandardLineAdjacentTriSol::GradSolIntPoints(ColMajorMatrix<su2double>          &matSolDOF,
                                                      vector<ColMajorMatrix<su2double> > &matGradSolInt) {

  /*--- Call the general functionality of gemmDOFs2Int with the appropriate
        arguments to compute the derivatives of the solution in the
        integration points of the face. ---*/
  gemmDOFs2Int->DOFs2Int(derLegBasisInt[0], matSolDOF.cols(), matSolDOF, matGradSolInt[0], nullptr);
  gemmDOFs2Int->DOFs2Int(derLegBasisInt[1], matSolDOF.cols(), matSolDOF, matGradSolInt[1], nullptr);
}

void CFEMStandardLineAdjacentTriSol::SolIntPoints(ColMajorMatrix<su2double> &matSolDOF,
                                                  ColMajorMatrix<su2double> &matSolInt) {

  /*--- Call the general functionality of gemmDOFs2Int with the appropriate
        arguments to compute the solution in the integration points
        of the face. ---*/
  gemmDOFs2Int->DOFs2Int(legBasisInt, matSolDOF.cols(), matSolDOF, matSolInt, nullptr);
}
