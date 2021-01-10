/*!
 * \file CFEMStandardQuadAdjacentPrismGrid.cpp
 * \brief Functions for the class CFEMStandardQuadAdjacentPrismGrid.
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

#include "../../include/fem/CFEMStandardQuadAdjacentPrismGrid.hpp"
#include "../../include/toolboxes/CSquareMatrixCM.hpp"

/*----------------------------------------------------------------------------------*/
/*            Public member functions of CFEMStandardQuadAdjacentPrismGrid.         */
/*----------------------------------------------------------------------------------*/

CFEMStandardQuadAdjacentPrismGrid::CFEMStandardQuadAdjacentPrismGrid(const unsigned short val_nPoly,
                                                                     const unsigned short val_orderExact,
                                                                     const unsigned short val_faceID_Elem,
                                                                     const unsigned short val_orientation,
                                                                     const bool           val_useLGL,
                                                                     CGemmBase           *val_gemm)
  : CFEMStandardPrismBase(),
    CFEMStandardQuadBase(val_nPoly, val_orderExact) {

  /*--- Convert the pointer for the gemm functionality. ---*/
  gemmDOFs2Int = dynamic_cast<CGemmStandard *> (val_gemm);
  if( !gemmDOFs2Int )
    SU2_MPI::Error(string("Dynamic cast failure. This should not happen"), CURRENT_FUNCTION);

  /*--- Determine the location of the grid DOFs. ---*/
  vector<passivedouble> rTriangleDOFs, sTriangleDOFs, rLineDOFs;
  if( val_useLGL ) {
    Location1DGridDOFsLGL(nPoly, rLineDOFs);
    LocationTriangleGridDOFsLGL(nPoly, rTriangleDOFs, sTriangleDOFs);
  }
  else {
    Location1DGridDOFsEquidistant(nPoly, rLineDOFs);
    LocationTriangleGridDOFsEquidistant(nPoly, rTriangleDOFs, sTriangleDOFs);
  }

  /*--- Determine the parametric coordinates of all DOFs of the prism. ---*/
  vector<passivedouble> rDOFs, sDOFs, tDOFs;
  LocationAllPointsPrism(rTriangleDOFs, sTriangleDOFs, rLineDOFs, rDOFs, sDOFs, tDOFs);

  /*--- Determine the inverse of the Vandermonde matrix of the DOFs. ---*/
  CSquareMatrixCM VInv(rDOFs.size());
  VandermondePrism(nPoly, rDOFs, sDOFs, tDOFs, VInv.GetMat());
  VInv.Invert();

  /*--- Convert the 2D parametric coordinates of the integration points of the
        quadrilateral face to the 3D parametric coordinates of the adjacent prism. ---*/
  vector<passivedouble> rInt, sInt, tInt;
  ConvertCoor2DQuadFaceTo3DPrism(rLineInt, val_faceID_Elem, val_orientation, rInt, sInt, tInt);

  /*--- Determine the Vandermonde matrix and its derivatives of the integration points.
        Make sure to allocate the number of rows to nIntTotPad and initialize them to zero. ---*/
  ColMajorMatrix<passivedouble> V(nIntegrationPad,  rDOFs.size()), VDr(nIntegrationPad,rDOFs.size()),
                                VDs(nIntegrationPad,rDOFs.size()), VDt(nIntegrationPad,rDOFs.size());
  V.setConstant(0.0);
  VDr.setConstant(0.0);
  VDs.setConstant(0.0);
  VDt.setConstant(0.0);

  VandermondePrism(nPoly, rInt, sInt, tInt, V);
  GradVandermondePrism(nPoly, rInt, sInt, tInt, VDr, VDs, VDt);

  /*--- The Lagrangian basis functions and its gradients can be obtained by
        multiplying V, VDr, VDs, VDt and VInv. ---*/
  derLagBasisInt.resize(3);

  VInv.MatMatMult('R', V,   lagBasisInt);
  VInv.MatMatMult('R', VDr, derLagBasisInt[0]);
  VInv.MatMatMult('R', VDs, derLagBasisInt[1]);
  VInv.MatMatMult('R', VDt, derLagBasisInt[2]);

  /*--- Check if the sum of the elements of the relevant rows of lagBasisInt
        is 1 and derLagBasisInt is 0. ---*/
  CheckRowSum(nIntegration, rDOFs.size(), 1.0, lagBasisInt);
  CheckRowSum(nIntegration, rDOFs.size(), 0.0, derLagBasisInt[0]);
  CheckRowSum(nIntegration, rDOFs.size(), 0.0, derLagBasisInt[1]);
  CheckRowSum(nIntegration, rDOFs.size(), 0.0, derLagBasisInt[2]);
}

void CFEMStandardQuadAdjacentPrismGrid::CoorIntPoints(const bool                notUsed,
                                                      ColMajorMatrix<su2double> &matCoorDOF,
                                                      ColMajorMatrix<su2double> &matCoorInt) {

  /*--- Call the general functionality of gemmDOFs2Int with the appropriate
        arguments to compute the coordinates in the integration points
        of the face. ---*/
  gemmDOFs2Int->DOFs2Int(lagBasisInt, 3, matCoorDOF, matCoorInt, nullptr);
}
