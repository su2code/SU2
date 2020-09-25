/*!
 * \file CFEMStandardVolumeTriGrid.cpp
 * \brief Functions for the class CFEMStandardVolumeTriGrid.
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

#include "../../include/fem/CFEMStandardVolumeTriGrid.hpp"
#include "../../include/toolboxes/CGeneralSquareMatrixCM.hpp"

/*----------------------------------------------------------------------------------*/
/*          Public member functions of CFEMStandardVolumeTriGrid.                   */
/*----------------------------------------------------------------------------------*/

CFEMStandardVolumeTriGrid::CFEMStandardVolumeTriGrid(const unsigned short val_nPoly,
                                                     const unsigned short val_orderExact)
  : CFEMStandardTri(val_nPoly, val_orderExact) {

  /*--- Compute the values of the Lagrangian basis functions in the integration
        points for both the equidistant and LGL point distribution. ---*/
  LagBasisIntPointsTriangle(rTriangleDOFsEqui, sTriangleDOFsEqui, lagBasisIntEqui);
  LagBasisIntPointsTriangle(rTriangleDOFsLGL,  sTriangleDOFsLGL,  lagBasisIntLGL);

  /*--- Compute the values of the derivatives of the Lagrangian basis functions in
        the integration points for both the equidistant and LGL point distribution. ---*/
  DerLagBasisIntPointsTriangle(rTriangleDOFsEqui, sTriangleDOFsEqui, derLagBasisIntEqui);
  DerLagBasisIntPointsTriangle(rTriangleDOFsLGL,  sTriangleDOFsLGL,  derLagBasisIntLGL);

  /*--- Set up the jitted gemm call, if supported. For this particular standard
        element the derivative of the coordinates are computed, which is 2. ---*/
  SetUpJittedGEMM(nIntegrationPad, 2, nDOFs);
}

void CFEMStandardVolumeTriGrid::DerivativesCoorVolumeIntPoints(const bool                         LGLDistribution,
                                                               ColMajorMatrix<su2double>          &matCoor,
                                                               vector<ColMajorMatrix<su2double> > &matDerCoor) {

  /*--- Check for which point distribution the derivatives must be computed. ---*/
  if( LGLDistribution ) {

    /*--- LGL distribution. Call the function OwnGemm 2 times to compute the derivatives
          of the Cartesian coordinates w.r.t. the three parametric coordinates. ---*/
    OwnGemm(nIntegrationPad, 2, nDOFs, derLagBasisIntLGL[0], matCoor, matDerCoor[0], nullptr);
    OwnGemm(nIntegrationPad, 2, nDOFs, derLagBasisIntLGL[1], matCoor, matDerCoor[1], nullptr);
  }
  else {

    /*--- Equidistant distribution. Call the function OwnGemm 2 times to compute the derivatives
          of the Cartesian coordinates w.r.t. the three parametric coordinates. ---*/
    OwnGemm(nIntegrationPad, 2, nDOFs, derLagBasisIntEqui[0], matCoor, matDerCoor[0], nullptr);
    OwnGemm(nIntegrationPad, 2, nDOFs, derLagBasisIntEqui[1], matCoor, matDerCoor[1], nullptr);
  }
}

void CFEMStandardVolumeTriGrid::DerLagBasisIntPointsTriangle(const vector<passivedouble>            &rDOFs,
                                                             const vector<passivedouble>            &sDOFs,
                                                             vector<ColMajorMatrix<passivedouble> > &derLag) {

  /*--- Determine the padded number of the total number of integration points. ---*/
  const unsigned short nIntTot    = rTriangleInt.size();
  const unsigned short nIntTotPad = ((nIntTot+vecLen-1)/vecLen)*vecLen;

  /*--- Determine the inverse of the Vandermonde matrix of the DOFs. ---*/
  CGeneralSquareMatrixCM VInv(rDOFs.size());
  VandermondeTriangle(rDOFs, sDOFs, VInv.GetMat());
  VInv.Invert();

  /*--- Determine the gradient of the Vandermonde matrix of the integration points. Make
        sure to allocate the number of rows to nIntTotPad and initialize them to zero. ---*/
  ColMajorMatrix<passivedouble> VDr(nIntTotPad,rDOFs.size()),
                                VDs(nIntTotPad,rDOFs.size());
  VDr.setConstant(0.0);
  VDs.setConstant(0.0);

  GradVandermondeTriangle(rTriangleInt, sTriangleInt, VDr, VDs);

  /*--- The gradients of the Lagrangian basis functions can be obtained by
        multiplying VDr, VDs and VInv. ---*/
  derLag.resize(2);
  VInv.MatMatMult('R', VDr, derLag[0]);
  VInv.MatMatMult('R', VDs, derLag[1]);

  /*--- Check if the sum of the elements of the relevant rows of derLag is 0. ---*/
  for(unsigned short i=0; i<nIntTot; ++i) {
    passivedouble rowSumDr = 0.0, rowSumDs = 0.0;
    for(unsigned short j=0; j<rDOFs.size(); ++j) {
      rowSumDr += derLag[0](i,j);
      rowSumDs += derLag[1](i,j);
    }

    assert(fabs(rowSumDr) < 1.e-6);
    assert(fabs(rowSumDs) < 1.e-6);
  }
}

void CFEMStandardVolumeTriGrid::LagBasisIntPointsTriangle(const vector<passivedouble>   &rDOFs,
                                                          const vector<passivedouble>   &sDOFs,
                                                          ColMajorMatrix<passivedouble> &lag) {

  /*--- Determine the padded number of the total number of integration points. ---*/
  const unsigned short nIntTot    = rTriangleInt.size();
  const unsigned short nIntTotPad = ((nIntTot+vecLen-1)/vecLen)*vecLen;

  /*--- Determine the inverse of the Vandermonde matrix of the DOFs. ---*/
  CGeneralSquareMatrixCM VInv(rDOFs.size());
  VandermondeTriangle(rDOFs, sDOFs, VInv.GetMat());
  VInv.Invert();

  /*--- Determine the Vandermonde matrix of the integration points. Make sure to
        allocate the number of rows to nIntTotPad and initialize them to zero. ---*/ 
  ColMajorMatrix<passivedouble> V(nIntTotPad,rDOFs.size());
  V.setConstant(0.0);
  VandermondeTriangle(rTriangleInt, sTriangleInt, V);

  /*--- The Lagrangian basis functions can be obtained by multiplying
        V and VInv. ---*/
  VInv.MatMatMult('R', V, lag);

  /*--- Check if the sum of the elements of the relevant rows of lag is 1. ---*/
  for(unsigned short i=0; i<nIntTot; ++i) {
    passivedouble rowSum = -1.0;
    for(unsigned short j=0; j<rDOFs.size(); ++j) rowSum += lag(i,j);
    assert(fabs(rowSum) < 1.e-6);
  }
}

void CFEMStandardVolumeTriGrid::GradVandermondeTriangle(const vector<passivedouble>   &r,
                                                        const vector<passivedouble>   &s,
                                                        ColMajorMatrix<passivedouble> &VDr,
                                                        ColMajorMatrix<passivedouble> &VDs) {

}

void CFEMStandardVolumeTriGrid::VandermondeTriangle(const vector<passivedouble>   &r,
                                                    const vector<passivedouble>   &s,
                                                    ColMajorMatrix<passivedouble> &V) {

}
