/*!
 * \file CFEMStandardTriGrid.cpp
 * \brief Functions for the class CFEMStandardTriGrid.
 * \author E. van der Weide
 * \version 7.0.7 "Blackbird"
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

#include "../../include/fem/CFEMStandardTriGrid.hpp"
#include "../../include/toolboxes/CGeneralSquareMatrixCM.hpp"

/*----------------------------------------------------------------------------------*/
/*             Public member functions of CFEMStandardTriGrid.                      */
/*----------------------------------------------------------------------------------*/

CFEMStandardTriGrid::CFEMStandardTriGrid(const unsigned short val_nPoly,
                                         const unsigned short val_orderExact,
                                         const bool           val_surfElement)
  : CFEMStandardTri(val_nPoly, val_orderExact) {

  /*--- Determine the number of space dimensions, which is 3 if this standard
        element is a surface element and 2 when it is a value element. ---*/
  nDim = val_surfElement ? 3 : 2;

  /*--- Compute the values of the Lagrangian basis functions in the integration
        points for both the equidistant and LGL point distribution. ---*/
  LagBasisIntPointsTriangle(rTriangleDOFsEqui, sTriangleDOFsEqui, lagBasisIntEqui);
  LagBasisIntPointsTriangle(rTriangleDOFsLGL,  sTriangleDOFsLGL,  lagBasisIntLGL);

  /*--- Compute the values of the derivatives of the Lagrangian basis functions in
        the integration points for both the equidistant and LGL point distribution. ---*/
  DerLagBasisIntPointsTriangle(rTriangleDOFsEqui, sTriangleDOFsEqui, derLagBasisIntEqui);
  DerLagBasisIntPointsTriangle(rTriangleDOFsLGL,  sTriangleDOFsLGL,  derLagBasisIntLGL);

  /*--- Create the local grid connectivities of the faces of the volume element.
        Only needed if this is a volume element. ---*/
  if( !val_surfElement ) LocalGridConnFaces();

  /*--- Set up the jitted gemm call, if supported. For this particular standard
        element the derivatives of the coordinates are computed, which is nDim. ---*/
  SetUpJittedGEMM(nIntegrationPad, nDim, nDOFs);
}

void CFEMStandardTriGrid::DerivativesCoorIntPoints(const bool                         LGLDistribution,
                                                   ColMajorMatrix<su2double>          &matCoor,
                                                   vector<ColMajorMatrix<su2double> > &matDerCoor) {

  /*--- Check for which point distribution the derivatives must be computed. ---*/
  if( LGLDistribution ) {

    /*--- LGL distribution. Call the function OwnGemm 2 times to compute the derivatives
          of the Cartesian coordinates w.r.t. the two parametric coordinates. The second
          argument in the function call is nDim, which corresponds to the number of Cartesian
          coordinates (3 for a surface element and 2 for a volume element). ---*/
    OwnGemm(nIntegrationPad, nDim, nDOFs, derLagBasisIntLGL[0], matCoor, matDerCoor[0], nullptr);
    OwnGemm(nIntegrationPad, nDim, nDOFs, derLagBasisIntLGL[1], matCoor, matDerCoor[1], nullptr);
  }
  else {

    /*--- LGL distribution. Call the function OwnGemm 2 times to compute the derivatives
          of the Cartesian coordinates w.r.t. the two parametric coordinates. The second
          argument in the function call is nDim, which corresponds to the number of Cartesian
          coordinates (3 for a surface element and 2 for a volume element). ---*/
    OwnGemm(nIntegrationPad, nDim, nDOFs, derLagBasisIntEqui[0], matCoor, matDerCoor[0], nullptr);
    OwnGemm(nIntegrationPad, nDim, nDOFs, derLagBasisIntEqui[1], matCoor, matDerCoor[1], nullptr);
  }
}

void CFEMStandardTriGrid::DerLagBasisIntPointsTriangle(const vector<passivedouble>            &rDOFs,
                                                       const vector<passivedouble>            &sDOFs,
                                                       vector<ColMajorMatrix<passivedouble> > &derLag) {

  /*--- Determine the padded number of the total number of integration points. ---*/
  const unsigned short nIntTot    = rTriangleInt.size();
  const unsigned short nIntTotPad = ((nIntTot+baseVectorLen-1)/baseVectorLen)*baseVectorLen;

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

void CFEMStandardTriGrid::LagBasisIntPointsTriangle(const vector<passivedouble>   &rDOFs,
                                                    const vector<passivedouble>   &sDOFs,
                                                    ColMajorMatrix<passivedouble> &lag) {

  /*--- Determine the padded number of the total number of integration points. ---*/
  const unsigned short nIntTot    = rTriangleInt.size();
  const unsigned short nIntTotPad = ((nIntTot+baseVectorLen-1)/baseVectorLen)*baseVectorLen;

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

void CFEMStandardTriGrid::GradVandermondeTriangle(const vector<passivedouble>   &r,
                                                  const vector<passivedouble>   &s,
                                                  ColMajorMatrix<passivedouble> &VDr,
                                                  ColMajorMatrix<passivedouble> &VDs) {

  /*--- For a triangle the orthogonal basis for the reference element is obtained
        by a combination of a Jacobi polynomial and a Legendre polynomial. This
        is the result of the orthonormalization of the monomial basis.
        Note that the sequence of the i and j loop must be identical to
        the evaluation of the Vandermonde matrix itself. ---*/
  unsigned short ii = 0;
  for(unsigned short i=0; i<=nPoly; ++i) {
    for(unsigned short j=0; j<=(nPoly-i); ++j, ++ii) {
      for(unsigned k=0; k<r.size(); ++k) {

        /*--- Determine the coefficients a and b. ---*/
        passivedouble a;
        if(fabs(s[k]-1.0) < 1.e-8) a = -1.0;
        else a = 2.0*(1.0+r[k])/(1.0-s[k]) - 1.0;

        const passivedouble b = s[k];

        /*--- Determine the value of the two 1D contributions to the 2D
              basis functions as well as the gradients of these basis
              functions w.r.t. to their arguments. ---*/
        const passivedouble fa  = NormJacobi(i,0,    0,a);
        const passivedouble gb  = NormJacobi(j,2*i+1,0,b);
        const passivedouble dfa = GradNormJacobi(i,0,    0,a);
        const passivedouble dgb = GradNormJacobi(j,2*i+1,0,b);

        /*--- Determine the gradients of the basis functions w.r.t. the
              coordinates r and s. The product rule must be used in order
              to change the derivative of a to the derivative of r and s. ---*/
        VDr(k,ii) = sqrt(2.0)*dfa*gb;
        VDs(k,ii) = VDr(k,ii);
        if(i > 0)
        {
          passivedouble tmp = 1.0;
          if( i-1 ) tmp = pow((1.0-b), (i-1));

          VDr(k,ii) = 2.0*tmp*VDr(k,ii);
          VDs(k,ii) = (a+1.0)*tmp*VDs(k,ii) - i*tmp*sqrt(2.0)*fa*gb;
        }

        passivedouble tmp = 1.0;
        if( i ) tmp = pow((1.0-b), i);
        VDs(k,ii) += sqrt(2.0)*fa*dgb*tmp;
      }
    }
  }
}

void CFEMStandardTriGrid::VandermondeTriangle(const vector<passivedouble>   &r,
                                              const vector<passivedouble>   &s,
                                              ColMajorMatrix<passivedouble> &V) {

  /*--- For a triangle the orthogonal basis for the reference element is obtained
        by a combination of a Jacobi polynomial and a Legendre polynomial. This
        is the result of the orthonormalization of the monomial basis. ---*/
  unsigned short ii = 0;
  for(unsigned short i=0; i<=nPoly; ++i) {
    for(unsigned short j=0; j<=(nPoly-i); ++j, ++ii) {
      for(unsigned short k=0; k<r.size(); ++k) {

        /*--- Determine the coefficients a and b. ---*/
        passivedouble a;
        if(fabs(s[k]-1.0) < 1.e-8) a = -1.0;
        else a = 2.0*(1.0+r[k])/(1.0-s[k]) - 1.0;

        const passivedouble b = s[k];

        /*--- Determine the value of the current basis function in this point. ---*/
        passivedouble tmp = 1.0;
        if( i ) tmp = pow((1.0-b),i);
        V(k,ii) = sqrt(2.0)*tmp*NormJacobi(i,0,0,a)*NormJacobi(j,2*i+1,0,b);
      }
    }
  }
}

void CFEMStandardTriGrid::LocalGridConnFaces(void) {

  /*--- Allocate the first index of gridConnFaces, which is equal to the number
        of faces of the triangle, which is 3. Reserve memory for the second
        index afterwards. ---*/
  gridConnFaces.resize(3);

  gridConnFaces[0].reserve(nPoly+1);
  gridConnFaces[1].reserve(nPoly+1);
  gridConnFaces[2].reserve(nPoly+1);

  /*--- For a triangular element the faces are lines. Loop over the nodes
        of the lines to set the connectivity. Make sure that the element
        is to the left of the faces. ---*/
  for(signed short i=0; i<=nPoly; ++i) gridConnFaces[0].push_back(i);
  for(signed short i=0; i<=nPoly; ++i) gridConnFaces[1].push_back((i+1)*(nPoly+1) - i*(i+1)/2 -1);
  for(signed short i=nPoly; i>=0; --i) gridConnFaces[2].push_back(i*(nPoly+1) - i*(i-1)/2);
}
