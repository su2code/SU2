/*!
 * \file CFEMStandardPrismGrid.cpp
 * \brief Functions for the class CFEMStandardPrismGrid.
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

#include "../../include/fem/CFEMStandardPrismGrid.hpp"
#include "../../include/toolboxes/CGeneralSquareMatrixCM.hpp"

/*----------------------------------------------------------------------------------*/
/*             Public member functions of CFEMStandardPrismGrid.                    */
/*----------------------------------------------------------------------------------*/

CFEMStandardPrismGrid::CFEMStandardPrismGrid(const unsigned short val_nPoly,
                                             const unsigned short val_orderExact)
  : CFEMStandardPrism(val_nPoly, val_orderExact) {

  /*--- Compute the values of the Lagrangian basis functions in the integration
        points for both the equidistant and LGL point distribution. ---*/
  LagBasisIntPointsPrism(rTriangleDOFsEqui, sTriangleDOFsEqui, rLineDOFsEqui, lagBasisIntEqui);
  LagBasisIntPointsPrism(rTriangleDOFsLGL,  sTriangleDOFsLGL,  rLineDOFsLGL,  lagBasisIntLGL);

  /*--- Compute the values of the derivatives of the Lagrangian basis functions in
        the integration points for both the equidistant and LGL point distribution. ---*/
  DerLagBasisIntPointsPrism(rTriangleDOFsEqui, sTriangleDOFsEqui, rLineDOFsEqui, derLagBasisIntEqui);
  DerLagBasisIntPointsPrism(rTriangleDOFsLGL,  sTriangleDOFsLGL,  rLineDOFsLGL,  derLagBasisIntLGL);

  /*--- Create the local grid connectivities of the faces of the volume element. ---*/
  LocalGridConnFaces();

  /*--- Set up the jitted gemm call, if supported. For this particular standard
        element the derivative of the coordinates are computed, which is 3. ---*/
  SetUpJittedGEMM(nIntegrationPad, 3, nDOFs);
}

void CFEMStandardPrismGrid::DerivativesCoorIntPoints(const bool                         LGLDistribution,
                                                     ColMajorMatrix<su2double>          &matCoor,
                                                     vector<ColMajorMatrix<su2double> > &matDerCoor) {

  /*--- Check for which point distribution the derivatives must be computed. ---*/
  if( LGLDistribution ) {

    /*--- LGL distribution. Call the function OwnGemm 3 times to compute the derivatives
          of the Cartesian coordinates w.r.t. the three parametric coordinates. ---*/
    OwnGemm(nIntegrationPad, 3, nDOFs, derLagBasisIntLGL[0], matCoor, matDerCoor[0], nullptr);
    OwnGemm(nIntegrationPad, 3, nDOFs, derLagBasisIntLGL[1], matCoor, matDerCoor[1], nullptr);
    OwnGemm(nIntegrationPad, 3, nDOFs, derLagBasisIntLGL[2], matCoor, matDerCoor[2], nullptr);
  }
  else {

    /*--- Equidistant distribution. Call the function OwnGemm 3 times to compute the derivatives
          of the Cartesian coordinates w.r.t. the three parametric coordinates. ---*/
    OwnGemm(nIntegrationPad, 3, nDOFs, derLagBasisIntEqui[0], matCoor, matDerCoor[0], nullptr);
    OwnGemm(nIntegrationPad, 3, nDOFs, derLagBasisIntEqui[1], matCoor, matDerCoor[1], nullptr);
    OwnGemm(nIntegrationPad, 3, nDOFs, derLagBasisIntEqui[2], matCoor, matDerCoor[2], nullptr);
  }
}

/*----------------------------------------------------------------------------------*/
/*             Private member functions of CFEMStandardPrismGrid.                   */
/*----------------------------------------------------------------------------------*/

void CFEMStandardPrismGrid::DerLagBasisIntPointsPrism(const vector<passivedouble>            &rTriangleDOFs,
                                                      const vector<passivedouble>            &sTriangleDOFs,
                                                      const vector<passivedouble>            &rLineDOFs,
                                                      vector<ColMajorMatrix<passivedouble> > &derLag) {

  /*--- Determine the parametric coordinates of all DOFs of the prism. ---*/
  vector<passivedouble> rDOFs, sDOFs, tDOFs;
  LocationAllDOFsPrism(rTriangleDOFs, sTriangleDOFs, rLineDOFs, rDOFs, sDOFs, tDOFs);

  /*--- Determine the parametric coordinates of all integration points of the prism. ---*/
  vector<passivedouble> rInt, sInt, tInt;
  LocationAllIntegrationPoints(rInt, sInt, tInt);

  /*--- Determine the padded number of the total number of integration points. ---*/
  const unsigned short nIntTot    = rInt.size();
  const unsigned short nIntTotPad = ((nIntTot+baseVectorLen-1)/baseVectorLen)*baseVectorLen;

  /*--- Determine the inverse of the Vandermonde matrix of the DOFs. ---*/
  CGeneralSquareMatrixCM VInv(rDOFs.size());
  VandermondePrism(rDOFs, sDOFs, tDOFs, VInv.GetMat());
  VInv.Invert();

  /*--- Determine the gradient of the Vandermonde matrix of the integration points. Make
        sure to allocate the number of rows to nIntTotPad and initialize them to zero. ---*/
  ColMajorMatrix<passivedouble> VDr(nIntTotPad,rDOFs.size()),
                                VDs(nIntTotPad,rDOFs.size()),
                                VDt(nIntTotPad,rDOFs.size());
  VDr.setConstant(0.0);
  VDs.setConstant(0.0);
  VDt.setConstant(0.0);

  GradVandermondePrism(rInt, sInt, tInt, VDr, VDs, VDt);

  /*--- The gradients of the Lagrangian basis functions can be obtained by
        multiplying VDr, VDs, VDt and VInv. ---*/
  derLag.resize(3);
  VInv.MatMatMult('R', VDr, derLag[0]);
  VInv.MatMatMult('R', VDs, derLag[1]);
  VInv.MatMatMult('R', VDt, derLag[2]);

  /*--- Check if the sum of the elements of the relevant rows of derLag is 0. ---*/
  for(unsigned short i=0; i<nIntTot; ++i) {
    passivedouble rowSumDr = 0.0, rowSumDs = 0.0, rowSumDt = 0.0;
    for(unsigned short j=0; j<rDOFs.size(); ++j) {
      rowSumDr += derLag[0](i,j);
      rowSumDs += derLag[1](i,j);
      rowSumDt += derLag[2](i,j);
    }

    assert(fabs(rowSumDr) < 1.e-6);
    assert(fabs(rowSumDs) < 1.e-6);
    assert(fabs(rowSumDt) < 1.e-6);
  }
}

void CFEMStandardPrismGrid::LagBasisIntPointsPrism(const vector<passivedouble>   &rTriangleDOFs,
                                                   const vector<passivedouble>   &sTriangleDOFs,
                                                   const vector<passivedouble>   &rLineDOFs,
                                                   ColMajorMatrix<passivedouble> &lag) {

  /*--- Determine the parametric coordinates of all DOFs of the prism. ---*/
  vector<passivedouble> rDOFs, sDOFs, tDOFs;
  LocationAllDOFsPrism(rTriangleDOFs, sTriangleDOFs, rLineDOFs, rDOFs, sDOFs, tDOFs);

  /*--- Determine the parametric coordinates of all integration points of the prism. ---*/
  vector<passivedouble> rInt, sInt, tInt;
  LocationAllIntegrationPoints(rInt, sInt, tInt);

  /*--- Determine the padded number of the total number of integration points. ---*/
  const unsigned short nIntTot    = rInt.size();
  const unsigned short nIntTotPad = ((nIntTot+baseVectorLen-1)/baseVectorLen)*baseVectorLen;

  /*--- Determine the inverse of the Vandermonde matrix of the DOFs. ---*/
  CGeneralSquareMatrixCM VInv(rDOFs.size());
  VandermondePrism(rDOFs, sDOFs, tDOFs, VInv.GetMat());
  VInv.Invert();

  /*--- Determine the Vandermonde matrix of the integration points. Make sure to
        allocate the number of rows to nIntTotPad and initialize them to zero. ---*/ 
  ColMajorMatrix<passivedouble> V(nIntTotPad,rDOFs.size());
  V.setConstant(0.0);
  VandermondePrism(rInt, sInt, tInt, V);

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

void CFEMStandardPrismGrid::LocationAllDOFsPrism(const vector<passivedouble>   &rTriangleDOFs,
                                                 const vector<passivedouble>   &sTriangleDOFs,
                                                 const vector<passivedouble>   &rLineDOFs,
                                                 vector<passivedouble>         &rDOFs,
                                                 vector<passivedouble>         &sDOFs,
                                                 vector<passivedouble>         &tDOFs) {

  /*--- Determine the number of DOFs of the triangle and line. ---*/
  const unsigned short nTri  = rTriangleDOFs.size();
  const unsigned short nLine = rLineDOFs.size();

  /*--- Determine the total number of DOFs and allocate the memory. ---*/
  const unsigned short nDOFsTot = nTri*nLine;
  rDOFs.resize(nDOFsTot);
  sDOFs.resize(nDOFsTot);
  tDOFs.resize(nDOFsTot);

  /*--- Determine the location of all the DOFs. ---*/
  unsigned short ii = 0;
  for(unsigned short k=0; k<nLine; ++k) {
    for(unsigned short j=0; j<nTri; ++j, ++ii) {
      rDOFs[ii] = rTriangleDOFs[j];
      sDOFs[ii] = sTriangleDOFs[j];
      tDOFs[ii] = rLineDOFs[k];
    }
  }
}

void CFEMStandardPrismGrid::GradVandermondePrism(const vector<passivedouble>   &r,
                                                 const vector<passivedouble>   &s,
                                                 const vector<passivedouble>   &t,
                                                 ColMajorMatrix<passivedouble> &VDr,
                                                 ColMajorMatrix<passivedouble> &VDs,
                                                 ColMajorMatrix<passivedouble> &VDt) {

  /*--- For a prism the orthogonal basis for the reference element is a tensor
        product of the 1D basis functions in the structured direction of the prism
        and the basis functions of a triangle. For that triangle the orthogonal
        basis is obtained by a combination of a Jacobi polynomial and a Legendre
        polynomial. This is the result of the orthonormalization of the
        monomial basis. Note that the sequence of the i, j and k loop must be
        identical to the evaluation of the Vandermonde matrix itself.  ---*/
  unsigned short ii = 0;
  for(unsigned short i=0; i<=nPoly; ++i) {
    for(unsigned short j=0; j<=(nPoly-i); ++j) {
      for(unsigned short k=0; k<=nPoly; ++k, ++ii) {
        for(unsigned short l=0; l<r.size(); ++l) {

          /*--- Determine the coefficients a and b. ---*/
          passivedouble a;
          if(fabs(s[l]-1.0) < 1.e-8) a = -1.0;
          else a = 2.0*(1.0+r[l])/(1.0-s[l]) - 1.0;

          const passivedouble b = s[l];

          /*--- Determine the value of the two 1D contributions to the 2D
                basis functions of the triangle as well as the gradients of
                these basis functions w.r.t. to its argument. ---*/
          const passivedouble fa  = NormJacobi(i,0,    0,a);
          const passivedouble gb  = NormJacobi(j,2*i+1,0,b);
          const passivedouble dfa = GradNormJacobi(i,0,    0,a);
          const passivedouble dgb = GradNormJacobi(j,2*i+1,0,b);

          /*--- Determine the gradients of the basis functions w.r.t. the
                coordinates r and s. The product rule must be used in order
                to change the derivative of a to the derivative of r and s. ---*/
          VDr(l,ii) = sqrt(2.0)*dfa*gb;
          VDs(l,ii) = VDr(l,ii);
          if(i > 0)
          {
            passivedouble tmp = 1.0;
            if(i > 1) tmp = pow((1.0-b), (i-1));

            VDr(l,ii) = 2.0*tmp*VDr(l,ii);
            VDs(l,ii) = (a+1.0)*tmp*VDs(l,ii) - i*tmp*sqrt(2.0)*fa*gb;
          }

          passivedouble tmp = 1.0;
          if(i > 0) tmp = pow((1.0-b), i);

          VDs(l,ii) += sqrt(2.0)*fa*dgb*tmp;

          /*--- Multiply VDr and VDs with the contribution from the structured
                direction of the prism. ---*/
          VDr(l,ii) *= NormJacobi(k,0,0,t[l]);
          VDs(l,ii) *= NormJacobi(k,0,0,t[l]);

          /*--- Compute the derivative of the basis function in the t-direction,
                which is the structured direction. ---*/
          VDt(l,ii) = sqrt(2.0)*tmp*fa*gb*GradNormJacobi(k,0,0,t[l]);
        }
      }
    }
  }
}

void CFEMStandardPrismGrid::VandermondePrism(const vector<passivedouble>   &r,
                                             const vector<passivedouble>   &s,
                                             const vector<passivedouble>   &t,
                                             ColMajorMatrix<passivedouble> &V) {

  /*--- For a prism the orthogonal basis for the reference element is a tensor
        product of the 1D basis functions in the structured direction of the prism
        and the basis functions of a triangle. For that triangle the orthogonal
        basis is obtained by a combination of a Jacobi polynomial and a Legendre
        polynomial. This is the result of the orthonormalization of the
        monomial basis. ---*/
  unsigned short ii = 0;
  for(unsigned short i=0; i<=nPoly; ++i) {
    for(unsigned short j=0; j<=(nPoly-i); ++j) {
      for(unsigned short k=0; k<=nPoly; ++k, ++ii) {
        for(unsigned short l=0; l<r.size(); ++l) {

          /*--- Determine the coefficients a and b. ---*/
          passivedouble a;
          if(fabs(s[l]-1.0) < 1.e-8) a = -1.0;
          else a = 2.0*(1.0+r[l])/(1.0-s[l]) - 1.0;

          const passivedouble b = s[l];

          /*--- Determine the value of the current basis function in this point. ---*/
          passivedouble tmp = 1.0;
          if( i ) tmp = pow((1.0-b),i);

          V(l,ii) = sqrt(2.0)*tmp*NormJacobi(i,0,0,a)*NormJacobi(j,2*i+1,0,b)
                  * NormJacobi(k,0,0,t[l]);
        }
      }
    }
  }
}

void CFEMStandardPrismGrid::LocalGridConnFaces(void) {

  /*--- Allocate the first index of gridConnFaces, which is equal to the number
        of faces of the prism, which is 5. Reserve memory for the second
        index afterwards. ---*/
  const unsigned short nDOFsQuad     = (nPoly+1)*(nPoly+1);
  const unsigned short nDOFsTriangle = (nPoly+1)*(nPoly+2)/2;
  gridConnFaces.resize(5);

  gridConnFaces[0].reserve(nDOFsTriangle);
  gridConnFaces[1].reserve(nDOFsTriangle);
  gridConnFaces[2].reserve(nDOFsQuad);
  gridConnFaces[3].reserve(nDOFsQuad);
  gridConnFaces[4].reserve(nDOFsQuad);

  /*--- Loop over all the nodes of the prism and pick the correct
        ones for the faces. ---*/
  unsigned int ii = 0;
  for(unsigned short k=0; k<=nPoly; ++k) {
    for(unsigned short j=0; j<=nPoly; ++j) {
      const unsigned short uppBoundI = nPoly - j;
      for(unsigned short i=0; i<=uppBoundI; ++i, ++ii) {
        if(k == 0)         gridConnFaces[0].push_back(ii);
        if(k == nPoly)     gridConnFaces[1].push_back(ii);
        if(j == 0)         gridConnFaces[2].push_back(ii);
        if(i == 0)         gridConnFaces[3].push_back(ii);
        if((i+j) == nPoly) gridConnFaces[4].push_back(ii);
      }
    }
  }

  /*--- Make sure that the element is to the left of the faces. ---*/
  const unsigned short n0 = 0;
  const unsigned short n1 = nPoly;
  const unsigned short n2 = nDOFsTriangle -1;
  const unsigned short n3 = n0 + nDOFsTriangle*nPoly;
  const unsigned short n4 = n1 + nDOFsTriangle*nPoly;
  const unsigned short n5 = n2 + nDOFsTriangle*nPoly;

  ChangeDirectionTriangleConn(gridConnFaces[0], n0, n1, n2);
  ChangeDirectionTriangleConn(gridConnFaces[1], n3, n5, n4);
  ChangeDirectionQuadConn(gridConnFaces[2], n0, n3, n4, n1);
  ChangeDirectionQuadConn(gridConnFaces[3], n0, n2, n5, n3);
  ChangeDirectionQuadConn(gridConnFaces[4], n1, n4, n5, n2);
}
