/*!
 * \file CFEMStandardPrismBase.cpp
 * \brief Functions for the class CFEMStandardPrismBase.
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

#include "../../include/fem/CFEMStandardPrismBase.hpp"
#include "../../include/fem/fem_gauss_jacobi_quadrature.hpp"
#include "../../include/toolboxes/CSquareMatrixCM.hpp"

/*----------------------------------------------------------------------------------*/
/*               Public member functions of CFEMStandardPrismBase.                  */
/*----------------------------------------------------------------------------------*/

CFEMStandardPrismBase::CFEMStandardPrismBase(const unsigned short val_nPoly,
                                             const unsigned short val_orderExact)
  : CFEMStandardQuadBase(),
    CFEMStandardTriBase() {

  /*--- Store the command line arguments. ---*/
  VTK_Type   = PRISM;
  nPoly      = val_nPoly;
  orderExact = val_orderExact;

  /*--- Determine the number of DOFs in 1D, the number of DOFs of the base triangle
        and the total number of DOFs.  Also determine the padded value of the 
        total number of DOFs. ---*/
  nDOFs1D       = nPoly + 1;
  nDOFsTriangle = nDOFs1D*(nDOFs1D+1)/2;
  nDOFs         = nDOFs1D*nDOFsTriangle;
  nDOFsPad      = PaddedValue(nDOFs);

  /*--- Determine the 1D integration points of a line, which corresponds to the
        direction normal to the base triangle of the prism. ---*/
  nInt1D = orderExact/2 + 1;
  rLineInt.resize(nInt1D);
  wLineInt.resize(nInt1D);

  CGaussJacobiQuadrature GaussJacobi;
  GaussJacobi.GetQuadraturePoints(0.0, 0.0, -1.0, 1.0, rLineInt, wLineInt);

  /*--- Determine the parametric location and weights of the
        integration rule of the base triangle. ---*/
  IntegrationPointsTriangle();
  nIntTriangle = rTriangleInt.size();

  /*--- The 3D quadrature rule is a tensor product of the 1D Gauss-Legendre
        quadrature rule and the integration rule of the triangle. Determine
        the total number of integration points and its padded value. ---*/
  nIntegration    = nInt1D*nIntTriangle;
  nIntegrationPad = PaddedValue(nIntegration);

  /*--- Allocate the memory for the padded number of integration points and
        initialize the weights to zero. This is done such that the padded
        values are initialized appropriately. ---*/
  wIntegration.resize(nIntegrationPad);
  wIntegration.setConstant(0.0);

  /*--- Determine the integration weights of the prism. ---*/
  unsigned short ii = 0;
  for(unsigned short j=0; j<nInt1D; ++j)
    for(unsigned short i=0; i<nIntTriangle; ++i, ++ii)
        wIntegration(ii) = wTriangleInt[i]*wLineInt[j];
}

/*----------------------------------------------------------------------------------*/
/*             Protected member functions of CFEMStandardPrismBase.                 */
/*----------------------------------------------------------------------------------*/

void CFEMStandardPrismBase::ConvertCoor2DQuadFaceTo3DPrism(const vector<passivedouble> &rLine,
                                                           const unsigned short        faceID_Elem,
                                                           const unsigned short        orientation,
                                                           vector<passivedouble>       &rPrism,
                                                           vector<passivedouble>       &sPrism,
                                                           vector<passivedouble>       &tPrism) {

  /*--- Determine the number of points on the quadrilateral face, which is a tensor
        product of rLine in two dimensions. Afterwards, allocate the memory for
        rPrism, sPrism and tPrism. ---*/
  const unsigned short nPoints1D = rLine.size();
  const unsigned short nPoints   = nPoints1D*nPoints1D;

  rPrism.resize(nPoints);
  sPrism.resize(nPoints);
  tPrism.resize(nPoints);

  /*--- Create the parametric coordinates in the frame of the face via a tensor product. ---*/
  vector<passivedouble> rF(nPoints), sF(nPoints);

  unsigned short k = 0;
  for(unsigned short j=0; j<nPoints1D; ++j) {
    for(unsigned short i=0; i<nPoints1D; ++i, ++k) {
      rF[k] = rLine[i];
      sF[k] = rLine[j];
    }
  }

  /*--- Make a distinction between face ID of the prism on which the current face resides. ---*/
  switch( faceID_Elem ) {

    case 2: {

      /*--- Face 2 of the prism, which corresponds to s == -1. Set this value. ---*/
      for(k=0; k<nPoints; ++k) sPrism[k] = -1.0;

      /*--- The r- and t-coordinates depend on the orientation of the face. ---*/
      switch( orientation ) {
        case 0: for(k=0; k<nPoints; ++k){rPrism[k] =  sF[k]; tPrism[k] =  rF[k];} break;
        case 1: for(k=0; k<nPoints; ++k){rPrism[k] =  rF[k]; tPrism[k] =  sF[k];} break;
        case 2: for(k=0; k<nPoints; ++k){rPrism[k] =  sF[k]; tPrism[k] = -rF[k];} break;
        case 3: for(k=0; k<nPoints; ++k){rPrism[k] = -rF[k]; tPrism[k] = -sF[k];} break;
        case 4: for(k=0; k<nPoints; ++k){rPrism[k] = -sF[k]; tPrism[k] =  rF[k];} break;
        default:
          SU2_MPI::Error(string("Invalid orientation for face 2. This should not happen."),
                         CURRENT_FUNCTION);
      }
      break;
    }

    case 3: {

      /*--- Face 3 of the prism, which corresponds to r == -1. Set this value. ---*/
      for(k=0; k<nPoints; ++k) rPrism[k] = -1.0;

      /*--- The s- and t-coordinates depend on the orientation of the face. ---*/
      switch( orientation ) {
        case 0: for(k=0; k<nPoints; ++k){sPrism[k] =  rF[k]; tPrism[k] =  sF[k];} break;
        case 1: for(k=0; k<nPoints; ++k){sPrism[k] =  sF[k]; tPrism[k] =  rF[k];} break;
        case 2: for(k=0; k<nPoints; ++k){sPrism[k] = -rF[k]; tPrism[k] =  sF[k];} break;
        case 3: for(k=0; k<nPoints; ++k){sPrism[k] = -sF[k]; tPrism[k] = -rF[k];} break;
        case 4: for(k=0; k<nPoints; ++k){sPrism[k] =  rF[k]; tPrism[k] = -sF[k];} break;
        default:
          SU2_MPI::Error(string("Invalid orientation for face 3. This should not happen."),
                         CURRENT_FUNCTION);
      }
      break;
    }

    case 4: {

      /*--- Face 4 of the prism, which corresponds to r+s == 0.
            Set the parametric coordinates s and t, depending on the
            orientation of the face. ---*/
      switch( orientation ) {
        case 0: for(k=0; k<nPoints; ++k){sPrism[k] =  sF[k]; tPrism[k] =  rF[k];} break;
        case 1: for(k=0; k<nPoints; ++k){sPrism[k] =  rF[k]; tPrism[k] =  sF[k];} break;
        case 2: for(k=0; k<nPoints; ++k){sPrism[k] =  sF[k]; tPrism[k] = -rF[k];} break;
        case 3: for(k=0; k<nPoints; ++k){sPrism[k] = -rF[k]; tPrism[k] = -sF[k];} break;
        case 4: for(k=0; k<nPoints; ++k){sPrism[k] = -sF[k]; tPrism[k] =  rF[k];} break;
        default:
          SU2_MPI::Error(string("Invalid orientation for face 4. This should not happen."),
                         CURRENT_FUNCTION);
      }

      /*--- Set the value of rPrism = -sPrism. ---*/
      for(k=0; k<nPoints; ++k) rPrism[k] = -sPrism[k];
      break;
    }

    default:
      SU2_MPI::Error(string("Invalid faceID. This should not happen."), CURRENT_FUNCTION);
  }
}

void CFEMStandardPrismBase::ConvertCoor2DTriFaceTo3DPrism(const vector<passivedouble> &rF,
                                                          const vector<passivedouble> &sF,
                                                          const unsigned short        faceID_Elem,
                                                          const unsigned short        orientation,
                                                          vector<passivedouble>       &rTrianglePrism,
                                                          vector<passivedouble>       &sTrianglePrism,
                                                          vector<passivedouble>       &rLinePrism) {

  /*--- Determine the value of rLinePrism, depending on the face of the prism. ---*/
  rLinePrism.resize(1);

  switch( faceID_Elem ) {
    case 0: rLinePrism[0] = -1.0; break;
    case 1: rLinePrism[0] =  1.0; break;
    default:
      SU2_MPI::Error(string("Invalid faceID. This should not happen."), CURRENT_FUNCTION);
  }

  /*--- Allocate the memory for rTrianglePrism and sTrianglePrism. ---*/
  const unsigned short nPoints = rF.size();
  rTrianglePrism.resize(nPoints);
  sTrianglePrism.resize(nPoints);

  /*--- Abbreviate rTrianglePrism and sTrianglePrism to avoid long lines below. ---*/
  vector<passivedouble> &r = rTrianglePrism, &s = sTrianglePrism;

  /*--- The values of rTrianglePrism and sTrianglePrism depend on both the face ID in the
        numbering of the prism as well as the orientation of the face w.r.t. the prism.
        Make this distinction and set the values accordingly. ---*/
  unsigned short i;
  if(faceID_Elem == 0) {
    switch( orientation ) {
      case 0: for(i=0; i<nPoints; ++i) {r[i]= rF[i];           s[i]= sF[i];} break;
      case 1: for(i=0; i<nPoints; ++i) {r[i]= sF[i];           s[i]= rF[i];} break;
      case 2: for(i=0; i<nPoints; ++i) {r[i]=-1.0-rF[i]-sF[i]; s[i]= sF[i];} break;
      case 3: for(i=0; i<nPoints; ++i) {r[i]= rF[i];           s[i]=-1.0-rF[i]-sF[i];} break;
      default:
        SU2_MPI::Error(string("Invalid orientation for face 0. This should not happen."),
                       CURRENT_FUNCTION);
    }
  }
  else if(faceID_Elem == 1) {
    switch( orientation ) {
      case 0: for(i=0; i<nPoints; ++i) {r[i]= sF[i];           s[i]= rF[i];} break;
      case 1: for(i=0; i<nPoints; ++i) {r[i]= rF[i];           s[i]= sF[i];} break;
      case 2: for(i=0; i<nPoints; ++i) {r[i]= sF[i];           s[i]=-1.0-rF[i]-sF[i];} break;
      case 3: for(i=0; i<nPoints; ++i) {r[i]=-1.0-rF[i]-sF[i]; s[i]= rF[i];} break;
      default:
        SU2_MPI::Error(string("Invalid orientation for face 1. This should not happen."),
                       CURRENT_FUNCTION);
    }
  }
}

void CFEMStandardPrismBase::DerLagBasisIntPointsPrism(const unsigned short                   mPoly,
                                                      const vector<passivedouble>            &rTriangleDOFs,
                                                      const vector<passivedouble>            &sTriangleDOFs,
                                                      const vector<passivedouble>            &rLineDOFs,
                                                      const vector<passivedouble>            &rTriangleInts,
                                                      const vector<passivedouble>            &sTrianglInts,
                                                      const vector<passivedouble>            &rLineInts,
                                                      vector<ColMajorMatrix<passivedouble> > &derLag) {

  /*--- Determine the parametric coordinates of all DOFs of the prism. ---*/
  vector<passivedouble> rDOFs, sDOFs, tDOFs;
  LocationAllPointsPrism(rTriangleDOFs, sTriangleDOFs, rLineDOFs, rDOFs, sDOFs, tDOFs);

  /*--- Determine the parametric coordinates of all integration points of the prism. ---*/
  vector<passivedouble> rInt, sInt, tInt;
  LocationAllPointsPrism(rTriangleInts, sTrianglInts, rLineInts, rInt, sInt, tInt);

  /*--- Determine the padded number of the total number of integration points. ---*/
  const unsigned short nIntTot    = rInt.size();
  const unsigned short nIntTotPad = PaddedValue(nIntTot);

  /*--- Determine the inverse of the Vandermonde matrix of the DOFs. ---*/
  CSquareMatrixCM VInv(rDOFs.size());
  VandermondePrism(mPoly, rDOFs, sDOFs, tDOFs, VInv.GetMat());
  VInv.Invert();

  /*--- Determine the gradient of the Vandermonde matrix of the integration points. Make
        sure to allocate the number of rows to nIntTotPad and initialize them to zero. ---*/
  ColMajorMatrix<passivedouble> VDr(nIntTotPad,rDOFs.size()),
                                VDs(nIntTotPad,rDOFs.size()),
                                VDt(nIntTotPad,rDOFs.size());
  VDr.setConstant(0.0);
  VDs.setConstant(0.0);
  VDt.setConstant(0.0);

  GradVandermondePrism(mPoly, rInt, sInt, tInt, VDr, VDs, VDt);

  /*--- The gradients of the Lagrangian basis functions can be obtained by
        multiplying VDr, VDs, VDt and VInv. ---*/
  derLag.resize(3);
  VInv.MatMatMult('R', VDr, derLag[0]);
  VInv.MatMatMult('R', VDs, derLag[1]);
  VInv.MatMatMult('R', VDt, derLag[2]);

  /*--- Check if the sum of the elements of the relevant rows of derLag is 0. ---*/
  CheckRowSum(nIntTot, rDOFs.size(), 0.0, derLag[0]);
  CheckRowSum(nIntTot, rDOFs.size(), 0.0, derLag[1]);
  CheckRowSum(nIntTot, rDOFs.size(), 0.0, derLag[2]);
}

void CFEMStandardPrismBase::HesLagBasisIntPointsPrism(const unsigned short                   mPoly,
                                                      const vector<passivedouble>            &rTriangleDOFs,
                                                      const vector<passivedouble>            &sTriangleDOFs,
                                                      const vector<passivedouble>            &rLineDOFs,
                                                      const vector<passivedouble>            &rTriangleInts,
                                                      const vector<passivedouble>            &sTrianglInts,
                                                      const vector<passivedouble>            &rLineInts,
                                                      vector<ColMajorMatrix<passivedouble> > &hesLag) {

  /*--- Determine the parametric coordinates of all DOFs of the prism. ---*/
  vector<passivedouble> rDOFs, sDOFs, tDOFs;
  LocationAllPointsPrism(rTriangleDOFs, sTriangleDOFs, rLineDOFs, rDOFs, sDOFs, tDOFs);

  /*--- Determine the parametric coordinates of all integration points of the prism. ---*/
  vector<passivedouble> rInt, sInt, tInt;
  LocationAllPointsPrism(rTriangleInts, sTrianglInts, rLineInts, rInt, sInt, tInt);

  /*--- Determine the padded number of the total number of integration points. ---*/
  const unsigned short nIntTot    = rInt.size();
  const unsigned short nIntTotPad = PaddedValue(nIntTot);

  /*--- Determine the inverse of the Vandermonde matrix of the DOFs. ---*/
  CSquareMatrixCM VInv(rDOFs.size());
  VandermondePrism(mPoly, rDOFs, sDOFs, tDOFs, VInv.GetMat());
  VInv.Invert();

  /*--- Determine the Hessian of the Vandermonde matrix of the integration points. Make
        sure to allocate the number of rows to nIntTotPad and initialize them to zero. ---*/
  ColMajorMatrix<passivedouble> VDr2(nIntTotPad,rDOFs.size()),
                                VDs2(nIntTotPad,rDOFs.size()),
                                VDt2(nIntTotPad,rDOFs.size()),
                                VDrs(nIntTotPad,rDOFs.size()),
                                VDrt(nIntTotPad,rDOFs.size()),
                                VDst(nIntTotPad,rDOFs.size());
  VDr2.setConstant(0.0);
  VDs2.setConstant(0.0);
  VDt2.setConstant(0.0);
  VDrs.setConstant(0.0);
  VDrt.setConstant(0.0);
  VDst.setConstant(0.0);

  HesVandermondePrism(mPoly, rInt, sInt, tInt, VDr2, VDs2, VDt2, VDrs, VDrt, VDst);

  /*--- The Hessian of the Lagrangian basis functions can be obtained by
        multiplying VDr2, VDs2, etc and VInv. ---*/
  hesLag.resize(6);
  VInv.MatMatMult('R', VDr2, hesLag[0]);
  VInv.MatMatMult('R', VDs2, hesLag[1]);
  VInv.MatMatMult('R', VDt2, hesLag[2]);
  VInv.MatMatMult('R', VDrs, hesLag[3]);
  VInv.MatMatMult('R', VDrt, hesLag[4]);
  VInv.MatMatMult('R', VDst, hesLag[5]);

  /*--- Check if the sum of the elements of the relevant rows of hesLagBasisInt is 0. ---*/
  CheckRowSum(nIntTot, rDOFs.size(), 0.0, hesLag[0]);
  CheckRowSum(nIntTot, rDOFs.size(), 0.0, hesLag[1]);
  CheckRowSum(nIntTot, rDOFs.size(), 0.0, hesLag[2]);
  CheckRowSum(nIntTot, rDOFs.size(), 0.0, hesLag[3]);
  CheckRowSum(nIntTot, rDOFs.size(), 0.0, hesLag[4]);
  CheckRowSum(nIntTot, rDOFs.size(), 0.0, hesLag[5]);
}

void CFEMStandardPrismBase::LagBasisIntPointsPrism(const unsigned short          mPoly,
                                                   const vector<passivedouble>   &rTriangleDOFs,
                                                   const vector<passivedouble>   &sTriangleDOFs,
                                                   const vector<passivedouble>   &rLineDOFs,
                                                   const vector<passivedouble>   &rTriangleInts,
                                                   const vector<passivedouble>   &sTrianglInts,
                                                   const vector<passivedouble>   &rLineInts,
                                                   ColMajorMatrix<passivedouble> &lag) {

  /*--- Determine the parametric coordinates of all DOFs of the prism. ---*/
  vector<passivedouble> rDOFs, sDOFs, tDOFs;
  LocationAllPointsPrism(rTriangleDOFs, sTriangleDOFs, rLineDOFs, rDOFs, sDOFs, tDOFs);

  /*--- Determine the parametric coordinates of all integration points of the prism. ---*/
  vector<passivedouble> rInt, sInt, tInt;
  LocationAllPointsPrism(rTriangleInts, sTrianglInts, rLineInts, rInt, sInt, tInt);

  /*--- Determine the padded number of the total number of integration points. ---*/
  const unsigned short nIntTot    = rInt.size();
  const unsigned short nIntTotPad = PaddedValue(nIntTot);

  /*--- Determine the inverse of the Vandermonde matrix of the DOFs. ---*/
  CSquareMatrixCM VInv(rDOFs.size());
  VandermondePrism(mPoly, rDOFs, sDOFs, tDOFs, VInv.GetMat());
  VInv.Invert();

  /*--- Determine the Vandermonde matrix of the integration points. Make sure to
        allocate the number of rows to nIntTotPad and initialize them to zero. ---*/
  ColMajorMatrix<passivedouble> V(nIntTotPad,rDOFs.size());
  V.setConstant(0.0);
  VandermondePrism(mPoly, rInt, sInt, tInt, V);

  /*--- The Lagrangian basis functions can be obtained by multiplying
        V and VInv. ---*/
  VInv.MatMatMult('R', V, lag);

  /*--- Check if the sum of the elements of the relevant rows of lag is 1. ---*/
  CheckRowSum(nIntTot, rDOFs.size(), 1.0, lag);
}

void CFEMStandardPrismBase::GradVandermondePrism(const unsigned short          mPoly,
                                                 const vector<passivedouble>   &r,
                                                 const vector<passivedouble>   &s,
                                                 const vector<passivedouble>   &t,
                                                 ColMajorMatrix<passivedouble> &VDr,
                                                 ColMajorMatrix<passivedouble> &VDs,
                                                 ColMajorMatrix<passivedouble> &VDt) {

  /*--- Abbreviate sqrt(2), which is the scaling factor for the orthonormal basis
        functions of a triangle/prism. ---*/
  const passivedouble sqrt2 = sqrt(2.0);

  /*--- For a prism the orthogonal basis for the reference element is a tensor
        product of the 1D basis functions in the structured direction of the prism
        and the basis functions of a triangle. For that triangle the orthogonal
        basis is obtained by a combination of a Jacobi polynomial and a Legendre
        polynomial. This is the result of the orthonormalization of the
        monomial basis. Note that the sequence of the i, j and k loop must be
        identical to the evaluation of the Vandermonde matrix itself. ---*/
  unsigned short ii = 0;
  for(unsigned short i=0; i<=mPoly; ++i) {
    for(unsigned short j=0; j<=(mPoly-i); ++j) {
      for(unsigned short k=0; k<=mPoly; ++k, ++ii) {
        for(unsigned short l=0; l<r.size(); ++l) {

          /*--- Determine the coefficients a and b. ---*/
          passivedouble a;
          if(fabs(s[l]-1.0) < 1.e-8) a = -1.0;
          else a = 2.0*(1.0+r[l])/(1.0-s[l]) - 1.0;

          const passivedouble b = s[l];

          /*--- Determine the value of the two 1D contributions to the 2D
                basis functions of the triangle as well as the gradients of
                these basis functions w.r.t. their arguments. ---*/
          const passivedouble fa  = NormJacobi(i,0,    0,a);
          const passivedouble gb  = NormJacobi(j,2*i+1,0,b);
          const passivedouble dfa = GradNormJacobi(i,0,    0,a);
          const passivedouble dgb = GradNormJacobi(j,2*i+1,0,b);

          /*--- Determine the value of the 1D basis function in the structured
                direction of the prism as well as the gradient of this basis
                function w.r.t. its argument. ---*/
          const passivedouble ht  = NormJacobi(k,0,0,t[l]);
          const passivedouble dht = GradNormJacobi(k,0,0,t[l]);

          /*--- Computation of the powers of (1-b) that occur in the expressions for
                the gradients. Note the safeguard to avoid division by zero. This is
                allowed, because the implementation is such that when this clipping
                is active, it is multiplied by zero. ---*/
          const passivedouble tmpbi   = i > 0 ? pow((1.0-b), i)   : 1.0;
          const passivedouble tmpbim1 = i > 1 ? pow((1.0-b), i-1) : 1.0;

          /*--- Compute the derivatives of the basis function. ---*/
          VDr(l,ii) = sqrt2*tmpbim1* 2.0*dfa*gb*ht;
          VDs(l,ii) = sqrt2*tmpbim1*((a+1)*dfa*gb - i*fa*gb)*ht
                    + sqrt2*tmpbi  * fa*dgb*ht;
          VDt(l,ii) = sqrt2*tmpbi  * fa*gb*dht;
        }
      }
    }
  }
}

void CFEMStandardPrismBase::HesVandermondePrism(const unsigned short          mPoly,
                                                const vector<passivedouble>   &r,
                                                const vector<passivedouble>   &s,
                                                const vector<passivedouble>   &t,
                                                ColMajorMatrix<passivedouble> &VDr2,
                                                ColMajorMatrix<passivedouble> &VDs2,
                                                ColMajorMatrix<passivedouble> &VDt2,
                                                ColMajorMatrix<passivedouble> &VDrs,
                                                ColMajorMatrix<passivedouble> &VDrt,
                                                ColMajorMatrix<passivedouble> &VDst) {

  /*--- Abbreviate sqrt(2), which is the scaling factor for the orthonormal basis
        functions of a triangle/prism. ---*/
  const passivedouble sqrt2 = sqrt(2.0);

  /*--- For a prism the orthogonal basis for the reference element is a tensor
        product of the 1D basis functions in the structured direction of the prism
        and the basis functions of a triangle. For that triangle the orthogonal
        basis is obtained by a combination of a Jacobi polynomial and a Legendre
        polynomial. This is the result of the orthonormalization of the
        monomial basis. Note that the sequence of the i, j and k loop must be
        identical to the evaluation of the Vandermonde matrix itself. ---*/
  unsigned short ii = 0;
  for(unsigned short i=0; i<=mPoly; ++i) {
    for(unsigned short j=0; j<=(mPoly-i); ++j) {
      for(unsigned short k=0; k<=mPoly; ++k, ++ii) {
        for(unsigned short l=0; l<r.size(); ++l) {

          /*--- Determine the coefficients a and b. ---*/
          passivedouble a;
          if(fabs(s[l]-1.0) < 1.e-8) a = -1.0;
          else a = 2.0*(1.0+r[l])/(1.0-s[l]) - 1.0;

          const passivedouble b = s[l];

          /*--- Determine the value of the two 1D contributions to the 2D
                basis functions of the triangle as well as the 1st and 2nd
                derivatives of these contributions w.r.t. to their arguments. ---*/
          const passivedouble fa   = NormJacobi(i,0,    0,a);
          const passivedouble gb   = NormJacobi(j,2*i+1,0,b);
          const passivedouble dfa  = GradNormJacobi(i,0,    0,a);
          const passivedouble dgb  = GradNormJacobi(j,2*i+1,0,b);
          const passivedouble d2fa = HesNormJacobi(i,0,    0,a);
          const passivedouble d2gb = HesNormJacobi(j,2*i+1,0,b);

          /*--- Determine the value of the 1D basis function in the structured
                direction of the prism as well as the 1st and 2nd derivative
                of this basis function w.r.t. its argument. ---*/
          const passivedouble ht   = NormJacobi(k,0,0,t[l]);
          const passivedouble dht  = GradNormJacobi(k,0,0,t[l]);
          const passivedouble d2ht = HesNormJacobi(k,0,0,t[l]);

          /*--- Computation of the powers of (1-b) that occur in the expressions for
                the Hessian. Note the safeguard to avoid division by zero. This is
                allowed, because the implementation is such that when this clipping
                is active, it is multiplied by zero. ---*/
          const passivedouble tmpbi   = i > 0 ? pow((1.0-b), i)   : 1.0;
          const passivedouble tmpbim1 = i > 1 ? pow((1.0-b), i-1) : 1.0;
          const passivedouble tmpbim2 = i > 2 ? pow((1.0-b), i-2) : 1.0;

          /*--- Compute the 2nd derivative w.r.t. to r. ---*/
          VDr2(l,ii) = tmpbim2*4.0*d2fa*gb*ht;

          /*--- Compute the 2nd derivative w.r.t. s. ---*/
          VDs2(l,ii) = tmpbim2*((a+1.0)*(a+1.0)*d2fa*gb
                     +          i*(i-1)*gb*fa
                     -          2.0*(i-1)*(a+1.0)*gb*dfa)*ht
                     + tmpbim1* 2.0*dgb*(dfa*(a+1.0) - i*fa)*ht
                     + tmpbi  * fa*d2gb*ht;

          /*--- Compute the 2nd derivative w.r.t. t. ---*/
          VDt2(l,ii) = tmpbi*fa*gb*d2ht;

          /*--- Compute the cross derivative w.r.t. r and s. ---*/
          VDrs(l,ii) = tmpbim2*2.0*gb*((a+1.0)*d2fa - (i-1)*dfa)*ht
                     + tmpbim1*2.0*dfa*dgb*ht;

          /*--- Compute the cross derivative w.r.t. r and t. ---*/
          VDrt(l,ii) = tmpbim1*2.0*dfa*gb*dht;

          /*--- Compute the cross derivative w.r.t. s and t. ---*/
          VDst(l,ii) = tmpbim1*((a+1)*dfa*gb - i*fa*gb)*dht
                     + tmpbi  * fa*dgb*dht;

          /*--- Multiply all 2nd derivatives with sqrt(2) to obtain the
                correct expressions. ---*/
          VDr2(l,ii) *= sqrt2;
          VDs2(l,ii) *= sqrt2;
          VDt2(l,ii) *= sqrt2;
          VDrs(l,ii) *= sqrt2;
          VDrt(l,ii) *= sqrt2;
          VDst(l,ii) *= sqrt2;
        }
      }
    }
  }
}

void CFEMStandardPrismBase::VandermondePrism(const unsigned short          mPoly,
                                             const vector<passivedouble>   &r,
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
  for(unsigned short i=0; i<=mPoly; ++i) {
    for(unsigned short j=0; j<=(mPoly-i); ++j) {
      for(unsigned short k=0; k<=mPoly; ++k, ++ii) {
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

void CFEMStandardPrismBase::SubConnLinearElements(void) {

  /*--- The prism is split into several linear prisms.
        Set the VTK sub-types accordingly. ---*/
  VTK_SubType1 = PRISM;
  VTK_SubType2 = NONE;

  /*--- Determine the number of DOFs for a triangle. This is the offset in
        k-direction, the structured direction of a prisms. ---*/
  const unsigned short nDOFTria = (nPoly+1)*(nPoly+2)/2;

  /*--- Loop in k-direction, which is the structured direction of the prism. ---*/
  for(unsigned short k=0; k<nPoly; ++k) {

    /*--- Initialize the counter jj to the ID of the first vertex for this
          k-value. jj contains the ID of the first vertex on the "0-2 edge"
          of the current bottom triangle. ---*/
    unsigned short jj = k*nDOFTria;

    /*--- Loop over subedges of the left boundary of the standard triangle. ---*/
    for(unsigned short j=0; j<nPoly; ++j) {

      /*--- Check if the "down" elements must be written. ---*/
      if( j ) {

        /*--- Determine the offset of the relevant DOF on the previous row. ---*/
        const unsigned short kk = jj - (nPoly + 1 - j);

        /*--- Loop over the edges of this row. ---*/
        for(unsigned short i=0; i<(nPoly-j); ++i) {

          /*--- Determine the local connectivity of this subelement and add
                it to subConn1ForPlotting. ---*/
          const unsigned short n0 = jj + i;
          const unsigned short n1 = kk + i;
          const unsigned short n2 = n0 + 1;
          const unsigned short n3 = n0 + nDOFTria;
          const unsigned short n4 = n1 + nDOFTria;
          const unsigned short n5 = n2 + nDOFTria;

          subConn1ForPlotting.push_back(n0);
          subConn1ForPlotting.push_back(n1);
          subConn1ForPlotting.push_back(n2);
          subConn1ForPlotting.push_back(n3);
          subConn1ForPlotting.push_back(n4);
          subConn1ForPlotting.push_back(n5);
        }
      }

      /*--- The "upp" elements must always be written.
            Determine the offset of the DOF on the next row. ---*/
      const unsigned short kk = jj + (nPoly + 1 - j);

      /*--- Loop over the edges of this row. ---*/
      for(unsigned short i=0; i<(nPoly-j); ++i) {

        /*--- Determine the local connectivity of this subelement and add
              it to subConn1ForPlotting.           ---*/
        const unsigned short n0 = jj + i;
        const unsigned short n1 = n0 + 1;
        const unsigned short n2 = kk + i;
        const unsigned short n3 = n0 + nDOFTria;
        const unsigned short n4 = n1 + nDOFTria;
        const unsigned short n5 = n2 + nDOFTria;

        subConn1ForPlotting.push_back(n0);
        subConn1ForPlotting.push_back(n1);
        subConn1ForPlotting.push_back(n2);
        subConn1ForPlotting.push_back(n3);
        subConn1ForPlotting.push_back(n4);
        subConn1ForPlotting.push_back(n5);
      }

      /*--- Set jj to kk for the next edge. ---*/
      jj = kk;
    }
  }
}

void CFEMStandardPrismBase::LocationAllPointsPrism(const vector<passivedouble> &rTriangle,
                                                   const vector<passivedouble> &sTriangle,
                                                   const vector<passivedouble> &rLine,
                                                   vector<passivedouble>       &rPrism,
                                                   vector<passivedouble>       &sPrism,
                                                   vector<passivedouble>       &tPrism) {

  /*--- Determine the total number of points. ---*/
  const unsigned short nPointTriangle = rTriangle.size();
  const unsigned short nPointLine     = rLine.size();
  const unsigned short nPoint         = nPointTriangle*nPointLine;

  /*--- Allocate the memory for rPrism, sPrism and tPrism. ---*/
  rPrism.resize(nPoint);
  sPrism.resize(nPoint);
  tPrism.resize(nPoint);

  /*--- Determine the location of all the integration points. ---*/
  unsigned short ii = 0;
  for(unsigned short k=0; k<nPointLine; ++k) {
    for(unsigned short j=0; j<nPointTriangle; ++j, ++ii) {
      rPrism[ii] = rTriangle[j];
      sPrism[ii] = sTriangle[j];
      tPrism[ii] = rLine[k];
    }
  }
}

void CFEMStandardPrismBase::LocalGridConnFaces(void) {

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
