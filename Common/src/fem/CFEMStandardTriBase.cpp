/*!
 * \file CFEMStandardTriBase.cpp
 * \brief Functions for the class CFEMStandardTriBase.
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

#include "../../include/fem/CFEMStandardTriBase.hpp"
#include "../../include/toolboxes/CSquareMatrixCM.hpp"

/*----------------------------------------------------------------------------------*/
/*               Public member functions of CFEMStandardTriBase.                    */
/*----------------------------------------------------------------------------------*/

CFEMStandardTriBase::CFEMStandardTriBase(const unsigned short val_nPoly,
                                         const unsigned short val_orderExact)
  : CFEMStandardLineBase() {

  /*--- Store the command line arguments. ---*/
  VTK_Type   = TRIANGLE;
  nPoly      = val_nPoly;
  orderExact = val_orderExact;

  /*--- Determine the total number of DOFs and its padded version. ---*/
  nDOFs    = (nPoly+1)*(nPoly+2)/2;
  nDOFsPad = PaddedValue(nDOFs);

  /*--- Determine the parametric location and weights of the
        integration rule of the triangle. ---*/
  IntegrationPointsTriangle();

  /*--- Determine the total number of integration points
        and its padded version. ---*/
  nIntegration    = rTriangleInt.size();
  nIntegrationPad = PaddedValue(nIntegration);

  /*--- Allocate the memory for the padded number of integration points and
        initialize the weights to zero. This is done such that the padded
        values are initialized appropriately. ---*/
  wIntegration.resize(nIntegrationPad);
  wIntegration.setConstant(0.0);

  /*--- Copy the values from wTriangleInt to wIntegration. ---*/
  for(unsigned short i=0; i<nIntegration; ++i)
    wIntegration(i) = wTriangleInt[i];
}

/*----------------------------------------------------------------------------------*/
/*                Protected member functions of CFEMStandardTriBase.                */
/*----------------------------------------------------------------------------------*/

void CFEMStandardTriBase::ChangeDirectionTriangleConn(vector<unsigned short> &connTriangle,
                                                      const unsigned short   vert0,
                                                      const unsigned short   vert1,
                                                      const unsigned short   vert2) {

  /*--- Determine the indices of the 3 corner vertices of the triangle. ---*/
  const unsigned short ind0 = 0;
  const unsigned short ind1 = nPoly;
  const unsigned short ind2 = (nPoly+1)*(nPoly+2)/2 -1;

  /*--- There exists a linear mapping from the indices of the numbering used in the
        connectivity of this face to the indices of the target numbering. This
        mapping is of the form ii = a + b*i + c*j and jj = d + e*i + f*j, where
        ii,jj are the indices of the target numbering and i,j the indices of the
        numbering used for this face. The values of the coefficients a,b,c,d,e,f
        depend on how the corner points coincide with each other. This is
        determined below. The bool verticesDontMatch is there to check if vertices
        do not match. This should not happen, but it is checked for security. ---*/
  signed short a=0, b=0, c=0, d=0, e=0, f=0;
  bool verticesDontMatch = false;

  if(vert0 == connTriangle[ind0]) {

    /*--- Vert0 coincides with the first vertex of the face connectivity.
          Check the situation for the neighboring vertices.  ---*/
    if(vert1 == connTriangle[ind1]) {

      /*--- The vertex numbering is the same for both faces. ---*/
      if(vert2 != connTriangle[ind2]) verticesDontMatch = true;

      a = 0; b = 1; c = 0; d = 0; e = 0; f = 1;
    }
    else if(vert1 == connTriangle[ind2]) {

      /*--- The ii-index corresponds to the j-index and the
            jj-index corresponds to i-index.    ---*/
      if(vert2 != connTriangle[ind1]) verticesDontMatch = true;

      a = 0; b = 0; c = 1; d = 0; e = 1; f = 0;
    }
    else {
      verticesDontMatch = true;
    }
  }
  else if(vert0 == connTriangle[ind1]) {

    /*--- Vert0 coincides with the second vertex of the face connectivity.
          Check the situation for the neighboring vertices.    ---*/
    if(vert1 == connTriangle[ind0]) {

      /*--- The ii-index corresponds to a combination of the i and j index
            and the jj-index corresponds to the j-index.   ---*/
      if(vert2 != connTriangle[ind2]) verticesDontMatch = true;

      a = nPoly; b = -1; c = -1; d = 0; e = 0; f = 1;
    }
    else if(vert1 == connTriangle[ind2]) {

      /*--- The jj-index corresponds to a combination of the i and j index
            and the ii-index corresponds to the j-index.  ---*/
      if(vert2 != connTriangle[ind0]) verticesDontMatch = true;

      a = 0; b = 0; c = 1; d = nPoly; e = -1; f = -1;
    }
    else {
      verticesDontMatch = true;
    }
  }
  else if(vert0 == connTriangle[ind2]) {

    /*--- Vert0 coincides with the third vertex of the face connectivity.
          Check the situation for the neighboring vertices.    ---*/
    if(vert1 == connTriangle[ind0]) {

      /*--- The ii-index corresponds to a combination of the i and j index
            and the jj-index corresponds with the i-index.  ---*/
      if(vert2 != connTriangle[ind1]) verticesDontMatch = true;

      a = nPoly; b = -1; c = -1; d = 0; e = 1; f = 0;
    }
    else if(vert1 == connTriangle[ind1]) {

      /*--- The jj-index corresponds to a combination of the i and j index
            and the ii-index corresponds with the i-index.   ---*/
      if(vert2 != connTriangle[ind0]) verticesDontMatch = true;

      a = 0; b = 1; c = 0; d = nPoly; e = -1; f = -1;
    }
    else {
      verticesDontMatch = true;
    }
  }
  else {
    verticesDontMatch = true;
  }

  /*--- If non-matching vertices have been found, terminate with an error message. ---*/
  if( verticesDontMatch )
    SU2_MPI::Error("Corner vertices do not match. This should not happen.",
                   CURRENT_FUNCTION);

  /*--- Copy the connectivity, such that things works out correctly when carrying
        out the renumbering.  ---*/
  vector<unsigned short> connTriangleOr = connTriangle;

  /*--- Loop over the vertices of the original face to copy the connectivity data. ---*/
  unsigned short ind = 0;
  for(unsigned short j=0; j<=nPoly; ++j) {
    for(unsigned short i=0; i<=(nPoly-j); ++i, ++ind) {

      /*--- Determine the ii and jj indices of the target, convert it to
            a 1D index and shore the modified index in connTriangle. ---*/
      const unsigned short ii = a + i*b + j*c;
      const unsigned short jj = d + i*e + j*f;

      const unsigned short iind = jj*(nPoly+1) + ii - jj*(jj-1)/2;

      connTriangle[iind] = connTriangleOr[ind];
    }
  }
}

void CFEMStandardTriBase::ConvertCoor1DFaceTo2DTriangle(const vector<passivedouble> rLine,
                                                        const unsigned short        faceID,
                                                        const unsigned short        orientation,
                                                        vector<passivedouble>       &rTriangle,
                                                        vector<passivedouble>       &sTriangle) {

  /*--- Easier storage of the number of points and allocate the memory for
        rTriangle and sTriangle. ---*/
  const unsigned short nPoints = rLine.size();
  rTriangle.resize(nPoints);
  sTriangle.resize(nPoints);

  /*--- Determine the face ID and set the values of rTriangle and sTriangle for an
        orientation of 0, i.e. assuming that the element is on side 0 of the face. ---*/
  switch( faceID ) {
    case 0: {
      for(unsigned short i=0; i<nPoints; ++i) {
        rTriangle[i] =  rLine[i];
        sTriangle[i] = -1.0;
      }
      break;
    }

    case 1: {
      for(unsigned short i=0; i<nPoints; ++i) {
        rTriangle[i] = -rLine[i];
        sTriangle[i] =  rLine[i];
      }
      break;
    }

    case 2: {
      for(unsigned short i=0; i<nPoints; ++i) {
        rTriangle[i] = -1.0;
        sTriangle[i] = -rLine[i];
      }
      break;
    }

    default:
      SU2_MPI::Error("This should not happen.", CURRENT_FUNCTION); 
  }

  /*--- If this is an element on side 1 of the face, indicated by
        orientation == 1, the sequence of rTriangle and sTriangle
        should be reversed. ---*/
  if(orientation == 1) {
    reverse(rTriangle.begin(), rTriangle.end());
    reverse(sTriangle.begin(), sTriangle.end());
  }
}

void CFEMStandardTriBase::DerLagBasisIntPointsTriangle(const unsigned short                   mPoly,
                                                       const vector<passivedouble>            &rDOFs,
                                                       const vector<passivedouble>            &sDOFs,
                                                       const vector<passivedouble>            &rInt,
                                                       const vector<passivedouble>            &sInt,
                                                       vector<ColMajorMatrix<passivedouble> > &derLag) {

  /*--- Determine the padded number of the total number of integration points. ---*/
  const unsigned short nIntTot    = rInt.size();
  const unsigned short nIntTotPad = PaddedValue(nIntTot);

  /*--- Determine the inverse of the Vandermonde matrix of the DOFs. ---*/
  CSquareMatrixCM VInv(rDOFs.size());
  VandermondeTriangle(mPoly, rDOFs, sDOFs, VInv.GetMat());
  VInv.Invert();

  /*--- Determine the gradient of the Vandermonde matrix of the integration points. Make
        sure to allocate the number of rows to nIntTotPad and initialize them to zero. ---*/
  ColMajorMatrix<passivedouble> VDr(nIntTotPad,rDOFs.size()),
                                VDs(nIntTotPad,rDOFs.size());
  VDr.setConstant(0.0);
  VDs.setConstant(0.0);

  GradVandermondeTriangle(mPoly, rInt, sInt, VDr, VDs);

  /*--- The gradients of the Lagrangian basis functions can be obtained by
        multiplying VDr, VDs and VInv. ---*/
  derLag.resize(2);
  VInv.MatMatMult('R', VDr, derLag[0]);
  VInv.MatMatMult('R', VDs, derLag[1]);

  /*--- Check if the sum of the elements of the relevant rows of derLag is 0. ---*/
  CheckRowSum(nIntTot, rDOFs.size(), 0.0, derLag[0]);
  CheckRowSum(nIntTot, rDOFs.size(), 0.0, derLag[1]);
}

void CFEMStandardTriBase::HesLagBasisIntPointsTriangle(const unsigned short                   mPoly,
                                                       const vector<passivedouble>            &rDOFs,
                                                       const vector<passivedouble>            &sDOFs,
                                                       const vector<passivedouble>            &rInt,
                                                       const vector<passivedouble>            &sInt,
                                                       vector<ColMajorMatrix<passivedouble> > &hesLag) {

  /*--- Determine the padded number of the total number of integration points. ---*/
  const unsigned short nIntTot    = rInt.size();
  const unsigned short nIntTotPad = PaddedValue(nIntTot);

  /*--- Determine the inverse of the Vandermonde matrix of the DOFs. ---*/
  CSquareMatrixCM VInv(rDOFs.size());
  VandermondeTriangle(mPoly, rDOFs, sDOFs, VInv.GetMat());
  VInv.Invert();

  /*--- Determine the Hessian of the Vandermonde matrix of the integration points. Make
        sure to allocate the number of rows to nIntTotPad and initialize them to zero. ---*/
  ColMajorMatrix<passivedouble> VDr2(nIntTotPad,rDOFs.size()),
                                VDs2(nIntTotPad,rDOFs.size()),
                                VDrs(nIntTotPad,rDOFs.size());
  VDr2.setConstant(0.0);
  VDs2.setConstant(0.0);
  VDrs.setConstant(0.0);

  HesVandermondeTriangle(mPoly, rInt, sInt, VDr2, VDs2, VDrs);

  /*--- The Hessian of the Lagrangian basis functions can be obtained by
        multiplying VDr2, VDs2, VDrs and VInv. ---*/
  hesLag.resize(3);
  VInv.MatMatMult('R', VDr2, hesLag[0]);
  VInv.MatMatMult('R', VDs2, hesLag[1]);
  VInv.MatMatMult('R', VDrs, hesLag[2]);

  /*--- Check if the sum of the elements of the relevant rows of hesLagBasisInt is 0. ---*/
  CheckRowSum(nIntTot, rDOFs.size(), 0.0, hesLag[0]);
  CheckRowSum(nIntTot, rDOFs.size(), 0.0, hesLag[1]);
  CheckRowSum(nIntTot, rDOFs.size(), 0.0, hesLag[2]);
}

void CFEMStandardTriBase::LagBasisIntPointsTriangle(const unsigned short          mPoly,
                                                    const vector<passivedouble>   &rDOFs,
                                                    const vector<passivedouble>   &sDOFs,
                                                    const vector<passivedouble>   &rInt,
                                                    const vector<passivedouble>   &sInt,
                                                    ColMajorMatrix<passivedouble> &lag) {

  /*--- Determine the padded number of the total number of integration points. ---*/
  const unsigned short nIntTot    = rInt.size();
  const unsigned short nIntTotPad = PaddedValue(nIntTot);

  /*--- Determine the inverse of the Vandermonde matrix of the DOFs. ---*/
  CSquareMatrixCM VInv(rDOFs.size());
  VandermondeTriangle(mPoly, rDOFs, sDOFs, VInv.GetMat());
  VInv.Invert();

  /*--- Determine the Vandermonde matrix of the integration points. Make sure to
        allocate the number of rows to nIntTotPad and initialize them to zero. ---*/
  ColMajorMatrix<passivedouble> V(nIntTotPad,rDOFs.size());
  V.setConstant(0.0);
  VandermondeTriangle(mPoly, rInt, sInt, V);

  /*--- The Lagrangian basis functions can be obtained by multiplying
        V and VInv. ---*/
  VInv.MatMatMult('R', V, lag);

  /*--- Check if the sum of the elements of the relevant rows of lag is 1. ---*/
  CheckRowSum(nIntTot, rDOFs.size(), 1.0, lag);
}

void CFEMStandardTriBase::GradVandermondeTriangle(const unsigned short          mPoly,
                                                  const vector<passivedouble>   &r,
                                                  const vector<passivedouble>   &s,
                                                  ColMajorMatrix<passivedouble> &VDr,
                                                  ColMajorMatrix<passivedouble> &VDs) {

  /*--- Abbreviate sqrt(2), which is the scaling factor for the orthonormal basis
        functions of a triangle. ---*/
  const passivedouble sqrt2 = sqrt(2.0);

  /*--- For a triangle the orthogonal basis for the reference element is obtained
        by a combination of a Jacobi polynomial and a Legendre polynomial. This
        is the result of the orthonormalization of the monomial basis.
        Note that the sequence of the i and j loop must be identical to
        the evaluation of the Vandermonde matrix itself. ---*/
  unsigned short ii = 0;
  for(unsigned short i=0; i<=mPoly; ++i) {
    for(unsigned short j=0; j<=(mPoly-i); ++j, ++ii) {
      for(unsigned k=0; k<r.size(); ++k) {

        /*--- Determine the coefficients a and b. ---*/
        passivedouble a;
        if(fabs(s[k]-1.0) < 1.e-8) a = -1.0;
        else a = 2.0*(1.0+r[k])/(1.0-s[k]) - 1.0;

        const passivedouble b = s[k];

        /*--- Determine the value of the two 1D contributions to the 2D
              basis functions as well as the gradients of these 
              contributions w.r.t. their arguments. ---*/
        const passivedouble fa  = NormJacobi(i,0,    0,a);
        const passivedouble gb  = NormJacobi(j,2*i+1,0,b);
        const passivedouble dfa = GradNormJacobi(i,0,    0,a);
        const passivedouble dgb = GradNormJacobi(j,2*i+1,0,b);

        /*--- Computation of the powers of (1-b) that occur in the expressions for
              the gradients. Note the safeguard to avoid division by zero. This is
              allowed, because the implementation is such that when this clipping
              is active, it is multiplied by zero. ---*/
        const passivedouble tmpbi   = i > 0 ? pow((1.0-b), i)   : 1.0;
        const passivedouble tmpbim1 = i > 1 ? pow((1.0-b), i-1) : 1.0;

        /*--- Compute the derivatives of the basis function. ---*/
        VDr(k,ii) = sqrt2*tmpbim1* 2.0*dfa*gb;
        VDs(k,ii) = sqrt2*tmpbim1*((a+1)*dfa*gb - i*fa*gb)
                  + sqrt2*tmpbi  * fa*dgb;
      }
    }
  }
}

void CFEMStandardTriBase::HesVandermondeTriangle(const unsigned short          mPoly,
                                                 const vector<passivedouble>   &r,
                                                 const vector<passivedouble>   &s,
                                                 ColMajorMatrix<passivedouble> &VDr2,
                                                 ColMajorMatrix<passivedouble> &VDs2,
                                                 ColMajorMatrix<passivedouble> &VDrs) {

  /*--- Abbreviate sqrt(2), which is the scaling factor for the orthonormal basis
        functions of a triangle. ---*/
  const passivedouble sqrt2 = sqrt(2.0);

  /*--- For a triangle the orthogonal basis for the reference element is obtained
        by a combination of a Jacobi polynomial and a Legendre polynomial. This
        is the result of the orthonormalization of the monomial basis.
        Note that the sequence of the i and j loop must be identical to
        the evaluation of the Vandermonde matrix itself. ---*/
  unsigned short ii = 0;
  for(unsigned short i=0; i<=mPoly; ++i) {
    for(unsigned short j=0; j<=(mPoly-i); ++j, ++ii) {
      for(unsigned short k=0; k<r.size(); ++k) {

        /*--- Determine the coefficients a and b. ---*/
        passivedouble a;
        if(fabs(s[k]-1.0) < 1.e-8) a = -1.0;
        else a = 2.0*(1.0+r[k])/(1.0-s[k]) - 1.0;

        const passivedouble b = s[k];

        /*--- Determine the value of the two 1D contributions to the 2D
              basis functions as well as the 1st and 2nd derivatives of
              these contributions w.r.t. their arguments. ---*/
        const passivedouble fa   = NormJacobi(i,0,    0,a);
        const passivedouble gb   = NormJacobi(j,2*i+1,0,b);
        const passivedouble dfa  = GradNormJacobi(i,0,    0,a);
        const passivedouble dgb  = GradNormJacobi(j,2*i+1,0,b);
        const passivedouble d2fa = HesNormJacobi(i,0,    0,a);
        const passivedouble d2gb = HesNormJacobi(j,2*i+1,0,b);

        /*--- Computation of the powers of (1-b) that occur in the expressions for
              the Hessian. Note the safeguard to avoid division by zero. This is
              allowed, because the implementation is such that when this clipping
              is active, it is multiplied by zero. ---*/
        const passivedouble tmpbi   = i > 0 ? pow((1.0-b), i)   : 1.0;
        const passivedouble tmpbim1 = i > 1 ? pow((1.0-b), i-1) : 1.0;
        const passivedouble tmpbim2 = i > 2 ? pow((1.0-b), i-2) : 1.0;

        /*--- Compute the 2nd derivative w.r.t. to r. ---*/
        VDr2(k,ii) = tmpbim2*4.0*d2fa*gb;

        /*--- Compute the 2nd derivative w.r.t. s. ---*/
        VDs2(k,ii) = tmpbim2*((a+1.0)*(a+1.0)*d2fa*gb
                   +          i*(i-1)*gb*fa
                   -          2.0*(i-1)*(a+1.0)*gb*dfa)
                   + tmpbim1* 2.0*dgb*(dfa*(a+1.0) - i*fa)
                   + tmpbi  * fa*d2gb;

        /*--- Compute the cross derivative w.r.t. r and s. ---*/
        VDrs(k,ii) = tmpbim2*2.0*gb*((a+1.0)*d2fa - (i-1)*dfa)
                   + tmpbim1*2.0*dfa*dgb;

        /*--- Multiply all 2nd derivatives with sqrt(2) to obtain the
              correct expressions. ---*/
        VDr2(k,ii) *= sqrt2;
        VDs2(k,ii) *= sqrt2;
        VDrs(k,ii) *= sqrt2;
      }
    }
  }
}

void CFEMStandardTriBase::VandermondeTriangle(const unsigned short          mPoly,
                                              const vector<passivedouble>   &r,
                                              const vector<passivedouble>   &s,
                                              ColMajorMatrix<passivedouble> &V) {

  /*--- Abbreviate sqrt(2), which is the scaling factor for the orthonormal basis
        functions of a triangle. ---*/
  const passivedouble sqrt2 = sqrt(2.0);

  /*--- For a triangle the orthogonal basis for the reference element is obtained
        by a combination of a Jacobi polynomial and a Legendre polynomial. This
        is the result of the orthonormalization of the monomial basis. ---*/
  unsigned short ii = 0;
  for(unsigned short i=0; i<=mPoly; ++i) {
    for(unsigned short j=0; j<=(mPoly-i); ++j, ++ii) {
      for(unsigned short k=0; k<r.size(); ++k) {

        /*--- Determine the coefficients a and b. ---*/
        passivedouble a;
        if(fabs(s[k]-1.0) < 1.e-8) a = -1.0;
        else a = 2.0*(1.0+r[k])/(1.0-s[k]) - 1.0;

        const passivedouble b = s[k];

        /*--- Determine the value of the current basis function in this point. ---*/
        passivedouble tmp = 1.0;
        if( i ) tmp = pow((1.0-b),i);
        V(k,ii) = sqrt2*tmp*NormJacobi(i,0,0,a)*NormJacobi(j,2*i+1,0,b);
      }
    }
  }
}

void CFEMStandardTriBase::LocationTriangleGridDOFsEquidistant(const unsigned short  mPoly,
                                                              vector<passivedouble> &r,
                                                              vector<passivedouble> &s) {

  /*--- Determine the number of DOFs of the triangle.
        Allocate the memory for r and s afterwards. ---*/
  const unsigned short nD = (mPoly+1)*(mPoly+2)/2;
  r.resize(nD);
  s.resize(nD);

  /*--- Determine the equidistant spacing in parametric space. ---*/
  const passivedouble dh = 2.0/mPoly;

  /*--- Double loop to determine the parametric coordinates. ---*/
  unsigned short ii = 0;
  for(unsigned short j=0; j<=mPoly; ++j) {
    const unsigned short uppBoundI = mPoly - j;
    for(unsigned short i=0; i<=uppBoundI; ++i, ++ii) {
      r[ii] = -1.0 + i*dh;
      s[ii] = -1.0 + j*dh;
    }
  }
}

void CFEMStandardTriBase::LocationTriangleGridDOFsLGL(const unsigned short  mPoly,
                                                      vector<passivedouble> &r,
                                                      vector<passivedouble> &s) {

  /*--- The code to determine the parametric coordinates of the DOFs of the
        triangle is a translation of the Matlab code belonging to the book
        Nodal Discontinuous Galerkin Methods, Algorithms, Analysis and Applications,
        written by Jan S. Hesthaven and Tim Warburton. ---*/

  /*--- Local parameters. ---*/
  const passivedouble sqrt3 = sqrt(3.0);
  const passivedouble alphaOpt[] = {0.0000, 0.0000, 0.0000, 1.4152, 0.1001, 0.2751,
                                    0.9800, 1.0999, 1.2832, 1.3648, 1.4773, 1.4959,
                                    1.5743, 1.5770, 1.6223, 1.6258};

  /*--- Determine the number of DOFs of the triangle. As this function is also
        called for prisms, the member variable nDOFs cannot be used for this. ---*/
  const unsigned short nD = (mPoly+1)*(mPoly+2)/2;

  /*--- Allocate the memory for the help variables, as well as r and s. ---*/
  vector<passivedouble> L1(nD), L2(nD), L3(nD), rout(nD), warpf1(nD), warpf2(nD), warpf3(nD);

  r.resize(nD);
  s.resize(nD);

  /*--- Determine the uniform spacing. ---*/
  /*--- Determine the equidistant spacing in parametric space. ---*/
  const passivedouble dh = 1.0/mPoly;

  /*-----------------------------------------------------------------------------*/
  /*--- Step 1: Create the equidistributed nodes on the equilateral triangle. ---*/
  /*-----------------------------------------------------------------------------*/

  unsigned short ii = 0;
  for(unsigned short j=0; j<=mPoly; ++j) {
    const unsigned short uppBoundI = mPoly - j;
    for(unsigned short i=0; i<=uppBoundI; ++i, ++ii) {
      L1[ii] = j*dh;
      L3[ii] = i*dh;
      L2[ii] = 1.0 - L1[ii] - L3[ii];

      r[ii] = L3[ii] - L2[ii];
      s[ii] = (2.0*L1[ii] - L2[ii] - L3[ii])/sqrt3;
    }
  }

  /*---------------------------------------------*/
  /*--- Step 2: Modify the node distribution. ---*/
  /*---------------------------------------------*/

  /*--- Set the optimal value of alp for this case. ---*/
  const passivedouble alp = (mPoly < 16) ? alphaOpt[mPoly] : 5.0/3.0;

  /*--- Compute the warp factors for the DOFs for each edge. ---*/
  for(unsigned short i=0; i<nD; ++i) rout[i] = L3[i] - L2[i];
  WarpFactor(mPoly, rout, warpf1);

  for(unsigned short i=0; i<nD; ++i) rout[i] = L1[i] - L3[i];
  WarpFactor(mPoly, rout, warpf2);

  for(unsigned short i=0; i<nD; ++i) rout[i] = L2[i] - L1[i];
  WarpFactor(mPoly, rout, warpf3);

  /*--- Loop over the DOFs to correct the coordinates. ---*/
  for(unsigned short i=0; i<nD; ++i) {

    /*--- Compute the blending functions for each edge. ---*/
    const passivedouble blend1 = 4.0*L2[i]*L3[i];
    const passivedouble blend2 = 4.0*L1[i]*L3[i];
    const passivedouble blend3 = 4.0*L1[i]*L2[i];

    /*--- Compute the combined blending and warp factor. ---*/
    const passivedouble warp1 = blend1*warpf1[i]*(1.0 + alp*alp*L1[i]*L1[i]);
    const passivedouble warp2 = blend2*warpf2[i]*(1.0 + alp*alp*L2[i]*L2[i]);
    const passivedouble warp3 = blend3*warpf3[i]*(1.0 + alp*alp*L3[i]*L3[i]);

    /*--- Compute the new coordinates. The multiplication factors in front
          of the warp factors corresponds to the cosine and sine of 0,
          2*pi/3 and 4*pi/3 respectively. ---*/
    r[i] = r[i] + warp1 - 0.5*(warp2 + warp3);
    s[i] = s[i] +   0.5*sqrt3*(warp2 - warp3);
  }
}

void CFEMStandardTriBase::SubConnLinearElements(void) {

  /*--- The triangle is split into several linear triangles.
        Set the VTK sub-types accordingly. ---*/
  VTK_SubType1 = TRIANGLE;
  VTK_SubType2 = NONE;

  /*--- Initialize the counter for the edges jj. ---*/
  unsigned short jj = 0;

  /*--- Loop over subedges of the left boundary of the standard triangle. ---*/
  for(unsigned short j=0; j<nPoly; ++j) {

    /*--- Check if the "down" elements must be written. ---*/
    if( j ) {

      /*--- Offset of the relevant DOF on the previous row. ---*/
      const unsigned short kk = jj - (nPoly + 1 - j);

      /*--- Loop over the sub-elements of this edge. ---*/
      for(unsigned short i=0; i<(nPoly-j); ++i) {
        const unsigned short n0 = jj + i;
        const unsigned short n1 = kk + i;
        const unsigned short n2 = n0 + 1;

        subConn1ForPlotting.push_back(n0);
        subConn1ForPlotting.push_back(n1);
        subConn1ForPlotting.push_back(n2);
      }
    }

    /*--- The "upp" elements must always be written.
          Determine the offset of the DOF on the next row. ---*/
    const unsigned short kk = jj + (nPoly + 1 - j);

    /*--- Loop over the sub-elements of this edge. ---*/
    for(unsigned short i=0; i<(nPoly-j); ++i) {
      const unsigned short n0 = jj + i;
      const unsigned short n1 = n0 + 1;
      const unsigned short n2 = kk + i;

      subConn1ForPlotting.push_back(n0);
      subConn1ForPlotting.push_back(n1);
      subConn1ForPlotting.push_back(n2);
    }

    /*--- Set jj to kk for the next edge. ---*/
    jj = kk;
  }
}

void CFEMStandardTriBase::SubConnLinearElementsFace(int val_faceID_Elem) {

  /*--- This is an added functionality specifically for surface output.
        The utility is similar to SubConnLinearElements, but store the 
        nodes (available in gridConnFaces) using the face ID w.r.t the 
        volume when this base class is considered a face. ---*/
        
  /*--- The triangle is split into several linear triangles.
        Set the VTK sub-types accordingly. ---*/
  VTK_SubType1 = TRIANGLE;
  VTK_SubType2 = NONE;

  /*--- Initialize the counter for the edges jj. ---*/
  unsigned short jj = 0;

  /*--- Loop over subedges of the left boundary of the standard triangle. ---*/
  for(unsigned short j=0; j<nPoly; ++j) {

    /*--- Check if the "down" elements must be written. ---*/
    if( j ) {

      /*--- Offset of the relevant DOF on the previous row. ---*/
      const unsigned short kk = jj - (nPoly + 1 - j);

      /*--- Loop over the sub-elements of this edge. ---*/
      for(unsigned short i=0; i<(nPoly-j); ++i) {
        const unsigned short n0 = jj + i;
        const unsigned short n1 = kk + i;
        const unsigned short n2 = n0 + 1;

        subConn1ForPlotting.push_back(gridConnFaces[val_faceID_Elem][n0]);
        subConn1ForPlotting.push_back(gridConnFaces[val_faceID_Elem][n1]);
        subConn1ForPlotting.push_back(gridConnFaces[val_faceID_Elem][n2]);
      }
    }

    /*--- The "upp" elements must always be written.
          Determine the offset of the DOF on the next row. ---*/
    const unsigned short kk = jj + (nPoly + 1 - j);

    /*--- Loop over the sub-elements of this edge. ---*/
    for(unsigned short i=0; i<(nPoly-j); ++i) {
      const unsigned short n0 = jj + i;
      const unsigned short n1 = n0 + 1;
      const unsigned short n2 = kk + i;

      subConn1ForPlotting.push_back(gridConnFaces[val_faceID_Elem][n0]);
      subConn1ForPlotting.push_back(gridConnFaces[val_faceID_Elem][n1]);
      subConn1ForPlotting.push_back(gridConnFaces[val_faceID_Elem][n2]);
    }

    /*--- Set jj to kk for the next edge. ---*/
    jj = kk;
  }
}

/*----------------------------------------------------------------------------------*/
/*                Private member functions of CFEMStandardTriBase.                  */
/*----------------------------------------------------------------------------------*/

void CFEMStandardTriBase::WarpFactor(const unsigned short        mPoly,
                                     const vector<passivedouble> &rout,
                                     vector<passivedouble>       &warp) {
    /*--- The code to determine the warp factor is a translation of the Matlab code
        belonging to the book
        Nodal Discontinuous Galerkin Methods, Algorithms, Analysis and Applications,
        written by Jan S. Hesthaven and Tim Warburton. ---*/

  /*--- Determine the number of DOFs for which the warp factor
        must be computed. ---*/
  const unsigned short nD = rout.size();

  /*--- Compute the 1D Gauss-Lobato and an equidistant node distribution. ---*/
  vector<passivedouble> r1DEqui, r1DLGL;
  Location1DGridDOFsEquidistant(mPoly, r1DEqui);
  Location1DGridDOFsLGL(mPoly, r1DLGL);

  /*--- Determine the 1D Vandermonde matrix based on the equidistant node
        distribution. Invert the matrix afterwards. ---*/
  CSquareMatrixCM V1D(mPoly+1);
  Vandermonde1D(mPoly, r1DEqui, V1D.GetMat());
  V1D.Invert();

  /*--- Determine the Lagrange polynomials at rout. This is accomplished by
        first computing the Legendre polynomials at rout and multiply
        this value by V1D. The result is stored in Lmat. ---*/
  ColMajorMatrix<passivedouble> Pmat(nD,mPoly+1), Lmat;
  Vandermonde1D(mPoly, rout, Pmat);
  V1D.MatMatMult('R', Pmat, Lmat);

  /*--- Determine the difference between the 1D LGL points and
        the equidistant points. The result is stored in r1DLGL. ---*/
  for(unsigned short i=0; i<=mPoly; ++i) r1DLGL[i] -= r1DEqui[i];

  /*--- Compute the unscaled warp factors, which is Lmat*r1DLGL. ---*/
  for(unsigned short j=0; j<nD; ++j) {
    warp[j] = 0.0;
    for(unsigned short i=0; i<=mPoly; ++i)
      warp[j] += Lmat(j,i)*r1DLGL[i];
  }

  /*--- Scale the warp factor by 1-r*r. ---*/
  for(unsigned short i=0; i<nD; ++i) {

    /*--- Take care of the exceptional case that r = 1. ---*/
    passivedouble zerof = 0.0;
    if(fabs(rout[i]) < (1.0-1.e-10)) zerof = 1.0;

    /*--- Compute 1-r*r. ---*/
    const passivedouble sf = 1.0 - zerof*rout[i]*rout[i];

    /*--- Compute the scaled warp factor. ---*/
    warp[i] = warp[i]/sf + warp[i]*(zerof-1.0);
    if(fabs(warp[i]) < 1.e-10) warp[i] = 0.0;
  }
}

void CFEMStandardTriBase::LocalGridConnFaces(void) {

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