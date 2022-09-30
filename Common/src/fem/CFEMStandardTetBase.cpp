/*!
 * \file CFEMStandardTetBase.cpp
 * \brief Functions for the class CFEMStandardTetBase.
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

#include "../../include/fem/CFEMStandardTetBase.hpp"
#include "../../include/toolboxes/CSquareMatrixCM.hpp"

/*----------------------------------------------------------------------------------*/
/*               Public member functions of CFEMStandardTetBase.                    */
/*----------------------------------------------------------------------------------*/

CFEMStandardTetBase::CFEMStandardTetBase(const unsigned short val_nPoly,
                                         const unsigned short val_orderExact)
  : CFEMStandardTriBase() {

  /*--- Store the command line arguments. ---*/
  VTK_Type   = TETRAHEDRON;
  nPoly      = val_nPoly;
  orderExact = val_orderExact;

  /*--- Determine the total number of DOFs and its padded version. ---*/
  nDOFs    = (nPoly+1)*(nPoly+2)*(nPoly+3)/6;
  nDOFsPad = PaddedValue(nDOFs);

  /*--- Determine the parametric location and weights of the
        integration rule of the tetrahedron. ---*/
  IntegrationPointsTetrahedron();

  /*--- Determine the total number of integration points
        and its padded version. ---*/
  nIntegration    = rTetInt.size();
  nIntegrationPad = PaddedValue(nIntegration);

  /*--- Allocate the memory for the padded number of integration points and
        initialize the weights to zero. This is done such that the padded
        values are initialized appropriately. ---*/
  wIntegration.resize(nIntegrationPad);
  wIntegration.setConstant(0.0);

  /*--- Copy the values from wTetInt to wIntegration. ---*/
  for(unsigned short i=0; i<nIntegration; ++i)
    wIntegration(i) = wTetInt[i];
}

/*----------------------------------------------------------------------------------*/
/*                Protected member functions of CFEMStandardTetBase.                */
/*----------------------------------------------------------------------------------*/

void CFEMStandardTetBase::ConvertCoor2DTriFaceTo3DTet(const vector<passivedouble> &rF,
                                                      const vector<passivedouble> &sF,
                                                      const unsigned short        faceID_Elem,
                                                      const unsigned short        orientation,
                                                      vector<passivedouble>       &rTet,
                                                      vector<passivedouble>       &sTet,
                                                      vector<passivedouble>       &tTet) {

  /*--- Determine the number of points on the triangular face. Afterwards, allocate
        the memory for rTet, sTet and tTet. ---*/
  const unsigned short nP = rF.size();

  rTet.resize(nP);
  sTet.resize(nP);
  tTet.resize(nP);

  /*--- Abbreviate rTet, sTet and tTet, such that the different cases
        can be put on one line. ---*/
  vector<passivedouble> &r = rTet, &s = sTet, &t = tTet;

  /*--- The values of rTet, sTet and tTet depend on both the face ID in the numbering
        of the tetrahedron as well as the orientation of the face w.r.t. the tetrahedron.
        Make this distinction and set the values accordingly. ---*/
  unsigned short k;
  switch( faceID_Elem ) {

    case 0: {
      switch( orientation ) {
        case 0: for(k=0; k<nP; ++k) {r[k]= rF[k]; s[k]= sF[k]; t[k]=-1.0;} break; 
        case 1: for(k=0; k<nP; ++k) {r[k]= sF[k]; s[k]= rF[k]; t[k]=-1.0;} break;
        case 2: for(k=0; k<nP; ++k) {s[k]= sF[k]; t[k]=-1.0;   r[k]=-1.0-rF[k]-sF[k];} break;
        case 3: for(k=0; k<nP; ++k) {r[k]= rF[k]; t[k]=-1.0;   s[k]=-1.0-rF[k]-sF[k];} break;
        default:
          SU2_MPI::Error(string("Invalid orientation for face 0. This should not happen."),
                         CURRENT_FUNCTION);
      }
      break;
    }

    case 1: {
      switch( orientation ) {
        case 0: for(k=0; k<nP; ++k) {r[k]= sF[k]; s[k]=-1.0;   t[k]= rF[k];} break;
        case 1: for(k=0; k<nP; ++k) {r[k]= rF[k]; s[k]=-1.0;   t[k]= sF[k];} break;
        case 2: for(k=0; k<nP; ++k) {r[k]= sF[k]; s[k]=-1.0;   t[k]=-1.0-rF[k]-sF[k];} break;
        case 3: for(k=0; k<nP; ++k) {s[k]=-1.0;   t[k]= rF[k]; r[k]=-1.0-rF[k]-sF[k];} break;
        default:
          SU2_MPI::Error(string("Invalid orientation for face 1. This should not happen."),
                         CURRENT_FUNCTION);
      }
      break;
    }

    case 2: {
      switch( orientation ) {
        case 0: for(k=0; k<nP; ++k) {r[k]=-1.0; s[k]= rF[k]; t[k]= sF[k];} break;
        case 1: for(k=0; k<nP; ++k) {r[k]=-1.0; s[k]= sF[k]; t[k]= rF[k];} break;
        case 2: for(k=0; k<nP; ++k) {r[k]=-1.0; t[k]= sF[k]; s[k]=-1.0-rF[k]-sF[k];} break;
        case 3: for(k=0; k<nP; ++k) {r[k]=-1.0; s[k]= rF[k]; t[k]=-1.0-rF[k]-sF[k];} break;
        default:
          SU2_MPI::Error(string("Invalid orientation for face 2. This should not happen."),
                         CURRENT_FUNCTION);
      }
      break;
    }

    case 3: {
      switch( orientation ) {
        case 0: for(k=0; k<nP; ++k) {s[k]=sF[k]; t[k]=rF[k]; r[k]=-1.0-rF[k]-sF[k];} break;
        case 1: for(k=0; k<nP; ++k) {s[k]=rF[k]; t[k]=sF[k]; r[k]=-1.0-rF[k]-sF[k];} break;
        case 2: for(k=0; k<nP; ++k) {r[k]=rF[k]; s[k]=sF[k]; t[k]=-1.0-rF[k]-sF[k];} break;
        case 3: for(k=0; k<nP; ++k) {r[k]=sF[k]; t[k]=rF[k]; s[k]=-1.0-rF[k]-sF[k];} break;
        default:
          SU2_MPI::Error(string("Invalid orientation for face 3. This should not happen."),
                         CURRENT_FUNCTION);
      }
      break;
    }

    default:
      SU2_MPI::Error(string("Invalid faceID. This should not happen."), CURRENT_FUNCTION);
  }
}

void CFEMStandardTetBase::DerLagBasisIntPointsTet(const unsigned short                   mPoly,
                                                  const vector<passivedouble>            &rDOFs,
                                                  const vector<passivedouble>            &sDOFs,
                                                  const vector<passivedouble>            &tDOFs,
                                                  const vector<passivedouble>            &rInt,
                                                  const vector<passivedouble>            &sInt,
                                                  const vector<passivedouble>            &tInt,
                                                  vector<ColMajorMatrix<passivedouble> > &derLag) {

  /*--- Determine the padded number of the total number of integration points. ---*/
  const unsigned short nIntTot    = rInt.size();
  const unsigned short nIntTotPad = PaddedValue(nIntTot);

  /*--- Determine the inverse of the Vandermonde matrix of the DOFs. ---*/
  CSquareMatrixCM VInv(rDOFs.size());
  VandermondeTetrahedron(mPoly, rDOFs, sDOFs, tDOFs, VInv.GetMat());
  VInv.Invert();

  /*--- Determine the gradient of the Vandermonde matrix of the integration points. Make
        sure to allocate the number of rows to nIntTotPad and initialize them to zero. ---*/
  ColMajorMatrix<passivedouble> VDr(nIntTotPad,rDOFs.size()),
                                VDs(nIntTotPad,rDOFs.size()),
                                VDt(nIntTotPad,rDOFs.size());
  VDr.setConstant(0.0);
  VDs.setConstant(0.0);
  VDt.setConstant(0.0);

  GradVandermondeTetrahedron(mPoly, rInt, sInt, tInt, VDr, VDs, VDt);

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

void CFEMStandardTetBase::HesLagBasisIntPointsTet(const unsigned short                   mPoly,
                                                  const vector<passivedouble>            &rDOFs,
                                                  const vector<passivedouble>            &sDOFs,
                                                  const vector<passivedouble>            &tDOFs,
                                                  const vector<passivedouble>            &rInt,
                                                  const vector<passivedouble>            &sInt,
                                                  const vector<passivedouble>            &tInt,
                                                  vector<ColMajorMatrix<passivedouble> > &hesLag) {

  /*--- Determine the padded number of the total number of integration points. ---*/
  const unsigned short nIntTot    = rInt.size();
  const unsigned short nIntTotPad = PaddedValue(nIntTot);

  /*--- Determine the inverse of the Vandermonde matrix of the DOFs. ---*/
  CSquareMatrixCM VInv(rDOFs.size());
  VandermondeTetrahedron(mPoly, rDOFs, sDOFs, tDOFs, VInv.GetMat());
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

  HesVandermondeTetrahedron(mPoly, rInt, sInt, tInt, VDr2, VDs2, VDt2, VDrs, VDrt, VDst);

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

void CFEMStandardTetBase::LagBasisIntPointsTet(const unsigned short          mPoly,
                                               const vector<passivedouble>   &rDOFs,
                                               const vector<passivedouble>   &sDOFs,
                                               const vector<passivedouble>   &tDOFs,
                                               const vector<passivedouble>   &rInt,
                                               const vector<passivedouble>   &sInt,
                                               const vector<passivedouble>   &tInt,
                                               ColMajorMatrix<passivedouble> &lag) {

  /*--- Determine the padded number of the total number of integration points. ---*/
  const unsigned short nIntTot    = rInt.size();
  const unsigned short nIntTotPad = PaddedValue(nIntTot);

  /*--- Determine the inverse of the Vandermonde matrix of the DOFs. ---*/
  CSquareMatrixCM VInv(rDOFs.size());
  VandermondeTetrahedron(mPoly, rDOFs, sDOFs, tDOFs, VInv.GetMat());
  VInv.Invert();

  /*--- Determine the Vandermonde matrix of the integration points. Make sure to
        allocate the number of rows to nIntTotPad and initialize them to zero. ---*/
  ColMajorMatrix<passivedouble> V(nIntTotPad,rDOFs.size());
  V.setConstant(0.0);
  VandermondeTetrahedron(mPoly, rInt, sInt, tInt, V);

  /*--- The Lagrangian basis functions can be obtained by multiplying
        V and VInv. ---*/
  VInv.MatMatMult('R', V, lag);

  /*--- Check if the sum of the elements of the relevant rows of lag is 1. ---*/
  CheckRowSum(nIntTot, rDOFs.size(), 1.0, lag);
}

void CFEMStandardTetBase::GradVandermondeTetrahedron(const unsigned short          mPoly,
                                                     const vector<passivedouble>   &r,
                                                     const vector<passivedouble>   &s,
                                                     const vector<passivedouble>   &t,
                                                     ColMajorMatrix<passivedouble> &VDr,
                                                     ColMajorMatrix<passivedouble> &VDs,
                                                     ColMajorMatrix<passivedouble> &VDt) {


  /*--- Abbreviate sqrt(8), which is the scaling factor for the orthonormal basis
        functions of a tetrahedron. ---*/
  const passivedouble sqrt8 = sqrt(8.0);

  /*--- For a tetrahedron the orthogonal basis for the reference element is obtained by a
        combination of Jacobi polynomials (of which the Legendre polynomials is a special
        case). This is the result of the orthonormalization of the monomial basis.
        Note that the sequence of the i, j and k loop must be identical to
        the evaluation of the Vandermonde matrix itself. ---*/
  unsigned short ii = 0;
  for(unsigned short i=0; i<=mPoly; ++i) {
    for(unsigned short j=0; j<=(mPoly-i); ++j) {
      for(unsigned short k=0; k<=(mPoly-i-j); ++k, ++ii) {
        for(unsigned short l=0; l<r.size(); ++l) {

          /*--- Determine the coefficients a, b and c. ---*/
          passivedouble a, b;
          passivedouble tmp = s[l] + t[l];
          if(fabs(tmp) < 1.e-8) a = -1.0;
          else                  a = -1.0 - 2.0*(1.0+r[l])/tmp;

          tmp = 1.0 - t[l];
          if(fabs(tmp) < 1.e-8) b = -1.0;
          else                  b = -1.0 + 2.0*(1.0+s[l])/tmp;

          const passivedouble c = t[l];

          /*--- Determine the value of the three 1D contributions to the 3D basis functions
                as well as the gradients of these contributions w.r.t. their arguments. ---*/
          const passivedouble fa  = NormJacobi(i,0,    0,a);
          const passivedouble gb  = NormJacobi(j,2*i+1,0,b);
          const passivedouble hc  = NormJacobi(k,2*(i+j+1),0,c);
          const passivedouble dfa = GradNormJacobi(i,0,    0,a);
          const passivedouble dgb = GradNormJacobi(j,2*i+1,0,b);
          const passivedouble dhc = GradNormJacobi(k,2*(i+j+1),0,c);

          /*--- Computation of the powers of (1-b) and (1-c) that occur in the
                expressions for the Hessian. Note the safeguard to avoid division
                by zero. This is allowed, because the implementation is such that
                when this clipping is active, it is multiplied by zero. ---*/
          const passivedouble tmpbi   = i > 0 ? pow((1.0-b), i)   : 1.0;
          const passivedouble tmpbim1 = i > 1 ? pow((1.0-b), i-1) : 1.0;

          const passivedouble tmpcij   = (i+j) > 0 ? pow((1.0-c), i+j)   : 1.0;
          const passivedouble tmpcijm1 = (i+j) > 1 ? pow((1.0-c), i+j-1) : 1.0;

          /*--- Compute the derivatives of the basis function. ---*/
          VDr(l,ii) = tmpbim1*tmpcijm1*4.0*dfa*gb*hc;
          VDs(l,ii) = tmpbim1*tmpcijm1*2.0*gb*hc*((a+1.0)*dfa - i*fa)
                    + tmpbi  *tmpcijm1*2.0*fa*dgb*hc;
          VDt(l,ii) = tmpbim1*tmpcijm1*gb*hc*(2.0*(a+1.0)*dfa - (b+1.0)*i*fa)
                    + tmpbi  *tmpcijm1*fa*hc*((b+1.0)*dgb - (i+j)*gb)
                    + tmpbi  *tmpcij  *fa*gb*dhc;

          /*--- Multiply all derivatives with sqrt(8) to obtain the
                correct expressions. ---*/
          VDr(l,ii) *= sqrt8;
          VDs(l,ii) *= sqrt8;
          VDt(l,ii) *= sqrt8;
        }
      }
    }
  }
}

void CFEMStandardTetBase::HesVandermondeTetrahedron(const unsigned short          mPoly,
                                                    const vector<passivedouble>   &r,
                                                    const vector<passivedouble>   &s,
                                                    const vector<passivedouble>   &t,
                                                    ColMajorMatrix<passivedouble> &VDr2,
                                                    ColMajorMatrix<passivedouble> &VDs2,
                                                    ColMajorMatrix<passivedouble> &VDt2,
                                                    ColMajorMatrix<passivedouble> &VDrs,
                                                    ColMajorMatrix<passivedouble> &VDrt,
                                                    ColMajorMatrix<passivedouble> &VDst) {

  /*--- Abbreviate sqrt(8), which is the scaling factor for the orthonormal basis
        functions of a tetrahedron. ---*/
  const passivedouble sqrt8 = sqrt(8.0);

  /*--- For a tetrahedron the orthogonal basis for the reference element is obtained by a
        combination of Jacobi polynomials (of which the Legendre polynomials is a special
        case). This is the result of the orthonormalization of the monomial basis.
        Note that the sequence of the i, j and k loop must be identical to
        the evaluation of the Vandermonde matrix itself. ---*/
  unsigned short ii = 0;
  for(unsigned short i=0; i<=mPoly; ++i) {
    for(unsigned short j=0; j<=(mPoly-i); ++j) {
      for(unsigned short k=0; k<=(mPoly-i-j); ++k, ++ii) {
        for(unsigned short l=0; l<r.size(); ++l) {

          /*--- Determine the coefficients a, b and c. ---*/
          passivedouble a, b;
          passivedouble tmp = s[l] + t[l];
          if(fabs(tmp) < 1.e-8) a = -1.0;
          else                  a = -1.0 - 2.0*(1.0+r[l])/tmp;

          tmp = 1.0 - t[l];
          if(fabs(tmp) < 1.e-8) b = -1.0;
          else                  b = -1.0 + 2.0*(1.0+s[l])/tmp;

          const passivedouble c = t[l];

          /*--- Determine the value of the three 1D contributions to the 3D
                basis functions as well as the 1st and 2nd derivatives of
                these contributions w.r.t. their arguments. ---*/
          const passivedouble fa   = NormJacobi(i,0,    0,a);
          const passivedouble gb   = NormJacobi(j,2*i+1,0,b);
          const passivedouble hc   = NormJacobi(k,2*(i+j+1),0,c);
          const passivedouble dfa  = GradNormJacobi(i,0,    0,a);
          const passivedouble dgb  = GradNormJacobi(j,2*i+1,0,b);
          const passivedouble dhc  = GradNormJacobi(k,2*(i+j+1),0,c);
          const passivedouble d2fa = HesNormJacobi(i,0,    0,a);
          const passivedouble d2gb = HesNormJacobi(j,2*i+1,0,b);
          const passivedouble d2hc = HesNormJacobi(k,2*(i+j+1),0,c);

          /*--- Computation of the powers of (1-b) and (1-c) that occur in the
                expressions for the Hessian. Note the safeguard to avoid division
                by zero. This is allowed, because the implementation is such that
                when this clipping is active, it is multiplied by zero. ---*/
          const passivedouble tmpbi   = i > 0 ? pow((1.0-b), i)   : 1.0;
          const passivedouble tmpbim1 = i > 1 ? pow((1.0-b), i-1) : 1.0;
          const passivedouble tmpbim2 = i > 2 ? pow((1.0-b), i-2) : 1.0;

          const passivedouble tmpcij   = (i+j) > 0 ? pow((1.0-c), i+j)   : 1.0;
          const passivedouble tmpcijm1 = (i+j) > 1 ? pow((1.0-c), i+j-1) : 1.0;
          const passivedouble tmpcijm2 = (i+j) > 2 ? pow((1.0-c), i+j-2) : 1.0;

          /*--- Compute the 2nd derivative w.r.t. to r. ---*/
          VDr2(l,ii) = tmpbim2*tmpcijm2*16.0*d2fa*gb*hc;

          /*--- Compute the 2nd derivative w.r.t. to s. ---*/
          VDs2(l,ii) = tmpbim2*tmpcijm2*(4.0*(a+1.0)*(a+1.0)*d2fa*gb*hc
                     -                   8.0*(i-1)*(a+1.0)*dfa*gb*hc
                     +                   4.0*i*(i-1)*fa*gb*hc)
                     + tmpbim1*tmpcijm2*(8.0*(a+1.0)*dfa*dgb*hc
                     -                   8.0*i*fa*dgb*hc)
                     + tmpbi  *tmpcijm2* 4.0*fa*d2gb*hc;

          /*--- Compute the 2nd derivative w.r.t. to t. ---*/
          VDt2(l,ii) = tmpbim2*tmpcijm2*(4.0*(a+1.0)*(a+1.0)*d2fa*gb*hc
                     -                   2.0*(a+1.0)*(b+3.0)*(i-1)*dfa*gb*hc
                     +                   i*(i-1)*(b+1.0)*(b+1.0)*fa*gb*hc)
                     + tmpbim1*tmpcijm2*(4.0*(a+1.0)*(b+1.0)*dfa*dgb*hc
                     -                   2.0*(b+1.0)*(b+1.0)*i*fa*dgb*hc
                     -                   2.0*(a+1.0)*(i+2*j-1)*dfa*gb*hc
                     +                   2.0*(b+1.0)*i*(i+j-1)*fa*gb*hc)
                     + tmpbi  *tmpcijm2*((b+1.0)*(b+1.0)*fa*d2gb*hc
                     -                   2.0*(b+1.0)*(i+j-1)*fa*dgb*hc
                     +                   (i+j)*(i+j-1)*fa*gb*hc)
                     + tmpbim1*tmpcijm1*(4.0*(a+1.0)*dfa*gb*dhc
                     -                   2.0*(b+1.0)*i*fa*gb*dhc)
                     + tmpbi  *tmpcijm1*(2.0*(b+1.0)*fa*dgb*dhc
                     -                   2.0*(i+j)*fa*gb*dhc)
                     + tmpbi  *tmpcij  * fa*gb*d2hc;

          /*--- Compute the cross derivative w.r.t. r and s. ---*/
          VDrs(l,ii) = tmpbim2*tmpcijm2*8.0*gb*hc*((a+1.0)*d2fa - (i-1)*dfa)
                     + tmpbim1*tmpcijm2*8.0*dfa*dgb*hc;

          /*--- Compute the cross derivative w.r.t. r and t. ---*/
          VDrt(l,ii) = tmpbim2*tmpcijm2*(8.0*(a+1.0)*d2fa*gb*hc
                     -                   4.0*(b+1.0)*(i-1)*dfa*gb*hc)
                     + tmpbim1*tmpcijm2*(4.0*(b+1.0)*dfa*dgb*hc
                     -                   4.0*(i+j-1)*dfa*gb*hc)
                     + tmpbim1*tmpcijm1* 4.0*dfa*gb*dhc;

          /*--- Compute the cross derivative w.r.t. s and t. ---*/
          VDst(l,ii) = tmpbim2*tmpcijm2*(4.0*(a+1.0)*(a+1.0)*d2fa*gb*hc
                     -                   2.0*(a+1.0)*(b+1.0)*(i-1)*dfa*gb*hc
                     -                   4.0*(a+1.0)*(i-1)*dfa*gb*hc
                     +                   2.0*(b+1.0)*i*(i-1)*fa*gb*hc)
                     + tmpbim1*tmpcijm2*(2.0*(a+1.0)*(b+1.0)*dfa*dgb*hc
                     -                   2.0*(a+1.0)*(i+j-1)*dfa*gb*hc
                     +                   4.0*(a+1.0)*dfa*dgb*hc
                     -                   4.0*(b+1.0)*i*fa*dgb*hc
                     +                   2.0*i*(i+j-1)*fa*gb*hc)
                     + tmpbim1*tmpcijm1*(2.0*(a+1.0)*dfa*gb*dhc
                     -                   2.0*i*fa*gb*dhc)
                     + tmpbi  *tmpcijm2*(2.0*(b+1.0)*fa*d2gb*hc
                     -                   2.0*(i+j-1)*fa*dgb*hc)
                     + tmpbi  *tmpcijm1* 2.0*fa*dgb*dhc;

          /*--- Multiply all 2nd derivatives with sqrt(8) to obtain the
                correct expressions. ---*/
          VDr2(l,ii) *= sqrt8;
          VDs2(l,ii) *= sqrt8;
          VDt2(l,ii) *= sqrt8;
          VDrs(l,ii) *= sqrt8;
          VDrt(l,ii) *= sqrt8;
          VDst(l,ii) *= sqrt8;
        }
      }
    }
  }
}

void CFEMStandardTetBase::VandermondeTetrahedron(const unsigned short          mPoly,
                                                 const vector<passivedouble>   &r,
                                                 const vector<passivedouble>   &s,
                                                 const vector<passivedouble>   &t,
                                                 ColMajorMatrix<passivedouble> &V) {

  /*--- For a tetrahedron the orthogonal basis for the reference element is obtained by a
        combination of Jacobi polynomials (of which the Legendre polynomials is a special
        case). This is the result of the orthonormalization of the monomial basis. ---*/
  unsigned short ii = 0;
  for(unsigned short i=0; i<=mPoly; ++i) {
    for(unsigned short j=0; j<=(mPoly-i); ++j) {
      for(unsigned short k=0; k<=(mPoly-i-j); ++k, ++ii) {
        for(unsigned short l=0; l<r.size(); ++l) {

          /*--- Determine the coefficients a, b and c. ---*/
          passivedouble a, b;
          passivedouble tmp = s[l] + t[l];
          if(fabs(tmp) < 1.e-8) a = -1.0;
          else                  a = -1.0 - 2.0*(1.0+r[l])/tmp;

          tmp = 1.0 - t[l];
          if(fabs(tmp) < 1.e-8) b = -1.0;
          else                  b = -1.0 + 2.0*(1.0+s[l])/tmp;

          const passivedouble c = t[l];

          /*--- Determine the value of the current basis function in this point. ---*/
          passivedouble tmpb = 1.0;
          if( i )   tmpb = pow((1.0-b),i);

          passivedouble tmpc = 1.0;
          if( i+j ) tmpc = pow((1.0-c),i+j);

          V(l,ii) = sqrt(8.0)*tmpb*tmpc*NormJacobi(i,0,0,a)*NormJacobi(j,2*i+1,0,b)
                  * NormJacobi(k,2*(i+j+1),0,c);
        }
      }
    }
  }
}

void CFEMStandardTetBase::LocationTetGridDOFsEquidistant(const unsigned short  mPoly,
                                                         vector<passivedouble> &rDOFs,
                                                         vector<passivedouble> &sDOFs,
                                                         vector<passivedouble> &tDOFs) {

  /*--- For a tetrahedron it is not possible to apply
        a tensor product and therefore all DOFs are
        simply stored. Allocate the memory. ---*/
  const unsigned short nD = (mPoly+1)*(mPoly+2)*(mPoly+3)/6;
  rDOFs.resize(nD);
  sDOFs.resize(nD);
  tDOFs.resize(nD);

  /*--- Determine the equidistant spacing along an edge. ---*/
  const passivedouble dh = 2.0/mPoly;

  /*--- Triple loop to compute the location of the grid DOFs. ---*/
  unsigned short ii = 0;
  for(unsigned short k=0; k<=mPoly; ++k) {
    const unsigned short uppBoundJ = mPoly - k;
    for(unsigned short j=0; j<=uppBoundJ; ++j) {
      const unsigned short uppBoundI = mPoly - k - j;
      for(unsigned short i=0; i<=uppBoundI; ++i, ++ii) {
        rDOFs[ii] = -1.0 + i*dh;
        sDOFs[ii] = -1.0 + j*dh;
        tDOFs[ii] = -1.0 + k*dh;
      }
    }
  }
}

void CFEMStandardTetBase::LocationTetGridDOFsLGL(const unsigned short  mPoly,
                                                 vector<passivedouble> &rDOFs,
                                                 vector<passivedouble> &sDOFs,
                                                 vector<passivedouble> &tDOFs) {

  /*--- The code to determine the parametric coordinates of the DOFs of the
        tetrahedron is a translation of the Matlab code belonging to the book
        Nodal Discontinuous Galerkin Methods, Algorithms, Analysis and Applications,
        written by Jan S. Hesthaven and Tim Warburton. ---*/

  /*--- Local parameters. ---*/
  const passivedouble tol   = 1.e-10;
  const passivedouble sqrt3 = sqrt(3.0);
  const passivedouble sqrt6 = sqrt(6.0);
  const passivedouble alphaOpt[] = {0.0000, 0.0000, 0.0000, 0.0000, 0.1002, 1.1332,
                                    1.5608, 1.3413, 1.2577, 1.1603, 1.10153,
                                    0.6080, 0.4523, 0.8856, 0.8717, 0.9655};

  /*--- Determine the number of DOFs based on the current value of mPoly. ---*/
  const unsigned short nD = (mPoly+1)*(mPoly+2)*(mPoly+3)/6;

  /*--- Create the equidistributed nodes. ---*/
  vector<passivedouble> r, s, t;
  LocationTetGridDOFsEquidistant(mPoly, r, s, t);

  /*--- Create the barycentric coordinates. ---*/
  vector<passivedouble> L1(nD), L2(nD), L3(nD), L4(nD);

  for(unsigned short i=0; i<nD; ++i) {
    L1[i] =  0.5*(1.0 + t[i]);
    L2[i] =  0.5*(1.0 + s[i]);
    L3[i] = -0.5*(1.0 + r[i] + s[i] + t[i]);
    L4[i] =  0.5*(1.0 + r[i]);
  }

  /*--- Set the coordinates of the vertices of the tetrahedron. ---*/
  const passivedouble v1[] = {-1.0, -1.0/sqrt3, -1.0/sqrt6};
  const passivedouble v2[] = { 1.0, -1.0/sqrt3, -1.0/sqrt6};
  const passivedouble v3[] = { 0.0,  2.0/sqrt3, -1.0/sqrt6};
  const passivedouble v4[] = { 0.0,  0.0,        3.0/sqrt6};

  /*--- Create for the four faces the orthogonal axis which
        are tangent to the face. ---*/
  passivedouble t1[4][3], t2[4][3];
  for(unsigned short i=0; i<3; ++i) {
    t1[0][i] = v2[i]-v1[i]; t2[0][i] = v3[i] - 0.5*(v1[i]+v2[i]);
    t1[1][i] = v2[i]-v1[i]; t2[1][i] = v4[i] - 0.5*(v1[i]+v2[i]);
    t1[2][i] = v3[i]-v2[i]; t2[2][i] = v4[i] - 0.5*(v2[i]+v3[i]);
    t1[3][i] = v3[i]-v1[i]; t2[3][i] = v4[i] - 0.5*(v1[i]+v3[i]);
  }

  /*--- Normalize the rows of t1 and t2. ---*/
  for(unsigned short i=0; i<4; ++i) {
    passivedouble norm;
    norm = sqrt(t1[i][0]*t1[i][0] + t1[i][1]*t1[i][1] + t1[i][2]*t1[i][2]);
    t1[i][0] /= norm; t1[i][1] /= norm; t1[i][2] /= norm;

    norm = sqrt(t2[i][0]*t2[i][0] + t2[i][1]*t2[i][1] + t2[i][2]*t2[i][2]);
    t2[i][0] /= norm; t2[i][1] /= norm; t2[i][2] /= norm;
  }

  /*--- Set the optimal value of alp for this case. ---*/
  const passivedouble alp = (mPoly < 16) ? alphaOpt[mPoly] : 1.0;

  /*--- Allocate the memory for the parametric coordinates. ---*/
  rDOFs.resize(nD);
  sDOFs.resize(nD);
  tDOFs.resize(nD);

  /*--- Initialize the parametric coordinates to the uniform distribution
        on the equilateral tetrahedron. ---*/
  for(unsigned short i=0; i<nD; ++i) {
    rDOFs[i] = L3[i]*v1[0] + L4[i]*v2[0] + L2[i]*v3[0] + L1[i]*v4[0];
    sDOFs[i] = L3[i]*v1[1] + L4[i]*v2[1] + L2[i]*v3[1] + L1[i]*v4[1];
    tDOFs[i] = L3[i]*v1[2] + L4[i]*v2[2] + L2[i]*v3[2] + L1[i]*v4[2];
  }

  /*--- Initialize the corrections, which are stored in r, s and t, to zero. ---*/
  r.assign(nD, 0.0);
  s.assign(nD, 0.0);
  t.assign(nD, 0.0);

  /*--- Loop over the four faces of the tetrahedron. ---*/
  for(unsigned short face=0; face<4; ++face) {

    /*--- In order to use the same approach for each face, the barycentric
          coordinates are different for the faces. Store these La, Lb,
          Lc and Ld. ---*/
    vector<passivedouble> La, Lb, Lc, Ld;
    if(face == 0)      {La = L1; Lb = L2; Lc = L3; Ld = L4;}
    else if(face == 1) {La = L2; Lb = L1; Lc = L3; Ld = L4;}
    else if(face == 2) {La = L3; Lb = L1; Lc = L4; Ld = L2;}
    else               {La = L4; Lb = L1; Lc = L3; Ld = L2;}

    /*--- Compute the warping factors tangential to this face. ---*/
    vector<passivedouble> warp1, warp2;
    EvalShift(mPoly, alp, Lb, Lc, Ld, warp1, warp2);

    /*--- Compute the volume blending. ---*/
    vector<passivedouble> blend(nD);
    for(unsigned short i=0; i<nD; ++i) {
      blend[i] = Lb[i]*Lc[i]*Ld[i];

      passivedouble denom = (Lb[i]+0.5*La[i])*(Lc[i]+0.5*La[i])*(Lc[i]+0.5*La[i]);
      denom = max(denom, tol);
      blend[i] *= (1.0 + alp*alp*La[i]*La[i])/denom;
    }

    /*--- Update the corrections. ---*/
    for(unsigned short i=0; i<nD; ++i) {
      const passivedouble abv1 = blend[i]*warp1[i], abv2 = blend[i]*warp2[i];
      r[i] += abv1*t1[face][0] + abv2*t2[face][0];
      s[i] += abv1*t1[face][1] + abv2*t2[face][1];
      t[i] += abv1*t1[face][2] + abv2*t2[face][2];
    }

    /*--- Fix the face warping. ---*/
    for(unsigned short i=0; i<nD; ++i) {
      unsigned short val = 0;
      if(Lb[i] > tol) ++val;
      if(Lc[i] > tol) ++val;
      if(Ld[i] > tol) ++val;
      if((La[i] < tol) && val < 3) {
        r[i] = warp1[i]*t1[face][0] + warp2[i]*t2[face][0];
        s[i] = warp1[i]*t1[face][1] + warp2[i]*t2[face][1];
        t[i] = warp1[i]*t1[face][2] + warp2[i]*t2[face][2];
      }
    }
  }

  /*--- Store the parametric coordinates of the equilateral tetrahedron in r,s,t. ---*/
  for(unsigned short i=0; i<nD; ++i) {
    r[i] += rDOFs[i];
    s[i] += sDOFs[i];
    t[i] += tDOFs[i];
  }

  /*--- Create the matrix needed to convert the coordinates to the standard tetrahedron. ---*/
  CSquareMatrixCM A(3);
  A(0,0) = v2[0]-v1[0]; A(1,0) = v2[1]-v1[1]; A(2,0) = v2[2]-v1[2];
  A(0,1) = v3[0]-v1[0]; A(1,1) = v3[1]-v1[1]; A(2,1) = v3[2]-v1[2];
  A(0,2) = v4[0]-v1[0]; A(1,2) = v4[1]-v1[1]; A(2,2) = v4[2]-v1[2];

  A.Invert();

  /*--- Convert the parametric coordinates from the equilateral tetrahedron to the
        parametric coordinates of the standard tetrahedron. ---*/
  for(unsigned short i=0; i<nD; ++i) {
    const passivedouble rr = 2.0*r[i] + v1[0] - v2[0] - v3[0] - v4[0];
    const passivedouble ss = 2.0*s[i] + v1[1] - v2[1] - v3[1] - v4[1];
    const passivedouble tt = 2.0*t[i] + v1[2] - v2[2] - v3[2] - v4[2];

    rDOFs[i] = A(0,0)*rr + A(0,1)*ss + A(0,2)*tt;
    sDOFs[i] = A(1,0)*rr + A(1,1)*ss + A(1,2)*tt;
    tDOFs[i] = A(2,0)*rr + A(2,1)*ss + A(2,2)*tt;
  }
}

void CFEMStandardTetBase::SubConnLinearElements(void) {

  /*--- The tetrahedron is split into several linear tetrahedron.
        Set the VTK sub-types accordingly. ---*/
  VTK_SubType1 = TETRAHEDRON;
  VTK_SubType2 = NONE;

  /*--- Initialize the number of DOFs for the current edges to the number of
        DOFs of the edges present in the tetrahedron. Also initialize the
        current k offset to zero.    ---*/
  unsigned short nDOFsCurrentEdges = nPoly + 1;
  unsigned short offCurrentK       = 0;

  /*--- Loop in the k-direction of the tetrahedron, which is along the edge
        from the first vertex to the last vertex of the tet. ---*/
  for(unsigned short k=0; k<nPoly; ++k) {

    /*--- Determine the offset for the next k. ---*/
    const unsigned short offNextK = offCurrentK
                                  + nDOFsCurrentEdges*(nDOFsCurrentEdges+1)/2;

    /*------------------------------------------------------------------------*/
    /*      Step 1: The tetrahedron at the end of the current i-edge.         */
    /*------------------------------------------------------------------------*/

    unsigned short n0 = offCurrentK + nDOFsCurrentEdges - 2;
    unsigned short n1 = n0 + 1;
    unsigned short n2 = n0 + nDOFsCurrentEdges;
    unsigned short n3 = offNextK + nDOFsCurrentEdges - 2;

    subConn1ForPlotting.push_back(n0);
    subConn1ForPlotting.push_back(n1);
    subConn1ForPlotting.push_back(n2);
    subConn1ForPlotting.push_back(n3);

    /*------------------------------------------------------------------------*/
    /* Step 2: The prisms that run from the end of the j-edge to the base of  */
    /*         tet just created. These prisms are subdivided into tetrahedra. */
    /*------------------------------------------------------------------------*/

    for(unsigned short i=0; i<(nDOFsCurrentEdges-2); ++i) {

      /*--- Determine the lowest j-index on the current i-line that contributes
            to the subprism. Convert that index to the local vertex number in
            the tetrahedron. ---*/
      unsigned short j = nDOFsCurrentEdges-2 -i;
      n0 = j*nDOFsCurrentEdges + i - j*(j-1)/2 + offCurrentK;

      /*--- Increment the j index to obtain the n1 vertex of the subprism. ---*/
      ++j;
      n1 = j*nDOFsCurrentEdges + i - j*(j-1)/2 + offCurrentK;

      /*--- The n2 vertex of the prism is located on the next k-level.
            The i-index remains the same, but the j index must be adapted, because
            the number of DOFs on the edges on the next k-level is one less. ---*/
      j  = nDOFsCurrentEdges-2 -i;
      n2 = j*(nDOFsCurrentEdges-1) + i - j*(j-1)/2 + offNextK;

      /*--- The n3 vertex is part of the upper triangle of the prism. Hence the
            i-index must be incremented by one, stored in ii. The j-index must
            be computed accordingly and converted to the 1D numbering. ---*/
      unsigned short ii = i+1;
      j  = nDOFsCurrentEdges-2 -ii;
      n3 = j*nDOFsCurrentEdges + ii - j*(j-1)/2 + offCurrentK;

      /*--- Increment the j index to obtain the n4 vertex of the subprism. ---*/
      ++j;
      unsigned short n4 = j*nDOFsCurrentEdges + ii - j*(j-1)/2 + offCurrentK;

      /*--- The n5 vertex of the prism is located on the next k-level.
            The i-index remains the same, but the j index must be adapted, because
            the number of DOFs on the edges on the next k-level is one less. ---*/
      j  = nDOFsCurrentEdges-2 -ii;
      unsigned short n5 = j*(nDOFsCurrentEdges-1) + ii - j*(j-1)/2 + offNextK;

      /*--- Divide the subprism into 3 subtetrahedra. ---*/
      subConn1ForPlotting.push_back(n0);
      subConn1ForPlotting.push_back(n1);
      subConn1ForPlotting.push_back(n2);
      subConn1ForPlotting.push_back(n4);

      subConn1ForPlotting.push_back(n3);
      subConn1ForPlotting.push_back(n0);
      subConn1ForPlotting.push_back(n2);
      subConn1ForPlotting.push_back(n4);

      subConn1ForPlotting.push_back(n3);
      subConn1ForPlotting.push_back(n2);
      subConn1ForPlotting.push_back(n5);
      subConn1ForPlotting.push_back(n4);
    }

    /*------------------------------------------------------------------------*/
    /* Step 3: The remaining subelements for this k-level, which are prisms   */
    /*         and hexas. These are subdivided into tetrahedra again.         */
    /*------------------------------------------------------------------------*/

    for(unsigned short i=0; i<(nDOFsCurrentEdges-2); ++i) {

      /*--- Define the variables n4 to n7, because these indices will be used
            again for the prism treated after the hexahedra. ---*/
      unsigned short n4, n5, n6, n7;

      /*--- Initialize n3, n2, n7 and n6 to the quad on the line j = 0. ---*/
      n3 = offCurrentK + i; n2 = n3 + 1;
      n7 = offNextK    + i; n6 = n7 + 1;

      /*--- Loop in the j-direction for the hexahedra on this i-row. ---*/
      for(unsigned short j=0; j<(nDOFsCurrentEdges-3-i); ++j) {

        /*--- Set the values of n0, n1, n4 and n5 from the previous hex. ---*/
        n0 = n3; n1 = n2; n4 = n7; n5 = n6;

        /*--- Nodes n3 and n7 are the j-neighbors of n0 and n4 respectively.
              Convert these (i,j,k) indices to the 1D index again. ---*/
        n3 = (j+1)*nDOFsCurrentEdges + i - j*(j+1)/2 + offCurrentK;
        n7 = (j+1)*(nDOFsCurrentEdges-1) + i - j*(j+1)/2 + offNextK;

        /*--- Nodes n2 and n6 are the i-neighbors of n3 and n7 respectively.
              Just add an offset of 1 in the 1D numbering. ---*/
        n2 = n3 + 1;
        n6 = n7 + 1;

        /*--- Divide the hexahedron in 6 tetrahedra and add their connectivity
              to the vector to store the subtetrahedra. ---*/
        subConn1ForPlotting.push_back(n0);
        subConn1ForPlotting.push_back(n3);
        subConn1ForPlotting.push_back(n4);
        subConn1ForPlotting.push_back(n1);

        subConn1ForPlotting.push_back(n3);
        subConn1ForPlotting.push_back(n7);
        subConn1ForPlotting.push_back(n4);
        subConn1ForPlotting.push_back(n1);

        subConn1ForPlotting.push_back(n4);
        subConn1ForPlotting.push_back(n7);
        subConn1ForPlotting.push_back(n5);
        subConn1ForPlotting.push_back(n1);

        subConn1ForPlotting.push_back(n1);
        subConn1ForPlotting.push_back(n2);
        subConn1ForPlotting.push_back(n3);
        subConn1ForPlotting.push_back(n7);

        subConn1ForPlotting.push_back(n2);
        subConn1ForPlotting.push_back(n5);
        subConn1ForPlotting.push_back(n6);
        subConn1ForPlotting.push_back(n7);

        subConn1ForPlotting.push_back(n1);
        subConn1ForPlotting.push_back(n5);
        subConn1ForPlotting.push_back(n2);
        subConn1ForPlotting.push_back(n7);
      }

      /*--- The prism that is in between the last hexahedron treated above
            and one of the prisms treated in step 2.
            The ID's of four of the vertices of the prism have already been
            computed for the last hexahedron above. However, they must be
            put in the correct numbering of the prism, which is done here. ---*/
      n0 = n3;
      n1 = n2;
      n3 = n7;
      n4 = n6;

      /*--- Determine the j-index of the two remaining vertices and convert
            the (i,j,k) indices to the 1D index. ---*/
      unsigned short j = nDOFsCurrentEdges-2-i;

      n2 = j*nDOFsCurrentEdges + i - j*(j-1)/2 + offCurrentK;
      n5 = j*(nDOFsCurrentEdges-1) + i - j*(j-1)/2 + offNextK;

      /*--- Divide the prism in 3 tetrahedra and add their connectivity
            to the vector to store the subtetrahedra.     ---*/
      subConn1ForPlotting.push_back(n3);
      subConn1ForPlotting.push_back(n5);
      subConn1ForPlotting.push_back(n4);
      subConn1ForPlotting.push_back(n1);

      subConn1ForPlotting.push_back(n0);
      subConn1ForPlotting.push_back(n2);
      subConn1ForPlotting.push_back(n3);
      subConn1ForPlotting.push_back(n1);

      subConn1ForPlotting.push_back(n2);
      subConn1ForPlotting.push_back(n5);
      subConn1ForPlotting.push_back(n3);
      subConn1ForPlotting.push_back(n1);
    }

    /*--- Set offCurrentK to offNextK for the next k-value and decrement
          the value of nDOFsCurrentEdges, also for the next k-value. ---*/
    offCurrentK = offNextK;
    --nDOFsCurrentEdges;
  }
}

/*----------------------------------------------------------------------------------*/
/*                Private member functions of CFEMStandardTetBase.                  */
/*----------------------------------------------------------------------------------*/

void CFEMStandardTetBase::EvalShift(const unsigned short        mPoly,
                                    const passivedouble         alpha,
                                    const vector<passivedouble> &L1,
                                    const vector<passivedouble> &L2,
                                    const vector<passivedouble> &L3,
                                    vector<passivedouble>       &dx,
                                    vector<passivedouble>       &dy) {

  /*--- The code for this function is a translation of the Matlab code belonging to
        the book Nodal Discontinuous Galerkin Methods, Algorithms, Analysis and
        Applications, written by Jan S. Hesthaven and Tim Warburton. ---*/

  /*--- Determine the number of DOFs based on the current value of mPoly. ---*/
  const unsigned short nD = (mPoly+1)*(mPoly+2)*(mPoly+3)/6;

  /*--- Compute for each edge the nodal blending functions. ---*/
  vector<passivedouble> blend1(nD), blend2(nD), blend3(nD);

  for(unsigned short i=0; i<nD; ++i) {
    blend1[i] = L2[i]*L3[i];
    blend2[i] = L1[i]*L3[i];
    blend3[i] = L1[i]*L2[i];
  }

  /*--- Compute the warp factors for the 3 edges. ---*/
  vector<passivedouble> warp1, warp2, warp3;
  vector<passivedouble> tmp(nD);

  for(unsigned short i=0; i<nD; ++i) tmp[i] = L3[i] - L2[i];
  EvalWarp(mPoly, tmp, warp1);

  for(unsigned short i=0; i<nD; ++i) tmp[i] = L1[i] - L3[i];
  EvalWarp(mPoly, tmp, warp2);

  for(unsigned short i=0; i<nD; ++i) tmp[i] = L2[i] - L1[i];
  EvalWarp(mPoly, tmp, warp3);

  /*--- Combine the blending and the warping. ---*/
  for(unsigned short i=0; i<nD; ++i) {
    warp1[i] *= blend1[i]*(1.0 + alpha*alpha*L1[i]*L1[i]);
    warp2[i] *= blend2[i]*(1.0 + alpha*alpha*L2[i]*L2[i]);
    warp3[i] *= blend3[i]*(1.0 + alpha*alpha*L3[i]*L3[i]);
  }

  /*--- Evaluate the shift in the equilateral triangle, which is the
        final result to be stored in dx and dy. ---*/
  dx.resize(nD);
  dy.resize(nD);

  const passivedouble cos1 = cos(SU2_TYPE::GetValue(TWO3*PI_NUMBER));
  const passivedouble cos2 = cos(SU2_TYPE::GetValue(FOUR3*PI_NUMBER));
  const passivedouble sin1 = sin(SU2_TYPE::GetValue(TWO3*PI_NUMBER));
  const passivedouble sin2 = sin(SU2_TYPE::GetValue(FOUR3*PI_NUMBER));

  for(unsigned short i=0; i<nD; ++i) {
    dx[i] = warp1[i] + cos1*warp2[i] + cos2*warp3[i];
    dy[i] =            sin1*warp2[i] + sin2*warp3[i];
  }
}

void CFEMStandardTetBase::EvalWarp(const unsigned short        mPoly,
                                   const vector<passivedouble> &xOut,
                                   vector<passivedouble>       &warp) {

  /*--- The code for this function is a translation of the Matlab code belonging to
        the book Nodal Discontinuous Galerkin Methods, Algorithms, Analysis and
        Applications, written by Jan S. Hesthaven and Tim Warburton. ---*/

  /*--- Determine the number of DOFs based on the current value of mPoly. ---*/
  const unsigned short nD = (mPoly+1)*(mPoly+2)*(mPoly+3)/6;

  /*--- Determine the 1D equidistributed and LGL nodes. Note that
        in this function the reversed definition is used. ---*/
  vector<passivedouble> xEq, xLGL;
  Location1DGridDOFsEquidistant(mPoly, xEq);
  Location1DGridDOFsLGL(mPoly, xLGL);

  reverse(xEq.begin(),  xEq.end());
  reverse(xLGL.begin(), xLGL.end());

  /*--- Initialize the warping vector to zero. ---*/
  warp.assign(nD, 0.0);

  /*--- Loop over the internal points of the edge. ---*/
  for(unsigned short i=1; i<mPoly; ++i) {

    /*--- Initialize the vector d to the distance between the
          equidistributed and LGL node i. ---*/
    vector<passivedouble> d(nD, xLGL[i]-xEq[i]);

    /*--- Loop over the other internal nodes of the edge to correct
          the entries of vector d. ---*/
    for(unsigned short j=1; j<mPoly; ++j) {
      if(j != i) {
        const passivedouble denom = 1.0/(xEq[i]-xEq[j]);
        for(unsigned short k=0; k<nD; ++k)
          d[k] *= (xOut[k]-xEq[j])*denom; 
      }
    }

    /*--- Carry out the scaling relative to the end points of the edge. ---*/
    passivedouble denom = -1.0/(xEq[i]-xEq[0]);
    for(unsigned short k=0; k<nD; ++k)
      d[k] *= denom;

    denom = 1.0/(xEq[i]-xEq[mPoly]);
    for(unsigned short k=0; k<nD; ++k)
      d[k] *= denom;

    /*--- Add d to warp for the DOFs. Note that the factor 4 is included here instead
          of in EvalShift as in the Matlab code of Hesthaven and Warburton. ---*/
    for(unsigned short k=0; k<nD; ++k)
      warp[k] += 4.0*d[k];
  }
}
