/*!
 * \file CFEMStandardTetGrid.cpp
 * \brief Functions for the class CFEMStandardTetGrid.
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

#include "../../include/fem/CFEMStandardTetGrid.hpp"
#include "../../include/toolboxes/CGeneralSquareMatrixCM.hpp"

/*----------------------------------------------------------------------------------*/
/*             Public member functions of CFEMStandardTetGrid.                      */
/*----------------------------------------------------------------------------------*/

CFEMStandardTetGrid::CFEMStandardTetGrid(const unsigned short val_nPoly,
                                         const unsigned short val_orderExact)
  : CFEMStandardTet(val_nPoly, val_orderExact) {

  /*--- Compute the values of the Lagrangian basis functions in the integration
        points for both the equidistant and LGL point distribution. ---*/
  LagBasisIntPointsTet(rTetDOFsEqui, sTetDOFsEqui, tTetDOFsEqui, lagBasisIntEqui);
  LagBasisIntPointsTet(rTetDOFsLGL,  sTetDOFsLGL,  tTetDOFsLGL,  lagBasisIntLGL);

  /*--- Compute the values of the derivatives of the Lagrangian basis functions in
        the integration points for both the equidistant and LGL point distribution. ---*/
  DerLagBasisIntPointsTet(rTetDOFsEqui, sTetDOFsEqui, tTetDOFsEqui, derLagBasisIntEqui);
  DerLagBasisIntPointsTet(rTetDOFsLGL,  sTetDOFsLGL,  tTetDOFsLGL,  derLagBasisIntLGL);

  /*--- Set up the jitted gemm call, if supported. For this particular standard
        element the derivative of the coordinates are computed, which is 3. ---*/
  SetUpJittedGEMM(nIntegrationPad, 3, nDOFs);
}

void CFEMStandardTetGrid::DerivativesCoorIntPoints(const bool                         LGLDistribution,
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
/*             Private member functions of CFEMStandardTetGrid.                     */
/*----------------------------------------------------------------------------------*/

void CFEMStandardTetGrid::DerLagBasisIntPointsTet(const vector<passivedouble>            &rDOFs,
                                                  const vector<passivedouble>            &sDOFs,
                                                  const vector<passivedouble>            &tDOFs,
                                                  vector<ColMajorMatrix<passivedouble> > &derLag) {

  /*--- Determine the padded number of the total number of integration points. ---*/
  const unsigned short nIntTot    = rTetInt.size();
  const unsigned short nIntTotPad = ((nIntTot+baseVectorLen-1)/baseVectorLen)*baseVectorLen;

  /*--- Determine the inverse of the Vandermonde matrix of the DOFs. ---*/
  CGeneralSquareMatrixCM VInv(rDOFs.size());
  VandermondeTetrahedron(rDOFs, sDOFs, tDOFs, VInv.GetMat());
  VInv.Invert();

  /*--- Determine the gradient of the Vandermonde matrix of the integration points. Make
        sure to allocate the number of rows to nIntTotPad and initialize them to zero. ---*/
  ColMajorMatrix<passivedouble> VDr(nIntTotPad,rDOFs.size()),
                                VDs(nIntTotPad,rDOFs.size()),
                                VDt(nIntTotPad,rDOFs.size());
  VDr.setConstant(0.0);
  VDs.setConstant(0.0);
  VDt.setConstant(0.0);

  GradVandermondeTetrahedron(rTetInt, sTetInt, tTetInt, VDr, VDs, VDt);

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

void CFEMStandardTetGrid::LagBasisIntPointsTet(const vector<passivedouble>   &rDOFs,
                                               const vector<passivedouble>   &sDOFs,
                                               const vector<passivedouble>   &tDOFs,
                                               ColMajorMatrix<passivedouble> &lag) {

  /*--- Determine the padded number of the total number of integration points. ---*/
  const unsigned short nIntTot    = rTetInt.size();
  const unsigned short nIntTotPad = ((nIntTot+baseVectorLen-1)/baseVectorLen)*baseVectorLen;

  /*--- Determine the inverse of the Vandermonde matrix of the DOFs. ---*/
  CGeneralSquareMatrixCM VInv(rDOFs.size());
  VandermondeTetrahedron(rDOFs, sDOFs, tDOFs, VInv.GetMat());
  VInv.Invert();

  /*--- Determine the Vandermonde matrix of the integration points. Make sure to
        allocate the number of rows to nIntTotPad and initialize them to zero. ---*/ 
  ColMajorMatrix<passivedouble> V(nIntTotPad,rDOFs.size());
  V.setConstant(0.0);
  VandermondeTetrahedron(rTetInt, sTetInt, tTetInt, V);

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

void CFEMStandardTetGrid::GradVandermondeTetrahedron(const vector<passivedouble>   &r,
                                                     const vector<passivedouble>   &s,
                                                     const vector<passivedouble>   &t,
                                                     ColMajorMatrix<passivedouble> &VDr,
                                                     ColMajorMatrix<passivedouble> &VDs,
                                                     ColMajorMatrix<passivedouble> &VDt) {

  /*--- For a tetrahedron the orthogonal basis for the reference element is obtained by a
        combination of Jacobi polynomials (of which the Legendre polynomials is a special
        case). This is the result of the orthonormalization of the monomial basis.
        Note that the sequence of the i, j and k loop must be identical to
        the evaluation of the Vandermonde matrix itself. ---*/
  unsigned short ii = 0;
  for(unsigned short i=0; i<=nPoly; ++i) {
    for(unsigned short j=0; j<=(nPoly-i); ++j) {
      for(unsigned short k=0; k<=(nPoly-i-j); ++k, ++ii) {
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

          /*--- Determine the value of the three 1D contributions to the 3D basis functions as
                well as the gradients of these basis functions w.r.t. to their arguments. ---*/
          const passivedouble fa  = NormJacobi(i,0,    0,a);
          const passivedouble gb  = NormJacobi(j,2*i+1,0,b);
          const passivedouble hc  = NormJacobi(k,2*(i+j+1),0,c);
          const passivedouble dfa = GradNormJacobi(i,0,    0,a);
          const passivedouble dgb = GradNormJacobi(j,2*i+1,0,b);
          const passivedouble dhc = GradNormJacobi(k,2*(i+j+1),0,c);

          /*--- Compute the derivative of the basis function w.r.t. r. As r is only present in
                the parameter a the derivative of the basis function w.r.t. a is multiplied by
                dadr. Note that the implementation is such that all possible singularities are
                divided out of the expression. ---*/
          VDr(l,ii) = sqrt(8.0)*dfa*gb*hc;
          if(i > 0) {
            VDr(l,ii) *= 4.0;
            if(i > 1 ) VDr(l,ii) *= pow((1.0-b), (i-1));
          }

          if(i+j > 1) VDr(l,ii) *= pow((1.0-c), (i+j-1));

          /*--- Compute the derivative of the basis function w.r.t. s. As s is present in both
                the parameters a and b, both variables must be taken into account when the
                derivative is computed. Note that the implementation is such that all possible
                singularities are divided out of the expression. The first part is the derivative
                of the basis function w.r.t. b multiplied by dbds. This value is stored, because
                it is needed later on to compute the derivative w.r.t. t. ---*/
          VDs(l,ii) = dgb;
          if( i ) VDs(l,ii) *= pow((1.0-b), i);

          if(i > 0) {
            tmp = i*gb;
            if(i > 1) tmp *= pow((1.0-b), (i-1));
            VDs(l,ii) -= tmp;
          }

          if(i+j > 0) {
            VDs(l,ii) *= 2.0*sqrt(8.0)*fa*hc;
            if(i+j > 1) VDs(l,ii) *= pow((1.0-c), (i+j-1));
          }

          const passivedouble dPsidbXdbds = VDs(l,ii);

          /*--- Add the contribution from the derivative of the basis function
                w.r.t. a multiplied by dads. ---*/
          VDs(l,ii) += 0.5*(a+1.0)*VDr(l,ii);

          /*--- Compute the derivative of the basis function w.r.t. t. As t is present in a, b and c,
                all parameters must be taken into account when the derivative is computed. Note that
                the implementation is such that all possible singularities are divided out of the
                expression. The first part is the derivative of the basis function w.r.t. c,
                which is equal to t. ---*/
          VDt(l,ii) = dhc;
          if(i+j > 0) {
            VDt(l,ii) *= pow((1.0-c), (i+j));

            tmp = (i+j)*hc;
            if(i+j > 1) tmp *= pow((1.0-c), (i+j-1));

            VDt(l,ii) -= tmp;
          }

          VDt(l,ii) *= sqrt(8.0)*fa*gb;
          if( i ) VDt(l,ii) *= pow((1.0-b), i);

          /*--- Add the contribution from the derivative of the basis function w.r.t. a multiplied
                by dadt and the derivative w.r.t. b multiplied by dbdt. ---*/
          VDt(l,ii) += 0.5*(a+1.0)*VDr(l,ii) + 0.5*(b+1.0)*dPsidbXdbds;
        }
      }
    }
  }
}

void CFEMStandardTetGrid::VandermondeTetrahedron(const vector<passivedouble>   &r,
                                                 const vector<passivedouble>   &s,
                                                 const vector<passivedouble>   &t,
                                                 ColMajorMatrix<passivedouble> &V) {

  /*--- For a tetrahedron the orthogonal basis for the reference element is obtained by a
        combination of Jacobi polynomials (of which the Legendre polynomials is a special
        case). This is the result of the orthonormalization of the monomial basis. ---*/
  unsigned short ii = 0;
  for(unsigned short i=0; i<=nPoly; ++i) {
    for(unsigned short j=0; j<=(nPoly-i); ++j) {
      for(unsigned short k=0; k<=(nPoly-i-j); ++k, ++ii) {
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
