/*!
 * \file CFEMStandardLineBase.cpp
 * \brief Functions for the class CFEMStandardLineBase.
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

#include "../../include/fem/CFEMStandardLineBase.hpp"
#include "../../include/fem/fem_gauss_jacobi_quadrature.hpp"
#include "../../include/toolboxes/CSquareMatrixCM.hpp"

/*----------------------------------------------------------------------------------*/
/*              Public member functions of CFEMStandardLineBase.                    */
/*----------------------------------------------------------------------------------*/

CFEMStandardLineBase::CFEMStandardLineBase(const unsigned short val_nPoly,
                                           const unsigned short val_orderExact) {

  /*--- Store the command line arguments. ---*/
  VTK_Type   = LINE;
  nPoly      = val_nPoly;
  orderExact = val_orderExact;

  /*--- Determine the total number of DOFs and its padded version. ---*/
  nDOFs    = (nPoly+1);
  nDOFsPad = PaddedValue(nDOFs);

  /*--- Determine the total number of integration points
        and its padded version. ---*/
  nIntegration    = orderExact/2 + 1;
  nIntegrationPad = PaddedValue(nIntegration);

  /*--- Determine the location and the weights of the 1D integration points. ---*/
  rLineInt.resize(nIntegration);
  wLineInt.resize(nIntegration);

  CGaussJacobiQuadrature GaussJacobi;
  GaussJacobi.GetQuadraturePoints(0.0, 0.0, -1.0, 1.0, rLineInt, wLineInt);

  /*--- Allocate the memory for the padded number of integration points and
        initialize the weights to zero. This is done such that the padded
        values are initialized appropriately. ---*/
  wIntegration.resize(nIntegrationPad);
  wIntegration.setConstant(0.0);

  /*--- Copy the values from wLineInt to wIntegration. ---*/
  for(unsigned short i=0; i<nIntegration; ++i)
    wIntegration(i) = wLineInt[i];
}

void CFEMStandardLineBase::DerLagBasisIntPointsLine(const vector<passivedouble>   &rDOFs,
                                                    const vector<passivedouble>   &rInt,
                                                    const bool                    usePadding,
                                                    ColMajorMatrix<passivedouble> &derLag) {

  /*--- Determine the number of integration points along the line
        and its padded value, if necessary. ---*/
  const unsigned short nIntLine    = rInt.size();
  const unsigned short nIntPadLine = usePadding ? PaddedValue(nIntLine) : nIntLine;

  /*--- Determine the inverse of the Vandermonde matrix of rDOFs. ---*/
  CSquareMatrixCM VInv(rDOFs.size());
  Vandermonde1D(rDOFs.size()-1, rDOFs, VInv.GetMat());
  VInv.Invert();

  /*--- Determine the gradient of the Vandermonde matrix of rInt. Make sure to
        allocate the number of rows to nIntPadLine and initialize them to zero. ---*/
  ColMajorMatrix<passivedouble> gradV(nIntPadLine,rDOFs.size());
  gradV.setConstant(0.0);
  GradVandermonde1D(rDOFs.size()-1, rInt, gradV);

  /*--- The derivatives of the Lagrangian basis functions can be obtained
        by multiplying gradV and VInv. ---*/
  VInv.MatMatMult('R', gradV, derLag);

  /*--- Check if the sum of the elements of the relevant rows of derLag is 0. ---*/
  CheckRowSum(nIntLine, rDOFs.size(), 0.0, derLag);
}

void CFEMStandardLineBase::HesLagBasisIntPointsLine(const vector<passivedouble>   &rDOFs,
                                                    const vector<passivedouble>   &rInt,
                                                    const bool                    usePadding,
                                                    ColMajorMatrix<passivedouble> &hesLag) {

  /*--- Determine the number of integration points along the line
        and its padded value, if necessary. ---*/
  const unsigned short nIntLine    = rInt.size();
  const unsigned short nIntPadLine = usePadding ? PaddedValue(nIntLine) : nIntLine;

  /*--- Determine the inverse of the Vandermonde matrix of rDOFs. ---*/
  CSquareMatrixCM VInv(rDOFs.size());
  Vandermonde1D(rDOFs.size()-1, rDOFs, VInv.GetMat());
  VInv.Invert();

  /*--- Determine the Hessian of the Vandermonde matrix of rInt. Make sure to
        allocate the number of rows to nIntPadLine and initialize them to zero. ---*/
  ColMajorMatrix<passivedouble> hesV(nIntPadLine,rDOFs.size());
  hesV.setConstant(0.0);
  HesVandermonde1D(rDOFs.size()-1, rInt, hesV);

  /*--- The Hessian of the Lagrangian basis functions can be obtained
        by multiplying hesV and VInv. ---*/
  VInv.MatMatMult('R', hesV, hesLag);

  /*--- Check if the sum of the elements of the relevant rows of hesLag is 0. ---*/
  CheckRowSum(nIntLine, rDOFs.size(), 0.0, hesLag);
}

void CFEMStandardLineBase::LagBasisIntPointsLine(const vector<passivedouble>   &rDOFs,
                                                 const vector<passivedouble>   &rInt,
                                                 const bool                    usePadding,
                                                 ColMajorMatrix<passivedouble> &lag) {

  /*--- Determine the number of integration points along the line
        and its padded value, if necessary. ---*/
  const unsigned short nIntLine    = rInt.size();
  const unsigned short nIntPadLine = usePadding ? PaddedValue(nIntLine) : nIntLine;

  /*--- Determine the inverse of the Vandermonde matrix of rDOFs. ---*/
  CSquareMatrixCM VInv(rDOFs.size());
  Vandermonde1D(rDOFs.size()-1, rDOFs, VInv.GetMat());
  VInv.Invert();

  /*--- Determine the Vandermonde matrix of rInt. Make sure to allocate
        the number of rows to nIntPadLine and initialize them to zero. ---*/
  ColMajorMatrix<passivedouble> V(nIntPadLine,rDOFs.size());
  V.setConstant(0.0);
  Vandermonde1D(rDOFs.size()-1, rInt, V);

  /*--- The Lagrangian basis functions can be obtained by multiplying
        V and VInv. ---*/
  VInv.MatMatMult('R', V, lag);

  /*--- Check if the sum of the elements of the relevant rows of lag is 1. ---*/
  CheckRowSum(nIntLine, rDOFs.size(), 1.0, lag);
}

/*-----------------------------------------------------------------------------------*/
/*---                         Protected member functions.                         ---*/
/*-----------------------------------------------------------------------------------*/

void CFEMStandardLineBase::Location1DGridDOFsEquidistant(const unsigned short  mPoly,
                                                         vector<passivedouble> &r) {

  /*--- Allocate the memory and set the location of the DOFs using
        equidistant spacing. ---*/
  r.resize(mPoly+1);
  const passivedouble dh = 2.0/mPoly;

  for(unsigned short i=0; i<=mPoly; ++i)
    r[i] = -1.0 + i*dh;
}

void CFEMStandardLineBase::Location1DGridDOFsLGL(const unsigned short  mPoly,
                                                 vector<passivedouble> &r) {

  /*--- Allocate the memory. ---*/
  const unsigned short nPoints = mPoly+1;
  r.resize(nPoints);

  /*--- The distribution of points is symmetric. Hence only half the number
        of points must be computed. The first and last are at the end of
        the interval and for and add number the mid point is zero. ---*/
  const unsigned short nn = nPoints/2;

  r[0]     = -1.0;
  r[mPoly] =  1.0;

  if(2*nn < nPoints) r[nn] = 0.0;

  /*--- Constants used in the initial guess of the roots in the loop below. ---*/
  const passivedouble t1 = 1.0 - 3.0*(mPoly-1)/(8.0*mPoly*mPoly*mPoly);
  const passivedouble t2 = SU2_TYPE::GetValue(PI_NUMBER)/(4.0*mPoly + 1.0);

  /*--- The remaing points must be computed. These are the roots of P'_{n-1}(x),
        P_n is the classic Legendre polynomial of order n. Loop over roots to
        be computed. ---*/
  unsigned short ii = mPoly-1;
  for(unsigned short i=1; i<nn; ++i, --ii) {

    /*--- Initial guess of this root. ---*/
    passivedouble x = t1*cos(t2*(4*i+1));

    /*--- Determine the Legendre Polynomials P_{n-2} and P_{n-1} and
          the value f = P'_{n-1}(x). ---*/
    passivedouble Pnm1 = Legendre(mPoly,   x);
    passivedouble Pnm2 = Legendre(mPoly-1, x);
    passivedouble f    = mPoly*(Pnm2 - x*Pnm1)/(1.0-x*x);

    /*--- Solve the root using Halley's method.
          Loop until machine precision has been reached. ---*/
    for(;;) {

      /*--- Determine the value of the first and second derivative of f. ---*/
      const passivedouble df  = (2.0*x*f - nPoints*mPoly*Pnm1)/(1.0-x*x);
      const passivedouble d2f = (2.0*x*df - (nPoints*mPoly-2)*f)/(1.0-x*x);

      /*--- Compute the new value of the root. ---*/
      x = x - 2.0*f*df/(2.0*df*df - f*d2f);

      /*--- Determine the new value of the Legendre polynomials and
            compute the new value of f. Store the old value. ---*/
      const passivedouble fOld = f;

      Pnm1 = Legendre(mPoly,   x);
      Pnm2 = Legendre(mPoly-1, x);
      f    = mPoly*(Pnm2 - x*Pnm1)/(1.0-x*x);

      /*--- Convergence criterion. ---*/
      if(fabs(fOld) <= fabs(f)) break; 
    }

    /*--- Store the value as well as the symmetric equivalent. ---*/
    r[ii] =  x;
    r[i]  = -x;
  }
}

void CFEMStandardLineBase::GradVandermonde1D(const unsigned short          mPoly,
                                             const vector<passivedouble>   &r,
                                             ColMajorMatrix<passivedouble> &VDr) {

  /*--- Compute the gradient of the 1D Vandermonde matrix. ---*/
  for(unsigned short j=0; j<=mPoly; ++j)
    for(unsigned short i=0; i<r.size(); ++i)
      VDr(i,j) = GradNormJacobi(j, 0, 0, r[i]);
}

void CFEMStandardLineBase::HesVandermonde1D(const unsigned short          mPoly,
                                            const vector<passivedouble>   &r,
                                            ColMajorMatrix<passivedouble> &VD2r) {

  /*--- Compute the Hessian of the 1D Vandermonde matrix. ---*/
  for(unsigned short j=0; j<=mPoly; ++j)
    for(unsigned short i=0; i<r.size(); ++i)
      VD2r(i,j) = HesNormJacobi(j, 0, 0, r[i]);
}

void CFEMStandardLineBase::Vandermonde1D(const unsigned short          mPoly,
                                         const vector<passivedouble>   &r,
                                         ColMajorMatrix<passivedouble> &V) {

  /*--- Compute the 1D Vandermonde matrix. ---*/
  for(unsigned short j=0; j<=mPoly; ++j)
    for(unsigned short i=0; i<r.size(); ++i)
      V(i,j) = NormJacobi(j, 0, 0, r[i]);
}

void CFEMStandardLineBase::SubConnLinearElements(void) {

  /*--- The line is split into several linear lines.
        Set the VTK sub-types accordingly. ---*/
  VTK_SubType1 = LINE;
  VTK_SubType2 = NONE;

  /*--- Determine the local subconnectivity of the line element used for plotting
        purposes. This is rather trivial, because the line element is subdivided
        into nPoly linear line elements. ---*/
  unsigned short nnPoly = max(nPoly,(unsigned short) 1);
  for(unsigned short i=0; i<nnPoly; ++i) {
    subConn1ForPlotting.push_back(i);
    subConn1ForPlotting.push_back(i+1);
  }
}

void CFEMStandardLineBase::SubConnLinearElementsFace(int val_faceID_Elem) {

  /*--- This is an added functionality specifically for surface output.
        The utility is similar to SubConnLinearElements, but store the 
        nodes (available in gridConnFaces) using the face ID w.r.t the 
        volume when this base class is considered a face. ---*/
        
  /*--- The line is split into several linear lines.
        Set the VTK sub-types accordingly. ---*/
  VTK_SubType1 = LINE;
  VTK_SubType2 = NONE;

  /*--- Determine the local subconnectivity of the line element used for plotting
        purposes. This is rather trivial, because the line element is subdivided
        into nPoly linear line elements. ---*/
  unsigned short nnPoly = max(nPoly,(unsigned short) 1);
  for(unsigned short i=0; i<nnPoly; ++i) {
    subConn1ForPlotting.push_back(gridConnFaces[val_faceID_Elem][i]);
    subConn1ForPlotting.push_back(gridConnFaces[val_faceID_Elem][i+1]);
  }
}
/*-----------------------------------------------------------------------------------*/
/*---                          Private member functions.                          ---*/
/*-----------------------------------------------------------------------------------*/

passivedouble CFEMStandardLineBase::Legendre(unsigned short n,
                                             passivedouble  x) {

  /*--- Initialization of the polynomials Pnm1 and Pn. ---*/
  passivedouble Pnm1 = 1.0;
  passivedouble Pn   = x;

  /*--- Take care of the special situation of n == 0. ---*/
  if(n == 0) Pn = Pnm1;
  else {

    /*--- Recursive definition of Pn. ---*/
    for(unsigned short i=2; i<=n; ++i)
    {
      const passivedouble tmp = Pnm1;
      Pnm1 = Pn;

      Pn = ((2*i-1)*x*Pn - (i-1)*tmp)/i;
    }
  }

  /*--- Return Pn. ---*/
  return Pn;
}
