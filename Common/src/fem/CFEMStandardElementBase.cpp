/*!
 * \file CFEMStandardElementBase.cpp
 * \brief Functions for the base class CFEMStandardElementBase.
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

#include "../../include/fem/CFEMStandardElementBase.hpp"
#include "../../include/toolboxes/CGeneralSquareMatrixCM.hpp"
#include "../../include/omp_structure.hpp"

/*----------------------------------------------------------------------------------*/
/*          Constructor and destructor of CFEMStandardElementBase.                  */
/*----------------------------------------------------------------------------------*/

CFEMStandardElementBase::CFEMStandardElementBase() {

#if defined(PRIMAL_SOLVER) && defined(HAVE_MKL)
  jitter = nullptr;
#endif
}

CFEMStandardElementBase::~CFEMStandardElementBase() {

#if defined(PRIMAL_SOLVER) && defined(HAVE_MKL)
  if( jitter ) {
    mkl_jit_destroy(jitter);
    jitter = nullptr;
  }
#endif
}

/*----------------------------------------------------------------------------------*/
/*          Public member functions of CFEMStandardElementBase.                     */
/*----------------------------------------------------------------------------------*/

unsigned short CFEMStandardElementBase::GetNDOFsStatic(unsigned short VTK_Type,
                                                       unsigned short nPoly) {
  unsigned short nDOFsEdge = nPoly + 1;
  unsigned short nDOFs;

  switch(VTK_Type) {

    case LINE:
      nDOFs = nDOFsEdge;
      break;

    case TRIANGLE:
      nDOFs = nDOFsEdge*(nDOFsEdge+1)/2;
      break;

    case QUADRILATERAL:
      nDOFs = nDOFsEdge*nDOFsEdge;
      break;

    case TETRAHEDRON:
      nDOFs = nDOFsEdge*(nDOFsEdge+1)*(nDOFsEdge+2)/6;
      break;

    case HEXAHEDRON:
      nDOFs = nDOFsEdge*nDOFsEdge*nDOFsEdge;
      break;

    case PRISM:
      nDOFs = nDOFsEdge*nDOFsEdge*(nDOFsEdge+1)/2;
      break;

    case PYRAMID:
      nDOFs = nDOFsEdge*(nDOFsEdge+1)*(2*nDOFsEdge+1)/6;
      break;

    default:
      nDOFs = 0;  // Indicates an invalid element.
      break;
  }

  return nDOFs;
}

void CFEMStandardElementBase::MetricTermsVolumeIntPoints(const bool                         LGLDistribution,
                                                         ColMajorMatrix<su2double>          &matCoor,
                                                         vector<ColMajorMatrix<su2double> > &matMetricTerms,
                                                         su2activevector                    &Jacobians) {

  /*--- Compute the derivatives of the coordinates in the volume integration points.
        The matrix matMetricTerms is used as storage for this data. ---*/
  DerivativesCoorVolumeIntPoints(LGLDistribution, matCoor, matMetricTerms);

  /*--- Check the dimension of the problem. ---*/
  if(matMetricTerms.size() == 2) {

    /*--- 2D problem. Set the references to the derivatives of the Cartesian
          coordinates w.r.t. the parametric coordinates, which are currently
          stored in matMetricTerms. ---*/
    ColMajorMatrix<su2double> &dCoorDr = matMetricTerms[0];
    ColMajorMatrix<su2double> &dCoorDs = matMetricTerms[1];

    /*--- Set the references to store the derivatives of the parametric coordinates
          w.r.t. the Cartesian coordinates. ---*/
    ColMajorMatrix<su2double> &dParDx = matMetricTerms[0];
    ColMajorMatrix<su2double> &dParDy = matMetricTerms[1];

    /*--- Loop over the padded number of integration points
          to compute the 2D metric terms. ---*/
    SU2_OMP_SIMD
    for(unsigned short i=0; i<nIntegrationPad; ++i) {

      /*--- Easier storage of the derivatives of the Cartesian coordinates
            w.r.t. the parametric coordinates. ---*/
      const su2double dxdr = dCoorDr(i,0);
      const su2double dydr = dCoorDr(i,1);

      const su2double dxds = dCoorDs(i,0);
      const su2double dyds = dCoorDs(i,1);

      /*--- Compute the Jacobian. ---*/
      Jacobians(i) = dxdr*dyds - dxds*dydr;

      /*--- Compute the derivatives of the parametric coordinates w.r.t.
            to the Cartesian coordinates. They are multiplied by the Jacobian
            to avoid a division by this Jacobian at this stage. ---*/
      dParDx(i,0) =  dyds;   // J drdx
      dParDx(i,1) = -dydr;   // J dsdx

      dParDy(i,0) = -dxds;   // J drdy
      dParDy(i,1) =  dxdr;   // J dsdy
    }
  }
  else {

    /*--- 3D problem. Set the references to the derivatives of the Cartesian
          coordinates w.r.t. the parametric coordinates, which are currently
          stored in matMetricTerms. ---*/
    ColMajorMatrix<su2double> &dCoorDr = matMetricTerms[0];
    ColMajorMatrix<su2double> &dCoorDs = matMetricTerms[1];
    ColMajorMatrix<su2double> &dCoorDt = matMetricTerms[2];

    /*--- Set the references to store the derivatives of the parametric coordinates
          w.r.t. the Cartesian coordinates. ---*/
    ColMajorMatrix<su2double> &dParDx = matMetricTerms[0];
    ColMajorMatrix<su2double> &dParDy = matMetricTerms[1];
    ColMajorMatrix<su2double> &dParDz = matMetricTerms[2];

    /*--- Loop over the padded number of integration points
          to compute the 3D metric terms. ---*/
    SU2_OMP_SIMD
    for(unsigned short i=0; i<nIntegrationPad; ++i) {

      /*--- Easier storage of the derivatives of the Cartesian coordinates
            w.r.t. the parametric coordinates. ---*/
      const su2double dxdr = dCoorDr(i,0);
      const su2double dydr = dCoorDr(i,1);
      const su2double dzdr = dCoorDr(i,2);

      const su2double dxds = dCoorDs(i,0);
      const su2double dyds = dCoorDs(i,1);
      const su2double dzds = dCoorDs(i,2);

      const su2double dxdt = dCoorDt(i,0);
      const su2double dydt = dCoorDt(i,1);
      const su2double dzdt = dCoorDt(i,2);

      /*--- Compute the Jacobian. ---*/
      Jacobians(i) = dxdr*(dyds*dzdt - dzds*dydt)
                   - dxds*(dydr*dzdt - dzdr*dydt)
                   + dxdt*(dydr*dzds - dzdr*dyds);

      /*--- Compute the derivatives of the parametric coordinates w.r.t.
            to the Cartesian coordinates. They are multiplied by the Jacobian
            to avoid a division by this Jacobian at this stage. ---*/
      dParDx(i,0) = dyds*dzdt - dzds*dydt;  // J drdx
      dParDx(i,1) = dzdr*dydt - dydr*dzdt;  // J dsdx
      dParDx(i,2) = dydr*dzds - dzdr*dyds;  // J dtdx

      dParDy(i,0) = dzds*dxdt - dxds*dzdt;  // J drdy
      dParDy(i,1) = dxdr*dzdt - dzdr*dxdt;  // J dsdy
      dParDy(i,2) = dzdr*dxds - dxdr*dzds;  // J dtdy

      dParDz(i,0) = dxds*dydt - dyds*dxdt;  // J drdz
      dParDz(i,1) = dydr*dxdt - dxdr*dydt;  // J dsdz
      dParDz(i,2) = dxdr*dyds - dydr*dxds;  // J dtdz
    }
  }
}

void CFEMStandardElementBase::MinMaxJacobians(const bool                         LGLDistribution,
                                              ColMajorMatrix<su2double>          &matCoor,
                                              vector<ColMajorMatrix<su2double> > &matMetricTerms,
                                              su2activevector                    &Jacobians,
                                              su2double                          &jacMin,
                                              su2double                          &jacMax) {

  /*--- Compute the metric terms in the volume integration points. ---*/
  MetricTermsVolumeIntPoints(LGLDistribution, matCoor, matMetricTerms, Jacobians);

  /*--- Loop over the integration points and determine the minimum
        and maximum value of the Jacobian. ---*/
  jacMin = jacMax = Jacobians(0);
  for(unsigned short i=0; i<nIntegration; ++i) {
    jacMin = min(jacMin, Jacobians(i));
    jacMax = max(jacMax, Jacobians(i));
  }
}

/*----------------------------------------------------------------------------------*/
/*          Protected member functions of CFEMStandardElementBase.                  */
/*----------------------------------------------------------------------------------*/

void CFEMStandardElementBase::Location1DGridDOFsEquidistant(vector<passivedouble> &r) {

  /*--- Allocate the memory and set the location of the DOFs using
        equidistant spacing. ---*/
  r.resize(nPoly+1);
  const passivedouble dh = 2.0/nPoly;

  for(unsigned short i=0; i<=nPoly; ++i)
    r[i] = -1.0 + i*dh;
}

void CFEMStandardElementBase::Location1DGridDOFsLGL(vector<passivedouble> &r) {

  /*--- Allocate the memory. ---*/
  const unsigned short nPoints = nPoly+1;
  r.resize(nPoints);

  /*--- The distribution of points is symmetric. Hence only half the number
        of points must be computed. The first and last are at the end of
        the interval and for and add number the mid point is zero. ---*/
  const unsigned short nn = nPoints/2;

  r[0]     = -1.0;
  r[nPoly] =  1.0;

  if(2*nn < nPoints) r[nn] = 0.0;

  /*--- Constants used in the initial guess of the roots in the loop below. ---*/
  const passivedouble t1 = 1.0 - 3.0*(nPoly-1)/(8.0*nPoly*nPoly*nPoly);
  const passivedouble t2 = SU2_TYPE::GetValue(PI_NUMBER)/(4.0*nPoly + 1.0);

  /*--- The remaing points must be computed. These are the roots of P'_{n-1}(x),
        P_n is the classic Legendre polynomial of order n. Loop over roots to
        be computed. ---*/
  unsigned short ii = nPoly-1;
  for(unsigned short i=1; i<nn; ++i, --ii) {

    /*--- Initial guess of this root. ---*/
    passivedouble x = t1*cos(t2*(4*i+1));

    /*--- Determine the Legendre Polynomials P_{n-2} and P_{n-1} and
          the value f = P'_{n-1}(x). ---*/
    passivedouble Pnm1 = Legendre(nPoly,   x);
    passivedouble Pnm2 = Legendre(nPoly-1, x);
    passivedouble f    = nPoly*(Pnm2 - x*Pnm1)/(1.0-x*x);

    /*--- Solve the root using Halley's method.
          Loop until machine precision has been reached. ---*/
    for(;;) {

      /*--- Determine the value of the first and second derivative of f. ---*/
      const passivedouble df  = (2.0*x*f - nPoints*nPoly*Pnm1)/(1.0-x*x);
      const passivedouble d2f = (2.0*x*df - (nPoints*nPoly-2)*f)/(1.0-x*x);

      /*--- Compute the new value of the root. ---*/
      x = x - 2.0*f*df/(2.0*df*df - f*d2f);

      /*--- Determine the new value of the Legendre polynomials and
            compute the new value of f. Store the old value. ---*/
      const passivedouble fOld = f;

      Pnm1 = Legendre(nPoly,   x);
      Pnm2 = Legendre(nPoly-1, x);
      f    = nPoly*(Pnm2 - x*Pnm1)/(1.0-x*x);

      /*--- Convergence criterion. ---*/
      if(fabs(fOld) <= fabs(f)) break; 
    }

    /*--- Store the value as well as the symmetric equivalent. ---*/
    r[ii] =  x;
    r[i]  = -x;
  }
}

void CFEMStandardElementBase::LocationTriangleGridDOFsEquidistant(vector<passivedouble> &r,
                                                                  vector<passivedouble> &s) {

  /*--- Determine the number of DOFs of the triangle. As this function is also
        called for prisms, the member variable nDOFs cannot be used for this.
        Allocate the memory for r and s afterwards. ---*/
  const unsigned short nD = (nPoly+1)*(nPoly+2)/2;
  r.resize(nD);
  s.resize(nD);

  /*--- Determine the equidistant spacing in parametric space. ---*/
  const passivedouble dh = 2.0/nPoly;

  /*--- Double loop to determine the parametric coordinates. ---*/
  unsigned short ii = 0;
  for(unsigned short j=0; j<=nPoly; ++j) {
    const unsigned short uppBoundI = nPoly - j;
    for(unsigned short i=0; i<=uppBoundI; ++i, ++ii) {
      r[ii] = -1.0 + i*dh;
      s[ii] = -1.0 + j*dh;
    }
  }
}

void CFEMStandardElementBase::LocationTriangleGridDOFsLGL(vector<passivedouble> &r,
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
  const unsigned short nD = (nPoly+1)*(nPoly+2)/2;

  /*--- Allocate the memory for the help variables, as well as r and s. ---*/
  vector<passivedouble> L1(nD), L2(nD), L3(nD), rout(nD), warpf1(nD), warpf2(nD), warpf3(nD);

  r.resize(nD);
  s.resize(nD);

  /*--- Determine the uniform spacing. ---*/
  /*--- Determine the equidistant spacing in parametric space. ---*/
  const passivedouble dh = 1.0/nPoly;

  /*-----------------------------------------------------------------------------*/
  /*--- Step 1: Create the equidistributed nodes on the equilateral triangle. ---*/
  /*-----------------------------------------------------------------------------*/

  unsigned short ii = 0;
  for(unsigned short j=0; j<=nPoly; ++j) {
    const unsigned short uppBoundI = nPoly - j;
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
  const passivedouble alp = (nPoly < 16) ? alphaOpt[nPoly] : 5.0/3.0;

  /*--- Compute the warp factors for the DOFs for each edge. ---*/
  for(unsigned short i=0; i<nD; ++i) rout[i] = L3[i] - L2[i];
  WarpFactor(rout, warpf1);

  for(unsigned short i=0; i<nD; ++i) rout[i] = L1[i] - L3[i];
  WarpFactor(rout, warpf2);

  for(unsigned short i=0; i<nD; ++i) rout[i] = L2[i] - L1[i];
  WarpFactor(rout, warpf3);

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

void CFEMStandardElementBase::DerLagBasisIntPointsLine(const vector<passivedouble>   &rDOFs,
                                                       const vector<passivedouble>   &rInt,
                                                       ColMajorMatrix<passivedouble> &derLag) {

  /*--- Determine the number of integration points along the line
        and its padded value. ---*/
  const unsigned short nIntLine    = rInt.size();
  const unsigned short nIntPadLine = ((nIntLine+vecLen-1)/vecLen)*vecLen;

  /*--- Determine the inverse of the Vandermonde matrix of rDOFs. ---*/
  CGeneralSquareMatrixCM VInv(nPoly+1);
  Vandermonde1D(rDOFs, VInv.GetMat());
  VInv.Invert();

  /*--- Determine the gradient of the Vandermonde matrix of rInt. Make sure to
        allocate the number of rows to nIntPadLine and initialize them to zero. ---*/ 
  ColMajorMatrix<passivedouble> gradV(nIntPadLine,nPoly+1);
  gradV.setConstant(0.0);
  GradVandermonde1D(rInt, gradV);

  /*--- The derivatives of the Lagrangian basis functions can be obtained
        by multiplying gradV and VInv. ---*/
  VInv.MatMatMult('R', gradV, derLag);

  /*--- Check if the sum of the elements of the relevant rows of derLag is 0. ---*/
  for(unsigned short i=0; i<nIntLine; ++i) {
    passivedouble rowsum = 0.0;
    for(unsigned short j=0; j<=nPoly; ++j) rowsum += derLag(i,j);
    assert(fabs(rowsum) < 1.e-6);
  }
}

void CFEMStandardElementBase::LagBasisIntPointsLine(const vector<passivedouble>   &rDOFs,
                                                    const vector<passivedouble>   &rInt,
                                                    ColMajorMatrix<passivedouble> &lag) {

  /*--- Determine the number of integration points along the line
        and its padded value. ---*/
  const unsigned short nIntLine    = rInt.size();
  const unsigned short nIntPadLine = ((nIntLine+vecLen-1)/vecLen)*vecLen;

  /*--- Determine the inverse of the Vandermonde matrix of rDOFs. ---*/
  CGeneralSquareMatrixCM VInv(nPoly+1);
  Vandermonde1D(rDOFs, VInv.GetMat());
  VInv.Invert();

  /*--- Determine the Vandermonde matrix of rInt. Make sure to allocate
        the number of rows to nIntPadLine and initialize them to zero. ---*/ 
  ColMajorMatrix<passivedouble> V(nIntPadLine,nPoly+1);
  V.setConstant(0.0);
  Vandermonde1D(rInt, V);

  /*--- The Lagrangian basis functions can be obtained by multiplying
        V and VInv. ---*/
  VInv.MatMatMult('R', V, lag);

  /*--- Check if the sum of the elements of the relevant rows of lag is 1. ---*/
  for(unsigned short i=0; i<nIntLine; ++i) {
    passivedouble rowsum = -1.0;
    for(unsigned short j=0; j<=nPoly; ++j) rowsum += lag(i,j);
    assert(fabs(rowsum) < 1.e-6);
  }
}

passivedouble CFEMStandardElementBase::GradNormJacobi(unsigned short n,
                                                      unsigned short alpha,
                                                      unsigned short beta,
                                                      passivedouble  x) {

  /*--- Make a distinction for n == 0 and n > 0. For n == 0 the derivative is
        zero, because the polynomial itself is constant. ---*/
  passivedouble grad;
  if(n == 0) grad = 0.0;
  else
  {
    passivedouble tmp = n*(n+alpha+beta+1.0);
    grad              = sqrt(tmp)*NormJacobi(n-1, alpha+1, beta+1, x);
  }

  /*--- Return the gradient. ---*/
  return grad;
}

passivedouble CFEMStandardElementBase::NormJacobi(unsigned short n,
                                                  unsigned short alpha,
                                                  unsigned short beta,
                                                  passivedouble  x) {
  /*--- Some abbreviations. ---*/
  const passivedouble ap1   = alpha + 1;
  const passivedouble bp1   = beta  + 1;
  const passivedouble apb   = alpha + beta;
  const passivedouble apbp1 = apb + 1;
  const passivedouble apbp2 = apb + 2;
  const passivedouble apbp3 = apb + 3;
  const passivedouble b2ma2 = beta*beta - alpha*alpha;

  /*--- Determine the term, which involves the gamma function for P0. As the
        arguments are integers, this term can be computed easily, because
        Gamma(n+1) = n!. ---*/
  passivedouble Gamap1 = 1.0, Gambp1 = 1.0, Gamapbp2 = 1.0;
  for(unsigned short i=2; i<=alpha; ++i)          Gamap1   *= i;
  for(unsigned short i=2; i<=beta; ++i)           Gambp1   *= i;
  for(unsigned short i=2; i<=(alpha+beta+1); ++i) Gamapbp2 *= i;

  /*--- Initialize the normalized polynomials. ---*/
  passivedouble Pnm1 = sqrt(pow(0.5,apbp1)*Gamapbp2/(Gamap1*Gambp1));
  passivedouble Pn   = 0.5*Pnm1*(apbp2*x + alpha - beta)*sqrt(apbp3/(ap1*bp1));

  /*--- Take care of the special situation of n == 0. ---*/
  if(n == 0) Pn = Pnm1;
  else
  {
    /*--- The value of the normalized Jacobi polynomial must be obtained via recursion. ---*/
    for(unsigned short i=2; i<=n; ++i)
    {
      /*--- Compute the coefficients a for i and i-1 and the coefficient bi. ---*/
      unsigned short j = i-1;
      passivedouble tmp  = 2*j + apb;
      passivedouble aim1 = 2.0*sqrt(j*(j+apb)*(j+alpha)*(j+beta)/((tmp-1.0)*(tmp+1.0)))
                         / tmp;

      passivedouble bi = b2ma2/(tmp*(tmp+2.0));

      tmp              = 2*i + apb;
      passivedouble ai = 2.0*sqrt(i*(i+apb)*(i+alpha)*(i+beta)/((tmp-1.0)*(tmp+1.0)))
                       / tmp;

      /*--- Compute the new value of Pn and make sure to store Pnm1 correctly. ---*/
      tmp  = Pnm1;
      Pnm1 = Pn;

      Pn = ((x-bi)*Pn - aim1*tmp)/ai;
    }
  }

  /*--- Return Pn. ---*/
  return Pn;
}

void CFEMStandardElementBase::OwnGemm(const int                     M,
                                      const int                     N,
                                      const int                     K,
                                      ColMajorMatrix<passivedouble> &A,
                                      ColMajorMatrix<su2double>     &B,
                                      ColMajorMatrix<su2double>     &C,
                                      const CConfig                 *config) {

  /*--- Check if the jitted gemm call can be used. ---*/
#if defined(PRIMAL_SOLVER) && defined(HAVE_MKL)

  /*--- Carry out the timing, if desired and call the jitted gemm function. ---*/
#ifdef PROFILE
  double timeGemm;
  if( config ) config->GEMM_Tick(&timeGemm);
#endif

  my_dgemm(jitter, A.data(), B.data(), C.data());

#ifdef PROFILE
  if( config ) config->GEMM_Tock(timeGemm, M, N, K);
#endif

#else

  /*--- Use the interface to the more standard BLAS functionality. ---*/
  blasFunctions.gemm(M, N, K, A.data(), B.data(), C.data(), config);

#endif
}

void CFEMStandardElementBase::SetUpJittedGEMM(const int M,
                                              const int N,
                                              const int K) {

  /*--- Check if the jitted GEMM call can be created. ---*/
#if defined(PRIMAL_SOLVER) && defined(HAVE_MKL)

  /*--- Create the GEMM Kernel and check if it went okay. ---*/
  passivedouble alpha = 1.0;
  passivedouble beta  = 0.0;
  mkl_jit_status_t status = mkl_jit_create_dgemm(&jitter, MKL_COL_MAJOR, MKL_NOTRANS, MKL_NOTRANS,
                                                 M, N, K, alpha, M, K, beta, M);

  if(status == MKL_JIT_ERROR)
    SU2_MPI::Error(string("Jitted gemm kernel could not be created"), CURRENT_FUNCTION);

  /*--- Retrieve the function pointer to the DGEMM kernel
        void my_dgemm(void*, double*, double*, double*). ---*/
  my_dgemm = mkl_jit_get_dgemm_ptr(jitter);

#endif
}

/*----------------------------------------------------------------------------------*/
/*          Private member functions of CFEMStandardElementBase.                    */
/*----------------------------------------------------------------------------------*/

passivedouble CFEMStandardElementBase::Legendre(unsigned short n,
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

void CFEMStandardElementBase::GradVandermonde1D(const vector<passivedouble>   &r,
                                                ColMajorMatrix<passivedouble> &VDr) {

  /*--- Compute the gradient of the 1D Vandermonde matrix. ---*/
  for(unsigned short j=0; j<=nPoly; ++j)
    for(unsigned short i=0; i<r.size(); ++i)
      VDr(i,j) = GradNormJacobi(j, 0, 0, r[i]);
}

void CFEMStandardElementBase::Vandermonde1D(const vector<passivedouble>   &r,
                                            ColMajorMatrix<passivedouble> &V) {

  /*--- Compute the 1D Vandermonde matrix. ---*/
  for(unsigned short j=0; j<=nPoly; ++j)
    for(unsigned short i=0; i<r.size(); ++i)
      V(i,j) = NormJacobi(j, 0, 0, r[i]);
}

void CFEMStandardElementBase::WarpFactor(const vector<passivedouble> &rout,
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
  Location1DGridDOFsEquidistant(r1DEqui);
  Location1DGridDOFsLGL(r1DLGL);

  /*--- Determine the 1D Vandermonde matrix based on the equidistant node
        distribution. Invert the matrix afterwards. ---*/
  CGeneralSquareMatrixCM V1D(nPoly+1);
  Vandermonde1D(r1DEqui, V1D.GetMat());
  V1D.Invert();

  /*--- Determine the Lagrange polynomials at rout. This is accomplished by
        first computing the Legendre polynomials at rout and multiply
        this value by V1D. The result is stored in Lmat. ---*/
  ColMajorMatrix<passivedouble> Pmat(nD,nPoly+1), Lmat;
  Vandermonde1D(rout, Pmat);
  V1D.MatMatMult('R', Pmat, Lmat);

  /*--- Determine the difference between the 1D LGL points and
        the equidistant points. The result is stored in r1DLGL. ---*/
  for(unsigned short i=0; i<=nPoly; ++i) r1DLGL[i] -= r1DEqui[i];

  /*--- Compute the unscaled warp factors, which is Lmat*r1DLGL. ---*/
  for(unsigned short j=0; j<nD; ++j) {
    warp[j] = 0.0;
    for(unsigned short i=0; i<=nPoly; ++i)
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
