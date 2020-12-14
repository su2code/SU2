/*!
 * \file CFEMStandardTetGrid.cpp
 * \brief Functions for the class CFEMStandardTetGrid.
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

#include "../../include/fem/CFEMStandardTetGrid.hpp"
#include "../../include/toolboxes/CSquareMatrixCM.hpp"

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

  /*--- Create the local grid connectivities of the faces of the volume element. ---*/
  LocalGridConnFaces();

  /*--- Determine the local subconnectivity of this standard element when split
        in several linear elements. Used for a.o. plotting and searcing. ---*/
  SubConnLinearElements();

  /*--- Set up the jitted gemm call, if supported. For this particular standard
        element the derivative of the coordinates are computed, which is 3. ---*/
  SetUpJittedGEMM(nIntegrationPad, 3, nDOFs, jitter, my_dgemm);
}

CFEMStandardTetGrid::CFEMStandardTetGrid(const unsigned short val_nPolyGrid,
                                         const unsigned short val_nPolySol,
                                         const unsigned short val_orderExact,
                                         const unsigned short val_locGridDOFs)
  : CFEMStandardTet(val_nPolyGrid, val_orderExact) {

  /*--- Determine the location of the grid DOFs and build the appropriate
        Lagrangian basis functions. ---*/
  if(val_locGridDOFs == LGL) {

    /*--- LGL distribution. Compute the corresponding Lagrangian basis functions and
          its first and second derivatives in the integration points. ---*/
    LagBasisIntPointsTet(rTetDOFsLGL, sTetDOFsLGL, tTetDOFsLGL, lagBasisIntLGL);
    DerLagBasisIntPointsTet(rTetDOFsLGL, sTetDOFsLGL, tTetDOFsLGL, derLagBasisIntLGL);
    HesLagBasisIntPointsTet(rTetDOFsLGL, sTetDOFsLGL, tTetDOFsLGL);

    /*--- Determine the location of nodal solution DOFs. ---*/
    nPoly = val_nPolySol;
    LocationTetGridDOFsLGL(rTetSolDOFs, sTetSolDOFs, tTetSolDOFs);
    nPoly = val_nPolyGrid;

    /*--- Compute the Lagrangian basis functions and its derivatives
          in the nodal solution DOFs. ---*/
    LagBasisAndDerSolDOFsTet(rTetDOFsLGL, sTetDOFsLGL, tTetDOFsLGL);
  }
  else {

    /*--- Equidistant distribution. Compute the corresponding Lagrangian basis functions and
          its first and second derivatives in the integration points. ---*/
    LagBasisIntPointsTet(rTetDOFsEqui, sTetDOFsEqui, tTetDOFsEqui, lagBasisIntEqui);
    DerLagBasisIntPointsTet(rTetDOFsEqui, sTetDOFsEqui, tTetDOFsEqui, derLagBasisIntEqui);
    HesLagBasisIntPointsTet(rTetDOFsEqui, sTetDOFsEqui, tTetDOFsEqui);

    /*--- Determine the location of nodal solution DOFs. ---*/
    nPoly = val_nPolySol;
    LocationTetGridDOFsEquidistant(rTetSolDOFs, sTetSolDOFs, tTetSolDOFs);
    nPoly = val_nPolyGrid;

    /*--- Compute the Lagrangian basis functions and its derivatives
          in the nodal solution DOFs. ---*/
    LagBasisAndDerSolDOFsTet(rTetDOFsEqui, sTetDOFsEqui, tTetDOFsEqui);
  }

  /*--- Create the local grid connectivities of the faces of the volume element. ---*/
  LocalGridConnFaces();

  /*--- Determine the local subconnectivity of this standard element when split
        in several linear elements. Used for a.o. plotting and searcing. ---*/
  SubConnLinearElements();

  /*--- Set up the jitted gemm call, if supported. For this particular standard
        element the derivative of the coordinates are computed, which is 3. ---*/
  SetUpJittedGEMM(nIntegrationPad, 3, nDOFs, jitter, my_dgemm);

  /*--- Set up the jitted gemm call, if supported, to compute the data
        in the nodal solution DOFs. ---*/
  SetUpJittedGEMM(lagBasisSolDOFs.rows(), 3, nDOFs, jitterSolDOFs, dgemmSolDOFs);
}

CFEMStandardTetGrid::~CFEMStandardTetGrid() {

#if defined(PRIMAL_SOLVER) && defined(HAVE_MKL)
  if( jitterSolDOFs ) {
    mkl_jit_destroy(jitterSolDOFs);
    jitterSolDOFs = nullptr;
  }
#endif
}

void CFEMStandardTetGrid::CoorIntPoints(const bool                LGLDistribution,
                                        ColMajorMatrix<su2double> &matCoorDOF,
                                        ColMajorMatrix<su2double> &matCoorInt) {

  /*--- Check for which point distribution the derivatives must be computed. ---*/
  if( LGLDistribution ) {

    /*--- LGL distribution. Call the function OwnGemm to compute the
          Cartesian coordinates in the integration points. ---*/
    OwnGemm(my_dgemm, nIntegrationPad, 3, nDOFs, lagBasisIntLGL, matCoorDOF, matCoorInt, nullptr);
  }
  else {

    /*--- Equidistant distribution. Call the function OwnGemm to compute the
          Cartesian coordinates in the integration points. ---*/
    OwnGemm(my_dgemm, nIntegrationPad, 3, nDOFs, lagBasisIntEqui, matCoorDOF, matCoorInt, nullptr);
  }
}

void CFEMStandardTetGrid::DerivativesCoorIntPoints(const bool                         LGLDistribution,
                                                   ColMajorMatrix<su2double>          &matCoor,
                                                   vector<ColMajorMatrix<su2double> > &matDerCoor) {
  /*--- Check for which point distribution the derivatives must be computed. ---*/
  if( LGLDistribution ) {

    /*--- LGL distribution. Call the function OwnGemm 3 times to compute the derivatives
          of the Cartesian coordinates w.r.t. the three parametric coordinates. ---*/
    OwnGemm(my_dgemm, nIntegrationPad, 3, nDOFs, derLagBasisIntLGL[0], matCoor, matDerCoor[0], nullptr);
    OwnGemm(my_dgemm, nIntegrationPad, 3, nDOFs, derLagBasisIntLGL[1], matCoor, matDerCoor[1], nullptr);
    OwnGemm(my_dgemm, nIntegrationPad, 3, nDOFs, derLagBasisIntLGL[2], matCoor, matDerCoor[2], nullptr);
  }
  else {

    /*--- Equidistant distribution. Call the function OwnGemm 3 times to compute the derivatives
          of the Cartesian coordinates w.r.t. the three parametric coordinates. ---*/
    OwnGemm(my_dgemm, nIntegrationPad, 3, nDOFs, derLagBasisIntEqui[0], matCoor, matDerCoor[0], nullptr);
    OwnGemm(my_dgemm, nIntegrationPad, 3, nDOFs, derLagBasisIntEqui[1], matCoor, matDerCoor[1], nullptr);
    OwnGemm(my_dgemm, nIntegrationPad, 3, nDOFs, derLagBasisIntEqui[2], matCoor, matDerCoor[2], nullptr);
  }
}

void CFEMStandardTetGrid::DerivativesCoorSolDOFs(ColMajorMatrix<su2double>          &matCoor,
                                                 vector<ColMajorMatrix<su2double> > &matDerCoor) {

  /*--- Call OwnGemm 3 times to compute the derivatives of the Cartesian
        coordinates w.r.t. the three parametric coordinates. ---*/
  const unsigned short nSolDOFs = lagBasisSolDOFs.rows();
  OwnGemm(dgemmSolDOFs, nSolDOFs, 3, nDOFs, derLagBasisSolDOFs[0], matCoor, matDerCoor[0], nullptr);
  OwnGemm(dgemmSolDOFs, nSolDOFs, 3, nDOFs, derLagBasisSolDOFs[1], matCoor, matDerCoor[1], nullptr);
  OwnGemm(dgemmSolDOFs, nSolDOFs, 3, nDOFs, derLagBasisSolDOFs[2], matCoor, matDerCoor[2], nullptr);
}

passivedouble CFEMStandardTetGrid::WorkEstimateVolume(CConfig *config) {

  /*--- TEMPORARY IMPLEMENTATION. ---*/
  return nIntegration + 0.1*nDOFs;
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
  CSquareMatrixCM VInv(rDOFs.size());
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

void CFEMStandardTetGrid::HesLagBasisIntPointsTet(const vector<passivedouble> &rDOFs,
                                                  const vector<passivedouble> &sDOFs,
                                                  const vector<passivedouble> &tDOFs) {

  /*--- Determine the padded number of the total number of integration points. ---*/
  const unsigned short nIntTot    = rTetInt.size();
  const unsigned short nIntTotPad = ((nIntTot+baseVectorLen-1)/baseVectorLen)*baseVectorLen;

  /*--- Determine the inverse of the Vandermonde matrix of the DOFs. ---*/
  CSquareMatrixCM VInv(rDOFs.size());
  VandermondeTetrahedron(rDOFs, sDOFs, tDOFs, VInv.GetMat());
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

  HesVandermondeTetrahedron(rTetInt, sTetInt, tTetInt, VDr2, VDs2, VDt2,
                            VDrs, VDrt, VDst);

  /*--- The Hessian of the Lagrangian basis functions can be obtained by
        multiplying VDr2, VDs2, etc and VInv. ---*/
  hesLagBasisInt.resize(6);
  VInv.MatMatMult('R', VDr2, hesLagBasisInt[0]);
  VInv.MatMatMult('R', VDs2, hesLagBasisInt[1]);
  VInv.MatMatMult('R', VDt2, hesLagBasisInt[2]);
  VInv.MatMatMult('R', VDrs, hesLagBasisInt[3]);
  VInv.MatMatMult('R', VDrt, hesLagBasisInt[4]);
  VInv.MatMatMult('R', VDst, hesLagBasisInt[5]);

  /*--- Check if the sum of the elements of the relevant rows of hesLagBasisInt is 0. ---*/
  for(unsigned short i=0; i<nIntTot; ++i) {
    passivedouble rowSumDr2 = 0.0, rowSumDs2 = 0.0, rowSumDt2 = 0.0;
    passivedouble rowSumDrs = 0.0, rowSumDrt = 0.0, rowSumDst = 0.0;
    for(unsigned short j=0; j<rDOFs.size(); ++j) {
      rowSumDr2 += hesLagBasisInt[0](i,j);
      rowSumDs2 += hesLagBasisInt[1](i,j);
      rowSumDt2 += hesLagBasisInt[2](i,j);
      rowSumDrs += hesLagBasisInt[3](i,j);
      rowSumDrt += hesLagBasisInt[4](i,j);
      rowSumDst += hesLagBasisInt[5](i,j);
    }

    assert(fabs(rowSumDr2) < 1.e-6);
    assert(fabs(rowSumDs2) < 1.e-6);
    assert(fabs(rowSumDt2) < 1.e-6);
    assert(fabs(rowSumDrs) < 1.e-6);
    assert(fabs(rowSumDrt) < 1.e-6);
    assert(fabs(rowSumDst) < 1.e-6);
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
  CSquareMatrixCM VInv(rDOFs.size());
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

void CFEMStandardTetGrid::LagBasisAndDerSolDOFsTet(const vector<passivedouble>  &rDOFs,
                                                   const vector<passivedouble>  &sDOFs,
                                                   const vector<passivedouble>  &tDOFs) {

  /*--- Determine the inverse of the Vandermonde matrix of the DOFs. ---*/
  CSquareMatrixCM VInv(rDOFs.size());
  VandermondeTetrahedron(rDOFs, sDOFs, tDOFs, VInv.GetMat());
  VInv.Invert();

  /*--- Determine the Vandermonde matrix and its gradients of the solution DOFs. ---*/
  ColMajorMatrix<passivedouble> V(rTetSolDOFs.size(),rDOFs.size()),
                                VDr(rTetSolDOFs.size(),rDOFs.size()),
                                VDs(rTetSolDOFs.size(),rDOFs.size()),
                                VDt(rTetSolDOFs.size(),rDOFs.size());

  VandermondeTetrahedron(rTetSolDOFs, sTetSolDOFs, tTetSolDOFs, V);
  GradVandermondeTetrahedron(rTetSolDOFs, sTetSolDOFs, tTetSolDOFs, VDr, VDs, VDt);

  /*--- The Lagrangian basis functions and its gradients can be obtained by
        multiplying V, VDr, VDs, VDt and VInv. ---*/
  derLagBasisSolDOFs.resize(3);
  VInv.MatMatMult('R', V,   lagBasisSolDOFs);
  VInv.MatMatMult('R', VDr, derLagBasisSolDOFs[0]);
  VInv.MatMatMult('R', VDs, derLagBasisSolDOFs[1]);
  VInv.MatMatMult('R', VDt, derLagBasisSolDOFs[2]);

  /*--- Check if the sum of the elements of the relevant rows of
        lagBasisSolDOFs is 1 and of derLagBasisSolDOFs is 0. ---*/
  for(unsigned short i=0; i<rTetSolDOFs.size(); ++i) {
    passivedouble rowSum = -1.0, rowSumDr = 0.0, rowSumDs = 0.0, rowSumDt = 0.0;
    for(unsigned short j=0; j<rDOFs.size(); ++j) {
      rowSum   += lagBasisSolDOFs(i,j);
      rowSumDr += derLagBasisSolDOFs[0](i,j);
      rowSumDs += derLagBasisSolDOFs[1](i,j);
      rowSumDt += derLagBasisSolDOFs[2](i,j);
    }

    assert(fabs(rowSum)   < 1.e-6);
    assert(fabs(rowSumDr) < 1.e-6);
    assert(fabs(rowSumDs) < 1.e-6);
    assert(fabs(rowSumDt) < 1.e-6);
  }
}

void CFEMStandardTetGrid::GradVandermondeTetrahedron(const vector<passivedouble>   &r,
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

void CFEMStandardTetGrid::HesVandermondeTetrahedron(const vector<passivedouble>   &r,
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

void CFEMStandardTetGrid::LocalGridConnFaces(void) {

  /*--- Allocate the first index of gridConnFaces, which is equal to the number
        of faces of the tetrahedron, which is 4. Reserve memory for the second
        index afterwards. ---*/
  const unsigned short nDOFsTriangle = (nPoly+1)*(nPoly+2)/2;
  gridConnFaces.resize(4);

  gridConnFaces[0].reserve(nDOFsTriangle);
  gridConnFaces[1].reserve(nDOFsTriangle);
  gridConnFaces[2].reserve(nDOFsTriangle);
  gridConnFaces[3].reserve(nDOFsTriangle);

  /*--- Loop over all the nodes of the tetrahedron and pick the correct
        ones for the faces. ---*/
  unsigned int ii = 0;
  for(unsigned short k=0; k<=nPoly; ++k) {
    const unsigned short uppBoundJ = nPoly - k;
    for(unsigned short j=0; j<=uppBoundJ; ++j) {
      const unsigned short uppBoundI = nPoly - k - j;
      for(unsigned short i=0; i<=uppBoundI; ++i, ++ii) {
        if(k == 0)           gridConnFaces[0].push_back(ii);
        if(j == 0)           gridConnFaces[1].push_back(ii);
        if(i == 0)           gridConnFaces[2].push_back(ii);
        if((i+j+k) == nPoly) gridConnFaces[3].push_back(ii);
      }
    }
  }

  /*--- Make sure that the element is to the left of the faces. ---*/
  const unsigned short n0 = 0;
  const unsigned short n1 = nPoly;
  const unsigned short n2 = nDOFsTriangle -1;
  const unsigned short n3 = nDOFs -1;

  ChangeDirectionTriangleConn(gridConnFaces[0], n0, n1, n2);
  ChangeDirectionTriangleConn(gridConnFaces[1], n0, n3, n1);
  ChangeDirectionTriangleConn(gridConnFaces[2], n0, n2, n3);
  ChangeDirectionTriangleConn(gridConnFaces[3], n1, n3, n2);
}

void CFEMStandardTetGrid::SubConnLinearElements(void) {

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
