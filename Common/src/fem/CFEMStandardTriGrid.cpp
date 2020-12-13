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

  /*--- Determine the local subconnectivity of this standard element when split
        in several linear elements. Used for a.o. plotting and searcing. ---*/
  SubConnLinearElements();

  /*--- Set up the jitted gemm call, if supported. For this particular standard
        element the derivatives of the coordinates are computed, which is nDim. ---*/
  SetUpJittedGEMM(nIntegrationPad, nDim, nDOFs, jitter, my_dgemm);
}

CFEMStandardTriGrid::CFEMStandardTriGrid(const unsigned short val_nPolyGrid,
                                         const unsigned short val_nPolySol,
                                         const unsigned short val_orderExact,
                                         const unsigned short val_locGridDOFs)
  : CFEMStandardTri(val_nPolyGrid, val_orderExact) {

  /*--- Determine the number of space dimensions, which is 2 when this
        constructor is called. ---*/
  nDim = 2;

  /*--- Determine the location of the grid DOFs and build the appropriate
        Lagrangian basis functions. ---*/
  if(val_locGridDOFs == LGL) {

    /*--- LGL distribution. Compute the corresponding Lagrangian basis functions and
          its first and second derivatives in the integration points. ---*/
    LagBasisIntPointsTriangle(rTriangleDOFsLGL, sTriangleDOFsLGL, lagBasisIntLGL);
    DerLagBasisIntPointsTriangle(rTriangleDOFsLGL, sTriangleDOFsLGL, derLagBasisIntLGL);
    HesLagBasisIntPointsTriangle(rTriangleDOFsLGL, sTriangleDOFsLGL);

    /*--- Determine the location of nodal solution DOFs. ---*/
    nPoly = val_nPolySol;
    LocationTriangleGridDOFsLGL(rTriangleSolDOFs, sTriangleSolDOFs);
    nPoly = val_nPolyGrid;

    /*--- Compute the Lagrangian basis functions and its derivatives
          in the nodal solution DOFs. ---*/
    LagBasisAndDerSolDOFsTriangle(rTriangleDOFsLGL, sTriangleDOFsLGL);
  }
  else {

    /*--- Equidistant distribution. Compute the corresponding Lagrangian basis functions and
          its first and second derivatives in the integration points. ---*/
    LagBasisIntPointsTriangle(rTriangleDOFsEqui, sTriangleDOFsEqui, lagBasisIntEqui);
    DerLagBasisIntPointsTriangle(rTriangleDOFsEqui, sTriangleDOFsEqui, derLagBasisIntEqui);
    HesLagBasisIntPointsTriangle(rTriangleDOFsEqui, sTriangleDOFsEqui);

    /*--- Determine the location of nodal solution DOFs. ---*/
    nPoly = val_nPolySol;
    LocationTriangleGridDOFsEquidistant(rTriangleSolDOFs, sTriangleSolDOFs);
    nPoly = val_nPolyGrid;

    /*--- Compute the Lagrangian basis functions and its derivatives
          in the nodal solution DOFs. ---*/
    LagBasisAndDerSolDOFsTriangle(rTriangleDOFsEqui, sTriangleDOFsEqui);
  }

  /*--- Create the local grid connectivities of the faces of the volume element. ---*/
  LocalGridConnFaces();

  /*--- Determine the local subconnectivity of this standard element when split
        in several linear elements. Used for a.o. plotting and searcing. ---*/
  SubConnLinearElements();

  /*--- Set up the jitted gemm call, if supported. For this particular standard
        element the derivatives of the coordinates are computed, which is nDim. ---*/
  SetUpJittedGEMM(nIntegrationPad, nDim, nDOFs, jitter, my_dgemm);

  /*--- Set up the jitted gemm call, if supported, to compute the data
        in the nodal solution DOFs. ---*/
  SetUpJittedGEMM(lagBasisSolDOFs.rows(), nDim, nDOFs, jitterSolDOFs, dgemmSolDOFs);
}

CFEMStandardTriGrid::~CFEMStandardTriGrid() {

#if defined(PRIMAL_SOLVER) && defined(HAVE_MKL)
  if( jitterSolDOFs ) {
    mkl_jit_destroy(jitterSolDOFs);
    jitterSolDOFs = nullptr;
  }
#endif
}

void CFEMStandardTriGrid::CoorIntPoints(const bool                LGLDistribution,
                                        ColMajorMatrix<su2double> &matCoorDOF,
                                        ColMajorMatrix<su2double> &matCoorInt) {

  /*--- Check for which point distribution the derivatives must be computed. ---*/
  if( LGLDistribution ) {

    /*--- LGL distribution. Call the function OwnGemm to compute the Cartesian
          coordinates in the integration points. The second argument in the
          function call is nDim, which corresponds to the number of Cartesian
          coordinates (3 for a surface element and 2 for a volume element). ---*/
    OwnGemm(my_dgemm, nIntegrationPad, nDim, nDOFs, lagBasisIntLGL, matCoorDOF, matCoorInt, nullptr);
  }
  else {

    /*--- Equidistant distribution. Call the function OwnGemm to compute the Cartesian
          coordinates in the integration points. The second argument in the
          function call is nDim, which corresponds to the number of Cartesian
          coordinates (3 for a surface element and 2 for a volume element). ---*/
    OwnGemm(my_dgemm, nIntegrationPad, nDim, nDOFs, lagBasisIntEqui, matCoorDOF, matCoorInt, nullptr);
  }
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
    OwnGemm(my_dgemm, nIntegrationPad, nDim, nDOFs, derLagBasisIntLGL[0], matCoor, matDerCoor[0], nullptr);
    OwnGemm(my_dgemm, nIntegrationPad, nDim, nDOFs, derLagBasisIntLGL[1], matCoor, matDerCoor[1], nullptr);
  }
  else {

    /*--- Equidistant distribution. Call the function OwnGemm 2 times to compute the derivatives
          of the Cartesian coordinates w.r.t. the two parametric coordinates. The third
          argument in the function call is nDim, which corresponds to the number of Cartesian
          coordinates (3 for a surface element and 2 for a volume element). ---*/
    OwnGemm(my_dgemm, nIntegrationPad, nDim, nDOFs, derLagBasisIntEqui[0], matCoor, matDerCoor[0], nullptr);
    OwnGemm(my_dgemm, nIntegrationPad, nDim, nDOFs, derLagBasisIntEqui[1], matCoor, matDerCoor[1], nullptr);
  }
}

void CFEMStandardTriGrid::DerivativesCoorSolDOFs(ColMajorMatrix<su2double>          &matCoor,
                                                 vector<ColMajorMatrix<su2double> > &matDerCoor) {

  /*--- Call OwnGemm 2 times to compute the derivatives of the Cartesian
        coordinates w.r.t. the two parametric coordinates. The third
        argument in the function call is nDim, which corresponds to the number of Cartesian
        coordinates (3 for a surface element and 2 for a volume element). ---*/
  const unsigned short nSolDOFs = lagBasisSolDOFs.rows();
  OwnGemm(dgemmSolDOFs, nSolDOFs, nDim, nDOFs, derLagBasisSolDOFs[0], matCoor, matDerCoor[0], nullptr);
  OwnGemm(dgemmSolDOFs, nSolDOFs, nDim, nDOFs, derLagBasisSolDOFs[1], matCoor, matDerCoor[1], nullptr);
}

passivedouble CFEMStandardTriGrid::WorkEstimateBoundaryFace(CConfig              *config,
                                                            const unsigned short elemType) {

  /*--- Determine the number of DOFs of the neighboring element. ---*/
  const unsigned short nDOFsElem = GetNDOFsStatic(elemType, nPoly);

  /*--- TEMPORARY IMPLEMENTATION. ---*/
  return nIntegration + 0.05*nDOFsElem;
}

passivedouble CFEMStandardTriGrid::WorkEstimateInternalFace(CConfig              *config,
                                                            const unsigned short elemType0,
                                                            const unsigned short nPoly0,
                                                            const unsigned short elemType1,
                                                            const unsigned short nPoly1) {

  /*--- Determine the number of DOFs of the neighboring elements. ---*/
  const unsigned short nDOFsElem0 = GetNDOFsStatic(elemType0, nPoly0);
  const unsigned short nDOFsElem1 = GetNDOFsStatic(elemType1, nPoly1);

  /* TEMPORARY IMPLEMENTATION. */
  return 2.0*nIntegration + 0.05*(nDOFsElem0 + nDOFsElem1);
}

passivedouble CFEMStandardTriGrid::WorkEstimateVolume(CConfig *config) {

  /*--- TEMPORARY IMPLEMENTATION. ---*/
  return nIntegration + 0.1*nDOFs;
}

passivedouble CFEMStandardTriGrid::WorkEstimateWallFunctions(CConfig              *config,
                                                             const unsigned short nPointsWF,
                                                             const unsigned short elemType) {

  /*--- TEMPORARY IMPLEMENTATION. ---*/
  return 0.25*nIntegration*nPointsWF;
}

/*----------------------------------------------------------------------------------*/
/*             Private member functions of CFEMStandardTriGrid.                     */
/*----------------------------------------------------------------------------------*/

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

void CFEMStandardTriGrid::HesLagBasisIntPointsTriangle(const vector<passivedouble> &rDOFs,
                                                       const vector<passivedouble> &sDOFs) {

  /*--- Determine the padded number of the total number of integration points. ---*/
  const unsigned short nIntTot    = rTriangleInt.size();
  const unsigned short nIntTotPad = ((nIntTot+baseVectorLen-1)/baseVectorLen)*baseVectorLen;

  /*--- Determine the inverse of the Vandermonde matrix of the DOFs. ---*/
  CGeneralSquareMatrixCM VInv(rDOFs.size());
  VandermondeTriangle(rDOFs, sDOFs, VInv.GetMat());
  VInv.Invert();

  /*--- Determine the Hessian of the Vandermonde matrix of the integration points. Make
        sure to allocate the number of rows to nIntTotPad and initialize them to zero. ---*/
  ColMajorMatrix<passivedouble> VDr2(nIntTotPad,rDOFs.size()),
                                VDs2(nIntTotPad,rDOFs.size()),
                                VDrs(nIntTotPad,rDOFs.size());
  VDr2.setConstant(0.0);
  VDs2.setConstant(0.0);
  VDrs.setConstant(0.0);

  HesVandermondeTriangle(rTriangleInt, sTriangleInt, VDr2, VDs2, VDrs);

  /*--- The Hessian of the Lagrangian basis functions can be obtained by
        multiplying VDr2, VDs2, VDrs and VInv. ---*/
  hesLagBasisInt.resize(3);
  VInv.MatMatMult('R', VDr2, hesLagBasisInt[0]);
  VInv.MatMatMult('R', VDs2, hesLagBasisInt[1]);
  VInv.MatMatMult('R', VDrs, hesLagBasisInt[2]);

  /*--- Check if the sum of the elements of the relevant rows of hesLagBasisInt is 0. ---*/
  for(unsigned short i=0; i<nIntTot; ++i) {
    passivedouble rowSumDr2 = 0.0, rowSumDs2 = 0.0, rowSumDrs = 0.0;
    for(unsigned short j=0; j<rDOFs.size(); ++j) {
      rowSumDr2 += hesLagBasisInt[0](i,j);
      rowSumDs2 += hesLagBasisInt[1](i,j);
      rowSumDrs += hesLagBasisInt[2](i,j);
    }

    assert(fabs(rowSumDr2) < 1.e-6);
    assert(fabs(rowSumDs2) < 1.e-6);
    assert(fabs(rowSumDrs) < 1.e-6);
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

void CFEMStandardTriGrid::LagBasisAndDerSolDOFsTriangle(const vector<passivedouble> &rDOFs,
                                                        const vector<passivedouble> &sDOFs) {

  /*--- Determine the inverse of the Vandermonde matrix of the DOFs. ---*/
  CGeneralSquareMatrixCM VInv(rDOFs.size());
  VandermondeTriangle(rDOFs, sDOFs, VInv.GetMat());
  VInv.Invert();

  /*--- Determine the Vandermonde matrix and its gradients of the solution DOFs. ---*/
  ColMajorMatrix<passivedouble> V(rTriangleSolDOFs.size(),rDOFs.size()),
                                VDr(rTriangleSolDOFs.size(),rDOFs.size()),
                                VDs(rTriangleSolDOFs.size(),rDOFs.size());

  VandermondeTriangle(rTriangleSolDOFs, sTriangleSolDOFs, V);
  GradVandermondeTriangle(rTriangleSolDOFs, sTriangleSolDOFs, VDr, VDs);

  /*--- The Lagrangian basis functions and its gradients can be obtained by
        multiplying V, VDr, VDs and VInv. ---*/
  derLagBasisSolDOFs.resize(2);
  VInv.MatMatMult('R', V,   lagBasisSolDOFs);
  VInv.MatMatMult('R', VDr, derLagBasisSolDOFs[0]);
  VInv.MatMatMult('R', VDs, derLagBasisSolDOFs[1]);

  /*--- Check if the sum of the elements of the relevant rows of
        lagBasisSolDOFs is 1 and of derLagBasisSolDOFs is 0. ---*/
  for(unsigned short i=0; i<rTriangleSolDOFs.size(); ++i) {
    passivedouble rowSum = -1.0, rowSumDr = 0.0, rowSumDs = 0.0;
    for(unsigned short j=0; j<rDOFs.size(); ++j) {
      rowSum   += lagBasisSolDOFs(i,j);
      rowSumDr += derLagBasisSolDOFs[0](i,j);
      rowSumDs += derLagBasisSolDOFs[1](i,j);
    }

    assert(fabs(rowSum)   < 1.e-6);
    assert(fabs(rowSumDr) < 1.e-6);
    assert(fabs(rowSumDs) < 1.e-6);
  }
}

void CFEMStandardTriGrid::GradVandermondeTriangle(const vector<passivedouble>   &r,
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
  for(unsigned short i=0; i<=nPoly; ++i) {
    for(unsigned short j=0; j<=(nPoly-i); ++j, ++ii) {
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

void CFEMStandardTriGrid::HesVandermondeTriangle(const vector<passivedouble>   &r,
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
  for(unsigned short i=0; i<=nPoly; ++i) {
    for(unsigned short j=0; j<=(nPoly-i); ++j, ++ii) {
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

void CFEMStandardTriGrid::SubConnLinearElements(void) {

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
