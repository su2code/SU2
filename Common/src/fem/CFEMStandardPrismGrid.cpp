/*!
 * \file CFEMStandardPrismGrid.cpp
 * \brief Functions for the class CFEMStandardPrismGrid.
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

#include "../../include/fem/CFEMStandardPrismGrid.hpp"
#include "../../include/toolboxes/CSquareMatrixCM.hpp"

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

  /*--- Determine the local subconnectivity of this standard element when split
        in several linear elements. Used for a.o. plotting and searcing. ---*/
  SubConnLinearElements();

  /*--- Set up the jitted gemm call, if supported. For this particular standard
        element the derivative of the coordinates are computed, which is 3. ---*/
  SetUpJittedGEMM(nIntegrationPad, 3, nDOFs, jitter, my_dgemm);
}

CFEMStandardPrismGrid::CFEMStandardPrismGrid(const unsigned short val_nPolyGrid,
                                             const unsigned short val_nPolySol,
                                             const unsigned short val_orderExact,
                                             const unsigned short val_locGridDOFs)
  : CFEMStandardPrism(val_nPolyGrid, val_orderExact) {

  /*--- Determine the location of the grid DOFs and build the appropriate
        Lagrangian basis functions. ---*/
  if(val_locGridDOFs == LGL) {

    /*--- LGL distribution. Compute the corresponding Lagrangian basis functions and
          its first and second derivatives in the integration points. ---*/
    LagBasisIntPointsPrism(rTriangleDOFsLGL, sTriangleDOFsLGL, rLineDOFsLGL, lagBasisIntLGL);
    DerLagBasisIntPointsPrism(rTriangleDOFsLGL, sTriangleDOFsLGL, rLineDOFsLGL, derLagBasisIntLGL);
    HesLagBasisIntPointsPrism(rTriangleDOFsLGL, sTriangleDOFsLGL, rLineDOFsLGL);

    /*--- Determine the location of nodal solution DOFs. ---*/
    nPoly = val_nPolySol;
    Location1DGridDOFsLGL(rLineSolDOFs);
    LocationTriangleGridDOFsLGL(rTriangleSolDOFs, sTriangleSolDOFs);
    nPoly = val_nPolyGrid;

    /*--- Compute the Lagrangian basis functions and its derivatives
          in the nodal solution DOFs. ---*/
    LagBasisAndDerSolDOFsPrism(rTriangleDOFsLGL, sTriangleDOFsLGL, rLineDOFsLGL);
  }
  else {

    /*--- Equidistant distribution. Compute the corresponding Lagrangian basis functions and
          its first and second derivatives in the integration points. ---*/
    LagBasisIntPointsPrism(rTriangleDOFsEqui, sTriangleDOFsEqui, rLineDOFsEqui, lagBasisIntEqui);
    DerLagBasisIntPointsPrism(rTriangleDOFsEqui, sTriangleDOFsEqui, rLineDOFsEqui, derLagBasisIntEqui);
    HesLagBasisIntPointsPrism(rTriangleDOFsEqui, sTriangleDOFsEqui, rLineDOFsEqui);

    /*--- Determine the location of nodal solution DOFs. ---*/
    nPoly = val_nPolySol;
    Location1DGridDOFsEquidistant(rLineSolDOFs);
    LocationTriangleGridDOFsEquidistant(rTriangleSolDOFs, sTriangleSolDOFs);
    nPoly = val_nPolyGrid;

    /*--- Compute the Lagrangian basis functions and its derivatives
          in the nodal solution DOFs. ---*/
    LagBasisAndDerSolDOFsPrism(rTriangleDOFsEqui, sTriangleDOFsEqui, rLineDOFsEqui);
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

CFEMStandardPrismGrid::~CFEMStandardPrismGrid() {

#if defined(PRIMAL_SOLVER) && defined(HAVE_MKL)
  if( jitterSolDOFs ) {
    mkl_jit_destroy(jitterSolDOFs);
    jitterSolDOFs = nullptr;
  }
#endif
}

void CFEMStandardPrismGrid::CoorIntPoints(const bool                LGLDistribution,
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

void CFEMStandardPrismGrid::DerivativesCoorIntPoints(const bool                         LGLDistribution,
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

void CFEMStandardPrismGrid::Derivatives2ndCoorIntPoints(const bool                         LGLDistribution,
                                                        ColMajorMatrix<su2double>          &matCoor,
                                                        vector<ColMajorMatrix<su2double> > &matDer2ndCoor) {

  /*--- Call the function OwnGemm 6 times to compute the 2nd derivatives of
        the Cartesian coordinates w.r.t. the three parametric coordinates. ---*/
  OwnGemm(my_dgemm, nIntegrationPad, 3, nDOFs, hesLagBasisInt[0], matCoor, matDer2ndCoor[0], nullptr);
  OwnGemm(my_dgemm, nIntegrationPad, 3, nDOFs, hesLagBasisInt[1], matCoor, matDer2ndCoor[1], nullptr);
  OwnGemm(my_dgemm, nIntegrationPad, 3, nDOFs, hesLagBasisInt[2], matCoor, matDer2ndCoor[2], nullptr);
  OwnGemm(my_dgemm, nIntegrationPad, 3, nDOFs, hesLagBasisInt[3], matCoor, matDer2ndCoor[3], nullptr);
  OwnGemm(my_dgemm, nIntegrationPad, 3, nDOFs, hesLagBasisInt[4], matCoor, matDer2ndCoor[4], nullptr);
  OwnGemm(my_dgemm, nIntegrationPad, 3, nDOFs, hesLagBasisInt[5], matCoor, matDer2ndCoor[5], nullptr);
}

void CFEMStandardPrismGrid::DerivativesCoorSolDOFs(ColMajorMatrix<su2double>          &matCoor,
                                                   vector<ColMajorMatrix<su2double> > &matDerCoor) {

  /*--- Call OwnGemm 3 times to compute the derivatives of the Cartesian
        coordinates w.r.t. the three parametric coordinates. ---*/
  const unsigned short nSolDOFs = lagBasisSolDOFs.rows();
  OwnGemm(dgemmSolDOFs, nSolDOFs, 3, nDOFs, derLagBasisSolDOFs[0], matCoor, matDerCoor[0], nullptr);
  OwnGemm(dgemmSolDOFs, nSolDOFs, 3, nDOFs, derLagBasisSolDOFs[1], matCoor, matDerCoor[1], nullptr);
  OwnGemm(dgemmSolDOFs, nSolDOFs, 3, nDOFs, derLagBasisSolDOFs[2], matCoor, matDerCoor[2], nullptr);
}

passivedouble CFEMStandardPrismGrid::WorkEstimateVolume(CConfig *config) {

  /*--- TEMPORARY IMPLEMENTATION. ---*/
  return nIntegration + 0.1*nDOFs;
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
  CSquareMatrixCM VInv(rDOFs.size());
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

void CFEMStandardPrismGrid::HesLagBasisIntPointsPrism(const vector<passivedouble> &rTriangleDOFs,
                                                      const vector<passivedouble> &sTriangleDOFs,
                                                      const vector<passivedouble> &rLineDOFs) {
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
  CSquareMatrixCM VInv(rDOFs.size());
  VandermondePrism(rDOFs, sDOFs, tDOFs, VInv.GetMat());
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

  HesVandermondePrism(rInt, sInt, tInt, VDr2, VDs2, VDt2, VDrs, VDrt, VDst);

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
  CSquareMatrixCM VInv(rDOFs.size());
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

void CFEMStandardPrismGrid::LagBasisAndDerSolDOFsPrism(const vector<passivedouble> &rTriangleDOFs,
                                                       const vector<passivedouble> &sTriangleDOFs,
                                                       const vector<passivedouble> &rLineDOFs) {

  /*--- Determine the parametric coordinates of all DOFs of the prism. ---*/
  vector<passivedouble> rDOFs, sDOFs, tDOFs;
  LocationAllDOFsPrism(rTriangleDOFs, sTriangleDOFs, rLineDOFs, rDOFs, sDOFs, tDOFs);

  /*--- Determine the parametric coordinates of all solution DOFs of the prism. ---*/
  const unsigned short nTri     = rTriangleSolDOFs.size();
  const unsigned short nLine    = rLineSolDOFs.size();
  const unsigned short nSolDOFs = nTri*nLine;
  vector<passivedouble> rSol(nSolDOFs), sSol(nSolDOFs), tSol(nSolDOFs);

  unsigned short ii = 0;
  for(unsigned short k=0; k<nLine; ++k) {
    for(unsigned short j=0; j<nTri; ++j, ++ii) {
      rSol[ii] = rTriangleSolDOFs[j];
      sSol[ii] = sTriangleSolDOFs[j];
      tSol[ii] = rLineSolDOFs[k];
    }
  }

  /*--- Determine the inverse of the Vandermonde matrix of the DOFs. ---*/
  CSquareMatrixCM VInv(rDOFs.size());
  VandermondePrism(rDOFs, sDOFs, tDOFs, VInv.GetMat());
  VInv.Invert();

  /*--- Determine the Vandermonde matrix and its gradients of the solution DOFs. ---*/
  ColMajorMatrix<passivedouble> V(nSolDOFs,rDOFs.size()),   VDr(nSolDOFs,rDOFs.size()),
                                VDs(nSolDOFs,rDOFs.size()), VDt(nSolDOFs,rDOFs.size());

  VandermondePrism(rSol, sSol, tSol, V);
  GradVandermondePrism(rSol, sSol, tSol, VDr, VDs, VDt);

  /*--- The Lagrangian basis functions and its gradients can be obtained by
        multiplying V, VDr, VDs, VDt and VInv. ---*/
  derLagBasisSolDOFs.resize(3);
  VInv.MatMatMult('R', V,   lagBasisSolDOFs);
  VInv.MatMatMult('R', VDr, derLagBasisSolDOFs[0]);
  VInv.MatMatMult('R', VDs, derLagBasisSolDOFs[1]);
  VInv.MatMatMult('R', VDt, derLagBasisSolDOFs[2]);

  /*--- Check if the sum of the elements of the relevant rows of
        lagBasisSolDOFs is 1 and of derLagBasisSolDOFs is 0. ---*/
  for(unsigned short i=0; i<nSolDOFs; ++i) {
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

void CFEMStandardPrismGrid::HesVandermondePrism(const vector<passivedouble>   &r,
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

void CFEMStandardPrismGrid::SubConnLinearElements(void) {

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
