/*!
 * \file CFEMStandardPyraGrid.cpp
 * \brief Functions for the class CFEMStandardPyraGrid.
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

#include "../../include/fem/CFEMStandardPyraGrid.hpp"
#include "../../include/toolboxes/CSquareMatrixCM.hpp"

/*----------------------------------------------------------------------------------*/
/*             Public member functions of CFEMStandardPyraGrid.                     */
/*----------------------------------------------------------------------------------*/

CFEMStandardPyraGrid::CFEMStandardPyraGrid(const unsigned short val_nPoly,
                                           const unsigned short val_orderExact)
  : CFEMStandardPyra(val_nPoly, val_orderExact) {

  /*--- Compute the values of the Lagrangian basis functions in the integration
        points for both the equidistant and LGL point distribution. ---*/
  LagBasisIntPointsPyra(rPyraDOFsEqui, sPyraDOFsEqui, tPyraDOFsEqui, lagBasisIntEqui);
  LagBasisIntPointsPyra(rPyraDOFsLGL,  sPyraDOFsLGL,  tPyraDOFsLGL,  lagBasisIntLGL);

  /*--- Compute the values of the derivatives of the Lagrangian basis functions in
        the integration points for both the equidistant and LGL point distribution. ---*/
  DerLagBasisIntPointsPyra(rPyraDOFsEqui, sPyraDOFsEqui, tPyraDOFsEqui, derLagBasisIntEqui);
  DerLagBasisIntPointsPyra(rPyraDOFsLGL,  sPyraDOFsLGL,  tPyraDOFsLGL,  derLagBasisIntLGL);

  /*--- Create the local grid connectivities of the faces of the volume element. ---*/
  LocalGridConnFaces();

  /*--- Determine the local subconnectivity of this standard element when split
        in several linear elements. Used for a.o. plotting and searcing. ---*/
  SubConnLinearElements();

  /*--- Set up the jitted gemm call, if supported. For this particular standard
        element the derivative of the coordinates are computed, which is 3. ---*/
  SetUpJittedGEMM(nIntegrationPad, 3, nDOFs, jitter, my_dgemm);
}

CFEMStandardPyraGrid::CFEMStandardPyraGrid(const unsigned short val_nPolyGrid,
                                           const unsigned short val_nPolySol,
                                           const unsigned short val_orderExact,
                                           const unsigned short val_locGridDOFs)
  : CFEMStandardPyra(val_nPolyGrid, val_orderExact) {

  /*--- Determine the location of the grid DOFs and build the appropriate
        Lagrangian basis functions. ---*/
  if(val_locGridDOFs == LGL) {

    /*--- LGL distribution. Compute the corresponding Lagrangian basis functions and
          its first and second derivatives in the integration points. ---*/
    LagBasisIntPointsPyra(rPyraDOFsLGL, sPyraDOFsLGL, tPyraDOFsLGL, lagBasisIntLGL);
    DerLagBasisIntPointsPyra(rPyraDOFsLGL, sPyraDOFsLGL, tPyraDOFsLGL, derLagBasisIntLGL);
    HesLagBasisIntPointsPyra(rPyraDOFsLGL, sPyraDOFsLGL, tPyraDOFsLGL);

    /*--- Determine the location of nodal solution DOFs. ---*/
    nPoly = val_nPolySol;
    LocationPyramidGridDOFsLGL(rPyraSolDOFs, sPyraSolDOFs, tPyraSolDOFs);
    nPoly = val_nPolyGrid;

    /*--- Compute the Lagrangian basis functions and its derivatives
          in the nodal solution DOFs. ---*/
    LagBasisAndDerSolDOFsPyra(rPyraDOFsLGL, sPyraDOFsLGL, tPyraDOFsLGL);
  }
  else {

    /*--- Equidistant distribution. Compute the corresponding Lagrangian basis functions and
          its first and second derivatives in the integration points. ---*/
    LagBasisIntPointsPyra(rPyraDOFsEqui, sPyraDOFsEqui, tPyraDOFsEqui, lagBasisIntEqui);
    DerLagBasisIntPointsPyra(rPyraDOFsEqui, sPyraDOFsEqui, tPyraDOFsEqui, derLagBasisIntEqui);
    HesLagBasisIntPointsPyra(rPyraDOFsEqui, sPyraDOFsEqui, tPyraDOFsEqui);

    /*--- Determine the location of nodal solution DOFs. ---*/
    nPoly = val_nPolySol;
    LocationPyramidGridDOFsEquidistant(rPyraSolDOFs, sPyraSolDOFs, tPyraSolDOFs);
    nPoly = val_nPolyGrid;

    /*--- Compute the Lagrangian basis functions and its derivatives
          in the nodal solution DOFs. ---*/
    LagBasisAndDerSolDOFsPyra(rPyraDOFsEqui, sPyraDOFsEqui, tPyraDOFsEqui);
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

CFEMStandardPyraGrid::~CFEMStandardPyraGrid() {

#if defined(PRIMAL_SOLVER) && defined(HAVE_MKL)
  if( jitterSolDOFs ) {
    mkl_jit_destroy(jitterSolDOFs);
    jitterSolDOFs = nullptr;
  }
#endif
}

void CFEMStandardPyraGrid::CoorIntPoints(const bool                LGLDistribution,
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

void CFEMStandardPyraGrid::DerivativesCoorIntPoints(const bool                         LGLDistribution,
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

void CFEMStandardPyraGrid::Derivatives2ndCoorIntPoints(const bool                         LGLDistribution,
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

void CFEMStandardPyraGrid::DerivativesCoorSolDOFs(ColMajorMatrix<su2double>          &matCoor,
                                                  vector<ColMajorMatrix<su2double> > &matDerCoor) {

  /*--- Call OwnGemm 3 times to compute the derivatives of the Cartesian
        coordinates w.r.t. the three parametric coordinates. ---*/
  const unsigned short nSolDOFs = lagBasisSolDOFs.rows();
  OwnGemm(dgemmSolDOFs, nSolDOFs, 3, nDOFs, derLagBasisSolDOFs[0], matCoor, matDerCoor[0], nullptr);
  OwnGemm(dgemmSolDOFs, nSolDOFs, 3, nDOFs, derLagBasisSolDOFs[1], matCoor, matDerCoor[1], nullptr);
  OwnGemm(dgemmSolDOFs, nSolDOFs, 3, nDOFs, derLagBasisSolDOFs[2], matCoor, matDerCoor[2], nullptr);
}

passivedouble CFEMStandardPyraGrid::WorkEstimateVolume(CConfig *config) {

  /*--- TEMPORARY IMPLEMENTATION. ---*/
  return nIntegration + 0.1*nDOFs;
}

/*----------------------------------------------------------------------------------*/
/*             Private member functions of CFEMStandardPyraGrid.                    */
/*----------------------------------------------------------------------------------*/

void CFEMStandardPyraGrid::DerLagBasisIntPointsPyra(const vector<passivedouble>            &rDOFs,
                                                    const vector<passivedouble>            &sDOFs,
                                                    const vector<passivedouble>            &tDOFs,
                                                    vector<ColMajorMatrix<passivedouble> > &derLag) {

  /*--- Determine the parametric coordinates of all integration points
        of the pyramid. ---*/
  vector<passivedouble> rInt, sInt, tInt;
  LocationAllIntegrationPoints(rInt, sInt, tInt);

  /*--- Determine the padded number of the total number of integration points. ---*/
  const unsigned short nIntTot    = rInt.size();
  const unsigned short nIntTotPad = ((nIntTot+baseVectorLen-1)/baseVectorLen)*baseVectorLen;

  /*--- Determine the inverse of the Vandermonde matrix of the DOFs. ---*/
  CSquareMatrixCM VInv(rDOFs.size());
  VandermondePyramid(rDOFs, sDOFs, tDOFs, VInv.GetMat());
  VInv.Invert();

  /*--- Determine the gradient of the Vandermonde matrix of the integration points. Make
        sure to allocate the number of rows to nIntTotPad and initialize them to zero. ---*/
  ColMajorMatrix<passivedouble> VDr(nIntTotPad,rDOFs.size()),
                                VDs(nIntTotPad,rDOFs.size()),
                                VDt(nIntTotPad,rDOFs.size());
  VDr.setConstant(0.0);
  VDs.setConstant(0.0);
  VDt.setConstant(0.0);

  GradVandermondePyramid(rInt, sInt, tInt, VDr, VDs, VDt);

  /*--- TEST ----*/
/*const passivedouble dr = 1.e-4;
  vector<passivedouble> rMin = rInt, rPlus = rInt;
  vector<passivedouble> sMin = sInt, sPlus = sInt;
  vector<passivedouble> tMin = tInt, tPlus = tInt;

  for(unsigned short i=0; i<rInt.size(); ++i) {
    rMin[i] -= dr; rPlus[i] += dr; 
    sMin[i] -= dr; sPlus[i] += dr;
    tMin[i] -= dr; tPlus[i] += dr;
  }

  ColMajorMatrix<passivedouble> VMin(nIntTotPad,rDOFs.size()),
                                VPlus(nIntTotPad,rDOFs.size());

  VMin.setConstant(0.0);
  VPlus.setConstant(0.0);

  VandermondePyramid(rInt, sInt, tMin,  VMin);
  VandermondePyramid(rInt, sInt, tPlus, VPlus);

  passivedouble maxDiff = 0.0;
  for(unsigned i=0; i<rDOFs.size(); ++i) {
    for(unsigned short j=0; j<rInt.size(); ++j) {

      const passivedouble diff = (VPlus(j,i) - VMin(j,i))/(2.0*dr);
      cout << "t comparison: " << i << " " << j << ": "
           << VDt(j,i) << " " << diff << " " << fabs(VDt(j,i)-diff) << endl;
      maxDiff = max(maxDiff, fabs(VDt(j,i)-diff));
    }
  }

  cout << endl << "maxDiff: " << maxDiff << endl << endl;
  SU2_MPI::Error("Test", CURRENT_FUNCTION); */

  /*--- End TEST. ---*/

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

void CFEMStandardPyraGrid::HesLagBasisIntPointsPyra(const vector<passivedouble> &rDOFs,
                                                    const vector<passivedouble> &sDOFs,
                                                    const vector<passivedouble> &tDOFs) {

  /*--- Determine the parametric coordinates of all integration points
        of the pyramid. ---*/
  vector<passivedouble> rInt, sInt, tInt;
  LocationAllIntegrationPoints(rInt, sInt, tInt);

  /*--- Determine the padded number of the total number of integration points. ---*/
  const unsigned short nIntTot    = rInt.size();
  const unsigned short nIntTotPad = ((nIntTot+baseVectorLen-1)/baseVectorLen)*baseVectorLen;

  /*--- Determine the inverse of the Vandermonde matrix of the DOFs. ---*/
  CSquareMatrixCM VInv(rDOFs.size());
  VandermondePyramid(rDOFs, sDOFs, tDOFs, VInv.GetMat());
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

  HesVandermondePyramid(rInt, sInt, tInt, VDr2, VDs2, VDt2, VDrs, VDrt, VDst);

  /*--- TEST ----*/
/*const passivedouble dr = 1.e-4;
  vector<passivedouble> rMin = rInt, rPlus = rInt;
  vector<passivedouble> sMin = sInt, sPlus = sInt;
  vector<passivedouble> tMin = tInt, tPlus = tInt;

  for(unsigned short i=0; i<rInt.size(); ++i) {
    rMin[i] -= dr; rPlus[i] += dr; 
    sMin[i] -= dr; sPlus[i] += dr;
    tMin[i] -= dr; tPlus[i] += dr;
  } */

/*ColMajorMatrix<passivedouble> VMin(nIntTotPad,rDOFs.size()),
                                V(nIntTotPad,rDOFs.size()),
                                VPlus(nIntTotPad,rDOFs.size());

  VMin.setConstant(0.0);
  V.setConstant(0.0);
  VPlus.setConstant(0.0);

  VandermondePyramid(rInt, sInt, tMin,  VMin);
  VandermondePyramid(rInt, sInt, tInt,  V);
  VandermondePyramid(rInt, sInt, tPlus, VPlus);

  passivedouble maxDiff = 0.0;
  for(unsigned i=0; i<rDOFs.size(); ++i) {
    for(unsigned short j=0; j<rInt.size(); ++j) {

      const passivedouble d2t = (VMin(j,i) - 2.0*V(j,i) + VPlus(j,i))/(dr*dr);
      cout << "t2 comparison: " << i << " " << j << ": "
           << VDt2(j,i) << " " << d2t << " " << fabs(VDt2(j,i)-d2t) << endl;
      maxDiff = max(maxDiff, fabs(VDt2(j,i)-d2t));
    }
  }

  cout << endl << "maxDiff: " << maxDiff << endl << endl; */

/*ColMajorMatrix<passivedouble> VMinMin(nIntTotPad,rDOFs.size()),
                                VMinPlus(nIntTotPad,rDOFs.size()),
                                VPlusMin(nIntTotPad,rDOFs.size()),
                                VPlusPlus(nIntTotPad,rDOFs.size());

  VMinMin.setConstant(0.0);
  VMinPlus.setConstant(0.0);
  VPlusMin.setConstant(0.0);
  VPlusPlus.setConstant(0.0);

  VandermondePyramid(rInt, sMin,  tMin,  VMinMin);
  VandermondePyramid(rInt, sMin,  tPlus, VMinPlus);
  VandermondePyramid(rInt, sPlus, tMin,  VPlusMin);
  VandermondePyramid(rInt, sPlus, tPlus, VPlusPlus);

  passivedouble maxDiff = 0.0;
  for(unsigned i=0; i<rDOFs.size(); ++i) {
    for(unsigned short j=0; j<rInt.size(); ++j) {

      const passivedouble dst = (VPlusPlus(j,i) - VMinPlus(j,i) - VPlusMin(j,i) + VMinMin(j,i))/(4.0*dr*dr);
      cout << "st comparison: " << i << " " << j << ": "
           << VDst(j,i) << " " << dst << " " << fabs(VDst(j,i)-dst) << endl;
      maxDiff = max(maxDiff, fabs(VDst(j,i)-dst));
    }
  }

  cout << endl << "maxDiff: " << maxDiff << endl << endl; */
/*SU2_MPI::Error("Test", CURRENT_FUNCTION); */

  /*--- End TEST. ---*/

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

void CFEMStandardPyraGrid::LagBasisIntPointsPyra(const vector<passivedouble>   &rDOFs,
                                                 const vector<passivedouble>   &sDOFs,
                                                 const vector<passivedouble>   &tDOFs,
                                                 ColMajorMatrix<passivedouble> &lag) {

  /*--- Determine the parametric coordinates of all integration points
        of the pyramid. ---*/
  vector<passivedouble> rInt, sInt, tInt;
  LocationAllIntegrationPoints(rInt, sInt, tInt);

  /*--- Determine the padded number of the total number of integration points. ---*/
  const unsigned short nIntTot    = rInt.size();
  const unsigned short nIntTotPad = ((nIntTot+baseVectorLen-1)/baseVectorLen)*baseVectorLen;

  /*--- Determine the inverse of the Vandermonde matrix of the DOFs. ---*/
  CSquareMatrixCM VInv(rDOFs.size());
  VandermondePyramid(rDOFs, sDOFs, tDOFs, VInv.GetMat());
  VInv.Invert();

  /*--- Determine the Vandermonde matrix of the integration points. Make sure to
        allocate the number of rows to nIntTotPad and initialize them to zero. ---*/ 
  ColMajorMatrix<passivedouble> V(nIntTotPad,rDOFs.size());
  V.setConstant(0.0);
  VandermondePyramid(rInt, sInt, tInt, V);

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

void CFEMStandardPyraGrid::LagBasisAndDerSolDOFsPyra(const vector<passivedouble>  &rDOFs,
                                                     const vector<passivedouble>  &sDOFs,
                                                     const vector<passivedouble>  &tDOFs) {

  /*--- Determine the inverse of the Vandermonde matrix of the DOFs. ---*/
  CSquareMatrixCM VInv(rDOFs.size());
  VandermondePyramid(rDOFs, sDOFs, tDOFs, VInv.GetMat());
  VInv.Invert();

  /*--- Determine the Vandermonde matrix and its gradients of the solution DOFs. ---*/
  ColMajorMatrix<passivedouble> V(rPyraSolDOFs.size(),rDOFs.size()),
                                VDr(rPyraSolDOFs.size(),rDOFs.size()),
                                VDs(rPyraSolDOFs.size(),rDOFs.size()),
                                VDt(rPyraSolDOFs.size(),rDOFs.size());

  VandermondePyramid(rPyraSolDOFs, sPyraSolDOFs, tPyraSolDOFs, V);
  GradVandermondePyramid(rPyraSolDOFs, sPyraSolDOFs, tPyraSolDOFs, VDr, VDs, VDt);

  /*--- The Lagrangian basis functions and its gradients can be obtained by
        multiplying V, VDr, VDs, VDt and VInv. ---*/
  derLagBasisSolDOFs.resize(3);
  VInv.MatMatMult('R', V,   lagBasisSolDOFs);
  VInv.MatMatMult('R', VDr, derLagBasisSolDOFs[0]);
  VInv.MatMatMult('R', VDs, derLagBasisSolDOFs[1]);
  VInv.MatMatMult('R', VDt, derLagBasisSolDOFs[2]);

  /*--- Check if the sum of the elements of the relevant rows of
        lagBasisSolDOFs is 1 and of derLagBasisSolDOFs is 0. ---*/
  for(unsigned short i=0; i<rPyraSolDOFs.size(); ++i) {
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

void CFEMStandardPyraGrid::GradVandermondePyramid(const vector<passivedouble>   &r,
                                                  const vector<passivedouble>   &s,
                                                  const vector<passivedouble>   &t,
                                                  ColMajorMatrix<passivedouble> &VDr,
                                                  ColMajorMatrix<passivedouble> &VDs,
                                                  ColMajorMatrix<passivedouble> &VDt) {

  /*--- For a pyramid the orthogonal basis for the reference element is
        obtained by a combination of Jacobi polynomials (of which the Legendre
        polynomials is a special case). This is the result of the
        orthonormalization of the monomial basis.
        Note that the sequence of the i, j and k loop must be identical to
        the evaluation of the Vandermonde matrix itself.  ---*/
  unsigned short ii = 0;
  for(unsigned short i=0; i<=nPoly; ++i) {
    for(unsigned short j=0; j<=nPoly; ++j) {
      unsigned short muij = max(i,j);
      for(unsigned short k=0; k<=(nPoly-muij); ++k, ++ii) {
        for(unsigned short l=0; l<r.size(); ++l) {

          /*--- Determine the coefficients a, b, c and d. ---*/
          passivedouble a, b;
          const passivedouble d = 1.0-t[l];
          if(fabs(d) < 1.e-8) a = b = 0.0;
          else {
            a = 2.0*r[l]/d;
            b = 2.0*s[l]/d;
          }

          const passivedouble c = t[l];

          /*--- Determine the value of the three 1D contributions to the 
                3D basis functions as well as the gradients of these
                contributions functions w.r.t. to their arguments. ---*/
          const passivedouble fa  = NormJacobi(i,0,         0,a);
          const passivedouble gb  = NormJacobi(j,0,         0,b);
          const passivedouble hc  = NormJacobi(k,2*(muij+1),0,c);
          const passivedouble dfa = GradNormJacobi(i,0,         0,a);
          const passivedouble dgb = GradNormJacobi(j,0,         0,b);
          const passivedouble dhc = GradNormJacobi(k,2*(muij+1),0,c);

          /*--- Computation of the powers of d that occur in the expressions for
                the gradients. Note the safeguard to avoid division by zero. This is
                allowed, because the implementation is such that when this clipping
                is active, it is multiplied by zero. ---*/
          const passivedouble tmpdmu   = muij > 0 ? pow(d, muij)   : 1.0;
          const passivedouble tmpdmum1 = muij > 1 ? pow(d, muij-1) : 1.0;

          /*--- Compute the derivatives of the basis function. ---*/
          VDr(l,ii) = 2.0*tmpdmum1* 2.0*dfa*gb*hc;
          VDs(l,ii) = 2.0*tmpdmum1* 2.0*fa*dgb*hc;
          VDt(l,ii) = 2.0*tmpdmum1*(a*dfa*gb*hc + b*fa*dgb*hc - muij*fa*gb*hc)
                    + 2.0*tmpdmu  * fa*gb*dhc;
        }
      }
    }
  }
}

void CFEMStandardPyraGrid::HesVandermondePyramid(const vector<passivedouble>   &r,
                                                 const vector<passivedouble>   &s,
                                                 const vector<passivedouble>   &t,
                                                 ColMajorMatrix<passivedouble> &VDr2,
                                                 ColMajorMatrix<passivedouble> &VDs2,
                                                 ColMajorMatrix<passivedouble> &VDt2,
                                                 ColMajorMatrix<passivedouble> &VDrs,
                                                 ColMajorMatrix<passivedouble> &VDrt,
                                                 ColMajorMatrix<passivedouble> &VDst) {

  /*--- For a pyramid the orthogonal basis for the reference element is
        obtained by a combination of Jacobi polynomials (of which the Legendre
        polynomials is a special case). This is the result of the
        orthonormalization of the monomial basis.
        Note that the sequence of the i, j and k loop must be identical to
        the evaluation of the Vandermonde matrix itself.  ---*/
  unsigned short ii = 0;
  for(unsigned short i=0; i<=nPoly; ++i) {
    for(unsigned short j=0; j<=nPoly; ++j) {
      unsigned short muij = max(i,j);
      for(unsigned short k=0; k<=(nPoly-muij); ++k, ++ii) {
        for(unsigned short l=0; l<r.size(); ++l) {

          /*--- Determine the coefficients a, b, c and d. ---*/
          const passivedouble d = 1.0-t[l];
          const passivedouble dInv = (fabs(d) < 1.e-8) ? 1.e+8 : 1.0/d; 
          const passivedouble a = 2.0*r[l]*dInv;
          const passivedouble b = 2.0*s[l]*dInv;
          const passivedouble c = t[l];

          /*--- Determine the value of the three 1D contributions to the 
                3D basis functions as well as the 1st and 2nd derivative of
                these contributions functions w.r.t. to their arguments. ---*/
          const passivedouble fa   = NormJacobi(i,0,         0,a);
          const passivedouble gb   = NormJacobi(j,0,         0,b);
          const passivedouble hc   = NormJacobi(k,2*(muij+1),0,c);
          const passivedouble dfa  = GradNormJacobi(i,0,         0,a);
          const passivedouble dgb  = GradNormJacobi(j,0,         0,b);
          const passivedouble dhc  = GradNormJacobi(k,2*(muij+1),0,c);
          const passivedouble d2fa = HesNormJacobi(i,0,         0,a);
          const passivedouble d2gb = HesNormJacobi(j,0,         0,b);
          const passivedouble d2hc = HesNormJacobi(k,2*(muij+1),0,c);

          /*--- Computation of the powers of d that occur in the expressions for
                the gradients. Note the safeguard to avoid division by zero. This is
                allowed, because the implementation is such that when this clipping
                is active, it is multiplied by zero. However, for the pyramid with
                its rational basis functions there are some terms in which a division
                by d cannot be avoided. This is handled via dInv. ---*/
          const passivedouble tmpdmu   = muij > 0 ? pow(d, muij)   : 1.0;
          const passivedouble tmpdmum1 = muij > 1 ? pow(d, muij-1) : 1.0;
          const passivedouble tmpdmum2 = muij > 2 ? pow(d, muij-2) : 1.0;

          /*--- Compute the 2nd derivative w.r.t. to r. ---*/
          VDr2(l,ii) = tmpdmum2*4.0*d2fa*gb*hc;

          /*--- Compute the 2nd derivative w.r.t. to s. ---*/
          VDs2(l,ii) = tmpdmum2*4.0*fa*d2gb*hc;

          /*--- Compute the 2nd derivative w.r.t. to r. ---*/
          VDt2(l,ii) = tmpdmum2*(muij*(muij-1)*fa*gb*hc
                     -           2.0*a*(muij-1)*dfa*gb*hc
                     -           2.0*b*(muij-1)*fa*dgb*hc
                     +           a*a*d2fa*gb*hc
                     +           b*b*fa*d2gb*hc)
                     + tmpdmum1*(2.0*a*dfa*gb*dhc
                     +           2.0*b*fa*dgb*dhc
                     -           2.0*muij*fa*gb*dhc)
                     + tmpdmu  *fa*gb*d2hc
                     + tmpdmum1*2.0*a*b*dfa*dgb*hc*dInv;

          /*--- Compute the cross derivative w.r.t. r and s. ---*/
          VDrs(l,ii) = tmpdmum1*4.0*dfa*dgb*hc*dInv;

          /*--- Compute the cross derivative w.r.t. r and t. ---*/
          VDrt(l,ii) = tmpdmum2*2.0*(a*d2fa*gb*hc - (muij-1)*dfa*gb*hc)
                     + tmpdmum1*2.0*(dfa*gb*dhc + b*dfa*dgb*hc*dInv);

          /*--- Compute the cross derivative w.r.t. s and t. ---*/
          VDst(l,ii) = tmpdmum2*2.0*(b*fa*d2gb*hc - (muij-1)*fa*dgb*hc)
                     + tmpdmum1*2.0*(fa*dgb*dhc + a*dfa*dgb*hc*dInv);

          /*--- Multiply all 2nd derivatives with 2 to obtain the
                correct expressions. ---*/
          VDr2(l,ii) *= 2.0;
          VDs2(l,ii) *= 2.0;
          VDt2(l,ii) *= 2.0;
          VDrs(l,ii) *= 2.0;
          VDrt(l,ii) *= 2.0;
          VDst(l,ii) *= 2.0;
        }
      }
    }
  }
}

void CFEMStandardPyraGrid::VandermondePyramid(const vector<passivedouble>   &r,
                                              const vector<passivedouble>   &s,
                                              const vector<passivedouble>   &t,
                                              ColMajorMatrix<passivedouble> &V) {

  /*--- For a pyramid the orthogonal basis for the reference element is
        obtained by a combination of Jacobi polynomials (of which the Legendre
        polynomials is a special case). This is the result of the
        orthonormalization of the monomial basis. ---*/
  unsigned short ii = 0;
  for(unsigned short i=0; i<=nPoly; ++i) {
    for(unsigned short j=0; j<=nPoly; ++j) {
      unsigned short muij = max(i,j);
      for(unsigned short k=0; k<=(nPoly-muij); ++k, ++ii) {
        for(unsigned short l=0; l<r.size(); ++l) {

          /*--- Determine the coefficients a, b, c and d. ---*/
          passivedouble a, b;
          const passivedouble d = 1.0-t[l];
          if(fabs(d) < 1.e-8) a = b = 0.0;
          else {
            a = 2.0*r[l]/d;
            b = 2.0*s[l]/d;
          }

          const passivedouble c = t[l];

          /*--- Determine the value of the current basis function in this point. ---*/
          V(l,ii) = 2.0*pow(d,muij)*NormJacobi(i,0,0,a)*NormJacobi(j,0,0,b)
                  * NormJacobi(k,2*(muij+1),0,c);
        }
      }
    }
  }
}

void CFEMStandardPyraGrid::LocalGridConnFaces(void) {

  /*--- Allocate the first index of gridConnFaces, which is equal to the number
        of faces of the pyramid, which is 5. Reserve memory for the second
        index afterwards. ---*/
  const unsigned short nDOFsQuad     = (nPoly+1)*(nPoly+1);
  const unsigned short nDOFsTriangle = (nPoly+1)*(nPoly+2)/2;
  gridConnFaces.resize(5);

  gridConnFaces[0].reserve(nDOFsQuad);
  gridConnFaces[1].reserve(nDOFsTriangle);
  gridConnFaces[2].reserve(nDOFsTriangle);
  gridConnFaces[3].reserve(nDOFsTriangle);
  gridConnFaces[4].reserve(nDOFsTriangle);

  /*--- Loop over all the nodes of the pyramid and pick the correct
        ones for the faces. ---*/
  unsigned short mPoly = nPoly;
  unsigned int ii = 0;
  for(unsigned short k=0; k<=nPoly; ++k, --mPoly) {
    for(unsigned short j=0; j<=mPoly; ++j) {
      for(unsigned short i=0; i<=mPoly; ++i, ++ii) {
        if(k == 0)     gridConnFaces[0].push_back(ii);
        if(j == 0)     gridConnFaces[1].push_back(ii);
        if(j == mPoly) gridConnFaces[2].push_back(ii);
        if(i == 0)     gridConnFaces[3].push_back(ii);
        if(i == mPoly) gridConnFaces[4].push_back(ii);
      }
    }
  }

  /*--- Make sure that the element is to the left of the faces. ---*/
  const unsigned short n0 = 0;
  const unsigned short n1 = nPoly;
  const unsigned short n2 = nDOFsQuad -1;
  const unsigned short n3 = n2 - nPoly;
  const unsigned short n4 = nDOFs -1;

  ChangeDirectionQuadConn(gridConnFaces[0], n0, n1, n2, n3);
  ChangeDirectionTriangleConn(gridConnFaces[1], n0, n4, n1);
  ChangeDirectionTriangleConn(gridConnFaces[2], n3, n2, n4);
  ChangeDirectionTriangleConn(gridConnFaces[3], n0, n3, n4);
  ChangeDirectionTriangleConn(gridConnFaces[4], n1, n4, n2);
}

void CFEMStandardPyraGrid::SubConnLinearElements(void) {

  /*--- The pyramid is split into several linear pyramids and tets.
        Set the VTK sub-types accordingly. ---*/
  VTK_SubType1 = PYRAMID;
  VTK_SubType2 = TETRAHEDRON;

  /*--- Initialize the number of DOFs for the current edges to the number of
        DOFs of the edges on the base of the pyramid. Also initialize the
        current k offset to zero.     ---*/
  unsigned short nDOFsCurrentEdges = nPoly + 1;
  unsigned short offCurrentK       = 0;

  /*--- Loop in the k-direction of the pyramid. ---*/
  for(unsigned short k=0; k<nPoly; ++k) {

    /*------------------------------------------------------------------------*/
    /*       Sub-pyramids in the same direction as the original pyramid.      */
    /*------------------------------------------------------------------------*/

    /*--- Determine the index of the first vertex of the quadrilateral of the
          next k value. ---*/
    unsigned short kk = offCurrentK + nDOFsCurrentEdges*nDOFsCurrentEdges;

    /*--- Loop in j-direction of the current quad. ---*/
    for(unsigned short j=0; j<(nDOFsCurrentEdges-1); ++j) {

      /*--- Index of the first vertex along the j-row of the current quad. ---*/
      const unsigned short jj = offCurrentK + j*nDOFsCurrentEdges;

      /*--- Loop in i-direction of the current quad. ---*/
      for(unsigned short i=0; i<(nDOFsCurrentEdges-1); ++i) {

        /*--- Determine the local indices of the corners of the quadrilateral
              of this subpyramid as well as the top of the subpyramid.
              Store the connectivity in subConn1ForPlotting.  ---*/
        const unsigned short n0 = jj + i;
        const unsigned short n1 = n0 + 1;
        const unsigned short n2 = n1 + nDOFsCurrentEdges;
        const unsigned short n3 = n0 + nDOFsCurrentEdges;
        const unsigned short n4 = kk + i;

        subConn1ForPlotting.push_back(n0);
        subConn1ForPlotting.push_back(n1);
        subConn1ForPlotting.push_back(n2);
        subConn1ForPlotting.push_back(n3);
        subConn1ForPlotting.push_back(n4);
      }

      /*--- Update kk for the next j-row. ---*/
      kk += nDOFsCurrentEdges - 1;
    }

    /*------------------------------------------------------------------------*/
    /*    Sub-pyramids in the opposite direction as the original pyramid.     */
    /*------------------------------------------------------------------------*/

    /*--- Reset the value of kk to the index of the first vertex of the
          quadrilateral of the next k-plane. ---*/
    kk = offCurrentK + nDOFsCurrentEdges*nDOFsCurrentEdges;

    /*--- Loop in j-direction of the current quad. Note that the starting index
          of this loop is 1. ---*/
    for(unsigned short j=1; j<(nDOFsCurrentEdges-1); ++j) {

      /*--- Index of the first vertex along the j-row of the current quad. ---*/
      const unsigned short jj = offCurrentK + j*nDOFsCurrentEdges;

      /*--- Loop in the i-direction of this quad. Again the starting index is 1. ---*/
      for(unsigned short i=1; i<(nDOFsCurrentEdges-1); ++i) {

        /*--- Determine the local indices of the corners of the quadrilateral
              of this subpyramid as well as the top of the subpyramid.  ---*/
        const unsigned short n0 = kk + i - 1;
        const unsigned short n1 = n0 + 1;
        const unsigned short n2 = n1 + nDOFsCurrentEdges-1;
        const unsigned short n3 = n0 + nDOFsCurrentEdges-1;
        const unsigned short n4 = jj + i;

        /*--- Store the connectivity of this subpyramid. Note that n1 and n3 are
              swapped, such that a positive volume is obtained according to the
              right hand rule. ---*/
        subConn1ForPlotting.push_back(n0);
        subConn1ForPlotting.push_back(n3);
        subConn1ForPlotting.push_back(n2);
        subConn1ForPlotting.push_back(n1);
        subConn1ForPlotting.push_back(n4);
      }

      /*--- Update kk for the next j-row. ---*/
      kk += nDOFsCurrentEdges - 1;
    }

    /*------------------------------------------------------------------------*/
    /*                   Sub-tetrahedra in the j-direction.                   */
    /*------------------------------------------------------------------------*/

    /*--- Reset the value of kk again. ---*/
    kk = offCurrentK + nDOFsCurrentEdges*nDOFsCurrentEdges;

    /*--- Loop in the i-direction of the current quad. Note that the starting
          index must be 1. ---*/
    for(unsigned short i=1; i<(nDOFsCurrentEdges-1); ++i) {

      /*--- Loop in the j-direction of the current quad. This loop starts at 0. ---*/
      for(unsigned short j=0; j<(nDOFsCurrentEdges-1); ++j) {

        /*--- Determine the local indices of the 4 corner points of this tet.
              Its connectivity is stored in subConn2ForPlotting, because in
              subConn1ForPlotting the pyramids are stored. ---*/
        const unsigned short n0 = kk + i-1 + j*(nDOFsCurrentEdges-1);
        const unsigned short n1 = n0 + 1;
        const unsigned short n2 = offCurrentK + i + j*nDOFsCurrentEdges;
        const unsigned short n3 = n2 + nDOFsCurrentEdges;

        subConn2ForPlotting.push_back(n0);
        subConn2ForPlotting.push_back(n1);
        subConn2ForPlotting.push_back(n2);
        subConn2ForPlotting.push_back(n3);
      }
    }

    /*------------------------------------------------------------------------*/
    /*                Sub-tetrahedra in the i-direction.                      */
    /*------------------------------------------------------------------------*/

    /*--- Loop in the j-direction of the current quad. Note that the starting
          index must be 1. ---*/
    for(unsigned short j=1; j<(nDOFsCurrentEdges-1); ++j) {

      /*--- Index of the first vertex along the j-row of the current quad. ---*/
      const unsigned short jj = offCurrentK + j*nDOFsCurrentEdges;

      /*--- Loop in the i-direction of the current quad. This loop starts at 0. ---*/
      for(unsigned short i=0; i<(nDOFsCurrentEdges-1); ++i) {

        /*--- Determine the local indices of the 4 corner points of this tet
              and store its connectivity in subConn2ForPlotting.   ---*/
        const unsigned short n0 = kk + i;
        const unsigned short n1 = jj + i;
        const unsigned short n2 = n0 + nDOFsCurrentEdges-1;
        const unsigned short n3 = n1 + 1;

        subConn2ForPlotting.push_back(n0);
        subConn2ForPlotting.push_back(n1);
        subConn2ForPlotting.push_back(n2);
        subConn2ForPlotting.push_back(n3);
      }

      /*--- Update kk for the next j-row. ---*/
      kk += nDOFsCurrentEdges - 1;
    }

    /*--- Update the value of offCurrentK with the amounts of DOFs present in the
          current quadrilateral plane and decrement nDOFsCurrentEdges, such that it
          contains the number of DOFs along an edge of the next quadrilateral plane. ---*/
    offCurrentK += nDOFsCurrentEdges*nDOFsCurrentEdges;
    --nDOFsCurrentEdges;
  }
}
