/*!
 * \file CFEMStandardTetGrid.cpp
 * \brief Functions for the class CFEMStandardTetGrid.
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

#include "../../include/fem/CFEMStandardTetGrid.hpp"

/*----------------------------------------------------------------------------------*/
/*             Public member functions of CFEMStandardTetGrid.                      */
/*----------------------------------------------------------------------------------*/

CFEMStandardTetGrid::CFEMStandardTetGrid(const unsigned short val_nPolyGrid,
                                         const unsigned short val_nPolySol,
                                         const unsigned short val_orderExact,
                                         const unsigned short val_locGridDOFs)
  : CFEMStandardTetBase(val_nPolyGrid, val_orderExact) {

  /*--- Determine the location of the grid and nodal solution DOFs. ---*/
  if(val_locGridDOFs == LGL) {
    LocationTetGridDOFsLGL(nPoly, rTetDOFs, sTetDOFs, tTetDOFs);
    LocationTetGridDOFsLGL(val_nPolySol, rTetSolDOFs, sTetSolDOFs, tTetSolDOFs);
  }
  else {
    LocationTetGridDOFsEquidistant(nPoly, rTetDOFs, sTetDOFs, tTetDOFs);
    LocationTetGridDOFsEquidistant(val_nPolySol, rTetSolDOFs, sTetSolDOFs, tTetSolDOFs);
  }

  /*--- Compute the corresponding Lagrangian basis functions and
        its first and second derivatives in the integration points. ---*/
  LagBasisIntPointsTet(nPoly, rTetDOFs, sTetDOFs, tTetDOFs,
                       rTetInt, sTetInt, tTetInt, lagBasisInt);
  DerLagBasisIntPointsTet(nPoly, rTetDOFs, sTetDOFs, tTetDOFs, 
                          rTetInt, sTetInt, tTetInt, derLagBasisInt);
  HesLagBasisIntPointsTet(nPoly, rTetDOFs, sTetDOFs, tTetDOFs,
                          rTetInt, sTetInt, tTetInt, hesLagBasisInt);

  /*--- Compute the Lagrangian basis functions and its derivatives
        in the nodal solution DOFs. Use LagBasisIntPointsTet and
        DerLagBasisIntPointsTet with different arguments. ---*/
  LagBasisIntPointsTet(nPoly, rTetDOFs, sTetDOFs, tTetDOFs, rTetSolDOFs,
                       sTetSolDOFs, tTetSolDOFs, lagBasisSolDOFs);
  DerLagBasisIntPointsTet(nPoly, rTetDOFs, sTetDOFs, tTetDOFs, rTetSolDOFs,
                          sTetSolDOFs, tTetSolDOFs, derLagBasisSolDOFs);

  /*--- Determine the local subconnectivity of this standard element when split
        in several linear elements. Used for a.o. plotting and searcing. ---*/
  CFEMStandardTetBase::SubConnLinearElements();

  /*--- Set up the jitted gemm calls, if supported. ---*/
  SetUpJittedGEMM(nIntegrationPad, 3, nDOFs, nIntegrationPad, nDOFs,
                  nIntegrationPad, jitterDOFs2Int, gemmDOFs2Int);
  SetUpJittedGEMM(lagBasisSolDOFs.rows(), 3, nDOFs, lagBasisSolDOFs.rows(), nDOFs,
                  lagBasisSolDOFs.rows(), jitterDOFs2SolDOFs, gemmDOFs2SolDOFs);
}

CFEMStandardTetGrid::~CFEMStandardTetGrid() {

#if defined(PRIMAL_SOLVER) && defined(HAVE_MKL)
  if( jitterDOFs2Int ) {
    mkl_jit_destroy(jitterDOFs2Int);
    jitterDOFs2Int = nullptr;
  }

  if( jitterDOFs2SolDOFs ) {
    mkl_jit_destroy(jitterDOFs2SolDOFs);
    jitterDOFs2SolDOFs = nullptr;
  }  
#endif
}

void CFEMStandardTetGrid::CoorIntPoints(const bool                notUsed,
                                        ColMajorMatrix<su2double> &matCoorDOF,
                                        ColMajorMatrix<su2double> &matCoorInt) {

  /*--- Call OwnGemm with the appropriate arguments to compute the data. ---*/
  OwnGemm(gemmDOFs2Int, jitterDOFs2Int, nIntegrationPad, 3, nDOFs,
          nIntegrationPad, nDOFs, nIntegrationPad,
          lagBasisInt, matCoorDOF, matCoorInt, nullptr);
}

void CFEMStandardTetGrid::DerivativesCoorIntPoints(const bool                         notUsed,
                                                   ColMajorMatrix<su2double>          &matCoor,
                                                   vector<ColMajorMatrix<su2double> > &matDerCoor) {


  /*--- Call OwnGemm with the appropriate arguments to compute the data. ---*/
  for(unsigned short nn=0; nn<3; ++nn)
    OwnGemm(gemmDOFs2Int, jitterDOFs2Int, nIntegrationPad, 3, nDOFs,
            nIntegrationPad, nDOFs, nIntegrationPad,
            derLagBasisInt[nn], matCoor, matDerCoor[nn], nullptr);
}

void CFEMStandardTetGrid::Derivatives2ndCoorIntPoints(ColMajorMatrix<su2double>          &matCoor,
                                                      vector<ColMajorMatrix<su2double> > &matDer2ndCoor) {

  /*--- Call the function OwnGemm 6 times to compute the 2nd derivatives of
        the Cartesian coordinates w.r.t. the three parametric coordinates. ---*/
  for(unsigned short nn=0; nn<6; ++nn)
    OwnGemm(gemmDOFs2Int, jitterDOFs2Int, nIntegrationPad, 3, nDOFs,
            nIntegrationPad, nDOFs, nIntegrationPad,
            hesLagBasisInt[nn], matCoor, matDer2ndCoor[nn], nullptr);
}

void CFEMStandardTetGrid::CoorSolDOFs(ColMajorMatrix<su2double> &matCoorDOF,
                                      ColMajorMatrix<su2double> &matCoorSolDOF) {

  /*--- Call OwnGemm with the appropriate arguments to compute the data. ---*/
  const unsigned short nSolDOFs = lagBasisSolDOFs.rows();
  OwnGemm(gemmDOFs2SolDOFs, jitterDOFs2SolDOFs, nSolDOFs, 3, nDOFs,
          nSolDOFs, nDOFs, nSolDOFs,
          lagBasisSolDOFs, matCoorDOF, matCoorSolDOF, nullptr);
}

void CFEMStandardTetGrid::DerivativesCoorSolDOFs(ColMajorMatrix<su2double>          &matCoor,
                                                 vector<ColMajorMatrix<su2double> > &matDerCoor) {

  /*--- Call OwnGemm with the appropriate arguments to compute the data. ---*/
  const unsigned short nSolDOFs = lagBasisSolDOFs.rows();
  for(unsigned short nn=0; nn<3; ++nn)
    OwnGemm(gemmDOFs2SolDOFs, jitterDOFs2SolDOFs, nSolDOFs, 3, nDOFs,
            nSolDOFs, nDOFs, nSolDOFs,
            derLagBasisSolDOFs[nn], matCoor, matDerCoor[nn], nullptr);
}

void CFEMStandardTetGrid::EvalCoorAndGradCoor(ColMajorMatrix<su2double> &matCoor,
                                              const su2double           *par,
                                              su2double                 x[3],
                                              su2double                 dxdpar[3][3]) {

  /*--- Convert the parametric coordinates to a passivedouble and
        store them in vectors of size 1. ---*/
  vector<passivedouble> rPar(1, SU2_TYPE::GetValue(par[0]));
  vector<passivedouble> sPar(1, SU2_TYPE::GetValue(par[1]));
  vector<passivedouble> tPar(1, SU2_TYPE::GetValue(par[2]));

  /*--- Compute the Lagrangian basis functions and its first
        derivatives in this parametric point. ---*/
  ColMajorMatrix<passivedouble> lag;
  vector<ColMajorMatrix<passivedouble> > derLag;

  LagBasisIntPointsTet(nPoly, rTetDOFs, sTetDOFs, tTetDOFs, rPar, sPar, tPar, lag);
  DerLagBasisIntPointsTet(nPoly, rTetDOFs, sTetDOFs, tTetDOFs, rPar, sPar, tPar, derLag);

  /*--- Initialize the coordinates and its derivatives. ---*/
  for(unsigned short j=0; j<3; ++j) {
    x[j] = 0.0;
    for(unsigned short i=0; i<3; ++i)
      dxdpar[j][i] = 0.0;
  }

  /*--- Loop to compute the coordinates and its derivatives. ---*/
  for(unsigned short l=0; l<3; ++l) {
    for(unsigned short i=0; i<nDOFs; ++i) {
      x[l]         += matCoor(i,l)*lag(0,i);
      dxdpar[l][0] += matCoor(i,l)*derLag[0](0,i);
      dxdpar[l][1] += matCoor(i,l)*derLag[1](0,i);
      dxdpar[l][2] += matCoor(i,l)*derLag[2](0,i);
    }
  }
}

void CFEMStandardTetGrid::InterpolCoorSubElem(const unsigned short subElem,
                                              const su2double      *weights,
                                              su2double            *parCoor) {

  /*--- Easier storage of the connectivity of the sub-element. ---*/
  const unsigned short *conn = subConn1ForPlotting.data() + 4*subElem;
  const unsigned short n0    = conn[0];
  const unsigned short n1    = conn[1];
  const unsigned short n2    = conn[2];
  const unsigned short n3    = conn[3];

  /*--- Compute the parametric coordinates. ---*/
  parCoor[0] = weights[0]*rTetDOFs[n0] + weights[1]*rTetDOFs[n1]
             + weights[2]*rTetDOFs[n2] + weights[3]*rTetDOFs[n3];
  parCoor[1] = weights[0]*sTetDOFs[n0] + weights[1]*sTetDOFs[n1]
             + weights[2]*sTetDOFs[n2] + weights[3]*sTetDOFs[n3];
  parCoor[2] = weights[0]*tTetDOFs[n0] + weights[1]*tTetDOFs[n1]
             + weights[2]*tTetDOFs[n2] + weights[3]*tTetDOFs[n3];
}
