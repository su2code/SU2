/*!
 * \file CFEMStandardTriGrid.cpp
 * \brief Functions for the class CFEMStandardTriGrid.
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

#include "../../include/fem/CFEMStandardTriGrid.hpp"

/*----------------------------------------------------------------------------------*/
/*             Public member functions of CFEMStandardTriGrid.                      */
/*----------------------------------------------------------------------------------*/

CFEMStandardTriGrid::CFEMStandardTriGrid(const unsigned short val_nPolyGrid,
                                         const unsigned short val_nPolySol,
                                         const unsigned short val_orderExact,
                                         const unsigned short val_locGridDOFs)
  : CFEMStandardTriBase(val_nPolyGrid, val_orderExact) {

  /*--- Determine the location of the grid and nodal solution DOFs. ---*/
  if(val_locGridDOFs == LGL) {
    LocationTriangleGridDOFsLGL(nPoly, rTriangleDOFs, sTriangleDOFs);    
    LocationTriangleGridDOFsLGL(val_nPolySol, rTriangleSolDOFs, sTriangleSolDOFs);
  }
  else {
    LocationTriangleGridDOFsEquidistant(nPoly, rTriangleDOFs, sTriangleDOFs);
    LocationTriangleGridDOFsEquidistant(val_nPolySol, rTriangleSolDOFs, sTriangleSolDOFs);
  }

  /*--- Compute the corresponding Lagrangian basis functions and
        its first and second derivatives in the integration points. ---*/
  LagBasisIntPointsTriangle(nPoly, rTriangleDOFs,  sTriangleDOFs,
                            rTriangleInt, sTriangleInt, lagBasisInt);
  DerLagBasisIntPointsTriangle(nPoly, rTriangleDOFs, sTriangleDOFs,
                               rTriangleInt, sTriangleInt, derLagBasisInt); 
  HesLagBasisIntPointsTriangle(nPoly, rTriangleDOFs, sTriangleDOFs,
                               rTriangleInt, sTriangleInt, hesLagBasisInt);

  /*--- Compute the Lagrangian basis functions and its derivatives
        in the nodal solution DOFs. Use LagBasisIntPointsTriangle and
        DerLagBasisIntPointsTriangle with different arguments. ---*/
  LagBasisIntPointsTriangle(nPoly, rTriangleDOFs,  sTriangleDOFs,
                            rTriangleSolDOFs, sTriangleSolDOFs, lagBasisSolDOFs);
  DerLagBasisIntPointsTriangle(nPoly, rTriangleDOFs, sTriangleDOFs,
                               rTriangleSolDOFs, sTriangleSolDOFs, derLagBasisSolDOFs);

  /*--- Determine the local subconnectivity of this standard element when split
        in several linear elements. Used for a.o. plotting and searcing. ---*/
  CFEMStandardTriBase::SubConnLinearElements();

  /*--- Set up the jitted gemm calls, if supported. ---*/
  SetUpJittedGEMM(nIntegrationPad, 2, nDOFs, nIntegrationPad, nDOFs,
                  nIntegrationPad, jitterDOFs2Int, gemmDOFs2Int);
  SetUpJittedGEMM(lagBasisSolDOFs.rows(), 2, nDOFs, lagBasisSolDOFs.rows(), nDOFs,
                  lagBasisSolDOFs.rows(), jitterDOFs2SolDOFs, gemmDOFs2SolDOFs);
}

CFEMStandardTriGrid::~CFEMStandardTriGrid() {

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

void CFEMStandardTriGrid::CoorIntPoints(const bool                notUsed,
                                        ColMajorMatrix<su2double> &matCoorDOF,
                                        ColMajorMatrix<su2double> &matCoorInt) {

  /*--- Call OwnGemm with the appropriate arguments to compute the data. ---*/
  OwnGemm(gemmDOFs2Int, jitterDOFs2Int, nIntegrationPad, 2, nDOFs,
          nIntegrationPad, nDOFs, nIntegrationPad,
          lagBasisInt, matCoorDOF, matCoorInt, nullptr);
}

void CFEMStandardTriGrid::DerivativesCoorIntPoints(const bool                         notUsed,
                                                   ColMajorMatrix<su2double>          &matCoor,
                                                   vector<ColMajorMatrix<su2double> > &matDerCoor) {

  /*--- Call OwnGemm with the appropriate arguments to compute the data. ---*/
  for(unsigned short nn=0; nn<2; ++nn)
    OwnGemm(gemmDOFs2Int, jitterDOFs2Int, nIntegrationPad, 2, nDOFs,
            nIntegrationPad, nDOFs, nIntegrationPad,
            derLagBasisInt[nn], matCoor, matDerCoor[nn], nullptr);
}

void CFEMStandardTriGrid::Derivatives2ndCoorIntPoints(ColMajorMatrix<su2double>          &matCoor,
                                                      vector<ColMajorMatrix<su2double> > &matDer2ndCoor) {

  /*--- Call OwnGemm with the appropriate arguments to compute the data. ---*/
  for(unsigned short nn=0; nn<3; ++nn)
    OwnGemm(gemmDOFs2Int, jitterDOFs2Int, nIntegrationPad, 2, nDOFs,
            nIntegrationPad, nDOFs, nIntegrationPad,
            hesLagBasisInt[nn], matCoor, matDer2ndCoor[nn], nullptr);
}

void CFEMStandardTriGrid::CoorSolDOFs(ColMajorMatrix<su2double> &matCoorDOF,
                                      ColMajorMatrix<su2double> &matCoorSolDOF) {

  /*--- Call OwnGemm with the appropriate arguments to compute the data. ---*/
  const unsigned short nSolDOFs = lagBasisSolDOFs.rows();
  OwnGemm(gemmDOFs2SolDOFs, jitterDOFs2SolDOFs, nSolDOFs, 2, nDOFs,
          nSolDOFs, nDOFs, nSolDOFs,
          lagBasisSolDOFs, matCoorDOF, matCoorSolDOF, nullptr);
}

void CFEMStandardTriGrid::DerivativesCoorSolDOFs(ColMajorMatrix<su2double>          &matCoor,
                                                 vector<ColMajorMatrix<su2double> > &matDerCoor) {

  /*--- Call OwnGemm with the appropriate arguments to compute the data. ---*/
  const unsigned short nSolDOFs = lagBasisSolDOFs.rows();
  for(unsigned short nn=0; nn<2; ++nn)
    OwnGemm(gemmDOFs2SolDOFs, jitterDOFs2SolDOFs, nSolDOFs, 2, nDOFs,
            nSolDOFs, nDOFs, nSolDOFs,
            derLagBasisSolDOFs[nn], matCoor, matDerCoor[nn], nullptr);
}

void CFEMStandardTriGrid::EvalCoorAndGradCoor(ColMajorMatrix<su2double> &matCoor,
                                              const su2double           *par,
                                              su2double                 x[3],
                                              su2double                 dxdpar[3][3]) {

  /*--- Convert the parametric coordinates to a passivedouble and
        store them in vectors of size 1. ---*/
  vector<passivedouble> rPar(1, SU2_TYPE::GetValue(par[0]));
  vector<passivedouble> sPar(1, SU2_TYPE::GetValue(par[1]));

  /*--- Compute the Lagrangian basis functions and its first
        derivatives in this parametric point. ---*/
  ColMajorMatrix<passivedouble> lag;
  vector<ColMajorMatrix<passivedouble> > derLag;

  LagBasisIntPointsTriangle(nPoly, rTriangleDOFs, sTriangleDOFs, rPar, sPar, lag);
  DerLagBasisIntPointsTriangle(nPoly, rTriangleDOFs, sTriangleDOFs, rPar, sPar, derLag);

  /*--- Initialize the coordinates and its derivatives. ---*/
  for(unsigned short j=0; j<2; ++j) {
    x[j] = 0.0;
    for(unsigned short i=0; i<2; ++i)
      dxdpar[j][i] = 0.0;
  }

  /*--- Loop to compute the coordinates and its derivatives. ---*/
  for(unsigned short l=0; l<2; ++l) {
    for(unsigned short i=0; i<nDOFs; ++i) {
      x[l]         += matCoor(i,l)*lag(0,i);
      dxdpar[l][0] += matCoor(i,l)*derLag[0](0,i);
      dxdpar[l][1] += matCoor(i,l)*derLag[1](0,i);
    }
  }
}

void CFEMStandardTriGrid::InterpolCoorSubElem(const unsigned short subElem,
                                              const su2double      *weights,
                                              su2double            *parCoor) {

  /*--- Easier storage of the connectivity of the sub-element. ---*/
  const unsigned short *conn = subConn1ForPlotting.data() + 3*subElem;
  const unsigned short n0    = conn[0];
  const unsigned short n1    = conn[1];
  const unsigned short n2    = conn[2];
  
  /*--- Compute the parametric coordinates. ---*/
  parCoor[0] = weights[0]*rTriangleDOFs[n0] + weights[1]*rTriangleDOFs[n1] + weights[2]*rTriangleDOFs[n2];
  parCoor[1] = weights[0]*sTriangleDOFs[n0] + weights[1]*sTriangleDOFs[n1] + weights[2]*sTriangleDOFs[n2];
}
