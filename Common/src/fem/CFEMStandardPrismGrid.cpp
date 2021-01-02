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

/*----------------------------------------------------------------------------------*/
/*             Public member functions of CFEMStandardPrismGrid.                    */
/*----------------------------------------------------------------------------------*/

CFEMStandardPrismGrid::CFEMStandardPrismGrid(const unsigned short val_nPolyGrid,
                                             const unsigned short val_nPolySol,
                                             const unsigned short val_orderExact,
                                             const unsigned short val_locGridDOFs)
  : CFEMStandardPrismBase(val_nPolyGrid, val_orderExact) {

  /*--- Determine the location of the grid and nodal solution DOFs. ---*/
  if(val_locGridDOFs == LGL) {
    Location1DGridDOFsLGL(nPoly, rLineDOFs);
    Location1DGridDOFsLGL(val_nPolySol, rLineSolDOFs);

    LocationTriangleGridDOFsLGL(nPoly, rTriangleDOFs, sTriangleDOFs);
    LocationTriangleGridDOFsLGL(val_nPolySol, rTriangleSolDOFs, sTriangleSolDOFs);
  }
  else {
    Location1DGridDOFsEquidistant(nPoly, rLineDOFs);
    Location1DGridDOFsEquidistant(val_nPolySol, rLineSolDOFs);

    LocationTriangleGridDOFsEquidistant(nPoly, rTriangleDOFs, sTriangleDOFs);
    LocationTriangleGridDOFsEquidistant(val_nPolySol, rTriangleSolDOFs, sTriangleSolDOFs);
  }

  /*--- Compute the corresponding Lagrangian basis functions and
        its first and second derivatives in the integration points. ---*/
  LagBasisIntPointsPrism(nPoly, rTriangleDOFs, sTriangleDOFs, rLineDOFs,
                         rTriangleInt, sTriangleInt, rLineInt, lagBasisInt);
  DerLagBasisIntPointsPrism(nPoly, rTriangleDOFs, sTriangleDOFs, rLineDOFs,
                            rTriangleInt, sTriangleInt, rLineInt, derLagBasisInt);
  HesLagBasisIntPointsPrism(nPoly, rTriangleDOFs, sTriangleDOFs, rLineDOFs,
                            rTriangleInt, sTriangleInt, rLineInt, hesLagBasisInt);

  /*--- Compute the Lagrangian basis functions and its derivatives
        in the nodal solution DOFs. Use LagBasisIntPointsPrism and
        DerLagBasisIntPointsPrism with different arguments. ---*/
  LagBasisIntPointsPrism(nPoly, rTriangleDOFs, sTriangleDOFs, rLineDOFs,
                         rTriangleSolDOFs, sTriangleSolDOFs, rLineSolDOFs,
                         lagBasisSolDOFs);
  DerLagBasisIntPointsPrism(nPoly, rTriangleDOFs, sTriangleDOFs, rLineDOFs,
                            rTriangleSolDOFs, sTriangleSolDOFs, rLineSolDOFs,
                            derLagBasisSolDOFs);

  /*--- Determine the local subconnectivity of this standard element when split
        in several linear elements. Used for a.o. plotting and searcing. ---*/
  CFEMStandardPrismBase::SubConnLinearElements();

  /*--- Set up the jitted gemm calls, if supported. ---*/
  SetUpJittedGEMM(nIntegrationPad, 3, nDOFs, nIntegrationPad, nDOFs,
                  nIntegrationPad, jitterDOFs2Int, gemmDOFs2Int);
  SetUpJittedGEMM(lagBasisSolDOFs.rows(), 3, nDOFs, lagBasisSolDOFs.rows(), nDOFs,
                  lagBasisSolDOFs.rows(), jitterDOFs2SolDOFs, gemmDOFs2SolDOFs);
}

CFEMStandardPrismGrid::~CFEMStandardPrismGrid() {

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

void CFEMStandardPrismGrid::CoorIntPoints(const bool                notUsed,
                                          ColMajorMatrix<su2double> &matCoorDOF,
                                          ColMajorMatrix<su2double> &matCoorInt) {

  /*--- Call OwnGemm with the appropriate arguments to compute the data. ---*/
  OwnGemm(gemmDOFs2Int, jitterDOFs2Int, nIntegrationPad, 3, nDOFs,
          nIntegrationPad, nDOFs, nIntegrationPad,
          lagBasisInt, matCoorDOF, matCoorInt, nullptr);
}

void CFEMStandardPrismGrid::DerivativesCoorIntPoints(const bool                         notUsed,
                                                     ColMajorMatrix<su2double>          &matCoor,
                                                     vector<ColMajorMatrix<su2double> > &matDerCoor) {

  /*--- Call OwnGemm with the appropriate arguments to compute the data. ---*/
  for(unsigned short nn=0; nn<3; ++nn)
    OwnGemm(gemmDOFs2Int, jitterDOFs2Int, nIntegrationPad, 3, nDOFs,
            nIntegrationPad, nDOFs, nIntegrationPad,
            derLagBasisInt[nn], matCoor, matDerCoor[nn], nullptr);
}

void CFEMStandardPrismGrid::Derivatives2ndCoorIntPoints(ColMajorMatrix<su2double>          &matCoor,
                                                        vector<ColMajorMatrix<su2double> > &matDer2ndCoor) {

  /*--- Call the function OwnGemm 6 times to compute the 2nd derivatives of
        the Cartesian coordinates w.r.t. the three parametric coordinates. ---*/
  for(unsigned short nn=0; nn<6; ++nn)
    OwnGemm(gemmDOFs2Int, jitterDOFs2Int, nIntegrationPad, 3, nDOFs,
            nIntegrationPad, nDOFs, nIntegrationPad,
            hesLagBasisInt[nn], matCoor, matDer2ndCoor[nn], nullptr);
}

void CFEMStandardPrismGrid::DerivativesCoorSolDOFs(ColMajorMatrix<su2double>          &matCoor,
                                                   vector<ColMajorMatrix<su2double> > &matDerCoor) {

  /*--- Call OwnGemm with the appropriate arguments to compute the data. ---*/
  const unsigned short nSolDOFs = lagBasisSolDOFs.rows();
  for(unsigned short nn=0; nn<3; ++nn)
    OwnGemm(gemmDOFs2SolDOFs, jitterDOFs2SolDOFs, nSolDOFs, 3, nDOFs,
            nSolDOFs, nDOFs, nSolDOFs,
            derLagBasisSolDOFs[nn], matCoor, matDerCoor[nn], nullptr);
}
