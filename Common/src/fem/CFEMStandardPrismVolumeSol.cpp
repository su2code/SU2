/*!
 * \file CFEMStandardPrismVolumeSol.cpp
 * \brief Functions for the class CFEMStandardPrismVolumeSol.
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

#include "../../include/fem/CFEMStandardPrismVolumeSol.hpp"

/*----------------------------------------------------------------------------------*/
/*             Public member functions of CFEMStandardPrismVolumeSol.               */
/*----------------------------------------------------------------------------------*/

CFEMStandardPrismVolumeSol::CFEMStandardPrismVolumeSol(const unsigned short val_nPoly,
                                                       const unsigned short val_orderExact,
                                                       const unsigned short val_locGridDOFs,
                                                       const unsigned short val_nVar)
  : CFEMStandardPrismBase(val_nPoly, val_orderExact) {

  /*--- Determine the location of the nodal solution DOFs. ---*/
  if(val_locGridDOFs == LGL) {
    Location1DGridDOFsLGL(nPoly, rLineSolDOFs);
    LocationTriangleGridDOFsLGL(nPoly, rTriangleSolDOFs, sTriangleSolDOFs);
  }
  else {
    Location1DGridDOFsEquidistant(nPoly, rLineSolDOFs);
    LocationTriangleGridDOFsEquidistant(nPoly, rTriangleSolDOFs, sTriangleSolDOFs);
  }

  /*--- Determine the location of all solution DOFs. ---*/
  vector<passivedouble> rDOFs, sDOFs, tDOFs;
  LocationAllPointsPrism(rTriangleSolDOFs, sTriangleSolDOFs, rLineSolDOFs,
                         rDOFs, sDOFs, tDOFs);

  /*--- Determine the location of all integration points. ---*/
  vector<passivedouble> rInt, sInt, tInt;
  LocationAllPointsPrism(rTriangleInt, sTriangleInt, rLineInt, rInt, sInt, tInt);

  /*--- Allocate the memory for the Legendre basis functions and its
        1st and 2nd derivatives in the integration points. ---*/
  legBasisInt.resize(nIntegrationPad, nDOFs); legBasisInt.setConstant(0.0);

  derLegBasisInt.resize(3);
  derLegBasisInt[0].resize(nIntegrationPad, nDOFs); derLegBasisInt[0].setConstant(0.0);
  derLegBasisInt[1].resize(nIntegrationPad, nDOFs); derLegBasisInt[1].setConstant(0.0);
  derLegBasisInt[2].resize(nIntegrationPad, nDOFs); derLegBasisInt[2].setConstant(0.0);

  hesLegBasisInt.resize(6);
  hesLegBasisInt[0].resize(nIntegrationPad, nDOFs); hesLegBasisInt[0].setConstant(0.0);
  hesLegBasisInt[1].resize(nIntegrationPad, nDOFs); hesLegBasisInt[1].setConstant(0.0);
  hesLegBasisInt[2].resize(nIntegrationPad, nDOFs); hesLegBasisInt[2].setConstant(0.0);
  hesLegBasisInt[3].resize(nIntegrationPad, nDOFs); hesLegBasisInt[3].setConstant(0.0);
  hesLegBasisInt[4].resize(nIntegrationPad, nDOFs); hesLegBasisInt[4].setConstant(0.0);
  hesLegBasisInt[5].resize(nIntegrationPad, nDOFs); hesLegBasisInt[5].setConstant(0.0);

  /*--- Compute the Legendre basis functions and its first
        and second derivatives in the integration points. ---*/
  VandermondePrism(nPoly, rInt, sInt, tInt, legBasisInt);
  GradVandermondePrism(nPoly, rInt, sInt, tInt, derLegBasisInt[0],
                       derLegBasisInt[1], derLegBasisInt[2]);
  HesVandermondePrism(nPoly, rInt, sInt, tInt, hesLegBasisInt[0],
                      hesLegBasisInt[1], hesLegBasisInt[2], hesLegBasisInt[3],
                      hesLegBasisInt[4], hesLegBasisInt[5]);

  /*--- Allocate the memory for the Legendre basis functions and its
        1st derivatives in the solution DOFs. ---*/
  legBasisSolDOFs.resize(nDOFs, nDOFs); legBasisSolDOFs.setConstant(0.0);

  derLegBasisSolDOFs.resize(3);
  derLegBasisSolDOFs[0].resize(nDOFs, nDOFs); derLegBasisSolDOFs[0].setConstant(0.0);
  derLegBasisSolDOFs[1].resize(nDOFs, nDOFs); derLegBasisSolDOFs[1].setConstant(0.0);
  derLegBasisSolDOFs[2].resize(nDOFs, nDOFs); derLegBasisSolDOFs[2].setConstant(0.0);

  /*--- Compute the Legendre basis functions and its first
        derivatives in the solution DOFs. ---*/
  VandermondePrism(nPoly, rDOFs, sDOFs, tDOFs, legBasisSolDOFs);
  GradVandermondePrism(nPoly, rDOFs, sDOFs, tDOFs, derLegBasisSolDOFs[0],
                       derLegBasisSolDOFs[1], derLegBasisSolDOFs[2]);

  /*--- Determine the local subconnectivity of this standard element when split
        in several linear elements. Used for a.o. plotting and searching. ---*/
  CFEMStandardPrismBase::SubConnLinearElements();

  /*--- Set up the jitted gemm calls, if supported. ---*/
  SetUpJittedGEMM(nIntegrationPad, val_nVar, nDOFs, jitterDOFs2Int, gemmDOFs2Int);
  SetUpJittedGEMM(nDOFs, val_nVar, nDOFs, jitterDOFs2SolDOFs, gemmDOFs2SolDOFs);
}

CFEMStandardPrismVolumeSol::~CFEMStandardPrismVolumeSol() {

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
