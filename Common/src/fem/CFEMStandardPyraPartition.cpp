/*!
 * \file CFEMStandardPyraPartition.cpp
 * \brief Functions for the class CFEMStandardPyraPartition.
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

#include "../../include/fem/CFEMStandardPyraPartition.hpp"

/*----------------------------------------------------------------------------------*/
/*             Public member functions of CFEMStandardPyraPartition.                */
/*----------------------------------------------------------------------------------*/

CFEMStandardPyraPartition::CFEMStandardPyraPartition(const unsigned short val_nPoly,
                                                     const unsigned short val_orderExact)
  : CFEMStandardPyraBase(val_nPoly, val_orderExact) {

  /*--- Determine the parametric locations of the grid DOFs of the pyramid. ---*/
  LocationPyramidGridDOFsEquidistant(nPoly, rPyraDOFsEqui, sPyraDOFsEqui, tPyraDOFsEqui);
  LocationPyramidGridDOFsLGL(nPoly, rPyraDOFsLGL, sPyraDOFsLGL, tPyraDOFsLGL);

  /*--- Determine the location of all integration points of the pyramid. ---*/
  vector<passivedouble> rInt, sInt, tInt;
  LocationAllIntegrationPoints(rInt, sInt, tInt);

  /*--- Compute the values of the Lagrangian basis functions in the integration
        points for both the equidistant and LGL point distribution. ---*/
  LagBasisIntPointsPyra(nPoly, rPyraDOFsEqui, sPyraDOFsEqui, tPyraDOFsEqui,
                        rInt, sInt, tInt, lagBasisIntEqui);
  LagBasisIntPointsPyra(nPoly, rPyraDOFsLGL,  sPyraDOFsLGL,  tPyraDOFsLGL,
                        rInt, sInt, tInt, lagBasisIntLGL);

  /*--- Compute the values of the derivatives of the Lagrangian basis functions in
        the integration points for both the equidistant and LGL point distribution. ---*/
  DerLagBasisIntPointsPyra(nPoly, rPyraDOFsEqui, sPyraDOFsEqui, tPyraDOFsEqui,
                           rInt, sInt, tInt, derLagBasisIntEqui);
  DerLagBasisIntPointsPyra(nPoly, rPyraDOFsLGL,  sPyraDOFsLGL,  tPyraDOFsLGL,
                           rInt, sInt, tInt, derLagBasisIntLGL);

  /*--- Create the local grid connectivities of the faces of the volume element. ---*/
  LocalGridConnFaces();

  /*--- Determine the local subconnectivity of this standard element when split
        in several linear elements. Used for a.o. plotting and searcing. ---*/
  CFEMStandardPyraBase::SubConnLinearElements();

  /*--- Set up the jitted gemm call, if supported. For this particular standard
        element the derivative of the coordinates are computed, which is 3. ---*/
  SetUpJittedGEMM(nIntegrationPad, 3, nDOFs, nIntegrationPad, nDOFs,
                  nIntegrationPad, true, jitterDOFs2Int, gemmDOFs2Int);
}

CFEMStandardPyraPartition::~CFEMStandardPyraPartition() {

#if defined(PRIMAL_SOLVER) && defined(HAVE_MKL)
  if( jitterDOFs2Int ) {
    mkl_jit_destroy(jitterDOFs2Int);
    jitterDOFs2Int = nullptr;
  }
#endif
}

void CFEMStandardPyraPartition::CoorIntPoints(const bool                LGLDistribution,
                                              ColMajorMatrix<su2double> &matCoorDOF,
                                              ColMajorMatrix<su2double> &matCoorInt) {

  /*--- Check for which point distribution the derivatives must be computed. ---*/
  if( LGLDistribution ) {

    /*--- LGL distribution. Call the function OwnGemm to compute the
          Cartesian coordinates in the integration points. ---*/
    OwnGemm(gemmDOFs2Int, jitterDOFs2Int, nIntegrationPad, 3, nDOFs,
            nIntegrationPad, nDOFs, nIntegrationPad, true,
            lagBasisIntLGL, matCoorDOF, matCoorInt, nullptr);
  }
  else {

    /*--- Equidistant distribution. Call the function OwnGemm to compute the
          Cartesian coordinates in the integration points. ---*/
    OwnGemm(gemmDOFs2Int, jitterDOFs2Int, nIntegrationPad, 3, nDOFs,
            nIntegrationPad, nDOFs, nIntegrationPad, true,
            lagBasisIntEqui, matCoorDOF, matCoorInt, nullptr);
  }
}

void CFEMStandardPyraPartition::DerivativesCoorIntPoints(const bool                         LGLDistribution,
                                                         ColMajorMatrix<su2double>          &matCoor,
                                                         vector<ColMajorMatrix<su2double> > &matDerCoor) {

  /*--- Check for which point distribution the derivatives must be computed. ---*/
  if( LGLDistribution ) {

    /*--- LGL distribution. Call the function OwnGemm 3 times to compute the derivatives
          of the Cartesian coordinates w.r.t. the three parametric coordinates. ---*/
    for(unsigned short nn=0; nn<3; ++nn)
      OwnGemm(gemmDOFs2Int, jitterDOFs2Int, nIntegrationPad, 3, nDOFs,
              nIntegrationPad, nDOFs, nIntegrationPad, true,
              derLagBasisIntLGL[nn], matCoor, matDerCoor[nn], nullptr);
  }
  else {

    /*--- Equidistant distribution. Call the function OwnGemm 3 times to compute the derivatives
          of the Cartesian coordinates w.r.t. the three parametric coordinates. ---*/
    for(unsigned short nn=0; nn<3; ++nn)
      OwnGemm(gemmDOFs2Int, jitterDOFs2Int, nIntegrationPad, 3, nDOFs,
              nIntegrationPad, nDOFs, nIntegrationPad, true,
              derLagBasisIntEqui[nn], matCoor, matDerCoor[nn], nullptr);
  }
}

passivedouble CFEMStandardPyraPartition::WorkEstimateVolume(CConfig *config) {

  /*--- TEMPORARY IMPLEMENTATION. ---*/
  return nIntegration + 0.1*nDOFs;
}