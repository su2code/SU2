/*!
 * \file CFEMStandardPyraPartition.cpp
 * \brief Functions for the class CFEMStandardPyraPartition.
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
  SetUpJittedGEMM(nIntegrationPad, 3, nDOFs, jitterDOFs2Int, gemmDOFs2Int);
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
    OwnGemm(gemmDOFs2Int, nIntegrationPad, 3, nDOFs, lagBasisIntLGL, matCoorDOF, matCoorInt, nullptr);
  }
  else {

    /*--- Equidistant distribution. Call the function OwnGemm to compute the
          Cartesian coordinates in the integration points. ---*/
    OwnGemm(gemmDOFs2Int, nIntegrationPad, 3, nDOFs, lagBasisIntEqui, matCoorDOF, matCoorInt, nullptr);
  }
}

void CFEMStandardPyraPartition::DerivativesCoorIntPoints(const bool                         LGLDistribution,
                                                         ColMajorMatrix<su2double>          &matCoor,
                                                         vector<ColMajorMatrix<su2double> > &matDerCoor) {

  /*--- Check for which point distribution the derivatives must be computed. ---*/
  if( LGLDistribution ) {

    /*--- LGL distribution. Call the function OwnGemm 3 times to compute the derivatives
          of the Cartesian coordinates w.r.t. the three parametric coordinates. ---*/
    OwnGemm(gemmDOFs2Int, nIntegrationPad, 3, nDOFs, derLagBasisIntLGL[0], matCoor, matDerCoor[0], nullptr);
    OwnGemm(gemmDOFs2Int, nIntegrationPad, 3, nDOFs, derLagBasisIntLGL[1], matCoor, matDerCoor[1], nullptr);
    OwnGemm(gemmDOFs2Int, nIntegrationPad, 3, nDOFs, derLagBasisIntLGL[2], matCoor, matDerCoor[2], nullptr);
  }
  else {

    /*--- Equidistant distribution. Call the function OwnGemm 3 times to compute the derivatives
          of the Cartesian coordinates w.r.t. the three parametric coordinates. ---*/
    OwnGemm(gemmDOFs2Int, nIntegrationPad, 3, nDOFs, derLagBasisIntEqui[0], matCoor, matDerCoor[0], nullptr);
    OwnGemm(gemmDOFs2Int, nIntegrationPad, 3, nDOFs, derLagBasisIntEqui[1], matCoor, matDerCoor[1], nullptr);
    OwnGemm(gemmDOFs2Int, nIntegrationPad, 3, nDOFs, derLagBasisIntEqui[2], matCoor, matDerCoor[2], nullptr);
  }
}

passivedouble CFEMStandardPyraPartition::WorkEstimateVolume(CConfig *config) {

  /*--- TEMPORARY IMPLEMENTATION. ---*/
  return nIntegration + 0.1*nDOFs;
}

/*----------------------------------------------------------------------------------*/
/*             Private member functions of CFEMStandardPyraPartition.               */
/*----------------------------------------------------------------------------------*/

void CFEMStandardPyraPartition::LocalGridConnFaces(void) {

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
