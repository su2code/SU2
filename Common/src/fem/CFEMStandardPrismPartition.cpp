/*!
 * \file CFEMStandardPrismPartition.cpp
 * \brief Functions for the class CFEMStandardPrismPartition.
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

#include "../../include/fem/CFEMStandardPrismPartition.hpp"

/*----------------------------------------------------------------------------------*/
/*             Public member functions of CFEMStandardPrismPartition.               */
/*----------------------------------------------------------------------------------*/

CFEMStandardPrismPartition::CFEMStandardPrismPartition(const unsigned short val_nPoly,
                                                       const unsigned short val_orderExact)
  : CFEMStandardPrismBase(val_nPoly, val_orderExact) {

  /*--- Determine the 1D parametric locations of the grid DOFs. These are needed
        as the 3D grid DOFs are obtained by taking the tensor product of the
        1D grid DOFs and the grid DOFs of the base triangle. ---*/
  Location1DGridDOFsEquidistant(nPoly, rLineDOFsEqui);
  Location1DGridDOFsLGL(nPoly, rLineDOFsLGL);

  /*--- Determine the parametric locations of the grid DOFs of the triangle. ---*/
  LocationTriangleGridDOFsEquidistant(nPoly, rTriangleDOFsEqui, sTriangleDOFsEqui);
  LocationTriangleGridDOFsLGL(nPoly, rTriangleDOFsLGL, sTriangleDOFsLGL);

  /*--- Compute the values of the Lagrangian basis functions in the integration
        points for both the equidistant and LGL point distribution. ---*/
  LagBasisIntPointsPrism(nPoly, rTriangleDOFsEqui, sTriangleDOFsEqui, rLineDOFsEqui,
                         rTriangleInt, sTriangleInt, rLineInt, lagBasisIntEqui);
  LagBasisIntPointsPrism(nPoly, rTriangleDOFsLGL,  sTriangleDOFsLGL,  rLineDOFsLGL,
                         rTriangleInt, sTriangleInt, rLineInt, lagBasisIntLGL);

  /*--- Compute the values of the derivatives of the Lagrangian basis functions in
        the integration points for both the equidistant and LGL point distribution. ---*/
  DerLagBasisIntPointsPrism(nPoly, rTriangleDOFsEqui, sTriangleDOFsEqui, rLineDOFsEqui,
                            rTriangleInt, sTriangleInt, rLineInt, derLagBasisIntEqui);
  DerLagBasisIntPointsPrism(nPoly, rTriangleDOFsLGL,  sTriangleDOFsLGL,  rLineDOFsLGL,
                            rTriangleInt, sTriangleInt, rLineInt, derLagBasisIntLGL);

  /*--- Create the local grid connectivities of the faces of the volume element. ---*/
  LocalGridConnFaces();

  /*--- Determine the local subconnectivity of this standard element when split
        in several linear elements. Used for a.o. plotting and searcing. ---*/
  CFEMStandardPrismBase::SubConnLinearElements();

  /*--- Set up the jitted gemm call, if supported. For this particular standard
        element the derivative of the coordinates are computed, which is 3. ---*/
  SetUpJittedGEMM(nIntegrationPad, 3, nDOFs, nIntegrationPad, nDOFs,
                  nIntegrationPad, jitterDOFs2Int, gemmDOFs2Int);
}

CFEMStandardPrismPartition::~CFEMStandardPrismPartition() {

#if defined(PRIMAL_SOLVER) && defined(HAVE_MKL)
  if( jitterDOFs2Int ) {
    mkl_jit_destroy(jitterDOFs2Int);
    jitterDOFs2Int = nullptr;
  }
#endif
}

void CFEMStandardPrismPartition::CoorIntPoints(const bool                LGLDistribution,
                                               ColMajorMatrix<su2double> &matCoorDOF,
                                               ColMajorMatrix<su2double> &matCoorInt) {

  /*--- Check for which point distribution the derivatives must be computed. ---*/
  if( LGLDistribution ) {

    /*--- LGL distribution. Call the function OwnGemm to compute the
          Cartesian coordinates in the integration points. ---*/
    OwnGemm(gemmDOFs2Int, jitterDOFs2Int, nIntegrationPad, 3, nDOFs,
            nIntegrationPad, nDOFs, nIntegrationPad,
            lagBasisIntLGL, matCoorDOF, matCoorInt, nullptr);
  }
  else {

    /*--- Equidistant distribution. Call the function OwnGemm to compute the
          Cartesian coordinates in the integration points. ---*/
    OwnGemm(gemmDOFs2Int, jitterDOFs2Int, nIntegrationPad, 3, nDOFs,
            nIntegrationPad, nDOFs, nIntegrationPad,
            lagBasisIntEqui, matCoorDOF, matCoorInt, nullptr);
  }
}

void CFEMStandardPrismPartition::DerivativesCoorIntPoints(const bool                         LGLDistribution,
                                                          ColMajorMatrix<su2double>          &matCoor,
                                                          vector<ColMajorMatrix<su2double> > &matDerCoor) {

  /*--- Check for which point distribution the derivatives must be computed. ---*/
  if( LGLDistribution ) {

    /*--- LGL distribution. Call the function OwnGemm 3 times to compute the derivatives
          of the Cartesian coordinates w.r.t. the three parametric coordinates. ---*/
    for(unsigned short nn=0; nn<3; ++nn)
      OwnGemm(gemmDOFs2Int, jitterDOFs2Int, nIntegrationPad, 3, nDOFs,
              nIntegrationPad, nDOFs, nIntegrationPad,
              derLagBasisIntLGL[nn], matCoor, matDerCoor[nn], nullptr);
  }
  else {

    /*--- Equidistant distribution. Call the function OwnGemm 3 times to compute the derivatives
          of the Cartesian coordinates w.r.t. the three parametric coordinates. ---*/
    for(unsigned short nn=0; nn<3; ++nn)
    OwnGemm(gemmDOFs2Int, jitterDOFs2Int, nIntegrationPad, 3, nDOFs,
            nIntegrationPad, nDOFs, nIntegrationPad,
            derLagBasisIntEqui[nn], matCoor, matDerCoor[nn], nullptr);
  }
}

passivedouble CFEMStandardPrismPartition::WorkEstimateVolume(CConfig *config) {

  /*--- TEMPORARY IMPLEMENTATION. ---*/
  return nIntegration + 0.1*nDOFs;
}

/*----------------------------------------------------------------------------------*/
/*             Private member functions of CFEMStandardPrismPartition.              */
/*----------------------------------------------------------------------------------*/

void CFEMStandardPrismPartition::LocalGridConnFaces(void) {

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
