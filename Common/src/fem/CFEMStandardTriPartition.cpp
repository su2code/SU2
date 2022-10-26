/*!
 * \file CFEMStandardTriPartition.cpp
 * \brief Functions for the class CFEMStandardTriPartition.
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

#include "../../include/fem/CFEMStandardTriPartition.hpp"

/*----------------------------------------------------------------------------------*/
/*             Public member functions of CFEMStandardTriPartition.                 */
/*----------------------------------------------------------------------------------*/

CFEMStandardTriPartition::CFEMStandardTriPartition(const unsigned short val_nPoly,
                                                   const unsigned short val_orderExact,
                                                   const bool           val_surfElement)
  : CFEMStandardTriBase(val_nPoly, val_orderExact) {

  /*--- Determine the number of space dimensions, which is 3 if this standard
        element is a surface element and 2 when it is a value element. ---*/
  nDim = val_surfElement ? 3 : 2;

  /*--- Determine the parametric locations of the grid DOFs of the triangle. ---*/
  LocationTriangleGridDOFsEquidistant(nPoly, rTriangleDOFsEqui, sTriangleDOFsEqui);
  LocationTriangleGridDOFsLGL(nPoly, rTriangleDOFsLGL, sTriangleDOFsLGL);

  /*--- Compute the values of the Lagrangian basis functions in the integration
        points for both the equidistant and LGL point distribution. ---*/
  LagBasisIntPointsTriangle(nPoly, rTriangleDOFsEqui, sTriangleDOFsEqui,
                            rTriangleInt, sTriangleInt, lagBasisIntEqui);
  LagBasisIntPointsTriangle(nPoly, rTriangleDOFsLGL,  sTriangleDOFsLGL,
                            rTriangleInt, sTriangleInt, lagBasisIntLGL);

  /*--- Compute the values of the derivatives of the Lagrangian basis functions in
        the integration points for both the equidistant and LGL point distribution. ---*/
  DerLagBasisIntPointsTriangle(nPoly, rTriangleDOFsEqui, sTriangleDOFsEqui,
                               rTriangleInt, sTriangleInt, derLagBasisIntEqui);
  DerLagBasisIntPointsTriangle(nPoly, rTriangleDOFsLGL,  sTriangleDOFsLGL,
                               rTriangleInt, sTriangleInt, derLagBasisIntLGL);

  /*--- Create the local grid connectivities of the faces of the volume element.
        Only needed if this is a volume element. ---*/
  if( !val_surfElement ) LocalGridConnFaces();

  /*--- Determine the local subconnectivity of this standard element when split
        in several linear elements. Used for a.o. plotting and searcing. ---*/
  CFEMStandardTriBase::SubConnLinearElements();

  /*--- Set up the jitted gemm call, if supported. For this particular standard
        element the derivatives of the coordinates are computed, which is nDim. ---*/
  SetUpJittedGEMM(nIntegrationPad, nDim, nDOFs, nIntegrationPad, nDOFs,
                  nIntegrationPad, true, jitterDOFs2Int, gemmDOFs2Int);
}

CFEMStandardTriPartition::~CFEMStandardTriPartition() {

#if defined(PRIMAL_SOLVER) && defined(HAVE_MKL)
  if( jitterDOFs2Int ) {
    mkl_jit_destroy(jitterDOFs2Int);
    jitterDOFs2Int = nullptr;
  }
#endif
}

void CFEMStandardTriPartition::CoorIntPoints(const bool                LGLDistribution,
                                             ColMajorMatrix<su2double> &matCoorDOF,
                                             ColMajorMatrix<su2double> &matCoorInt) {

  /*--- Check for which point distribution the derivatives must be computed. ---*/
  if( LGLDistribution ) {

    /*--- LGL distribution. Call the function OwnGemm to compute the Cartesian
          coordinates in the integration points. The fourth argument in the
          function call is nDim, which corresponds to the number of Cartesian
          coordinates (3 for a surface element and 2 for a volume element). ---*/
    OwnGemm(gemmDOFs2Int, jitterDOFs2Int, nIntegrationPad, nDim, nDOFs,
            nIntegrationPad, nDOFs, nIntegrationPad, true,
            lagBasisIntLGL, matCoorDOF, matCoorInt, nullptr);
  }
  else {

    /*--- Equidistant distribution. Call the function OwnGemm to compute the Cartesian
          coordinates in the integration points. The fourth argument in the
          function call is nDim, which corresponds to the number of Cartesian
          coordinates (3 for a surface element and 2 for a volume element). ---*/
    OwnGemm(gemmDOFs2Int, jitterDOFs2Int, nIntegrationPad, nDim, nDOFs,
            nIntegrationPad, nDOFs, nIntegrationPad, true,
            lagBasisIntEqui, matCoorDOF, matCoorInt, nullptr);
  }
}

void CFEMStandardTriPartition::DerivativesCoorIntPoints(const bool                         LGLDistribution,
                                                        ColMajorMatrix<su2double>          &matCoor,
                                                        vector<ColMajorMatrix<su2double> > &matDerCoor) {

  /*--- Check for which point distribution the derivatives must be computed. ---*/
  if( LGLDistribution ) {

    /*--- LGL distribution. Call the function OwnGemm 2 times to compute the derivatives
          of the Cartesian coordinates w.r.t. the two parametric coordinates. The fourth
          argument in the function call is nDim, which corresponds to the number of Cartesian
          coordinates (3 for a surface element and 2 for a volume element). ---*/
    for(unsigned short nn=0; nn<2; ++nn)
      OwnGemm(gemmDOFs2Int, jitterDOFs2Int, nIntegrationPad, nDim, nDOFs,
              nIntegrationPad, nDOFs, nIntegrationPad, true,
              derLagBasisIntLGL[nn], matCoor, matDerCoor[nn], nullptr);
  }
  else {

    /*--- Equidistant distribution. Call the function OwnGemm 2 times to compute the derivatives
          of the Cartesian coordinates w.r.t. the two parametric coordinates. The fourth
          argument in the function call is nDim, which corresponds to the number of Cartesian
          coordinates (3 for a surface element and 2 for a volume element). ---*/
    for(unsigned short nn=0; nn<2; ++nn)
      OwnGemm(gemmDOFs2Int, jitterDOFs2Int, nIntegrationPad, nDim, nDOFs,
              nIntegrationPad, nDOFs, nIntegrationPad, true,
              derLagBasisIntEqui[nn], matCoor, matDerCoor[nn], nullptr);
  }
}

passivedouble CFEMStandardTriPartition::WorkEstimateBoundaryFace(CConfig              *config,
                                                                 const unsigned short elemType) {

  /*--- Determine the number of DOFs of the neighboring element. ---*/
  const unsigned short nDOFsElem = GetNDOFsStatic(elemType, nPoly);

  /*--- TEMPORARY IMPLEMENTATION. ---*/
  return nIntegration + 0.05*nDOFsElem;
}

passivedouble CFEMStandardTriPartition::WorkEstimateInternalFace(CConfig              *config,
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

passivedouble CFEMStandardTriPartition::WorkEstimateVolume(CConfig *config) {

  /*--- TEMPORARY IMPLEMENTATION. ---*/
  return nIntegration + 0.1*nDOFs;
}

passivedouble CFEMStandardTriPartition::WorkEstimateWallFunctions(CConfig              *config,
                                                                  const unsigned short nPointsWF,
                                                                  const unsigned short elemType) {

  /*--- TEMPORARY IMPLEMENTATION. ---*/
  return 0.25*nIntegration*nPointsWF;
}