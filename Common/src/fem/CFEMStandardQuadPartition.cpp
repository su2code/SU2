/*!
 * \file CFEMStandardQuadPartition.cpp
 * \brief Functions for the class CFEMStandardQuadPartition.
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

#include "../../include/fem/CFEMStandardQuadPartition.hpp"

/*----------------------------------------------------------------------------------*/
/*             Public member functions of CFEMStandardQuadPartition.                */
/*----------------------------------------------------------------------------------*/

CFEMStandardQuadPartition::CFEMStandardQuadPartition(const unsigned short val_nPoly,
                                                     const unsigned short val_orderExact,
                                                     const bool           val_surfElement)
  : CFEMStandardQuadBase(val_nPoly, val_orderExact) {

  /*--- Determine the number of space dimensions, which is 3 if this standard
        element is a surface element and 2 when it is a volume element. ---*/
  nDim = val_surfElement ? 3 : 2;

  /*--- Determine the 1D parametric locations of the grid DOFs. 1D is enough,
        because a tensor product is used to obtain the 2D coordinates. ---*/
  Location1DGridDOFsEquidistant(nPoly, rLineDOFsEqui);
  Location1DGridDOFsLGL(nPoly, rLineDOFsLGL);

  /*--- Compute the values of the 1D Lagrangian basis functions in the integration
        points for both the equidistant and LGL point distribution. ---*/
  LagBasisIntPointsLine(rLineDOFsEqui, rLineInt, true, lagBasisLineIntEqui);
  LagBasisIntPointsLine(rLineDOFsLGL,  rLineInt, true, lagBasisLineIntLGL);

  /*--- Compute the values of the derivatives of the 1D Lagrangian basis functions in
        the integration points for both the equidistant and LGL point distribution. ---*/
  DerLagBasisIntPointsLine(rLineDOFsEqui, rLineInt, true, derLagBasisLineIntEqui);
  DerLagBasisIntPointsLine(rLineDOFsLGL,  rLineInt, true, derLagBasisLineIntLGL);

  /*--- Create the local grid connectivities of the faces of the volume element.
        Only needed if this is a volume element. ---*/
  if( !val_surfElement ) LocalGridConnFaces();

  /*--- Determine the local subconnectivity of this standard element when split
        in several linear elements. Used for a.o. plotting and searching. ---*/
  CFEMStandardQuadBase::SubConnLinearElements();

  /*--- Set the function pointers for the tensor product multiplications to
        determine the data in the volume integration points. ---*/
  SetFunctionPointerVolumeDataQuad(nDOFs1D, nInt1D, TensorProductDataVolIntPoints);
}

void CFEMStandardQuadPartition::CoorIntPoints(const bool                LGLDistribution,
                                              ColMajorMatrix<su2double> &matCoorDOF,
                                              ColMajorMatrix<su2double> &matCoorInt) {

  /*--- Check for which point distribution the derivatives must be computed. ---*/
  if( LGLDistribution ) {

    /*--- LGL distribution. Call the function TensorProductVolumeDataQuad to compute the
          Cartesian coordinates in the integration points. The second argument in the function
          call is nDim, which corresponds to the number of Cartesian coordinates (3 for a
          surface element and 2 for a volume element). ---*/
    TensorProductVolumeDataQuad(TensorProductDataVolIntPoints, nDim, nDOFs1D, nInt1D,
		                    lagBasisLineIntLGL, lagBasisLineIntLGL, matCoorDOF,
                                matCoorInt, true, nullptr);
  }
  else {

    /*--- Equidistant distribution. Call the function TensorProductVolumeDataQuad to compute the
          Cartesian coordinates in the integration points. The second argument in the function
          call is nDim, which corresponds to the number of Cartesian coordinates (3 for a
          surface element and 2 for a volume element). ---*/
    TensorProductVolumeDataQuad(TensorProductDataVolIntPoints, nDim, nDOFs1D, nInt1D,
                                lagBasisLineIntEqui, lagBasisLineIntEqui, matCoorDOF,
                                matCoorInt, true, nullptr);
  }
}

void CFEMStandardQuadPartition::DerivativesCoorIntPoints(const bool                         LGLDistribution,
                                                         ColMajorMatrix<su2double>          &matCoor,
                                                         vector<ColMajorMatrix<su2double> > &matDerCoor) {

  /*--- Check for which point distribution the derivatives must be computed. ---*/
  if( LGLDistribution ) {

    /*--- LGL distribution. Call the function TensorProductVolumeDataQuad 2 times to compute the
          derivatives of the Cartesian coordinates w.r.t. the two parametric coordinates. The second
          argument in the function call is nDim, which corresponds to the number of Cartesian
          coordinates (3 for a surface element and 2 for a volume element). ---*/
    TensorProductVolumeDataQuad(TensorProductDataVolIntPoints,nDim, nDOFs1D, nInt1D,
                                derLagBasisLineIntLGL, lagBasisLineIntLGL, matCoor,
                                matDerCoor[0], true, nullptr);
    TensorProductVolumeDataQuad(TensorProductDataVolIntPoints, nDim, nDOFs1D, nInt1D,
                                lagBasisLineIntLGL, derLagBasisLineIntLGL, matCoor,
                                matDerCoor[1], true, nullptr);
  }
  else {

    /*--- Equidistant distribution. Call the function TensorProductVolumeDataQuad 2 times to compute the
          derivatives of the Cartesian coordinates w.r.t. the two parametric coordinates. The second
          argument in the function call is nDim, which corresponds to the number of Cartesian
          coordinates (3 for a surface element and 2 for a volume element). ---*/
    TensorProductVolumeDataQuad(TensorProductDataVolIntPoints, nDim, nDOFs1D, nInt1D,
                                derLagBasisLineIntEqui, lagBasisLineIntEqui, matCoor,
                                matDerCoor[0], true, nullptr);
    TensorProductVolumeDataQuad(TensorProductDataVolIntPoints, nDim, nDOFs1D, nInt1D,
                                lagBasisLineIntEqui, derLagBasisLineIntEqui, matCoor,
                                matDerCoor[1], true, nullptr);
  }
}

passivedouble CFEMStandardQuadPartition::WorkEstimateBoundaryFace(CConfig              *config,
                                                                  const unsigned short elemType) {

  /*--- Determine the number of DOFs of the neighboring element. ---*/
  const unsigned short nDOFsElem = GetNDOFsStatic(elemType, nPoly);

  /*--- TEMPORARY IMPLEMENTATION. ---*/
  return nIntegration + 0.05*nDOFsElem;
}

passivedouble CFEMStandardQuadPartition::WorkEstimateInternalFace(CConfig              *config,
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

passivedouble CFEMStandardQuadPartition::WorkEstimateVolume(CConfig *config) {

  /*--- TEMPORARY IMPLEMENTATION. ---*/
  return nIntegration + 0.1*nDOFs;
}

passivedouble CFEMStandardQuadPartition::WorkEstimateWallFunctions(CConfig              *config,
                                                                   const unsigned short nPointsWF,
                                                                   const unsigned short elemType) {

  /*--- TEMPORARY IMPLEMENTATION. ---*/
  return 0.25*nIntegration*nPointsWF;
}