/*!
 * \file CFEMStandardQuadGrid.cpp
 * \brief Functions for the class CFEMStandardQuadGrid.
 * \author E. van der Weide
 * \version 7.0.7 "Blackbird"
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

#include "../../include/fem/CFEMStandardQuadGrid.hpp"

/*----------------------------------------------------------------------------------*/
/*             Public member functions of CFEMStandardQuadGrid.                     */
/*----------------------------------------------------------------------------------*/

CFEMStandardQuadGrid::CFEMStandardQuadGrid(const unsigned short val_nPoly,
                                           const unsigned short val_orderExact,
                                           const bool           val_surfElement)
  : CFEMStandardQuad(val_nPoly, val_orderExact) {

  /*--- Determine the number of space dimensions, which is 3 if this standard
        element is a surface element and 2 when it is a value element. ---*/
  nDim = val_surfElement ? 3 : 2;

  /*--- Compute the values of the 1D Lagrangian basis functions in the integration
        points for both the equidistant and LGL point distribution. ---*/
  LagBasisIntPointsLine(rLineDOFsEqui, rLineInt, lagBasisLineIntEqui);
  LagBasisIntPointsLine(rLineDOFsLGL,  rLineInt, lagBasisLineIntLGL);

  /*--- Compute the values of the derivatives of the 1D Lagrangian basis functions in
        the integration points for both the equidistant and LGL point distribution. ---*/
  DerLagBasisIntPointsLine(rLineDOFsEqui, rLineInt, derLagBasisLineIntEqui);
  DerLagBasisIntPointsLine(rLineDOFsLGL,  rLineInt, derLagBasisLineIntLGL);

  /*--- Create the local grid connectivities of the faces of the volume element.
        Only needed if this is a volume element. ---*/
  if( !val_surfElement ) LocalGridConnFaces();

  /*--- Determine the local subconnectivity of this standard element when split
        in several linear elements. Used for a.o. plotting and searcing. ---*/
  SubConnLinearElements();
}

void CFEMStandardQuadGrid::CoorIntPoints(const bool                LGLDistribution,
                                         ColMajorMatrix<su2double> &matCoorDOF,
                                         ColMajorMatrix<su2double> &matCoorInt) {

  /*--- Check for which point distribution the derivatives must be computed. ---*/
  if( LGLDistribution ) {

    /*--- LGL distribution. Call the function TensorProductIntegrationPoints to compute the
          Cartesian coordinates in the integration points. The first argument in the function
          call is nDim, which corresponds to the number of Cartesian coordinates (3 for a
          surface element and 2 for a volume element). ---*/
    TensorProductIntegrationPoints(nDim, lagBasisLineIntLGL, lagBasisLineIntLGL,
                                   matCoorDOF, matCoorInt, nullptr);
  }
  else {

    /*--- Equidistant distribution. Call the function TensorProductIntegrationPoints to compute the
          Cartesian coordinates in the integration points. The first argument in the function
          call is nDim, which corresponds to the number of Cartesian coordinates (3 for a
          surface element and 2 for a volume element). ---*/
    TensorProductIntegrationPoints(nDim, lagBasisLineIntEqui, lagBasisLineIntEqui,
                                   matCoorDOF, matCoorInt, nullptr);
  }
}

void CFEMStandardQuadGrid::DerivativesCoorIntPoints(const bool                         LGLDistribution,
                                                    ColMajorMatrix<su2double>          &matCoor,
                                                    vector<ColMajorMatrix<su2double> > &matDerCoor) {

  /*--- Check for which point distribution the derivatives must be computed. ---*/
  if( LGLDistribution ) {

    /*--- LGL distribution. Call the function TensorProductIntegrationPoints 2 times to compute the
          derivatives of the Cartesian coordinates w.r.t. the two parametric coordinates. The first
          argument in the function call is nDim, which corresponds to the number of Cartesian
          coordinates (3 for a surface element and 2 for a volume element). ---*/
    TensorProductIntegrationPoints(nDim, derLagBasisLineIntLGL, lagBasisLineIntLGL,
                                   matCoor, matDerCoor[0], nullptr);
    TensorProductIntegrationPoints(nDim, lagBasisLineIntLGL, derLagBasisLineIntLGL,
                                   matCoor, matDerCoor[1], nullptr);
  }
  else {

    /*--- Equidistant distribution. Call the function TensorProductIntegrationPoints 3 times to compute the
          derivatives of the Cartesian coordinates w.r.t. the three parametric coordinates. The first
          argument in the function call is nDim, which corresponds to the number of Cartesian
          coordinates (3 for a surface element and 2 for a volume element). ---*/
    TensorProductIntegrationPoints(nDim, derLagBasisLineIntEqui, lagBasisLineIntEqui,
                                   matCoor, matDerCoor[0], nullptr);
    TensorProductIntegrationPoints(nDim, lagBasisLineIntEqui, derLagBasisLineIntEqui,
                                   matCoor, matDerCoor[1], nullptr);
  }
}

/*----------------------------------------------------------------------------------*/
/*             Private member functions of CFEMStandardQuadGrid.                    */
/*----------------------------------------------------------------------------------*/

void CFEMStandardQuadGrid::LocalGridConnFaces(void) {

  /*--- Allocate the first index of gridConnFaces, which is equal to the number
        of faces of the quadrilateral, which is 4. Reserve memory for the second
        index afterwards. ---*/
  gridConnFaces.resize(4);

  gridConnFaces[0].reserve(nPoly+1);
  gridConnFaces[1].reserve(nPoly+1);
  gridConnFaces[2].reserve(nPoly+1);
  gridConnFaces[3].reserve(nPoly+1);

  /*--- Define the corner vertices of the quadrilateral. ---*/
  const unsigned short n0 = 0, n1 = nPoly, n2 = nDOFs-1, n3 = nPoly*(nPoly+1);

  /*--- For a quad element the faces are lines. Loop over the nodes of the
        lines to set the connectivity. Make sure that the element
        is to the left of the faces. ---*/
  for(signed short i=n0; i<=n1; ++i)          gridConnFaces[0].push_back(i);
  for(signed short i=n1; i<=n2; i+=(nPoly+1)) gridConnFaces[1].push_back(i);
  for(signed short i=n2; i>=n3; --i)          gridConnFaces[2].push_back(i);
  for(signed short i=n3; i>=n0; i-=(nPoly+1)) gridConnFaces[3].push_back(i);
}

void CFEMStandardQuadGrid::SubConnLinearElements(void) {

  /*--- The quadrilateral is split into several linear quads.
        Set the VTK sub-types accordingly. ---*/
  VTK_SubType1 = QUADRILATERAL;
  VTK_SubType2 = NONE;

  /*--- Determine the local subconnectivity of the quadrilateral element used for
        plotting purposes. Note that the connectivity of the linear subelements
        obey the VTK connectivity rule of a quadrilateral, which is different
        from the connectivity for the high order quadrilateral. ---*/
  unsigned short nnPoly = max(nPoly,(unsigned short) 1);
  for(unsigned short j=0; j<nnPoly; ++j) {
    unsigned short jj = j*(nnPoly+1);
    for(unsigned short i=0; i<nnPoly; ++i) {
      const unsigned short n0 = jj + i;
      const unsigned short n1 = n0 + 1;
      const unsigned short n2 = n1 + nPoly+1;
      const unsigned short n3 = n2 - 1;

      subConn1ForPlotting.push_back(n0);
      subConn1ForPlotting.push_back(n1);
      subConn1ForPlotting.push_back(n2);
      subConn1ForPlotting.push_back(n3);
    }
  }
}
