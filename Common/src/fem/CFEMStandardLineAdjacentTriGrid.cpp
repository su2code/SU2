/*!
 * \file CFEMStandardLineAdjacentTriGrid.cpp
 * \brief Functions for the class CFEMStandardLineAdjacentTriGrid.
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

#include "../../include/fem/CFEMStandardLineAdjacentTriGrid.hpp"

/*----------------------------------------------------------------------------------*/
/*            Public member functions of CFEMStandardLineAdjacentTriGrid.           */
/*----------------------------------------------------------------------------------*/

CFEMStandardLineAdjacentTriGrid::CFEMStandardLineAdjacentTriGrid(const unsigned short val_nPoly,
                                                                 const unsigned short val_orderExact,
                                                                 const unsigned short val_faceID_Elem,
                                                                 const unsigned short val_orientation,
                                                                 const bool           val_useLGL,
                                                                 CGemmBase           *val_gemm)
  : CFEMStandardTriBase(),
    CFEMStandardLineBase(val_nPoly, val_orderExact) {

  /*--- Store the faceID of the element and the orientation. ---*/
  faceID_Elem = val_faceID_Elem;
  orientation = val_orientation;

  /*--- Convert the pointer for the gemm functionality. ---*/
  gemmDOFs2Int = dynamic_cast<CGemmStandard *> (val_gemm);
  if( !gemmDOFs2Int )
    SU2_MPI::Error(string("Dynamic cast failure. This should not happen"), CURRENT_FUNCTION);

  /*--- Determine the location of the grid DOFs. ---*/
  if( val_useLGL ) LocationTriangleGridDOFsLGL(nPoly, rTriangleDOFs, sTriangleDOFs);
  else             LocationTriangleGridDOFsEquidistant(nPoly, rTriangleDOFs, sTriangleDOFs);

  /*--- Convert the 1D parametric coordinates of the integration points of the
        face to the 2D parametric coordinates of the adjacent triangle. ---*/
  ConvertCoor1DFaceTo2DTriangle(rLineInt, faceID_Elem, orientation,
                                rTriangleInt, sTriangleInt);

  /*--- Compute the corresponding Lagrangian basis functions and
        its first derivatives in the integration points. ---*/
  LagBasisIntPointsTriangle(nPoly, rTriangleDOFs,  sTriangleDOFs,
                            rTriangleInt, sTriangleInt, lagBasisInt);
  DerLagBasisIntPointsTriangle(nPoly, rTriangleDOFs, sTriangleDOFs,
                               rTriangleInt, sTriangleInt, derLagBasisInt);

  /*--- Determine the local subconnectivity of this standard element when split
        in several linear elements. Used for a.o. plotting and searching. ---*/
  CFEMStandardLineBase::SubConnLinearElements();
}

void CFEMStandardLineAdjacentTriGrid::CoorIntPoints(const bool                notUsed,
                                                    ColMajorMatrix<su2double> &matCoorDOF,
                                                    ColMajorMatrix<su2double> &matCoorInt) {

  /*--- Call the general functionality of gemmDOFs2Int with the appropriate
        arguments to compute the coordinates in the integration points
        of the face. ---*/
  gemmDOFs2Int->gemm(lagBasisInt, 2, matCoorDOF, matCoorInt, nullptr);
}

void CFEMStandardLineAdjacentTriGrid::DerivativesCoorIntPoints(const bool                         notUsed,
                                                               ColMajorMatrix<su2double>          &matCoorDOF,
                                                               vector<ColMajorMatrix<su2double> > &matDerCoorInt) {
  /*--- Call the general functionality of gemmDOFs2Int with the appropriate
        arguments to compute the derivatives of the coordinates in the
        integration points of the face. ---*/
  gemmDOFs2Int->gemm(derLagBasisInt[0], 2, matCoorDOF, matDerCoorInt[0], nullptr);
  gemmDOFs2Int->gemm(derLagBasisInt[1], 2, matCoorDOF, matDerCoorInt[1], nullptr);
}

/*----------------------------------------------------------------------------------*/
/*            Private member functions of CFEMStandardLineAdjacentTriGrid.          */
/*----------------------------------------------------------------------------------*/

void CFEMStandardLineAdjacentTriGrid::ConvertVolumeToSurfaceGradients(vector<ColMajorMatrix<su2double> > &matDerVol,
                                                                      vector<ColMajorMatrix<su2double> > &matDerFace) {

  /*--- The conversion of the gradients only takes place for elements on side 0
        of the element, i.e. orientation == 0. Check this. ---*/
  assert(orientation == 0);

  /*--- Allocate the memory for matDerFace. ---*/
  const unsigned short nRows = matDerVol[0].rows();
  const unsigned short nCols = matDerVol[0].cols();

  matDerFace.resize(1);
  matDerFace[0].resize(nRows, nCols);

  /*--- Determine on which face of the element the current surface resides and
        set the surface gradient accordingly. ---*/
  switch( faceID_Elem ) {
    case 0: {
      for(unsigned short j=0; j<nCols; ++j) {
        SU2_OMP_SIMD
        for(unsigned short i=0; i<nRows; ++i)
          matDerFace[0](i,j) = matDerVol[0](i,j);
      }
      break;
    }

    case 1: {
      for(unsigned short j=0; j<nCols; ++j) {
        SU2_OMP_SIMD_IF_NOT_AD
        for(unsigned short i=0; i<nRows; ++i)
          matDerFace[0](i,j) = matDerVol[1](i,j) - matDerVol[0](i,j);
      }
      break;
    }

    case 2: {
      for(unsigned short j=0; j<nCols; ++j) {
        SU2_OMP_SIMD_IF_NOT_AD
        for(unsigned short i=0; i<nRows; ++i)
          matDerFace[0](i,j) = -matDerVol[1](i,j);
      }
      break;
    }
  }
}
