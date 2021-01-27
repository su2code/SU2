/*!
 * \file CFEMStandardQuadAdjacentPyraGrid.cpp
 * \brief Functions for the class CFEMStandardQuadAdjacentPyraGrid.
 * \author E. van der Weide
 * \version 7.1.0 "Blackbird"
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

#include "../../include/fem/CFEMStandardQuadAdjacentPyraGrid.hpp"

/*----------------------------------------------------------------------------------*/
/*            Public member functions of CFEMStandardQuadAdjacentPyraGrid.          */
/*----------------------------------------------------------------------------------*/

CFEMStandardQuadAdjacentPyraGrid::CFEMStandardQuadAdjacentPyraGrid(const unsigned short val_nPoly,
                                                                   const unsigned short val_orderExact,
                                                                   const unsigned short val_faceID_Elem,
                                                                   const unsigned short val_orientation,
                                                                   const bool           val_useLGL,
                                                                   CGemmBase           *val_gemm)
  : CFEMStandardPyraBase(),
    CFEMStandardQuadBase(val_nPoly, val_orderExact) {

  /*--- Store the faceID of the element and the orientation. ---*/
  faceID_Elem = val_faceID_Elem;
  orientation = val_orientation;

  /*--- Convert the pointer for the gemm functionality. ---*/
  gemmDOFs2Int = dynamic_cast<CGemmStandard *> (val_gemm);
  if( !gemmDOFs2Int )
    SU2_MPI::Error(string("Dynamic cast failure. This should not happen"), CURRENT_FUNCTION);

  /*--- Determine the location of the grid DOFs. ---*/
  vector<passivedouble> rPyraDOFs, sPyraDOFs, tPyraDOFs;
  if( val_useLGL ) LocationPyramidGridDOFsLGL(nPoly, rPyraDOFs, sPyraDOFs, tPyraDOFs);
  else             LocationPyramidGridDOFsEquidistant(nPoly, rPyraDOFs, sPyraDOFs, tPyraDOFs);

  /*--- Convert the 2D parametric coordinates of the integration points of the
        quadrilateral face to the 3D parametric coordinates of the adjacent pyramid. ---*/
  vector<passivedouble> rInt, sInt, tInt;
  ConvertCoor2DQuadFaceTo3DPyra(rLineInt, faceID_Elem, orientation, rInt, sInt, tInt);

  /*--- Compute the corresponding Lagrangian basis functions and
        its first derivatives in the integration points. ---*/
  LagBasisIntPointsPyra(nPoly, rPyraDOFs, sPyraDOFs, tPyraDOFs,
                        rInt, sInt, tInt, lagBasisInt);
  DerLagBasisIntPointsPyra(nPoly, rPyraDOFs, sPyraDOFs, tPyraDOFs,
                           rInt, sInt, tInt, derLagBasisInt);

  /*--- Determine the local subconnectivity of this standard element when split
        in several linear elements. Used for a.o. plotting and searching. ---*/
  CFEMStandardQuadBase::SubConnLinearElements();
}

void CFEMStandardQuadAdjacentPyraGrid::CoorIntPoints(const bool                notUsed,
                                                     ColMajorMatrix<su2double> &matCoorDOF,
                                                     ColMajorMatrix<su2double> &matCoorInt) {

  /*--- Call the general functionality of gemmDOFs2Int with the appropriate
        arguments to compute the coordinates in the integration points
        of the face. ---*/
  gemmDOFs2Int->DOFs2Int(lagBasisInt, 3, matCoorDOF, matCoorInt, nullptr);
}

void CFEMStandardQuadAdjacentPyraGrid::DerivativesCoorIntPoints(const bool                         notUsed,
                                                                ColMajorMatrix<su2double>          &matCoorDOF,
                                                                vector<ColMajorMatrix<su2double> > &matDerCoorInt) {
  /*--- Call the general functionality of gemmDOFs2Int with the appropriate
        arguments to compute the derivatives of the coordinates in the
        integration points of the face. ---*/
  gemmDOFs2Int->DOFs2Int(derLagBasisInt[0], 3, matCoorDOF, matDerCoorInt[0], nullptr);
  gemmDOFs2Int->DOFs2Int(derLagBasisInt[1], 3, matCoorDOF, matDerCoorInt[1], nullptr);
  gemmDOFs2Int->DOFs2Int(derLagBasisInt[2], 3, matCoorDOF, matDerCoorInt[2], nullptr);
}

/*----------------------------------------------------------------------------------*/
/*            Private member functions of CFEMStandardQuadAdjacentPyraGrid.         */
/*----------------------------------------------------------------------------------*/

void CFEMStandardQuadAdjacentPyraGrid::ConvertVolumeToSurfaceGradients(vector<ColMajorMatrix<su2double> > &matDerVol,
                                                                       vector<ColMajorMatrix<su2double> > &matDerFace) {

  /*--- The conversion of the gradients only takes place for elements on side 0
        of the element, i.e. orientation == 0. Furthermore, for a pyramid the
        face ID of the element should 0. Check this. ---*/
  assert(orientation == 0);
  assert(faceID_Elem == 0);

  /*--- Copy the appropriate data. ---*/
  matDerFace.resize(2);
  matDerFace[0] = matDerVol[0];
  matDerFace[1] = matDerVol[1];
}
