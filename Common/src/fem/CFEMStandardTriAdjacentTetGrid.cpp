/*!
 * \file CFEMStandardTriAdjacentTetGrid.cpp
 * \brief Functions for the class CFEMStandardTriAdjacentTetGrid.
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

#include "../../include/fem/CFEMStandardTriAdjacentTetGrid.hpp"

/*----------------------------------------------------------------------------------*/
/*            Public member functions of CFEMStandardTriAdjacentTetGrid.            */
/*----------------------------------------------------------------------------------*/

CFEMStandardTriAdjacentTetGrid::CFEMStandardTriAdjacentTetGrid(const unsigned short val_nPoly,
                                                               const unsigned short val_orderExact,
                                                               const unsigned short val_faceID_Elem,
                                                               const unsigned short val_orientation,
                                                               const bool           val_useLGL,
                                                               CGemmBase           *val_gemm)
  : CFEMStandardTetBase(),
    CFEMStandardTriBase(val_nPoly, val_orderExact) {

  /*--- Store the faceID of the element and the orientation. ---*/
  faceID_Elem = val_faceID_Elem;
  orientation = val_orientation;

  /*--- Convert the pointer for the gemm functionality. ---*/
  gemmDOFs2Int = dynamic_cast<CGemmStandard *> (val_gemm);
  if( !gemmDOFs2Int )
    SU2_MPI::Error(string("Dynamic cast failure. This should not happen"), CURRENT_FUNCTION);

  /*--- Determine the location of the grid DOFs. ---*/
  vector<passivedouble> rTetDOFs, sTetDOFs, tTetDOFs;
  if( val_useLGL) LocationTetGridDOFsLGL(nPoly, rTetDOFs, sTetDOFs, tTetDOFs);
  else            LocationTetGridDOFsEquidistant(nPoly, rTetDOFs, sTetDOFs, tTetDOFs);

  /*--- Convert the 2D parametric coordinates of the integration points of the
        triangular face to the 3D parametric coordinates of the adjacent tetrahedron. ---*/
  vector<passivedouble> rInt, sInt, tInt;
  ConvertCoor2DTriFaceTo3DTet(rTriangleInt, sTriangleInt, val_faceID_Elem,
                              val_orientation, rInt, sInt, tInt);

  /*--- Compute the corresponding Lagrangian basis functions and
        its first derivatives in the integration points. ---*/
  LagBasisIntPointsTet(nPoly, rTetDOFs, sTetDOFs, tTetDOFs,
                       rInt, sInt, tInt, lagBasisInt);
  DerLagBasisIntPointsTet(nPoly, rTetDOFs, sTetDOFs, tTetDOFs, 
                          rInt, sInt, tInt, derLagBasisInt);

  /*--- Determine the local subconnectivity of this standard element when split
        in several linear elements. Used for a.o. plotting and searching. ---*/
  CFEMStandardTriBase::SubConnLinearElements();
}

void CFEMStandardTriAdjacentTetGrid::CoorIntPoints(const bool                notUsed,
                                                   ColMajorMatrix<su2double> &matCoorDOF,
                                                   ColMajorMatrix<su2double> &matCoorInt) {

  /*--- Call the general functionality of gemmDOFs2Int with the appropriate
        arguments to compute the coordinates in the integration points
        of the face. ---*/
  gemmDOFs2Int->DOFs2Int(lagBasisInt, 3, matCoorDOF, matCoorInt, nullptr);
}

void CFEMStandardTriAdjacentTetGrid::DerivativesCoorIntPoints(const bool                         notUsed,
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
/*            Private member functions of CFEMStandardTriAdjacentTetGrid.           */
/*----------------------------------------------------------------------------------*/

void CFEMStandardTriAdjacentTetGrid::ConvertVolumeToSurfaceGradients(vector<ColMajorMatrix<su2double> > &matDerVol,
                                                                     vector<ColMajorMatrix<su2double> > &matDerFace) {

  /*--- The conversion of the gradients only takes place for elements on side 0
        of the element, i.e. orientation == 0. Check this. ---*/
  assert(orientation == 0);

  /*--- Set the indices of the volume gradients to copy to the surface gradients.
        Note that for face 3 this is not correct, but is used as a first step. ---*/
  unsigned short ind0,  ind1;

  switch( faceID_Elem ) {
    case 0: ind0 = 0; ind1 = 1; break;
    case 1: ind0 = 2; ind1 = 0; break;
    case 2: ind0 = 1; ind1 = 2; break;
    case 3: ind0 = 2; ind1 = 1; break;
  }

  /*--- Copy the surface gradients from the appropriate volume gradients. ---*/
  matDerFace.resize(2);
  matDerFace[0] = matDerVol[ind0];
  matDerFace[1] = matDerVol[ind1];

  /*--- Correct the surface gradients for face 3. ---*/
  if(faceID_Elem == 3) {
    const unsigned short nCols = matDerFace[0].cols();
    const unsigned short nRows = matDerFace[0].rows();
    for(unsigned short j=0; j<nCols; ++j) {
      SU2_OMP_SIMD_IF_NOT_AD
      for(unsigned short i=0; i<nRows; ++i) {
        matDerFace[0](i,j) -= matDerVol[0](i,j);
        matDerFace[1](i,j) -= matDerVol[0](i,j);
      }
    }
  }
}
