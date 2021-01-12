/*!
 * \file CFEMStandardTriAdjacentTetGrid.cpp
 * \brief Functions for the class CFEMStandardTriAdjacentTetGrid.
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
}

void CFEMStandardTriAdjacentTetGrid::CoorIntPoints(const bool                notUsed,
                                                   ColMajorMatrix<su2double> &matCoorDOF,
                                                   ColMajorMatrix<su2double> &matCoorInt) {

  /*--- Call the general functionality of gemmDOFs2Int with the appropriate
        arguments to compute the coordinates in the integration points
        of the face. ---*/
  gemmDOFs2Int->DOFs2Int(lagBasisInt, 3, matCoorDOF, matCoorInt, nullptr);
}
