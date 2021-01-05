/*!
 * \file CFEMStandardQuadAdjacentPyraGrid.cpp
 * \brief Functions for the class CFEMStandardQuadAdjacentPyraGrid.
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

  /*--- Store the pointer for the gemm functionality. ---*/
  gemmDOFs2Int = val_gemm;

  /*--- Determine the location of the grid DOFs. ---*/
  vector<passivedouble> rPyraDOFs, sPyraDOFs, tPyraDOFs;
  if( val_useLGL ) LocationPyramidGridDOFsLGL(nPoly, rPyraDOFs, sPyraDOFs, tPyraDOFs);
  else             LocationPyramidGridDOFsEquidistant(nPoly, rPyraDOFs, sPyraDOFs, tPyraDOFs);

  /*--- Convert the 2D parametric coordinates of the integration points of the
        quadrilateral face to the 3D parametric coordinates of the adjacent pyramid. ---*/
  vector<passivedouble> rInt, sInt, tInt;
  ConvertCoor2DQuadFaceTo3DPyra(rLineInt, val_faceID_Elem, val_orientation, rInt, sInt, tInt);

  /*--- Compute the corresponding Lagrangian basis functions and
        its first derivatives in the integration points. ---*/
  LagBasisIntPointsPyra(nPoly, rPyraDOFs, sPyraDOFs, tPyraDOFs,
                        rInt, sInt, tInt, lagBasisInt);
  DerLagBasisIntPointsPyra(nPoly, rPyraDOFs, sPyraDOFs, tPyraDOFs,
                           rInt, sInt, tInt, derLagBasisInt);
}
