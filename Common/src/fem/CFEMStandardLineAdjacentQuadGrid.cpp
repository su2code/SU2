/*!
 * \file CFEMStandardLineAdjacentQuadGrid.cpp
 * \brief Functions for the class CFEMStandardLineAdjacentQuadGrid.
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

#include "../../include/fem/CFEMStandardLineAdjacentQuadGrid.hpp"

/*----------------------------------------------------------------------------------*/
/*            Public member functions of CFEMStandardLineAdjacentQuadGrid.          */
/*----------------------------------------------------------------------------------*/

CFEMStandardLineAdjacentQuadGrid::CFEMStandardLineAdjacentQuadGrid(const unsigned short val_nPoly,
                                                                   const unsigned short val_orderExact,
                                                                   const unsigned short val_faceID_Elem,
                                                                   const unsigned short val_orientation,
                                                                   const bool           val_useLGL,
                                                                   CGemmBase           *val_gemm)
  : CFEMStandardQuadBase(),
    CFEMStandardLineBase(val_nPoly, val_orderExact) {

  /*--- Store the pointer for the gemm functionality. ---*/
  gemmDOFs2Int = val_gemm;

  /*--- Determine the 1D parametric locations of the grid DOFs. 1D is enough,
        because a tensor product is used to obtain the 2D coordinates. ---*/
  if( val_useLGL) Location1DGridDOFsLGL(nPoly, rLineDOFs);
  else            Location1DGridDOFsEquidistant(nPoly, rLineDOFs);

  /*--- Create the standard 1D Lagrangian basis functions. ---*/
  ColMajorMatrix<passivedouble> lagInt1D, derLagInt1D;
  ColMajorMatrix<passivedouble> lagM1, lagP1, derLagM1, derLagP1;

  LagBasisIntPointsLine(rLineDOFs, rLineInt, lagInt1D);
  DerLagBasisIntPointsLine(rLineDOFs, rLineInt, derLagInt1D);

  vector<passivedouble> rBound(1, -1.0);
  LagBasisIntPointsLine(rLineDOFs, rBound, lagM1);
  DerLagBasisIntPointsLine(rLineDOFs, rBound, derLagM1);

  rBound[0] = 1.0;
  LagBasisIntPointsLine(rLineDOFs, rBound, lagP1);
  DerLagBasisIntPointsLine(rLineDOFs, rBound, derLagP1);

  /*--- Create the tensor components of the Lagrangian basis functions for the
        current situation from the standard basis functions created above. ---*/
  CreateTensorContributionsLineAdjQuad(nIntegration, nPoly+1, val_faceID_Elem, val_orientation,
                                       lagInt1D, derLagInt1D, lagM1, lagP1, derLagM1, derLagP1,
                                       tensorSol, tensorDSolDr, tensorDSolDs);

}
