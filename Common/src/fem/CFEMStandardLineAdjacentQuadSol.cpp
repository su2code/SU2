/*!
 * \file CFEMStandardLineAdjacentQuadSol.cpp
 * \brief Functions for the class CFEMStandardLineAdjacentQuadSol.
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

#include "../../include/fem/CFEMStandardLineAdjacentQuadSol.hpp"

/*----------------------------------------------------------------------------------*/
/*            Public member functions of CFEMStandardLineAdjacentQuadSol.           */
/*----------------------------------------------------------------------------------*/

CFEMStandardLineAdjacentQuadSol::CFEMStandardLineAdjacentQuadSol(const unsigned short val_nPoly,
                                                                 const unsigned short val_orderExact,
                                                                 const unsigned short val_faceID_Elem,
                                                                 const unsigned short val_orientation,
                                                                 CGemmBase           *val_gemm_1,
                                                                 CGemmBase           *val_gemm_2)
  : CFEMStandardQuadBase(),
    CFEMStandardLineBase(val_nPoly, val_orderExact) {

  /*--- Store the faceID of the element and the orientation. ---*/
  faceID_Elem = val_faceID_Elem;
  orientation = val_orientation;

  /*--- Convert the pointers for the gemm functionalities. ---*/
  gemmDOFs2Int = dynamic_cast<CGemmFaceQuad *> (val_gemm_1);
  gemmInt2DOFs = dynamic_cast<CGemmFaceQuad *> (val_gemm_2);
  if(!gemmDOFs2Int || !gemmInt2DOFs)
    SU2_MPI::Error(string("Dynamic cast failure. This should not happen"), CURRENT_FUNCTION);

  /*--- Convert the 1D coordinates of the standard line to the parametric 
        coordinates in normal and tangential direction of the adjacent
        quadrilateral. ---*/
  vector<passivedouble> rNormal, rTangential;
  ConvertCoor1DFaceTo2DQuad(rLineInt, faceID_Elem, orientation, rNormal, rTangential);

  /*--- Create the 1D Legendre basis functions for the parametric
        coordinates in normal direction. ---*/
  ColMajorMatrix<passivedouble> legN(1, nPoly+1), derLegN(1, nPoly+1);

  Vandermonde1D(nPoly, rNormal, legN);
  GradVandermonde1D(nPoly, rNormal, derLegN);

  /*--- Create the 1D Legendre basis functions for the parametric
        coordinates in tangential direction. ---*/
  ColMajorMatrix<passivedouble> legT, derLegT;
  legT.resize(nIntegrationPad, nPoly+1);    legT.setConstant(0.0);
  derLegT.resize(nIntegrationPad, nPoly+1); derLegT.setConstant(0.0);

  Vandermonde1D(nPoly, rTangential, legT);
  GradVandermonde1D(nPoly, rTangential, derLegT);

  /*--- Allocate the memory for the first index of tensorSol, etc. As this
        function is only called for 2D simulations, there are two
        contributions to the tensor products. ---*/
  tensorSol.resize(2);
  tensorDSolDr.resize(2);
  tensorDSolDs.resize(2);

  /*--- Set the components of the tensors. For the derivatives this
        depend on the faceID in the element. ---*/
  tensorSol[0] = legN;
  tensorSol[1] = legT;

  if((faceID_Elem == 0) || (faceID_Elem == 2)) {
    tensorDSolDr[0] = legN;    tensorDSolDr[1] = derLegT;
    tensorDSolDs[0] = derLegN; tensorDSolDs[1] = legT;
  }
  else {
    tensorDSolDr[0] = derLegN; tensorDSolDr[1] = legT;
    tensorDSolDs[0] = legN;    tensorDSolDs[1] = derLegT;
  }
}
