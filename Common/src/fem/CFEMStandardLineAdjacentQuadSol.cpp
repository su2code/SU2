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

  /*--- Store the pointers for the gemm functionalities. ---*/
  gemmDOFs2Int = val_gemm_1;
  gemmInt2DOFs = val_gemm_2;

  /*--- Create the standard 1D Legendre basis functions in the integration points. ---*/
  ColMajorMatrix<passivedouble> legInt1D, derLegInt1D;
  legInt1D.resize(nIntegrationPad, nPoly+1);    legInt1D.setConstant(0.0);
  derLegInt1D.resize(nIntegrationPad, nPoly+1); derLegInt1D.setConstant(0.0);

  Vandermonde1D(nPoly, rLineInt, legInt1D);
  GradVandermonde1D(nPoly, rLineInt, derLegInt1D);

  /*--- Create the standard 1D Legendre basis function on the left end (r == -1). ---*/
  vector<passivedouble> rBound(1, -1.0);
  ColMajorMatrix<passivedouble> legM1(1, nPoly+1), derLegM1(1, nPoly+1);

  Vandermonde1D(nPoly, rBound, legM1);
  GradVandermonde1D(nPoly, rBound, derLegM1);

  /*--- Create the standard 1D Legendre basis function on the right end (r == 1). ---*/
  rBound[0] = 1.0;
  ColMajorMatrix<passivedouble> legP1(1, nPoly+1), derLegP1(1, nPoly+1);

  Vandermonde1D(nPoly, rBound, legP1);
  GradVandermonde1D(nPoly, rBound, derLegP1);

  /*--- Create the tensor components of the Legendre basis functions for the
        current situation from the standard basis functions created above. ---*/
  CreateTensorContributionsLineAdjQuad(nIntegration, nPoly+1, val_faceID_Elem, val_orientation,
                                       legInt1D, derLegInt1D, legM1, legP1, derLegM1, derLegP1,
                                       tensorSol, tensorDSolDr, tensorDSolDs);

}
