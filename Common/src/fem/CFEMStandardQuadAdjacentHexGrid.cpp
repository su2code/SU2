/*!
 * \file CFEMStandardQuadAdjacentHexGrid.cpp
 * \brief Functions for the class CFEMStandardQuadAdjacentHexGrid.
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

#include "../../include/fem/CFEMStandardQuadAdjacentHexGrid.hpp"

/*----------------------------------------------------------------------------------*/
/*            Public member functions of CFEMStandardQuadAdjacentHexGrid.           */
/*----------------------------------------------------------------------------------*/

CFEMStandardQuadAdjacentHexGrid::CFEMStandardQuadAdjacentHexGrid(const unsigned short val_nPoly,
                                                                 const unsigned short val_orderExact,
                                                                 const unsigned short val_faceID_Elem,
                                                                 const unsigned short val_orientation,
                                                                 const bool           val_useLGL,
                                                                 CGemmBase           *val_gemm)
  : CFEMStandardHexBase(),
    CFEMStandardQuadBase(val_nPoly, val_orderExact) {

  /*--- Store the faceID of the element and the orientation. ---*/
  faceID_Elem = val_faceID_Elem;
  orientation = val_orientation;

  /*--- Convert the pointer for the gemm functionality. ---*/
  gemmDOFs2Int = dynamic_cast<CGemmFaceHex *> (val_gemm);
  if( !gemmDOFs2Int )
    SU2_MPI::Error(string("Dynamic cast failure. This should not happen"), CURRENT_FUNCTION);


  /*--- Determine the 1D parametric locations of the grid DOFs. 1D is enough,
        because a tensor product is used to obtain the 3D coordinates. ---*/
  if( val_useLGL) Location1DGridDOFsLGL(nPoly, rLineDOFs);
  else            Location1DGridDOFsEquidistant(nPoly, rLineDOFs);

  /*--- Convert the 2D coordinates of the standard quadrilateral to the parametric 
        coordinates in normal and two tangential directions of the adjacent
        hexahedron. ---*/
  vector<passivedouble> rNorm, rTang0, rTang1;
  ConvertCoor2DQuadFaceTo3DHex(rLineInt, faceID_Elem, orientation, swapTangInTensor,
                               rNorm, rTang0, rTang1);

  /*--- Create the 1D Lagrangian basis functions for the parametric
        coordinates in normal direction. ---*/
  ColMajorMatrix<passivedouble> lagN, derLagN;

  LagBasisIntPointsLine(rLineDOFs, rNorm, false, lagN);
  DerLagBasisIntPointsLine(rLineDOFs, rNorm, false, derLagN);

  /*--- Create the 1D Lagrangian basis functions for the parametric
        coordinates in first tangential direction. ---*/
  ColMajorMatrix<passivedouble> lagT0, derLagT0;

  LagBasisIntPointsLine(rLineDOFs, rTang0, true, lagT0);
  DerLagBasisIntPointsLine(rLineDOFs, rTang0, true, derLagT0);

  /*--- Create the 1D Lagrangian basis functions for the parametric
        coordinates in second tangential direction. ---*/
  ColMajorMatrix<passivedouble> lagT1, derLagT1;

  LagBasisIntPointsLine(rLineDOFs, rTang1, true, lagT1);
  DerLagBasisIntPointsLine(rLineDOFs, rTang1, true, derLagT1);

  /*--- Allocate the memory for the first index of tensorSol, etc. As this
        function is only called for 3D simulations, there are three
        contributions to the tensor products. ---*/
  tensorSol.resize(3);
  tensorDSolDr.resize(3);
  tensorDSolDs.resize(3);
  tensorDSolDt.resize(3);

  /*--- Set the components of the tensors. For the derivatives this
        depends on the faceID in the element. ---*/
  tensorSol[0] = lagN;
  tensorSol[1] = lagT0;
  tensorSol[2] = lagT1;

  if((faceID_Elem == 0) || (faceID_Elem == 1)) {
    tensorDSolDr[0] = lagN;    tensorDSolDr[1] = derLagT0; tensorDSolDr[2] = lagT1;
    tensorDSolDs[0] = lagN;    tensorDSolDs[1] = lagT0;    tensorDSolDs[2] = derLagT1;
    tensorDSolDt[0] = derLagN; tensorDSolDt[1] = lagT0;    tensorDSolDt[2] = lagT1;
  }
  else if((faceID_Elem == 2) || (faceID_Elem == 3)) {
    tensorDSolDr[0] = lagN;    tensorDSolDr[1] = derLagT0; tensorDSolDr[2] = lagT1;
    tensorDSolDs[0] = derLagN; tensorDSolDs[1] = lagT0;    tensorDSolDs[2] = lagT1;
    tensorDSolDt[0] = lagN;    tensorDSolDt[1] = lagT0;    tensorDSolDt[2] = derLagT1;
  }
  else {
    tensorDSolDr[0] = derLagN; tensorDSolDr[1] = lagT0;    tensorDSolDr[2] = lagT1;
    tensorDSolDs[0] = lagN;    tensorDSolDs[1] = derLagT0; tensorDSolDs[2] = lagT1;
    tensorDSolDt[0] = lagN;    tensorDSolDt[1] = lagT0;    tensorDSolDt[2] = derLagT1;
  }
}

void CFEMStandardQuadAdjacentHexGrid::CoorIntPoints(const bool                notUsed,
                                                    ColMajorMatrix<su2double> &matCoorDOF,
                                                    ColMajorMatrix<su2double> &matCoorInt) {

  /*--- Call the general functionality of gemmDOFs2Int with the appropriate
        arguments to compute the coordinates in the integration points
        of the face. ---*/
  gemmDOFs2Int->DOFs2Int(tensorSol, faceID_Elem, swapTangInTensor, 3, matCoorDOF, matCoorInt);
}

void CFEMStandardQuadAdjacentHexGrid::DerivativesCoorIntPoints(const bool                         notUsed,
                                                               ColMajorMatrix<su2double>          &matCoorDOF,
                                                               vector<ColMajorMatrix<su2double> > &matDerCoorInt) {

  /*--- Call the general functionality of gemmDOFs2Int with the appropriate
        arguments to compute the derivatives of the coordinates in the
        integration points of the face. ---*/
  gemmDOFs2Int->DOFs2Int(tensorDSolDr, faceID_Elem, swapTangInTensor, 3, matCoorDOF, matDerCoorInt[0]);
  gemmDOFs2Int->DOFs2Int(tensorDSolDs, faceID_Elem, swapTangInTensor, 3, matCoorDOF, matDerCoorInt[1]);
  gemmDOFs2Int->DOFs2Int(tensorDSolDt, faceID_Elem, swapTangInTensor, 3, matCoorDOF, matDerCoorInt[2]);
}

/*----------------------------------------------------------------------------------*/
/*            Private member functions of CFEMStandardQuadAdjacentHexGrid.          */
/*----------------------------------------------------------------------------------*/

void CFEMStandardQuadAdjacentHexGrid::ConvertVolumeToSurfaceGradients(vector<ColMajorMatrix<su2double> > &matDerVol,
                                                                      vector<ColMajorMatrix<su2double> > &matDerFace) {

  /*--- The conversion of the gradients only takes place for elements on side 0
        of the element, i.e. orientation == 0. Check this. ---*/
  assert(orientation == 0);

  /*--- Set the indices of the volume gradients to copy to the surface gradients. ---*/
  unsigned short ind0,  ind1;

  switch( faceID_Elem ) {
    case 0: ind0 = 0; ind1 = 1; break;
    case 1: ind0 = 1; ind1 = 0; break;
    case 2: ind0 = 2; ind1 = 0; break;
    case 3: ind0 = 0; ind1 = 2; break;
    case 4: ind0 = 1; ind1 = 2; break;
    case 5: ind0 = 2; ind1 = 1; break;
  }

  /*--- Copy the surface gradients from the appropriate volume gradients. ---*/
  matDerFace.resize(2);
  matDerFace[0] = matDerVol[ind0];
  matDerFace[1] = matDerVol[ind1];
}
