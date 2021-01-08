/*!
 * \file CFEMStandardHexGrid.cpp
 * \brief Functions for the class CFEMStandardHexGrid.
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

#include "../../include/fem/CFEMStandardHexGrid.hpp"

/*----------------------------------------------------------------------------------*/
/*             Public member functions of CFEMStandardHexGrid.                      */
/*----------------------------------------------------------------------------------*/

CFEMStandardHexGrid::CFEMStandardHexGrid(const unsigned short val_nPolyGrid,
                                         const unsigned short val_nPolySol,
                                         const unsigned short val_orderExact,
                                         const unsigned short val_locGridDOFs)

  : CFEMStandardHexBase(val_nPolyGrid,val_orderExact) {

  /*--- Determine the 1D parametric locations of the grid DOFs. 1D is enough,
        because a tensor product is used to obtain the 3D coordinates. ---*/
  if(val_locGridDOFs == LGL) Location1DGridDOFsLGL(nPoly, rLineDOFs);
  else                       Location1DGridDOFsEquidistant(nPoly, rLineDOFs);

  /*--- Compute the 1D parametric coordinates of the solution DOFs. Only
        different from rLineDOFs when a different polynomial degree is
        used for the grid and solution. ---*/
  if(val_locGridDOFs == LGL) Location1DGridDOFsLGL(val_nPolySol, rLineSolDOFs);
  else                       Location1DGridDOFsEquidistant(val_nPolySol, rLineSolDOFs);

  /*--- Compute the 1D Lagrangian basis functions and its first
        and second derivatives in the integration points. ---*/
  LagBasisIntPointsLine(rLineDOFs, rLineInt, lagBasisLineInt);
  DerLagBasisIntPointsLine(rLineDOFs, rLineInt, derLagBasisLineInt);
  HesLagBasisIntPointsLine(rLineDOFs, rLineInt, hesLagBasisLineInt);

  /*--- Call LagBasisIntPointsLine and DerLagBasisIntPointsLine with the
        solution DOFs as argument to compute the Lagrangian basis functions
        and its derivatives in the solution DOFs. ---*/
  LagBasisIntPointsLine(rLineDOFs, rLineSolDOFs, lagBasisLineSolDOFs);
  DerLagBasisIntPointsLine(rLineDOFs, rLineSolDOFs, derLagBasisLineSolDOFs);

  /*--- Determine the local subconnectivity of this standard element when split
        in several linear elements. Used for a.o. plotting and searcing. ---*/
  CFEMStandardHexBase::SubConnLinearElements();

  /*--- Set the function pointers for the tensor product multiplications to
        determine the data in the volume integration points and volume
        nodal solution DOFs. ---*/
  SetFunctionPointerVolumeDataHex(nDOFs1D, nInt1D, TensorProductDataVolIntPoints);
  SetFunctionPointerVolumeDataHex(nDOFs1D, val_nPolySol+1, TensorProductDataVolSolDOFs);
}

void CFEMStandardHexGrid::CoorIntPoints(const bool                notUsed,
                                        ColMajorMatrix<su2double> &matCoorDOF,
                                        ColMajorMatrix<su2double> &matCoorInt) {

  /*--- Call the function TensorProductVolumeDataHex to compute the
        Cartesian coordinates in the integration points. ---*/
  TensorProductVolumeDataHex(TensorProductDataVolIntPoints, 3, nDOFs1D, nInt1D,
                             lagBasisLineInt, lagBasisLineInt, lagBasisLineInt,
                             matCoorDOF, matCoorInt, nullptr);
}

void CFEMStandardHexGrid::DerivativesCoorIntPoints(const bool                         notUsed,
                                                   ColMajorMatrix<su2double>          &matCoor,
                                                   vector<ColMajorMatrix<su2double> > &matDerCoor) {

  /*--- Call the function TensorProductVolumeDataHex 3 times to compute the derivatives
        of the Cartesian coordinates w.r.t. the three parametric coordinates. ---*/
  TensorProductVolumeDataHex(TensorProductDataVolIntPoints, 3, nDOFs1D, nInt1D,
                             derLagBasisLineInt, lagBasisLineInt, lagBasisLineInt,
                             matCoor, matDerCoor[0], nullptr);
  TensorProductVolumeDataHex(TensorProductDataVolIntPoints, 3, nDOFs1D, nInt1D,
                             lagBasisLineInt, derLagBasisLineInt, lagBasisLineInt,
                             matCoor, matDerCoor[1], nullptr);
  TensorProductVolumeDataHex(TensorProductDataVolIntPoints, 3, nDOFs1D, nInt1D,
                             lagBasisLineInt, lagBasisLineInt, derLagBasisLineInt,
                             matCoor, matDerCoor[2], nullptr);
}

void CFEMStandardHexGrid::Derivatives2ndCoorIntPoints(ColMajorMatrix<su2double>          &matCoor,
                                                      vector<ColMajorMatrix<su2double> > &matDer2ndCoor) {

  /*--- Call the function TensorProductVolumeDataHex 6 times to compute the 2nd derivatives
        of the Cartesian coordinates w.r.t. the three parametric coordinates. ---*/
  TensorProductVolumeDataHex(TensorProductDataVolIntPoints, 3, nDOFs1D, nInt1D,
                             hesLagBasisLineInt, lagBasisLineInt, lagBasisLineInt,
                             matCoor, matDer2ndCoor[0], nullptr);
  TensorProductVolumeDataHex(TensorProductDataVolIntPoints, 3, nDOFs1D, nInt1D,
                             lagBasisLineInt, hesLagBasisLineInt, lagBasisLineInt,
                             matCoor, matDer2ndCoor[1], nullptr);
  TensorProductVolumeDataHex(TensorProductDataVolIntPoints, 3, nDOFs1D, nInt1D,
                             lagBasisLineInt, lagBasisLineInt, hesLagBasisLineInt,
                             matCoor, matDer2ndCoor[2], nullptr);

  TensorProductVolumeDataHex(TensorProductDataVolIntPoints, 3, nDOFs1D, nInt1D,
                             derLagBasisLineInt, derLagBasisLineInt, lagBasisLineInt,
                             matCoor, matDer2ndCoor[3], nullptr);
  TensorProductVolumeDataHex(TensorProductDataVolIntPoints, 3, nDOFs1D, nInt1D,
                             derLagBasisLineInt, lagBasisLineInt, derLagBasisLineInt,
                             matCoor, matDer2ndCoor[4], nullptr);
  TensorProductVolumeDataHex(TensorProductDataVolIntPoints, 3, nDOFs1D, nInt1D,
                             lagBasisLineInt, derLagBasisLineInt, derLagBasisLineInt,
                             matCoor, matDer2ndCoor[5], nullptr);
}

void CFEMStandardHexGrid::CoorSolDOFs(ColMajorMatrix<su2double> &matCoorDOF,
                                      ColMajorMatrix<su2double> &matCoorSolDOF) {

  /*--- Call the function TensorProductVolumeDataHex to compute
        the Cartesian coordinates in the solution DOFs. ---*/
  TensorProductVolumeDataHex(TensorProductDataVolSolDOFs, 3, nDOFs1D, rLineSolDOFs.size(),
                             lagBasisLineSolDOFs, lagBasisLineSolDOFs, lagBasisLineSolDOFs,
                             matCoorDOF, matCoorSolDOF, nullptr);
}

void CFEMStandardHexGrid::DerivativesCoorSolDOFs(ColMajorMatrix<su2double>          &matCoor,
                                                 vector<ColMajorMatrix<su2double> > &matDerCoor) {

  /*--- Call the function TensorProductVolumeDataHex 3 times to compute the derivatives of
        the Cartesian coordinates w.r.t. the three parametric coordinates. ---*/
  TensorProductVolumeDataHex(TensorProductDataVolSolDOFs, 3, nDOFs1D, rLineSolDOFs.size(),
                             derLagBasisLineSolDOFs, lagBasisLineSolDOFs, lagBasisLineSolDOFs,
                             matCoor, matDerCoor[0], nullptr);
  TensorProductVolumeDataHex(TensorProductDataVolSolDOFs, 3, nDOFs1D, rLineSolDOFs.size(),
                             lagBasisLineSolDOFs, derLagBasisLineSolDOFs, lagBasisLineSolDOFs,
                             matCoor, matDerCoor[1], nullptr);
  TensorProductVolumeDataHex(TensorProductDataVolSolDOFs, 3, nDOFs1D, rLineSolDOFs.size(),
                             lagBasisLineSolDOFs, lagBasisLineSolDOFs, derLagBasisLineSolDOFs,
                             matCoor, matDerCoor[2], nullptr);
}
