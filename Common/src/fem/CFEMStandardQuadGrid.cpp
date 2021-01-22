/*!
 * \file CFEMStandardQuadGrid.cpp
 * \brief Functions for the class CFEMStandardQuadGrid.
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

#include "../../include/fem/CFEMStandardQuadGrid.hpp"

/*----------------------------------------------------------------------------------*/
/*             Public member functions of CFEMStandardQuadGrid.                     */
/*----------------------------------------------------------------------------------*/

CFEMStandardQuadGrid::CFEMStandardQuadGrid(const unsigned short val_nPolyGrid,
                                           const unsigned short val_nPolySol,
                                           const unsigned short val_orderExact,
                                           const unsigned short val_locGridDOFs)
  : CFEMStandardQuadBase(val_nPolyGrid, val_orderExact) {

  /*--- Determine the 1D parametric locations of the grid DOFs. 1D is enough,
        because a tensor product is used to obtain the 2D coordinates. ---*/
  if(val_locGridDOFs == LGL) Location1DGridDOFsLGL(nPoly, rLineDOFs);
  else                       Location1DGridDOFsEquidistant(nPoly, rLineDOFs);

  /*--- Compute the 1D parametric coordinates of the solution DOFs. Only
        different from rLineDOFs when a different polynomial degree is
        used for the grid and solution. ---*/
  if(val_locGridDOFs == LGL) Location1DGridDOFsLGL(val_nPolySol, rLineSolDOFs);
  else                       Location1DGridDOFsEquidistant(val_nPolySol, rLineSolDOFs);

  /*--- Compute the 1D Lagrangian basis functions and its first
        and second derivatives in the integration points. ---*/
  LagBasisIntPointsLine(rLineDOFs, rLineInt, true, lagBasisLineInt);
  DerLagBasisIntPointsLine(rLineDOFs, rLineInt, true, derLagBasisLineInt);
  HesLagBasisIntPointsLine(rLineDOFs, rLineInt, true, hesLagBasisLineInt);

  /*--- Call LagBasisIntPointsLine and DerLagBasisIntPointsLine with the
        solution DOFs as argument to compute the Lagrangian basis functions
        and its derivatives in the solution DOFs. ---*/
  LagBasisIntPointsLine(rLineDOFs, rLineSolDOFs, true, lagBasisLineSolDOFs);
  DerLagBasisIntPointsLine(rLineDOFs, rLineSolDOFs, true, derLagBasisLineSolDOFs);

  /*--- Determine the local subconnectivity of this standard element when split
        in several linear elements. Used for a.o. plotting and searching. ---*/
  CFEMStandardQuadBase::SubConnLinearElements();

  /*--- Set the function pointers for the tensor product multiplications to
        determine the data in the volume integration points and volume
        nodal solution DOFs. ---*/
  SetFunctionPointerVolumeDataQuad(nDOFs1D, nInt1D, TensorProductDataVolIntPoints);
  SetFunctionPointerVolumeDataQuad(nDOFs1D, val_nPolySol+1, TensorProductDataVolSolDOFs);
}

void CFEMStandardQuadGrid::CoorIntPoints(const bool                notUsed,
                                         ColMajorMatrix<su2double> &matCoorDOF,
                                         ColMajorMatrix<su2double> &matCoorInt) {

  /*--- Call the function TensorProductVolumeDataQuad to compute the
        Cartesian coordinates in the integration points. ---*/
  TensorProductVolumeDataQuad(TensorProductDataVolIntPoints, 2, nDOFs1D, nInt1D, lagBasisLineInt,
                              lagBasisLineInt, matCoorDOF, matCoorInt, nullptr);
}

void CFEMStandardQuadGrid::DerivativesCoorIntPoints(const bool                         notUsed,
                                                    ColMajorMatrix<su2double>          &matCoor,
                                                    vector<ColMajorMatrix<su2double> > &matDerCoor) {

  /*--- Call the function TensorProductVolumeDataQuad 2 times to compute the
        derivatives of the Cartesian coordinates w.r.t. the two parametric coordinates. ---*/
  TensorProductVolumeDataQuad(TensorProductDataVolIntPoints, 2, nDOFs1D, nInt1D, derLagBasisLineInt,
                              lagBasisLineInt, matCoor, matDerCoor[0], nullptr);
  TensorProductVolumeDataQuad(TensorProductDataVolIntPoints, 2, nDOFs1D, nInt1D, lagBasisLineInt,
                              derLagBasisLineInt, matCoor, matDerCoor[1], nullptr);
}

void CFEMStandardQuadGrid::Derivatives2ndCoorIntPoints(ColMajorMatrix<su2double>          &matCoor,
                                                       vector<ColMajorMatrix<su2double> > &matDer2ndCoor) {

  /*--- Call the function TensorProductVolumeDataQuad 3 times to compute the 2nd
        derivatives of the Cartesian coordinates w.r.t. the two parametric coordinates. ---*/
  TensorProductVolumeDataQuad(TensorProductDataVolIntPoints, 2, nDOFs1D, nInt1D, hesLagBasisLineInt,
                              lagBasisLineInt, matCoor, matDer2ndCoor[0], nullptr);
  TensorProductVolumeDataQuad(TensorProductDataVolIntPoints, 2, nDOFs1D, nInt1D, lagBasisLineInt,
                              hesLagBasisLineInt, matCoor, matDer2ndCoor[1], nullptr);
  TensorProductVolumeDataQuad(TensorProductDataVolIntPoints, 2, nDOFs1D, nInt1D, derLagBasisLineInt,
                              derLagBasisLineInt, matCoor, matDer2ndCoor[2], nullptr);
}

void CFEMStandardQuadGrid::CoorSolDOFs(ColMajorMatrix<su2double> &matCoorDOF,
                                       ColMajorMatrix<su2double> &matCoorSolDOF) {

  /*--- Call the function TensorProductVolumeDataQuad to compute
        the Cartesian coordinates in the solution DOFs. ---*/
  TensorProductVolumeDataQuad(TensorProductDataVolSolDOFs, 2, nDOFs1D, rLineSolDOFs.size(),
                              lagBasisLineSolDOFs, lagBasisLineSolDOFs, matCoorDOF, matCoorSolDOF, nullptr);
}

void CFEMStandardQuadGrid::DerivativesCoorSolDOFs(ColMajorMatrix<su2double>          &matCoor,
                                                  vector<ColMajorMatrix<su2double> > &matDerCoor) {

  /*--- Call the function TensorProductVolumeDataQuad 2 times to compute the derivatives
        of the Cartesian coordinates w.r.t. the two parametric coordinates. ---*/
  TensorProductVolumeDataQuad(TensorProductDataVolSolDOFs, 2, nDOFs1D, rLineSolDOFs.size(),
		              derLagBasisLineSolDOFs, lagBasisLineSolDOFs, matCoor, matDerCoor[0], nullptr);
  TensorProductVolumeDataQuad(TensorProductDataVolSolDOFs, 2, nDOFs1D, rLineSolDOFs.size(),
		              lagBasisLineSolDOFs, derLagBasisLineSolDOFs, matCoor, matDerCoor[1], nullptr);
}

void CFEMStandardQuadGrid::EvalCoorAndGradCoor(ColMajorMatrix<su2double> &matCoor,
                                               const su2double           *par,
                                               su2double                 x[3],
                                               su2double                 dxdpar[3][3]) {

  /*--- Convert the parametric coordinates to a passivedouble and
        store it in a vector. ---*/
  vector<passivedouble> rPar(2);
  rPar[0] = SU2_TYPE::GetValue(par[0]);
  rPar[1] = SU2_TYPE::GetValue(par[1]);

  /*--- Compute the 1D Lagrangian basis functions and its first
        derivatives in these parametric points. ---*/
  ColMajorMatrix<passivedouble> lag, derLag;

  LagBasisIntPointsLine(rLineDOFs, rPar, false, lag);
  DerLagBasisIntPointsLine(rLineDOFs, rPar, false, derLag);

  /*--- Initialize the coordinates and its derivatives. ---*/
  for(unsigned short j=0; j<2; ++j) {
    x[j] = 0.0;
    for(unsigned short i=0; i<2; ++i)
      dxdpar[j][i] = 0.0;
  }

  /*--- Double loop to compute the coordinates and its derivatives. ---*/
  unsigned short ii = 0;
  for(unsigned short j=0; j<nDOFs1D; ++j) {
    for(unsigned short i=0; i<nDOFs1D; ++i, ++ii) {

      const passivedouble lag3D    = lag(0,i)*lag(1,j);
      const passivedouble dLag3Ddr = derLag(0,i)*lag(1,j);
      const passivedouble dLag3Dds = lag(0,i)*derLag(1,j);

      for(unsigned short l=0; l<2; ++l) {
        x[l]         += matCoor(ii,l)*lag3D;
        dxdpar[l][0] += matCoor(ii,l)*dLag3Ddr;
        dxdpar[l][1] += matCoor(ii,l)*dLag3Dds;
      }
    }
  }
}

void CFEMStandardQuadGrid::InterpolCoorSubElem(const unsigned short subElem,
                                               const su2double      *weights,
                                               su2double            *parCoor) {

  /*--- Easier storage of the connectivity of the sub-element. ---*/
  const unsigned short *conn = subConn1ForPlotting.data() + 4*subElem;

  /*--- Determine the nodal offset in j-direction. ---*/
  const unsigned short jOff = nPoly+1;

  /*--- Initialize the parametric coordinates to zero. ---*/
  parCoor[0] = parCoor[1] = 0.0;

  /*--- Loop over the 4 nodes of the linear sub-element, which is a quad. ---*/
  for(unsigned short ind=0; ind<4; ++ind) {

    /*--- Determine the (i,j) indices of this vertex inside the parent element. ---*/
    const unsigned short j = conn[ind]/jOff;
    const unsigned short i = conn[ind] - j*jOff;

    /*--- Update the parametric coordinates. ---*/
    parCoor[0] = weights[ind]*rLineDOFs[i];
    parCoor[1] = weights[ind]*rLineDOFs[j];
  }
}
