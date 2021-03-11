/*!
 * \file CFEMStandardHexGrid.cpp
 * \brief Functions for the class CFEMStandardHexGrid.
 * \author E. van der Weide
 * \version 7.1.1 "Blackbird"
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

void CFEMStandardHexGrid::EvalCoorAndGradCoor(ColMajorMatrix<su2double> &matCoor,
                                              const su2double           *par,
                                              su2double                 x[3],
                                              su2double                 dxdpar[3][3]) {

  /*--- Convert the parametric coordinates to a passivedouble and
        store it in a vector. ---*/
  vector<passivedouble> rPar(3);
  rPar[0] = SU2_TYPE::GetValue(par[0]);
  rPar[1] = SU2_TYPE::GetValue(par[1]);
  rPar[2] = SU2_TYPE::GetValue(par[2]);

  /*--- Compute the 1D Lagrangian basis functions and its first
        derivatives in these parametric points. ---*/
  ColMajorMatrix<passivedouble> lag, derLag;

  LagBasisIntPointsLine(rLineDOFs, rPar, false, lag);
  DerLagBasisIntPointsLine(rLineDOFs, rPar, false, derLag);

  /*--- Initialize the coordinates and its derivatives. ---*/
  for(unsigned short j=0; j<3; ++j) {
    x[j] = 0.0;
    for(unsigned short i=0; i<3; ++i)
      dxdpar[j][i] = 0.0;
  }

  /*--- Triple loop to compute the coordinates and its derivatives. ---*/
  unsigned short ii = 0;
  for(unsigned short k=0; k<nDOFs1D; ++k) {
    for(unsigned short j=0; j<nDOFs1D; ++j) {
      for(unsigned short i=0; i<nDOFs1D; ++i, ++ii) {

        const passivedouble lag3D    = lag(0,i)*lag(1,j)*lag(2,k);
        const passivedouble dLag3Ddr = derLag(0,i)*lag(1,j)*lag(2,k);
        const passivedouble dLag3Dds = lag(0,i)*derLag(1,j)*lag(2,k);
        const passivedouble dLag3Ddt = lag(0,i)*lag(1,j)*derLag(2,k);

        for(unsigned short l=0; l<3; ++l) {
          x[l]         += matCoor(ii,l)*lag3D;
          dxdpar[l][0] += matCoor(ii,l)*dLag3Ddr;
          dxdpar[l][1] += matCoor(ii,l)*dLag3Dds;
          dxdpar[l][2] += matCoor(ii,l)*dLag3Ddt;
        }
      }
    }
  }
}

void CFEMStandardHexGrid::InterpolCoorSubElem(const unsigned short subElem,
                                              const su2double      *weights,
                                              su2double            *parCoor) {

  /*--- Easier storage of the connectivity of the sub-element. ---*/
  const unsigned short *conn = subConn1ForPlotting.data() + 8*subElem;

  /*--- Determine the nodal offset in j- and k-direction. ---*/
  const unsigned short jOff = nPoly+1;
  const unsigned short kOff = jOff*jOff;

  /*--- Initialize the parametric coordinates to zero. ---*/
  parCoor[0] = parCoor[1] = parCoor[2] = 0.0;

  /*--- Loop over the 8 nodes of the linear sub-element, which is a hex. ---*/
  for(unsigned short ind=0; ind<8; ++ind) {

    /*--- Determine the (i,j,k) indices of this vertex inside the parent element. ---*/
    const unsigned short k   = conn[ind]/kOff;
    const unsigned short tmp = conn[ind] - k*kOff;
    const unsigned short j   = tmp/jOff;
    const unsigned short i   = tmp - j*jOff;

    /*--- Update the parametric coordinates. ---*/
    parCoor[0] = weights[ind]*rLineDOFs[i];
    parCoor[1] = weights[ind]*rLineDOFs[j];
    parCoor[2] = weights[ind]*rLineDOFs[k];
  }
}
