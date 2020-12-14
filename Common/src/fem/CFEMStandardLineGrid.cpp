/*!
 * \file CFEMStandardLineGrid.cpp
 * \brief Functions for the class CFEMStandardLineGrid.
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

#include "../../include/fem/CFEMStandardLineGrid.hpp"

/*----------------------------------------------------------------------------------*/
/*            Public member functions of CFEMStandardLineGrid.                      */
/*----------------------------------------------------------------------------------*/

CFEMStandardLineGrid::CFEMStandardLineGrid(const unsigned short val_nPoly,
                                           const unsigned short val_orderExact)
  : CFEMStandardLine(val_nPoly, val_orderExact) {

  /*--- Compute the values of the 1D Lagrangian basis functions in the integration
        points for both the equidistant and LGL point distribution. ---*/
  LagBasisIntPointsLine(rLineDOFsEqui, rLineInt, lagBasisLineIntEqui);
  LagBasisIntPointsLine(rLineDOFsLGL,  rLineInt, lagBasisLineIntLGL);

  /*--- Compute the values of the derivatives of the 1D Lagrangian basis functions in
        the integration points for both the equidistant and LGL point distribution. ---*/
  DerLagBasisIntPointsLine(rLineDOFsEqui, rLineInt, derLagBasisLineIntEqui);
  DerLagBasisIntPointsLine(rLineDOFsLGL,  rLineInt, derLagBasisLineIntLGL);

  /*--- Determine the local subconnectivity of this standard element when split
        in several linear elements. Used for a.o. plotting and searcing. ---*/
  SubConnLinearElements();

  /*--- Set up the jitted gemm call, if supported. For this particular standard
        element the derivatives of the coordinates are computed, which is 2,
        because the line element is only present as a surface element in a
        2D simulation. ---*/
  SetUpJittedGEMM(nIntegrationPad, 2, nDOFs, jitter, my_dgemm);
}

void CFEMStandardLineGrid::CoorIntPoints(const bool                LGLDistribution,
                                         ColMajorMatrix<su2double> &matCoorDOF,
                                         ColMajorMatrix<su2double> &matCoorInt) {

  /*--- Check for which point distribution the derivatives must be computed. ---*/
  if( LGLDistribution ) {

    /*--- LGL distribution. Call the function OwnGemm to compute the
          Cartesian coordinates in the integration points. ---*/
    OwnGemm(my_dgemm, nIntegrationPad, 2, nDOFs, lagBasisLineIntLGL,
            matCoorDOF, matCoorInt, nullptr);
  }
  else {

    /*--- Equidistant distribution. Call the function OwnGemm to compute the
          Cartesian coordinates in the integration points. ---*/
    OwnGemm(my_dgemm, nIntegrationPad, 2, nDOFs, lagBasisLineIntEqui,
            matCoorDOF, matCoorInt, nullptr);
  }
}

void CFEMStandardLineGrid::DerivativesCoorIntPoints(const bool                         LGLDistribution,
                                                    ColMajorMatrix<su2double>          &matCoor,
                                                    vector<ColMajorMatrix<su2double> > &matDerCoor) {

  /*--- Check for which point distribution the derivatives must be computed. ---*/
  if( LGLDistribution ) {

    /*--- LGL distribution. Call the function OwnGemm to compute the derivatives
          of the Cartesian coordinates w.r.t. the parametric coordinate. ---*/
    OwnGemm(my_dgemm, nIntegrationPad, 2, nDOFs, derLagBasisLineIntLGL,
            matCoor, matDerCoor[0], nullptr);
  }
  else {

    /*--- Equidistant distribution. Call the function OwnGemm to compute the derivatives
          of the Cartesian coordinates w.r.t. the parametric coordinate. ---*/
    OwnGemm(my_dgemm, nIntegrationPad, 2, nDOFs, derLagBasisLineIntEqui,
            matCoor, matDerCoor[0], nullptr);
  }
}

passivedouble CFEMStandardLineGrid::WorkEstimateBoundaryFace(CConfig              *config,
                                                             const unsigned short elemType) {

  /*--- Determine the number of DOFs of the neighboring element. ---*/
  const unsigned short nDOFsElem = GetNDOFsStatic(elemType, nPoly);

  /*--- TEMPORARY IMPLEMENTATION. ---*/
  return nIntegration + 0.05*nDOFsElem;
}

passivedouble CFEMStandardLineGrid::WorkEstimateInternalFace(CConfig              *config,
                                                             const unsigned short elemType0,
                                                             const unsigned short nPoly0,
                                                             const unsigned short elemType1,
                                                             const unsigned short nPoly1) {

  /*--- Determine the number of DOFs of the neighboring elements. ---*/
  const unsigned short nDOFsElem0 = GetNDOFsStatic(elemType0, nPoly0);
  const unsigned short nDOFsElem1 = GetNDOFsStatic(elemType1, nPoly1);

  /* TEMPORARY IMPLEMENTATION. */
  return 2.0*nIntegration + 0.05*(nDOFsElem0 + nDOFsElem1);
}

passivedouble CFEMStandardLineGrid::WorkEstimateWallFunctions(CConfig              *config,
                                                              const unsigned short nPointsWF,
                                                              const unsigned short elemType) {

  /*--- TEMPORARY IMPLEMENTATION. ---*/
  return 0.25*nIntegration*nPointsWF;
}

/*----------------------------------------------------------------------------------*/
/*            Private member functions of CFEMStandardLineGrid.                     */
/*----------------------------------------------------------------------------------*/

void CFEMStandardLineGrid::SubConnLinearElements(void) {

  /*--- The line is split into several linear lines.
        Set the VTK sub-types accordingly. ---*/
  VTK_SubType1 = LINE;
  VTK_SubType2 = NONE;

  /*--- Determine the local subconnectivity of the line element used for plotting
        purposes. This is rather trivial, because the line element is subdivided
        into nPoly linear line elements. ---*/
  unsigned short nnPoly = max(nPoly,(unsigned short) 1);
  for(unsigned short i=0; i<nnPoly; ++i) {
    subConn1ForPlotting.push_back(i);
    subConn1ForPlotting.push_back(i+1);
  }
}
