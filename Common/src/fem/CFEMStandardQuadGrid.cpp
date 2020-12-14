/*!
 * \file CFEMStandardQuadGrid.cpp
 * \brief Functions for the class CFEMStandardQuadGrid.
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

#include "../../include/fem/CFEMStandardQuadGrid.hpp"

/*----------------------------------------------------------------------------------*/
/*             Public member functions of CFEMStandardQuadGrid.                     */
/*----------------------------------------------------------------------------------*/

CFEMStandardQuadGrid::CFEMStandardQuadGrid(const unsigned short val_nPoly,
                                           const unsigned short val_orderExact,
                                           const bool           val_surfElement)
  : CFEMStandardQuad(val_nPoly, val_orderExact) {

  /*--- Determine the number of space dimensions, which is 3 if this standard
        element is a surface element and 2 when it is a volume element. ---*/
  nDim = val_surfElement ? 3 : 2;

  /*--- Compute the values of the 1D Lagrangian basis functions in the integration
        points for both the equidistant and LGL point distribution. ---*/
  LagBasisIntPointsLine(rLineDOFsEqui, rLineInt, lagBasisLineIntEqui);
  LagBasisIntPointsLine(rLineDOFsLGL,  rLineInt, lagBasisLineIntLGL);

  /*--- Compute the values of the derivatives of the 1D Lagrangian basis functions in
        the integration points for both the equidistant and LGL point distribution. ---*/
  DerLagBasisIntPointsLine(rLineDOFsEqui, rLineInt, derLagBasisLineIntEqui);
  DerLagBasisIntPointsLine(rLineDOFsLGL,  rLineInt, derLagBasisLineIntLGL);

  /*--- Create the local grid connectivities of the faces of the volume element.
        Only needed if this is a volume element. ---*/
  if( !val_surfElement ) LocalGridConnFaces();

  /*--- Determine the local subconnectivity of this standard element when split
        in several linear elements. Used for a.o. plotting and searcing. ---*/
  SubConnLinearElements();
}

CFEMStandardQuadGrid::CFEMStandardQuadGrid(const unsigned short val_nPolyGrid,
                                           const unsigned short val_nPolySol,
                                           const unsigned short val_orderExact,
                                           const unsigned short val_locGridDOFs)
  : CFEMStandardQuad(val_nPolyGrid, val_orderExact) {

  /*--- Set the number of space dimensions, which is always 2 when this
        constructor is called for a quadrilateral standard element. ---*/
  nDim = 2;

  /*--- Determine the location of the grid DOFs and build the appropriate
        Lagrangian basis functions. ---*/
  if(val_locGridDOFs == LGL) {

    /*--- LGL distribution. Compute the 1D Lagrangian basis functions and
          its first and second derivatives in the integration points. ---*/
    LagBasisIntPointsLine(rLineDOFsLGL, rLineInt, lagBasisLineIntLGL);
    DerLagBasisIntPointsLine(rLineDOFsLGL, rLineInt, derLagBasisLineIntLGL);
    HesLagBasisIntPointsLine(rLineDOFsLGL, rLineInt, hesLagBasisLineInt);

    /*--- Compute the 1D parametric coordinates of the solution DOFs. Only
          different from rLineDOFsLGL when a different polynomial degree is
          used for the grid and solution. ---*/
    nPoly = val_nPolySol;
    Location1DGridDOFsLGL(rLineSolDOFs);
    nPoly = val_nPolyGrid;

    /*--- Call LagBasisIntPointsLine and DerLagBasisIntPointsLine with the
          solution DOFs as argument to compute the Lagrangian basis functions
          and its derivatives in the solution DOFs. ---*/
    LagBasisIntPointsLine(rLineDOFsLGL, rLineSolDOFs, lagBasisLineSolDOFs);
    DerLagBasisIntPointsLine(rLineDOFsLGL, rLineSolDOFs, derLagBasisLineSolDOFs);
  }
  else {

    /*--- Equidistant distribution. Compute the 1D Lagrangian basis functions and
          its first and second derivatives in the integration points. ---*/
    LagBasisIntPointsLine(rLineDOFsEqui, rLineInt, lagBasisLineIntEqui);
    DerLagBasisIntPointsLine(rLineDOFsEqui, rLineInt, derLagBasisLineIntEqui);
    HesLagBasisIntPointsLine(rLineDOFsEqui, rLineInt, hesLagBasisLineInt);

    /*--- Compute the 1D parametric coordinates of the solution DOFs. Only
          different from rLineDOFsEqui when a different polynomial degree is
          used for the grid and solution. ---*/
    nPoly = val_nPolySol;
    Location1DGridDOFsEquidistant(rLineSolDOFs);
    nPoly = val_nPolyGrid;

    /*--- Call LagBasisIntPointsLine and DerLagBasisIntPointsLine with the
          solution DOFs as argument to compute the Lagrangian basis functions
          and its derivatives in the solution DOFs. ---*/
    LagBasisIntPointsLine(rLineDOFsEqui, rLineSolDOFs, lagBasisLineSolDOFs);
    DerLagBasisIntPointsLine(rLineDOFsEqui, rLineSolDOFs, derLagBasisLineSolDOFs);
  }

  /*--- Create the local grid connectivities of the faces of the volume element. ---*/
  LocalGridConnFaces();

  /*--- Determine the local subconnectivity of this standard element when split
        in several linear elements. Used for a.o. plotting and searcing. ---*/
  SubConnLinearElements();

  /*--- Create the map with the function pointers to carry out the tensor product
        to compute the data in the 2D nodal solution DOFs of the quadrilateral. ---*/
  map<CUnsignedShort2T, TPI2D> mapFunctions;
  CreateMapTensorProductVolumeIntPoints2D(mapFunctions);

  /*--- Try to find the combination of the number of 1D DOFs and solution DOFs
        in mapFunctions. If not found, write a clear error message that this
        tensor product is not supported. ---*/
  CUnsignedShort2T nDOFsAndSolDOFs(nDOFs1D, val_nPolySol+1);
  auto MI = mapFunctions.find(nDOFsAndSolDOFs);
  if(MI == mapFunctions.end()) {
    std::ostringstream message;
    message << "The tensor product TensorProductVolumeIntPoints2D_" << nDOFs1D
            << "_" << val_nPolySol+1 << " not created by the automatic source code "
            << "generator. Modify this automatic source code creator";
    SU2_MPI::Error(message.str(), CURRENT_FUNCTION);
  }

  /*--- Set the function pointer to carry out tensor product. ---*/
  TensorProductDataVolSolDOFs = MI->second;
}

void CFEMStandardQuadGrid::CoorIntPoints(const bool                LGLDistribution,
                                         ColMajorMatrix<su2double> &matCoorDOF,
                                         ColMajorMatrix<su2double> &matCoorInt) {

  /*--- Check for which point distribution the derivatives must be computed. ---*/
  if( LGLDistribution ) {

    /*--- LGL distribution. Call the function TensorProductIntegrationPoints to compute the
          Cartesian coordinates in the integration points. The first argument in the function
          call is nDim, which corresponds to the number of Cartesian coordinates (3 for a
          surface element and 2 for a volume element). ---*/
    TensorProductIntegrationPoints(nDim, lagBasisLineIntLGL, lagBasisLineIntLGL,
                                   matCoorDOF, matCoorInt, nullptr);
  }
  else {

    /*--- Equidistant distribution. Call the function TensorProductIntegrationPoints to compute the
          Cartesian coordinates in the integration points. The first argument in the function
          call is nDim, which corresponds to the number of Cartesian coordinates (3 for a
          surface element and 2 for a volume element). ---*/
    TensorProductIntegrationPoints(nDim, lagBasisLineIntEqui, lagBasisLineIntEqui,
                                   matCoorDOF, matCoorInt, nullptr);
  }
}

void CFEMStandardQuadGrid::DerivativesCoorIntPoints(const bool                         LGLDistribution,
                                                    ColMajorMatrix<su2double>          &matCoor,
                                                    vector<ColMajorMatrix<su2double> > &matDerCoor) {

  /*--- Check for which point distribution the derivatives must be computed. ---*/
  if( LGLDistribution ) {

    /*--- LGL distribution. Call the function TensorProductIntegrationPoints 2 times to compute the
          derivatives of the Cartesian coordinates w.r.t. the two parametric coordinates. The first
          argument in the function call is nDim, which corresponds to the number of Cartesian
          coordinates (3 for a surface element and 2 for a volume element). ---*/
    TensorProductIntegrationPoints(nDim, derLagBasisLineIntLGL, lagBasisLineIntLGL,
                                   matCoor, matDerCoor[0], nullptr);
    TensorProductIntegrationPoints(nDim, lagBasisLineIntLGL, derLagBasisLineIntLGL,
                                   matCoor, matDerCoor[1], nullptr);
  }
  else {

    /*--- Equidistant distribution. Call the function TensorProductIntegrationPoints 2 times to compute the
          derivatives of the Cartesian coordinates w.r.t. the two parametric coordinates. The first
          argument in the function call is nDim, which corresponds to the number of Cartesian
          coordinates (3 for a surface element and 2 for a volume element). ---*/
    TensorProductIntegrationPoints(nDim, derLagBasisLineIntEqui, lagBasisLineIntEqui,
                                   matCoor, matDerCoor[0], nullptr);
    TensorProductIntegrationPoints(nDim, lagBasisLineIntEqui, derLagBasisLineIntEqui,
                                   matCoor, matDerCoor[1], nullptr);
  }
}

void CFEMStandardQuadGrid::DerivativesCoorSolDOFs(ColMajorMatrix<su2double>          &matCoor,
                                                  vector<ColMajorMatrix<su2double> > &matDerCoor) {

  /*--- Call the function TensorProductSolDOFs 2 times to compute the derivatives of
        the Cartesian coordinates w.r.t. the two parametric coordinates. The first
        argument in the function call is nDim, which corresponds to the number of Cartesian
        coordinates (3 for a surface element and 2 for a volume element). ---*/
  TensorProductSolDOFs(nDim, derLagBasisLineSolDOFs, lagBasisLineSolDOFs,
                       matCoor, matDerCoor[0], nullptr);
  TensorProductSolDOFs(nDim, lagBasisLineSolDOFs, derLagBasisLineSolDOFs,
                       matCoor, matDerCoor[1], nullptr);
}

passivedouble CFEMStandardQuadGrid::WorkEstimateBoundaryFace(CConfig              *config,
                                                             const unsigned short elemType) {

  /*--- Determine the number of DOFs of the neighboring element. ---*/
  const unsigned short nDOFsElem = GetNDOFsStatic(elemType, nPoly);

  /*--- TEMPORARY IMPLEMENTATION. ---*/
  return nIntegration + 0.05*nDOFsElem;
}

passivedouble CFEMStandardQuadGrid::WorkEstimateInternalFace(CConfig              *config,
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

passivedouble CFEMStandardQuadGrid::WorkEstimateVolume(CConfig *config) {

  /*--- TEMPORARY IMPLEMENTATION. ---*/
  return nIntegration + 0.1*nDOFs;
}

passivedouble CFEMStandardQuadGrid::WorkEstimateWallFunctions(CConfig              *config,
                                                              const unsigned short nPointsWF,
                                                              const unsigned short elemType) {

  /*--- TEMPORARY IMPLEMENTATION. ---*/
  return 0.25*nIntegration*nPointsWF;
}

/*----------------------------------------------------------------------------------*/
/*             Private member functions of CFEMStandardQuadGrid.                    */
/*----------------------------------------------------------------------------------*/

void CFEMStandardQuadGrid::LocalGridConnFaces(void) {

  /*--- Allocate the first index of gridConnFaces, which is equal to the number
        of faces of the quadrilateral, which is 4. Reserve memory for the second
        index afterwards. ---*/
  gridConnFaces.resize(4);

  gridConnFaces[0].reserve(nPoly+1);
  gridConnFaces[1].reserve(nPoly+1);
  gridConnFaces[2].reserve(nPoly+1);
  gridConnFaces[3].reserve(nPoly+1);

  /*--- Define the corner vertices of the quadrilateral. ---*/
  const unsigned short n0 = 0, n1 = nPoly, n2 = nDOFs-1, n3 = nPoly*(nPoly+1);

  /*--- For a quad element the faces are lines. Loop over the nodes of the
        lines to set the connectivity. Make sure that the element
        is to the left of the faces. ---*/
  for(signed short i=n0; i<=n1; ++i)          gridConnFaces[0].push_back(i);
  for(signed short i=n1; i<=n2; i+=(nPoly+1)) gridConnFaces[1].push_back(i);
  for(signed short i=n2; i>=n3; --i)          gridConnFaces[2].push_back(i);
  for(signed short i=n3; i>=n0; i-=(nPoly+1)) gridConnFaces[3].push_back(i);
}

void CFEMStandardQuadGrid::SubConnLinearElements(void) {

  /*--- The quadrilateral is split into several linear quads.
        Set the VTK sub-types accordingly. ---*/
  VTK_SubType1 = QUADRILATERAL;
  VTK_SubType2 = NONE;

  /*--- Determine the local subconnectivity of the quadrilateral element used for
        plotting purposes. Note that the connectivity of the linear subelements
        obey the VTK connectivity rule of a quadrilateral, which is different
        from the connectivity for the high order quadrilateral. ---*/
  unsigned short nnPoly = max(nPoly,(unsigned short) 1);
  for(unsigned short j=0; j<nnPoly; ++j) {
    unsigned short jj = j*(nnPoly+1);
    for(unsigned short i=0; i<nnPoly; ++i) {
      const unsigned short n0 = jj + i;
      const unsigned short n1 = n0 + 1;
      const unsigned short n2 = n1 + nPoly+1;
      const unsigned short n3 = n2 - 1;

      subConn1ForPlotting.push_back(n0);
      subConn1ForPlotting.push_back(n1);
      subConn1ForPlotting.push_back(n2);
      subConn1ForPlotting.push_back(n3);
    }
  }
}

void CFEMStandardQuadGrid::TensorProductSolDOFs(const int                           N,
                                                const ColMajorMatrix<passivedouble> &Ai,
                                                const ColMajorMatrix<passivedouble> &Aj,
                                                const ColMajorMatrix<su2double>     &B,
                                                ColMajorMatrix<su2double>           &C,
                                                const CConfig                       *config) {

  /*--- Call the function to which TensorProductDataVolSolDOFs points to carry out
        the actual tensor product. Perform the timing, if desired. ---*/
#ifdef PROFILE
  double timeGemm;
  if( config ) config->TensorProduct_Tick(&timeGemm);
#endif

  TensorProductDataVolSolDOFs(N, B.rows(), C.rows(), Ai.data(), Aj.data(), B.data(), C.data());

#ifdef PROFILE
  if( config ) config->TensorProduct_Tock(timeGemm, 2, N, nDOFs1D, rLineSolDOFs.size());
#endif
}
