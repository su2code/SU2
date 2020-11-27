/*!
 * \file CFEMStandardHexGrid.cpp
 * \brief Functions for the class CFEMStandardHexGrid.
 * \author E. van der Weide
 * \version 7.0.7 "Blackbird"
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

CFEMStandardHexGrid::CFEMStandardHexGrid(const unsigned short val_nPoly,
                                         const unsigned short val_orderExact)
  : CFEMStandardHex(val_nPoly,val_orderExact) {

  /*--- Compute the values of the 1D Lagrangian basis functions in the integration
        points for both the equidistant and LGL point distribution. ---*/
  LagBasisIntPointsLine(rLineDOFsEqui, rLineInt, lagBasisLineIntEqui);
  LagBasisIntPointsLine(rLineDOFsLGL,  rLineInt, lagBasisLineIntLGL);

  /*--- Compute the values of the derivatives of the 1D Lagrangian basis functions in
        the integration points for both the equidistant and LGL point distribution. ---*/
  DerLagBasisIntPointsLine(rLineDOFsEqui, rLineInt, derLagBasisLineIntEqui);
  DerLagBasisIntPointsLine(rLineDOFsLGL,  rLineInt, derLagBasisLineIntLGL);

  /*--- Create the local grid connectivities of the faces of the volume element. ---*/
  LocalGridConnFaces();

  /*--- Determine the local subconnectivity of this standard element when split
        in several linear elements. Used for a.o. plotting and searcing. ---*/
  SubConnLinearElements();
}

CFEMStandardHexGrid::CFEMStandardHexGrid(const unsigned short val_nPolyGrid,
                                         const unsigned short val_nPolySol,
                                         const unsigned short val_orderExact,
                                         const unsigned short val_locGridDOFs)

  : CFEMStandardHex(val_nPolyGrid,val_orderExact) {

  /*--- Determine the location of the grid DOFs and build the appropriate
        Lagrangian basis functions. ---*/
  if(val_locGridDOFs == LGL) {

    /*--- LGL distribution. Compute the 1D Lagrangian basis functions and
          its derivatives in the integration points. ---*/
    LagBasisIntPointsLine(rLineDOFsLGL, rLineInt, lagBasisLineIntLGL);
    DerLagBasisIntPointsLine(rLineDOFsLGL, rLineInt, derLagBasisLineIntLGL);

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
          its derivatives in the integration points. ---*/
    LagBasisIntPointsLine(rLineDOFsEqui, rLineInt, lagBasisLineIntEqui);
    DerLagBasisIntPointsLine(rLineDOFsEqui, rLineInt, derLagBasisLineIntEqui);

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
        to compute the data in the 3D nodal solution DOFs of the hexahedron. ---*/
  map<CUnsignedShort2T, TPI3D> mapFunctions;
  CreateMapTensorProductVolumeIntPoints3D(mapFunctions);

  /*--- Try to find the combination of the number of 1D DOFs and solution DOFs
        in mapFunctions. If not found, write a clear error message that this
        tensor product is not supported. ---*/
  CUnsignedShort2T nDOFsAndSolDOFs(nDOFs1D, val_nPolySol+1);
  auto MI = mapFunctions.find(nDOFsAndSolDOFs);
  if(MI == mapFunctions.end()) {
    std::ostringstream message;
    message << "The tensor product TensorProductVolumeIntPoints3D_" << nDOFs1D
            << "_" << val_nPolySol+1 << " not created by the automatic source code "
            << "generator. Modify this automatic source code creator";
    SU2_MPI::Error(message.str(), CURRENT_FUNCTION);
  }

  /*--- Set the function pointer to carry out tensor product. ---*/
  TensorProductDataVolSolDOFs = MI->second;
}

void CFEMStandardHexGrid::CoorIntPoints(const bool                LGLDistribution,
                                        ColMajorMatrix<su2double> &matCoorDOF,
                                        ColMajorMatrix<su2double> &matCoorInt) {

  /*--- Check for which point distribution the coordinates must be computed. ---*/
  if( LGLDistribution ) {

    /*--- LGL distribution. Call the function TensorProductIntegrationPoints to compute the
          Cartesian coordinates in the integration points. ---*/
    TensorProductIntegrationPoints(3, lagBasisLineIntLGL, lagBasisLineIntLGL, lagBasisLineIntLGL,
                                   matCoorDOF, matCoorInt, nullptr);
  }
  else {

    /*--- Equidistant distribution. Call the function TensorProductIntegrationPoints to compute the
          Cartesian coordinates in the integration points. ---*/
    TensorProductIntegrationPoints(3, lagBasisLineIntEqui, lagBasisLineIntEqui, lagBasisLineIntEqui,
                                   matCoorDOF, matCoorInt, nullptr);
  }
}

void CFEMStandardHexGrid::DerivativesCoorIntPoints(const bool                         LGLDistribution,
                                                   ColMajorMatrix<su2double>          &matCoor,
                                                   vector<ColMajorMatrix<su2double> > &matDerCoor) {

  /*--- Check for which point distribution the derivatives must be computed. ---*/
  if( LGLDistribution ) {

    /*--- LGL distribution. Call the function TensorProductIntegrationPoints 3 times to compute the
          derivatives of the Cartesian coordinates w.r.t. the three parametric coordinates. ---*/
    TensorProductIntegrationPoints(3, derLagBasisLineIntLGL, lagBasisLineIntLGL, lagBasisLineIntLGL,
                                   matCoor, matDerCoor[0], nullptr);
    TensorProductIntegrationPoints(3, lagBasisLineIntLGL, derLagBasisLineIntLGL, lagBasisLineIntLGL,
                                   matCoor, matDerCoor[1], nullptr);
    TensorProductIntegrationPoints(3, lagBasisLineIntLGL, lagBasisLineIntLGL, derLagBasisLineIntLGL,
                                   matCoor, matDerCoor[2], nullptr);
  }
  else {

    /*--- Equidistant distribution. Call the function TensorProductIntegrationPoints 3 times to compute the
          derivatives of the Cartesian coordinates w.r.t. the three parametric coordinates. ---*/
    TensorProductIntegrationPoints(3, derLagBasisLineIntEqui, lagBasisLineIntEqui, lagBasisLineIntEqui,
                                   matCoor, matDerCoor[0], nullptr);
    TensorProductIntegrationPoints(3, lagBasisLineIntEqui, derLagBasisLineIntEqui, lagBasisLineIntEqui,
                                   matCoor, matDerCoor[1], nullptr);
    TensorProductIntegrationPoints(3, lagBasisLineIntEqui, lagBasisLineIntEqui, derLagBasisLineIntEqui,
                                   matCoor, matDerCoor[2], nullptr);
  }
}

void CFEMStandardHexGrid::DerivativesCoorSolDOFs(ColMajorMatrix<su2double>          &matCoor,
                                                 vector<ColMajorMatrix<su2double> > &matDerCoor) {

  /*--- Call the function TensorProductSolDOFs 3 times to compute the derivatives of
        the Cartesian coordinates w.r.t. the three parametric coordinates. ---*/
  TensorProductSolDOFs(3, derLagBasisLineSolDOFs, lagBasisLineSolDOFs, lagBasisLineSolDOFs,
                       matCoor, matDerCoor[0], nullptr);
  TensorProductSolDOFs(3, lagBasisLineSolDOFs, derLagBasisLineSolDOFs, lagBasisLineSolDOFs,
                       matCoor, matDerCoor[1], nullptr);
  TensorProductSolDOFs(3, lagBasisLineSolDOFs, lagBasisLineSolDOFs, derLagBasisLineSolDOFs,
                       matCoor, matDerCoor[2], nullptr);
}

passivedouble CFEMStandardHexGrid::WorkEstimateVolume(CConfig *config) {

  /*--- TEMPORARY IMPLEMENTATION. ---*/
  return nIntegration + 0.1*nDOFs;
}

/*----------------------------------------------------------------------------------*/
/*             Private member functions of CFEMStandardHexGrid.                     */
/*----------------------------------------------------------------------------------*/

void CFEMStandardHexGrid::LocalGridConnFaces(void) {

  /*--- Allocate the first index of gridConnFaces, which is equal to the number
        of faces of the hexahedron, which is 6. Reserve memory for the second
        index afterwards. ---*/
  const unsigned short nDOFsQuad = (nPoly+1)*(nPoly+1);
  gridConnFaces.resize(6);

  for(unsigned short i=0; i<6; ++i) gridConnFaces[i].reserve(nDOFsQuad);

  /*--- Loop over all the nodes of the hexahedron and pick the correct
        ones for the faces. ---*/
  unsigned int ii = 0;
  for(unsigned short k=0; k<=nPoly; ++k) {
    for(unsigned short j=0; j<=nPoly; ++j) {
      for(unsigned short i=0; i<=nPoly; ++i, ++ii) {
        if(k == 0)     gridConnFaces[0].push_back(ii);
        if(k == nPoly) gridConnFaces[1].push_back(ii);
        if(j == 0)     gridConnFaces[2].push_back(ii);
        if(j == nPoly) gridConnFaces[3].push_back(ii);
        if(i == 0)     gridConnFaces[4].push_back(ii);
        if(i == nPoly) gridConnFaces[5].push_back(ii);
      }
    }
  }

  /*--- Make sure that the element is to the left of the faces. ---*/
  const unsigned short n0 = 0;
  const unsigned short n1 = nPoly;
  const unsigned short n2 = nDOFsQuad -1;
  const unsigned short n3 = n2 - nPoly;
  const unsigned short n4 = n0 + nDOFsQuad*nPoly;
  const unsigned short n5 = n1 + nDOFsQuad*nPoly;
  const unsigned short n6 = n2 + nDOFsQuad*nPoly;
  const unsigned short n7 = n3 + nDOFsQuad*nPoly;

  ChangeDirectionQuadConn(gridConnFaces[0], n0, n1, n2, n3);
  ChangeDirectionQuadConn(gridConnFaces[1], n4, n7, n6, n5);
  ChangeDirectionQuadConn(gridConnFaces[2], n0, n4, n5, n1);
  ChangeDirectionQuadConn(gridConnFaces[3], n3, n2, n6, n7);
  ChangeDirectionQuadConn(gridConnFaces[4], n0, n3, n7, n4);
  ChangeDirectionQuadConn(gridConnFaces[5], n1, n5, n6, n2);
}

void CFEMStandardHexGrid::SubConnLinearElements(void) {

  /*--- The hexahedron is split into several linear hexahedra.
        Set the VTK sub-types accordingly. ---*/
  VTK_SubType1 = HEXAHEDRON;
  VTK_SubType2 = NONE;

  /*--- Determine the nodal offset in j- and k-direction. ---*/
  const unsigned short jOff = nPoly+1;
  const unsigned short kOff = jOff*jOff;

  /*--- Loop over the subelements in k-direction. ---*/
  for(unsigned short k=0; k<nPoly; ++k) {

    /*--- Abbreviate the offset in k-direction used in the connectivity. ---*/
    const unsigned short kk = k*kOff;

    /*--- Loop over the subelements in j-direction. ---*/
    for(unsigned short j=0; j<nPoly; ++j) {

      /*--- Abbreviate the offset in j-direction used in the connectivity. ---*/
      const unsigned short jj = j*jOff;

      /*--- Loop over the subelements in i-direction. ---*/
      for(unsigned short i=0; i<nPoly; ++i) {

        /*--- Determine the 8 vertices of this subhexahedron and store
              them in subConn1ForPlotting. ---*/
        const unsigned short n0 = kk + jj + i;
        const unsigned short n1 = n0 + 1;
        const unsigned short n2 = n1 + jOff;
        const unsigned short n3 = n0 + jOff;
        const unsigned short n4 = n0 + kOff;
        const unsigned short n5 = n1 + kOff;
        const unsigned short n6 = n2 + kOff;
        const unsigned short n7 = n3 + kOff;

        subConn1ForPlotting.push_back(n0);
        subConn1ForPlotting.push_back(n1);
        subConn1ForPlotting.push_back(n2);
        subConn1ForPlotting.push_back(n3);
        subConn1ForPlotting.push_back(n4);
        subConn1ForPlotting.push_back(n5);
        subConn1ForPlotting.push_back(n6);
        subConn1ForPlotting.push_back(n7);
      }
    }
  }
}

void CFEMStandardHexGrid::TensorProductSolDOFs(const int                           N,
                                               const ColMajorMatrix<passivedouble> &Ai,
                                               const ColMajorMatrix<passivedouble> &Aj,
                                               const ColMajorMatrix<passivedouble> &Ak,
                                               const ColMajorMatrix<su2double>     &B,
                                               ColMajorMatrix<su2double>           &C,
                                               const CConfig                       *config) {

  /*--- Call the function to which TensorProductDataVolSolDOFs points to carry out
        the actual tensor product. Perform the timing, if desired. ---*/
#ifdef PROFILE
  double timeGemm;
  if( config ) config->TensorProduct_Tick(&timeGemm);
#endif

  TensorProductDataVolSolDOFs(N, B.rows(), C.rows(), Ai.data(), Aj.data(), Ak.data(),
                              B.data(), C.data());

#ifdef PROFILE
  if( config ) config->TensorProduct_Tock(timeGemm, 3, N, nDOFs1D, rLineSolDOFs.size());
#endif
}
