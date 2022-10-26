/*!
 * \file CFEMStandardHexBase.cpp
 * \brief Functions for the class CFEMStandardHexBase.
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

#include "../../include/fem/CFEMStandardHexBase.hpp"
#include "../../include/fem/fem_gauss_jacobi_quadrature.hpp"

/*----------------------------------------------------------------------------------*/
/*                  Public member functions of CFEMStandardHexBase.                 */
/*----------------------------------------------------------------------------------*/

CFEMStandardHexBase::CFEMStandardHexBase(const unsigned short val_nPoly,
                                         const unsigned short val_orderExact)
  : CFEMStandardQuadBase() {

  /*--- Store the command line arguments. ---*/
  VTK_Type   = HEXAHEDRON;
  nPoly      = val_nPoly;
  orderExact = val_orderExact;

  /*--- Determine the number of DOFs in 1D and the total number of DOFs.
        Also determine the padded value of the latter. ---*/
  nDOFs1D  = nPoly + 1;
  nDOFs    = nDOFs1D*nDOFs1D*nDOFs1D;
  nDOFsPad = PaddedValue(nDOFs);

  /*--- Determine the number of integration points in 1D as well as the total number
        of integration points. Also determine the padded value of the latter. ---*/
  nInt1D          = orderExact/2 + 1;
  nIntegration    = nInt1D*nInt1D*nInt1D;
  nIntegrationPad = PaddedValue(nIntegration);

  /*--- Determine the location and the weights of the 1D integration points.
        The 3D integration points are obtained via a tensor product. ---*/
  rLineInt.resize(nInt1D);
  wLineInt.resize(nInt1D);

  CGaussJacobiQuadrature GaussJacobi;
  GaussJacobi.GetQuadraturePoints(0.0, 0.0, -1.0, 1.0, rLineInt, wLineInt);

  /*--- Allocate the memory for the padded number of integration points and
        initialize the weights to zero. This is done such that the padded
        values are initialized appropriately. ---*/
  wIntegration.resize(nIntegrationPad);
  wIntegration.setConstant(0.0);

  /*--- Triple loop to set the integration weights in the 3D integration points. ---*/
  unsigned short ii = 0;
  for(unsigned short k=0; k<nInt1D; ++k)
    for(unsigned short j=0; j<nInt1D; ++j)
      for(unsigned short i=0; i<nInt1D; ++i, ++ii)
        wIntegration(ii) = wLineInt[i]*wLineInt[j]*wLineInt[k];
}

/*----------------------------------------------------------------------------------*/
/*                Protected member functions of CFEMStandardHex.                    */
/*----------------------------------------------------------------------------------*/

void CFEMStandardHexBase::ConvertCoor2DQuadFaceTo3DHex(const vector<passivedouble> rLine,
                                                       const unsigned short        faceID,
                                                       const unsigned short        orientation,
                                                       bool                        &swapTangInTensor,
                                                       vector<passivedouble>       &rNorm,
                                                       vector<passivedouble>       &rTang0,
                                                       vector<passivedouble>       &rTang1) {

  /*--- Determine the value of the coordinate in normal direction. This is
        either -1 or 1, depending on the face. ---*/
  rNorm.resize(1);
  switch( faceID ) {
    case 0: case 2: case 4: rNorm[0] = -1.0; break;
    case 1: case 3: case 5: rNorm[0] =  1.0; break;
    default:
      SU2_MPI::Error(string("Invalid faceID. This should not happen."), CURRENT_FUNCTION);
  }

  /*--- Copy rLine into rTang0 and rTang1 as initialization. ---*/
  rTang0 = rLine;
  rTang1 = rLine;

  /*--- Initialize swapTangInTensor, revT0 and revT1 to false. swapTangInTensor indicates
        whether or not the two tangential directions in the quad must be swapped when
        the tensor product is carried out, while revT0 and revT1 indicate whether or not
        the tangential directions in the reference frame of the quad are reversed
        compared to the frame of the hex. ---*/
  swapTangInTensor = false;
  bool revT0       = false;
  bool revT1       = false;

  /*--- Set the above booleans to true if this is the case for the given
        face ID and orientation. ---*/
  switch( faceID ) {
    case 0: case 3: case 4: {
      if((orientation == 1) || (orientation == 3)) swapTangInTensor = true;
      if((orientation == 2) || (orientation == 3)) revT0 = true;
      if((orientation == 3) || (orientation == 4)) revT1 = true;
      break;
    }

    case 1: case 2: case 5: {
      if((orientation == 0) || (orientation == 2) || (orientation == 4)) swapTangInTensor = true;
      if((orientation == 3) || (orientation == 4)) revT0 = true;
      if((orientation == 2) || (orientation == 3)) revT1 = true;
      break;
    }

    default:
      SU2_MPI::Error(string("Invalid faceID. This should not happen."), CURRENT_FUNCTION);
  }

  /*--- Reverse the tangential contributions, if needed. ---*/
  if( revT0 ) reverse(rTang0.begin(), rTang0.end());
  if( revT1 ) reverse(rTang1.begin(), rTang1.end());
}

void CFEMStandardHexBase::SubConnLinearElements(void) {

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

void CFEMStandardHexBase::SetFunctionPointerVolumeDataHex(const unsigned short K,
                                                          const unsigned short M,
                                                          TPI3D                &TPVolData) {

  /*--- Create the map with the function pointers to carry out the tensor product
        to compute the data in the volume points of the quadrilateral. ---*/
  map<CUnsignedShort2T, TPI3D> mapFunctions;
  CreateMapTensorProductVolumeIntPoints3D(mapFunctions);

  /*--- Try to find the combination of K and M in mapFunctions. If not found,
        write a clear error message that this tensor product is not supported. ---*/
  CUnsignedShort2T KM(K, M);
  auto MI = mapFunctions.find(KM);
  if(MI == mapFunctions.end()) {
    std::ostringstream message;
    message << "The tensor product TensorProductVolumeIntPoints3D_" << K
            << "_" << M << " not created by the automatic source code "
            << "generator. Modify this automatic source code creator";
    SU2_MPI::Error(message.str(), CURRENT_FUNCTION);
  }

  /*--- Set the function pointer to carry out tensor product. ---*/
  TPVolData = MI->second;
}

void CFEMStandardHexBase::TensorProductVolumeDataHex(TPI3D                               &TPVolData,
                                                     const int                           N,
                                                     const int                           K,
                                                     const int                           M,
                                                     const ColMajorMatrix<passivedouble> &Ai,
                                                     const ColMajorMatrix<passivedouble> &Aj,
                                                     const ColMajorMatrix<passivedouble> &Ak,
                                                     const ColMajorMatrix<su2double>     &B,
                                                     ColMajorMatrix<su2double>           &C,
                                                     const bool                          initZero,
                                                     const CConfig                       *config) {

  /*--- Call the function to which TPVolData points to carry out the
        actual tensor product. Perform the timing, if desired. ---*/
#ifdef PROFILE
  double timeGemm;
  if( config ) config->TensorProduct_Tick(&timeGemm);
#endif

  if( initZero ) C.setConstant(0.0);
  TPVolData(N, B.rows(), C.rows(), Ai.data(), Aj.data(), Ak.data(), B.data(), C.data());

#ifdef PROFILE
  if( config ) config->TensorProduct_Tock(timeGemm, 3, N, K, M);
#endif
}

void CFEMStandardHexBase::LocalGridConnFaces(void) {

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