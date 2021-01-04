/*!
 * \file CFEMStandardQuadBase.cpp
 * \brief Functions for the class CFEMStandardQuadBase.
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

#include "../../include/fem/CFEMStandardQuadBase.hpp"
#include "../../include/fem/fem_gauss_jacobi_quadrature.hpp"

/*----------------------------------------------------------------------------------*/
/*              Public member functions of CFEMStandardQuadBase.                    */
/*----------------------------------------------------------------------------------*/

CFEMStandardQuadBase::CFEMStandardQuadBase(const unsigned short val_nPoly,
                                           const unsigned short val_orderExact)
  : CFEMStandardLineBase() {

  /*--- Store the command line arguments. ---*/
  VTK_Type   = QUADRILATERAL;
  nPoly      = val_nPoly;
  orderExact = val_orderExact;

  /*--- Determine the number of DOFs in 1D and the total number of DOFs.
        Also determine the padded value of the latter. ---*/
  nDOFs1D  = nPoly + 1;
  nDOFs    = nDOFs1D*nDOFs1D;
  nDOFsPad = PaddedValue(nDOFs);

  /*--- Determine the number of integration points in 1D as well as the total number
        of integration points. Also determine the padded value of the latter. ---*/
  nInt1D          = orderExact/2 + 1;
  nIntegration    = nInt1D*nInt1D;
  nIntegrationPad = PaddedValue(nIntegration);

  /*--- Determine the location and the weights of the 1D integration points.
        The 2D integration points are obtained via a tensor product. ---*/
  rLineInt.resize(nInt1D);
  wLineInt.resize(nInt1D);

  CGaussJacobiQuadrature GaussJacobi;
  GaussJacobi.GetQuadraturePoints(0.0, 0.0, -1.0, 1.0, rLineInt, wLineInt);

  /*--- Allocate the memory for the padded number of integration points and
        initialize the weights to zero. This is done such that the padded
        values are initialized appropriately. ---*/
  wIntegration.resize(nIntegrationPad);
  wIntegration.setConstant(0.0);

  /*--- Double loop to set the integration weights in the 2D integration points. ---*/
  unsigned short ii = 0;
  for(unsigned short j=0; j<nInt1D; ++j)
    for(unsigned short i=0; i<nInt1D; ++i, ++ii)
      wIntegration(ii) = wLineInt[i]*wLineInt[j];
}

/*----------------------------------------------------------------------------------*/
/*                Protected member functions of CFEMStandardQuadBase.               */
/*----------------------------------------------------------------------------------*/

void CFEMStandardQuadBase::ChangeDirectionQuadConn(vector<unsigned short> &connQuad,
                                                   const unsigned short   vert0,
                                                   const unsigned short   vert1,
                                                   const unsigned short   vert2,
                                                   const unsigned short   vert3) {

  /*--- Determine the indices of the 4 corner vertices of the quad. ---*/
  const unsigned short ind0 = 0;
  const unsigned short ind1 = nPoly;
  const unsigned short ind2 = (nPoly+1)*(nPoly+1) -1;
  const unsigned short ind3 = ind2 - nPoly;

  /*--- There exists a linear mapping from the indices of the numbering used in the
        connectivity of this face to the indices of the target numbering. This
        mapping is of the form ii = a + b*i + c*j and jj = d + e*i + f*j, where
        ii,jj are the indices of the target numbering and i,j the indices of the
        numbering used for this face. The values of the coefficients a,b,c,d,e,f
        depend on how the corner points coincide with each other. This is
        determined below. The bool verticesDontMatch is there to check if vertices
        do not match. This should not happen, but it is checked for security. ---*/

  signed short a=0, b=0, c=0, d=0, e=0, f=0;  // Initialization to avoid a compiler warning.
  bool verticesDontMatch = false;

  if(vert0 == connQuad[ind0]) {

    /*--- Vert0 coincides with the first vertex of the face connectivity.
          Set the coefficients a and d accordingly.   ---*/
    a = d = 0;
    if(vert2 != connQuad[ind2]) verticesDontMatch = true;

    /*--- Check the situation for the neighboring vertices. ---*/
    if(vert1 == connQuad[ind1]) {

      /*--- The vertex numbering is the same for both faces. ---*/
      if(vert3 != connQuad[ind3]) verticesDontMatch = true;

      b = f = 1; c = e = 0;
    }
    else if(vert1 == connQuad[ind3]) {

      /*--- The i and j numbering are swapped. ---*/
      if(vert3 != connQuad[ind1]) verticesDontMatch = true;

      b = f = 0; c = e = 1;
    }
    else {
      verticesDontMatch = true;
    }
  }
  else if(vert1 == connQuad[ind0]) {

    /*--- Vert1 coincides with the first vertex of the face connectivity.
          Set the coefficients a and d accordingly.  ---*/
    a = nPoly; d = 0;
    if(vert3 != connQuad[ind2]) verticesDontMatch = true;

    /*--- Check the situation for the neighboring vertices. ---*/
    if(vert0 == connQuad[ind1]) {

      /*--- The i-direction is negated while the j-direction coincides. ---*/
      if(vert2 != connQuad[ind3]) verticesDontMatch = true;

      b = -1; f = 1; c = e = 0;
    }
    else if(vert0 == connQuad[ind3]) {

      /*--- The j-direction of the current face corresponds with the negative
            i-direction of the target, while the i-direction coincides with
            the j-direction of the target.     ---*/
      if(vert2 != connQuad[ind1]) verticesDontMatch = true;

      b = f = 0; c = -1; e = 1;
    }
    else {
      verticesDontMatch = true;
    }
  }
  else if(vert2 == connQuad[ind0]) {

    /*--- Vert2 coincides with the first vertex of the face connectivity.
          Set the coefficients a and d accordingly.  ---*/
    a = d = nPoly;
    if(vert0 != connQuad[ind2]) verticesDontMatch = true;

    /*--- Check the situation for the neighboring vertices. ---*/
    if(vert1 == connQuad[ind3]) {

      /*--- Both the i and j-direction are negated. ---*/
      if(vert3 != connQuad[ind1]) verticesDontMatch = true;

      b = f = -1; c = e = 0;
    }
    else if(vert1 == connQuad[ind1]) {

      /*--- The i and j-direction are negated and swapped. ---*/
      if(vert3 != connQuad[ind3]) verticesDontMatch = true;

      b = f = 0; c = e = -1;
    }
    else {
      verticesDontMatch = true;
    }
  }
  else if(vert3 == connQuad[ind0]) {

    /*--- Vert3 coincides with the first vertex of the face connectivity.
          Set the coefficients a and d accordingly.  ---*/
    a = 0; d = nPoly;
    if(vert1 != connQuad[ind2]) verticesDontMatch = true;

    /*--- Check the situation for the neighboring vertices. ---*/
    if(vert0 == connQuad[ind3]) {

      /*--- The i-direction coincides while the j-direction is negated. ---*/
      if(vert2 != connQuad[ind1]) verticesDontMatch = true;

      b = 1; f = -1; c = e = 0;
    }
    else if(vert0 == connQuad[ind1]) {

      /*--- The j-direction of the current face corresponds with the i-direction
            of the target, while the i-direction coincides with the negative
            j-direction of the target.    ---*/
      if(vert2 != connQuad[ind3]) verticesDontMatch = true;

      b = f = 0; c = 1; e = -1;
    }
    else {
      verticesDontMatch = true;
    }
  }
  else {
    verticesDontMatch = true;
  }

  /*--- If non-matching vertices have been found, terminate with an error message. ---*/
  if( verticesDontMatch )
    SU2_MPI::Error("Corner vertices do not match. This should not happen.",
                   CURRENT_FUNCTION);

  /*--- Copy the connectivity, such that things works out correctly when carrying
        out the renumbering.      ---*/
  vector<unsigned short> connQuadOr = connQuad;

  /*--- Loop over the vertices of the original face to copy the connectivity data. ---*/
  unsigned short ind = 0;
  for(unsigned short j=0; j<=nPoly; ++j) {
    for(unsigned short i=0; i<=nPoly; ++i, ++ind) {

      /*--- Determine the ii and jj indices of the target, convert it to
            a 1D index and shore the modified index in connQuad. ---*/
      const unsigned short ii = a + i*b + j*c;
      const unsigned short jj = d + i*e + j*f;

      const unsigned short iind = jj*(nPoly+1) + ii;

      connQuad[iind] = connQuadOr[ind];
    }
  }
}

void CFEMStandardQuadBase::CreateTensorContributionsLineAdjQuad(const unsigned short                   mPointsLine,
                                                                const unsigned short                   mDOFs1D,
                                                                const unsigned short                   faceID_Elem,
                                                                const unsigned short                   orientation,
                                                                const ColMajorMatrix<passivedouble>    &b1DPoints,
                                                                const ColMajorMatrix<passivedouble>    &derB1DPoints,
                                                                const ColMajorMatrix<passivedouble>    &b1DM1,
                                                                const ColMajorMatrix<passivedouble>    &b1DP1,
                                                                const ColMajorMatrix<passivedouble>    &derB1DM1,
                                                                const ColMajorMatrix<passivedouble>    &derB1DP1,
                                                                vector<ColMajorMatrix<passivedouble> > &tensorSol,
                                                                vector<ColMajorMatrix<passivedouble> > &tensorDSolDr,
                                                                vector<ColMajorMatrix<passivedouble> > &tensorDSolDs) {

  /*--- Allocate the memory for the first index of tensorSol, etc. As this
        function is only called for 2D simulations, there are two
        contributions to the tensor products. ---*/
  tensorSol.resize(2);
  tensorDSolDr.resize(2);
  tensorDSolDs.resize(2);

  /*--- Step 1. Set the contribution to the tensor product for orientation == 0,
                depending on the face ID of the element. Flag the tangential
                component, stored in the second entry, for reversing, but do
                not carry this out yet. It may overruled when orientation == 1. ---*/
  bool reverseTangential = false;
  switch( faceID_Elem ) {

    case 0: {

      /*--- Face ID 0 of the quadrilateral. Set the components of the
            tensor product accordingly. ---*/
      tensorSol[0]    = b1DM1;    tensorSol[1]    = b1DPoints;
      tensorDSolDr[0] = b1DM1;    tensorDSolDr[1] = derB1DPoints;
      tensorDSolDs[0] = derB1DM1; tensorDSolDs[1] = b1DPoints;
      break;
    }

    case 1: {

      /*--- Face ID 1 of the quadrilateral. Set the components of the
            tensor product accordingly. ---*/
      tensorSol[0]    = b1DP1;    tensorSol[1]    = b1DPoints;
      tensorDSolDr[0] = derB1DP1; tensorDSolDr[1] = b1DPoints;
      tensorDSolDs[0] = b1DP1;    tensorDSolDs[1] = derB1DPoints;
      break;
    }

    case 2: {

      /*--- Face ID 2 of the quadrilateral. Set the components of the
            tensor product accordingly. ---*/
      tensorSol[0]    = b1DP1;    tensorSol[1]    = b1DPoints;
      tensorDSolDr[0] = b1DP1;    tensorDSolDr[1] = derB1DPoints;
      tensorDSolDs[0] = derB1DP1; tensorDSolDs[1] = b1DPoints;
      reverseTangential = true;
      break;
    }

    case 3: {

      /*--- Face ID 3 of the quadrilateral. Set the components of the
            tensor product accordingly. ---*/
      tensorSol[0]    = b1DM1;    tensorSol[1]    = b1DPoints;
      tensorDSolDr[0] = derB1DM1; tensorDSolDr[1] = b1DPoints;
      tensorDSolDs[0] = b1DM1;    tensorDSolDs[1] = derB1DPoints;
      reverseTangential = true;
      break;
    }

    default:
      SU2_MPI::Error(string("Invalid faceID. This should not happen."), CURRENT_FUNCTION);
  }

  /*--- Step 2. Reverse the tangential component, if needed. As the direction of
                the tangential component is opposite when the orientation == 1,
                first negate reverseTangential when this is the case.
                Note that for the reversing mPointsLine and mDOFs1D must be
                used and not the dimension of the contribution to the tensor
                product, because these contributions may be padded. ---*/
  if(orientation == 1) reverseTangential = !reverseTangential;

  if( reverseTangential ) {
    const unsigned short mPointsHalf = mPointsLine/2;
    unsigned short jj = mPointsLine-1;
    for(unsigned short j=0; j<mPointsHalf; ++j, --jj) {
      for(unsigned short i=0; i<mDOFs1D; ++i) {
        swap(tensorSol[1](j,i),    tensorSol[1](jj,i)); 
        swap(tensorDSolDr[1](j,i), tensorDSolDr[1](jj,i));
        swap(tensorDSolDs[1](j,i), tensorDSolDs[1](jj,i));
      }
    }
  }
}

void CFEMStandardQuadBase::SubConnLinearElements(void) {

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

void CFEMStandardQuadBase::SetFunctionPointerVolumeDataQuad(const unsigned short K,
                                                            const unsigned short M,
                                                            TPI2D                &TPVolData) {

  /*--- Create the map with the function pointers to carry out the tensor product
        to compute the data in the volume points of the quadrilateral. ---*/
  map<CUnsignedShort2T, TPI2D> mapFunctions;
  CreateMapTensorProductVolumeIntPoints2D(mapFunctions);

  /*--- Try to find the combination of K and M in mapFunctions. If not found,
        write a clear error message that this tensor product is not supported. ---*/
  CUnsignedShort2T KM(K, M);
  auto MI = mapFunctions.find(KM);
  if(MI == mapFunctions.end()) {
    std::ostringstream message;
    message << "The tensor product TensorProductVolumeIntPoints2D_" << K
            << "_" << M << " not created by the automatic source code "
            << "generator. Modify this automatic source code creator";
    SU2_MPI::Error(message.str(), CURRENT_FUNCTION);
  }

  /*--- Set the function pointer to carry out tensor product. ---*/
  TPVolData = MI->second; 
}

void CFEMStandardQuadBase::TensorProductVolumeDataQuad(TPI2D                               &TPVolData,
                                                       const int                           N,
                                                       const int                           K,
                                                       const int                           M,
                                                       const ColMajorMatrix<passivedouble> &Ai,
                                                       const ColMajorMatrix<passivedouble> &Aj,
                                                       const ColMajorMatrix<su2double>     &B,
                                                       ColMajorMatrix<su2double>           &C,
                                                       const CConfig                       *config) {

  /*--- Call the function to which TPVolData points to carry out the
        actual tensor product. Perform the timing, if desired. ---*/
#ifdef PROFILE
  double timeGemm;
  if( config ) config->TensorProduct_Tick(&timeGemm);
#endif

  TPVolData(N, B.rows(), C.rows(), Ai.data(), Aj.data(), B.data(), C.data());

#ifdef PROFILE
  if( config ) config->TensorProduct_Tock(timeGemm, 2, N, K, M);
#endif
}
