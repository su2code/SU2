/*!
 * \file CFEMStandardPyraBase.cpp
 * \brief Functions for the class CFEMStandardPyraBase.
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

#include "../../include/fem/CFEMStandardPyraBase.hpp"
#include "../../include/fem/fem_gauss_jacobi_quadrature.hpp"
#include "../../include/toolboxes/CSquareMatrixCM.hpp"

/*----------------------------------------------------------------------------------*/
/*            Public member functions of CFEMStandardPyraBase.                      */
/*----------------------------------------------------------------------------------*/

CFEMStandardPyraBase::CFEMStandardPyraBase(const unsigned short val_nPoly,
                                           const unsigned short val_orderExact)
  : CFEMStandardQuadBase(),
    CFEMStandardTriBase() {

  /*--- Store the command line arguments. ---*/
  VTK_Type   = PYRAMID;
  nPoly      = val_nPoly;
  orderExact = val_orderExact;

  /*--- Determine the number of DOFs in 1D of the base quad and the total number of DOFs.
        Also determine the padded value of the latter. ---*/
  nDOFs1D  = nPoly + 1;
  nDOFs    = nDOFs1D*(nDOFs1D+1)*(2*nDOFs1D+1)/6;
  nDOFsPad = PaddedValue(nDOFs);

  /*--- The 3D quadrature rule for a pyramid is obtained by transforming the
        standard pyramid into a standard hexahedron by means of the Duffy
        transformation. Hence the integration rule is a tensor product,
        albeit a special one. ---*/
  nInt1D          = orderExact/2 + 1;
  nIntegration    = nInt1D*nInt1D*nInt1D;
  nIntegrationPad = PaddedValue(nIntegration);

  /*--- Determine the location and the weights of the 1D Gauss-Legendre
        integration points. ---*/
  rLineIntGL.resize(nInt1D);
  wLineIntGL.resize(nInt1D);

  CGaussJacobiQuadrature GaussJacobi;
  GaussJacobi.GetQuadraturePoints(0.0, 0.0, -1.0, 1.0, rLineIntGL, wLineIntGL);

  /*--- Determine the location and the weights of the 1D Gauss-Jacobi
        integration points for alpha = 2 and beta = 0. The alpha = 2
        comes from the transformation from the pyramid to the hexahedron.
        The determinant of the Jacobian of this transformation is 0.25*(1-t)^2,
        hence alpha = 2 in the Gauss-Jacobi integration rule. ---*/
  rLineIntGJ.resize(nInt1D);
  wLineIntGJ.resize(nInt1D);
  GaussJacobi.GetQuadraturePoints(2.0, 0.0, -1.0, 1.0, rLineIntGJ, wLineIntGJ);

  /*--- Allocate the memory for the padded number of integration points and
        initialize the weights to zero. This is done such that the padded
        values are initialized appropriately. ---*/
  wIntegration.resize(nIntegrationPad);
  wIntegration.setConstant(0.0);

  /*--- Triple loop to set the integration weights in the 3D integration points.
        In k-direction, which corresponds with the direction from bottom to top,
        Gauss-Jacobi must be used. Furthermore, a factor of 0.25 must be used,
        because the determinant of the transformation from pyramid to hexahedron
        is 0.25*(1-t)^2. The Gauss-Jacobi rule only accounts for (1-t)^2.---*/
  unsigned short ii = 0;
  for(unsigned short k=0; k<nInt1D; ++k)
    for(unsigned short j=0; j<nInt1D; ++j)
      for(unsigned short i=0; i<nInt1D; ++i, ++ii)
        wIntegration(ii) = 0.25*wLineIntGL[i]*wLineIntGL[j]*wLineIntGJ[k];
}

/*----------------------------------------------------------------------------------*/
/*            Protected member functions of CFEMStandardPyraBase.                   */
/*----------------------------------------------------------------------------------*/

void CFEMStandardPyraBase::ConvertCoor2DQuadFaceTo3DPyra(const vector<passivedouble> &rLine,
                                                         const unsigned short        faceID_Elem,
                                                         const unsigned short        orientation,
                                                         vector<passivedouble>       &rPyra,
                                                         vector<passivedouble>       &sPyra,
                                                         vector<passivedouble>       &tPyra) {

  /*--- A pyramid has one quadrilateral face, which should have a faceID == 0.
        Check this. ---*/
  if(faceID_Elem != 0)
    SU2_MPI::Error(string("Invalid faceID. This should not happen."), CURRENT_FUNCTION);

  /*--- Determine the number of points on the quadrilateral face, which is a tensor
        product of rLine in two dimensions. Afterwards, allocate the memory for
        rPyra, sPyra and tPyra. ---*/
  const unsigned short nPoints1D = rLine.size();
  const unsigned short nPoints   = nPoints1D*nPoints1D;

  rPyra.resize(nPoints);
  sPyra.resize(nPoints);
  tPyra.resize(nPoints);

  /*--- Create the parametric coordinates in the frame of the face via a tensor product. ---*/
  vector<passivedouble> rF(nPoints), sF(nPoints);

  unsigned short k = 0;
  for(unsigned short j=0; j<nPoints1D; ++j) {
    for(unsigned short i=0; i<nPoints1D; ++i, ++k) {
      rF[k] = rLine[i];
      sF[k] = rLine[j];
    }
  }

  /*--- The quadrilateral face corresponds to parametri coordinate t == -1.
        Set this value. ---*/
  for(k=0; k<nPoints; ++k) tPyra[k] = -1.0;

  /*--- The r- and s-coordinates depend on the orientation of the face. ---*/
  switch( orientation ) {

    case 0: for(k=0; k<nPoints; ++k){rPyra[k] =  rF[k]; sPyra[k] =  sF[k];} break;
    case 1: for(k=0; k<nPoints; ++k){rPyra[k] =  sF[k]; sPyra[k] =  rF[k];} break;
    case 2: for(k=0; k<nPoints; ++k){rPyra[k] = -rF[k]; sPyra[k] =  sF[k];} break;
    case 3: for(k=0; k<nPoints; ++k){rPyra[k] = -sF[k]; sPyra[k] = -rF[k];} break;
    case 4: for(k=0; k<nPoints; ++k){rPyra[k] =  rF[k]; sPyra[k] = -sF[k];} break;
    default:
      SU2_MPI::Error(string("Invalid orientation for face 0. This should not happen."),
                     CURRENT_FUNCTION);
  }
}

void CFEMStandardPyraBase::ConvertCoor2DTriFaceTo3DPyra(const vector<passivedouble> &rF,
                                                        const vector<passivedouble> &sF,
                                                        const unsigned short        faceID_Elem,
                                                        const unsigned short        orientation,
                                                        vector<passivedouble>       &rPyra,
                                                        vector<passivedouble>       &sPyra,
                                                        vector<passivedouble>       &tPyra) {

  /*--- Determine the number of points on the triangular face. Afterwards, allocate
        the memory for rPyra, sPyra and tPyra. ---*/
  const unsigned short nP = rF.size();

  rPyra.resize(nP);
  sPyra.resize(nP);
  tPyra.resize(nP);

  /*--- Abbreviate rPyra, sPyra and tPyra, such that the different cases
        can be put on one line. ---*/
  vector<passivedouble> &r = rPyra, &s = sPyra, &t = tPyra;

  /*--- The values of rPyra, sPyra and tPyra depend on both the face ID in the numbering
        of the pyramid as well as the orientation of the face w.r.t. the pyramid.
        Make this distinction and set the values accordingly. ---*/
  unsigned short k;
  switch( faceID_Elem ) {

    case 1: {
      switch( orientation ) {
        case 0: for(k=0; k<nP; ++k) {r[k]= 0.5*rF[k]+    sF[k]+0.5; s[k]= 0.5*rF[k]          -0.5; t[k]= rF[k];}           break;
        case 1: for(k=0; k<nP; ++k) {r[k]=     rF[k]+0.5*sF[k]-0.5; s[k]=           0.5*sF[k]-0.5; t[k]=       sF[k];}     break;
        case 2: for(k=0; k<nP; ++k) {r[k]=-0.5*rF[k]+0.5*sF[k];     s[k]=-0.5*rF[k]-0.5*sF[k]-1.0; t[k]=-rF[k]-sF[k]-1.0;} break;
        case 3: for(k=0; k<nP; ++k) {r[k]=-0.5*rF[k]-    sF[k]-0.5; s[k]= 0.5*rF[k]          -0.5; t[k]= rF[k];}           break;
        default:
          SU2_MPI::Error(string("Invalid orientation for face 1. This should not happen."),
                         CURRENT_FUNCTION);
      }
      break;
    }

    case 2: {
      switch( orientation ) {
        case 0: for(k=0; k<nP; ++k) {r[k]=     rF[k]+0.5*sF[k]+0.5; s[k]=          -0.5*sF[k]+0.5; t[k]=       sF[k];}     break;
        case 1: for(k=0; k<nP; ++k) {r[k]= 0.5*rF[k]+    sF[k]+0.5; s[k]=-0.5*rF[k]          +0.5; t[k]= rF[k];}           break;
        case 2: for(k=0; k<nP; ++k) {r[k]=    -rF[k]-0.5*sF[k]+0.5; s[k]=          -0.5*sF[k]+0.5; t[k]=       sF[k];}     break;
        case 3: for(k=0; k<nP; ++k) {r[k]= 0.5*rF[k]-0.5*sF[k];     s[k]= 0.5*rF[k]+0.5*sF[k]+1.0; t[k]=-rF[k]-sF[k]-1.0;} break;
        default:
          SU2_MPI::Error(string("Invalid orientation for face 2. This should not happen."),
                         CURRENT_FUNCTION);
      }
      break;
    }

    case 3: {
      switch( orientation ) {
        case 0: for(k=0; k<nP; ++k) {r[k]=           0.5*sF[k]-0.5; s[k]=     rF[k]+0.5*sF[k]+0.5; t[k]=       sF[k];}     break;
        case 1: for(k=0; k<nP; ++k) {r[k]= 0.5*rF[k]          -0.5; s[k]= 0.5*rF[k]+    sF[k]+0.5; t[k]= rF[k];}           break;
        case 2: for(k=0; k<nP; ++k) {r[k]=           0.5*sF[k]-0.5; s[k]=    -rF[k]-0.5*sF[k]-0.5; t[k]=       sF[k];}     break;
        case 3: for(k=0; k<nP; ++k) {r[k]=-0.5*rF[k]-0.5*sF[k]-1.0; s[k]= 0.5*rF[k]-0.5*sF[k];     t[k]=-rF[k]-sF[k]-1.0;} break;
        default:
          SU2_MPI::Error(string("Invalid orientation for face 3. This should not happen."),
                         CURRENT_FUNCTION);
      }
      break;
    }

    case 4: {
      switch( orientation ) {
        case 0: for(k=0; k<nP; ++k) {r[k]=-0.5*rF[k]          +0.5; s[k]= 0.5*rF[k]+    sF[k]+0.5; t[k]= rF[k];}           break;
        case 1: for(k=0; k<nP; ++k) {r[k]=          -0.5*sF[k]+0.5; s[k]=     rF[k]+0.5*sF[k]+0.5; t[k]=       sF[k];}     break;
        case 2: for(k=0; k<nP; ++k) {r[k]= 0.5*rF[k]+0.5*sF[k]+1.0; s[k]=-0.5*rF[k]+0.5*sF[k];     t[k]=-rF[k]-sF[k]-1.0;} break;
        case 3: for(k=0; k<nP; ++k) {r[k]=-0.5*rF[k]          +0.5; s[k]=-0.5*rF[k]-    sF[k]-0.5; t[k]= rF[k];}           break;
        default:
          SU2_MPI::Error(string("Invalid orientation for face 4. This should not happen."),
                         CURRENT_FUNCTION);
      }
      break;
    }

    default:
      SU2_MPI::Error(string("Invalid faceID. This should not happen."), CURRENT_FUNCTION);
  }
}

void CFEMStandardPyraBase::DerLagBasisIntPointsPyra(const unsigned short                   mPoly,
                                                    const vector<passivedouble>            &rDOFs,
                                                    const vector<passivedouble>            &sDOFs,
                                                    const vector<passivedouble>            &tDOFs,
                                                    const vector<passivedouble>            &rInt,
                                                    const vector<passivedouble>            &sInt,
                                                    const vector<passivedouble>            &tInt,
                                                    vector<ColMajorMatrix<passivedouble> > &derLag) {

  /*--- Determine the padded number of the total number of integration points. ---*/
  const unsigned short nIntTot    = rInt.size();
  const unsigned short nIntTotPad = PaddedValue(nIntTot);

  /*--- Determine the inverse of the Vandermonde matrix of the DOFs. ---*/
  CSquareMatrixCM VInv(rDOFs.size());
  VandermondePyramid(mPoly, rDOFs, sDOFs, tDOFs, VInv.GetMat());
  VInv.Invert();

  /*--- Determine the gradient of the Vandermonde matrix of the integration points. Make
        sure to allocate the number of rows to nIntTotPad and initialize them to zero. ---*/
  ColMajorMatrix<passivedouble> VDr(nIntTotPad,rDOFs.size()),
                                VDs(nIntTotPad,rDOFs.size()),
                                VDt(nIntTotPad,rDOFs.size());
  VDr.setConstant(0.0);
  VDs.setConstant(0.0);
  VDt.setConstant(0.0);

  GradVandermondePyramid(mPoly, rInt, sInt, tInt, VDr, VDs, VDt);

  /*--- TEST ----*/
/*const passivedouble dr = 1.e-4;
  vector<passivedouble> rMin = rInt, rPlus = rInt;
  vector<passivedouble> sMin = sInt, sPlus = sInt;
  vector<passivedouble> tMin = tInt, tPlus = tInt;

  for(unsigned short i=0; i<rInt.size(); ++i) {
    rMin[i] -= dr; rPlus[i] += dr; 
    sMin[i] -= dr; sPlus[i] += dr;
    tMin[i] -= dr; tPlus[i] += dr;
  }

  ColMajorMatrix<passivedouble> VMin(nIntTotPad,rDOFs.size()),
                                VPlus(nIntTotPad,rDOFs.size());

  VMin.setConstant(0.0);
  VPlus.setConstant(0.0);

  VandermondePyramid(mPoly, rInt, sInt, tMin,  VMin);
  VandermondePyramid(mPoly, rInt, sInt, tPlus, VPlus);

  passivedouble maxDiff = 0.0;
  for(unsigned i=0; i<rDOFs.size(); ++i) {
    for(unsigned short j=0; j<rInt.size(); ++j) {

      const passivedouble diff = (VPlus(j,i) - VMin(j,i))/(2.0*dr);
      cout << "t comparison: " << i << " " << j << ": "
           << VDt(j,i) << " " << diff << " " << fabs(VDt(j,i)-diff) << endl;
      maxDiff = max(maxDiff, fabs(VDt(j,i)-diff));
    }
  }

  cout << endl << "maxDiff: " << maxDiff << endl << endl;
  SU2_MPI::Error("Test", CURRENT_FUNCTION); */

  /*--- End TEST. ---*/

  /*--- The gradients of the Lagrangian basis functions can be obtained by
        multiplying VDr, VDs, VDt and VInv. ---*/
  derLag.resize(3);
  VInv.MatMatMult('R', VDr, derLag[0]);
  VInv.MatMatMult('R', VDs, derLag[1]);
  VInv.MatMatMult('R', VDt, derLag[2]);

  /*--- Check if the sum of the elements of the relevant rows of derLag is 0. ---*/
  CheckRowSum(nIntTot, rDOFs.size(), 0.0, derLag[0]);
  CheckRowSum(nIntTot, rDOFs.size(), 0.0, derLag[1]);
  CheckRowSum(nIntTot, rDOFs.size(), 0.0, derLag[2]);
}

void CFEMStandardPyraBase::HesLagBasisIntPointsPyra(const unsigned short                   mPoly,
                                                    const vector<passivedouble>            &rDOFs,
                                                    const vector<passivedouble>            &sDOFs,
                                                    const vector<passivedouble>            &tDOFs,
                                                    const vector<passivedouble>            &rInt,
                                                    const vector<passivedouble>            &sInt,
                                                    const vector<passivedouble>            &tInt,
                                                    vector<ColMajorMatrix<passivedouble> > &hesLag) {

  /*--- Determine the padded number of the total number of integration points. ---*/
  const unsigned short nIntTot    = rInt.size();
  const unsigned short nIntTotPad = PaddedValue(nIntTot);

  /*--- Determine the inverse of the Vandermonde matrix of the DOFs. ---*/
  CSquareMatrixCM VInv(rDOFs.size());
  VandermondePyramid(mPoly, rDOFs, sDOFs, tDOFs, VInv.GetMat());
  VInv.Invert();

  /*--- Determine the Hessian of the Vandermonde matrix of the integration points. Make
        sure to allocate the number of rows to nIntTotPad and initialize them to zero. ---*/
  ColMajorMatrix<passivedouble> VDr2(nIntTotPad,rDOFs.size()),
                                VDs2(nIntTotPad,rDOFs.size()),
                                VDt2(nIntTotPad,rDOFs.size()),
                                VDrs(nIntTotPad,rDOFs.size()),
                                VDrt(nIntTotPad,rDOFs.size()),
                                VDst(nIntTotPad,rDOFs.size());
  VDr2.setConstant(0.0);
  VDs2.setConstant(0.0);
  VDt2.setConstant(0.0);
  VDrs.setConstant(0.0);
  VDrt.setConstant(0.0);
  VDst.setConstant(0.0);

  HesVandermondePyramid(mPoly, rInt, sInt, tInt, VDr2, VDs2, VDt2, VDrs, VDrt, VDst);

  /*--- TEST ----*/
/*const passivedouble dr = 1.e-4;
  vector<passivedouble> rMin = rInt, rPlus = rInt;
  vector<passivedouble> sMin = sInt, sPlus = sInt;
  vector<passivedouble> tMin = tInt, tPlus = tInt;

  for(unsigned short i=0; i<rInt.size(); ++i) {
    rMin[i] -= dr; rPlus[i] += dr; 
    sMin[i] -= dr; sPlus[i] += dr;
    tMin[i] -= dr; tPlus[i] += dr;
  } */

/*ColMajorMatrix<passivedouble> VMin(nIntTotPad,rDOFs.size()),
                                V(nIntTotPad,rDOFs.size()),
                                VPlus(nIntTotPad,rDOFs.size());

  VMin.setConstant(0.0);
  V.setConstant(0.0);
  VPlus.setConstant(0.0);

  VandermondePyramid(mPoly, rInt, sInt, tMin,  VMin);
  VandermondePyramid(mPoly, rInt, sInt, tInt,  V);
  VandermondePyramid(mPoly, rInt, sInt, tPlus, VPlus);

  passivedouble maxDiff = 0.0;
  for(unsigned i=0; i<rDOFs.size(); ++i) {
    for(unsigned short j=0; j<rInt.size(); ++j) {

      const passivedouble d2t = (VMin(j,i) - 2.0*V(j,i) + VPlus(j,i))/(dr*dr);
      cout << "t2 comparison: " << i << " " << j << ": "
           << VDt2(j,i) << " " << d2t << " " << fabs(VDt2(j,i)-d2t) << endl;
      maxDiff = max(maxDiff, fabs(VDt2(j,i)-d2t));
    }
  }

  cout << endl << "maxDiff: " << maxDiff << endl << endl; */

/*ColMajorMatrix<passivedouble> VMinMin(nIntTotPad,rDOFs.size()),
                                VMinPlus(nIntTotPad,rDOFs.size()),
                                VPlusMin(nIntTotPad,rDOFs.size()),
                                VPlusPlus(nIntTotPad,rDOFs.size());

  VMinMin.setConstant(0.0);
  VMinPlus.setConstant(0.0);
  VPlusMin.setConstant(0.0);
  VPlusPlus.setConstant(0.0);

  VandermondePyramid(mPoly, rInt, sMin,  tMin,  VMinMin);
  VandermondePyramid(mPoly, rInt, sMin,  tPlus, VMinPlus);
  VandermondePyramid(mPoly, rInt, sPlus, tMin,  VPlusMin);
  VandermondePyramid(mPoly, rInt, sPlus, tPlus, VPlusPlus);

  passivedouble maxDiff = 0.0;
  for(unsigned i=0; i<rDOFs.size(); ++i) {
    for(unsigned short j=0; j<rInt.size(); ++j) {

      const passivedouble dst = (VPlusPlus(j,i) - VMinPlus(j,i) - VPlusMin(j,i) + VMinMin(j,i))/(4.0*dr*dr);
      cout << "st comparison: " << i << " " << j << ": "
           << VDst(j,i) << " " << dst << " " << fabs(VDst(j,i)-dst) << endl;
      maxDiff = max(maxDiff, fabs(VDst(j,i)-dst));
    }
  }

  cout << endl << "maxDiff: " << maxDiff << endl << endl; */
/*SU2_MPI::Error("Test", CURRENT_FUNCTION); */

  /*--- End TEST. ---*/

  /*--- The Hessian of the Lagrangian basis functions can be obtained by
        multiplying VDr2, VDs2, etc and VInv. ---*/
  hesLag.resize(6);
  VInv.MatMatMult('R', VDr2, hesLag[0]);
  VInv.MatMatMult('R', VDs2, hesLag[1]);
  VInv.MatMatMult('R', VDt2, hesLag[2]);
  VInv.MatMatMult('R', VDrs, hesLag[3]);
  VInv.MatMatMult('R', VDrt, hesLag[4]);
  VInv.MatMatMult('R', VDst, hesLag[5]);

  /*--- Check if the sum of the elements of the relevant rows of hesLag is 0. ---*/
  CheckRowSum(nIntTot, rDOFs.size(), 0.0, hesLag[0]);
  CheckRowSum(nIntTot, rDOFs.size(), 0.0, hesLag[1]);
  CheckRowSum(nIntTot, rDOFs.size(), 0.0, hesLag[2]);
  CheckRowSum(nIntTot, rDOFs.size(), 0.0, hesLag[3]);
  CheckRowSum(nIntTot, rDOFs.size(), 0.0, hesLag[4]);
  CheckRowSum(nIntTot, rDOFs.size(), 0.0, hesLag[5]);
}

void CFEMStandardPyraBase::LagBasisIntPointsPyra(const unsigned short          mPoly,
                                                 const vector<passivedouble>   &rDOFs,
                                                 const vector<passivedouble>   &sDOFs,
                                                 const vector<passivedouble>   &tDOFs,
                                                 const vector<passivedouble>   &rInt,
                                                 const vector<passivedouble>   &sInt,
                                                 const vector<passivedouble>   &tInt,
                                                 ColMajorMatrix<passivedouble> &lag) {

  /*--- Determine the padded number of the total number of integration points. ---*/
  const unsigned short nIntTot    = rInt.size();
  const unsigned short nIntTotPad = PaddedValue(nIntTot);

  /*--- Determine the inverse of the Vandermonde matrix of the DOFs. ---*/
  CSquareMatrixCM VInv(rDOFs.size());
  VandermondePyramid(mPoly, rDOFs, sDOFs, tDOFs, VInv.GetMat());
  VInv.Invert();

  /*--- Determine the Vandermonde matrix of the integration points. Make sure to
        allocate the number of rows to nIntTotPad and initialize them to zero. ---*/ 
  ColMajorMatrix<passivedouble> V(nIntTotPad,rDOFs.size());
  V.setConstant(0.0);
  VandermondePyramid(mPoly, rInt, sInt, tInt, V);

  /*--- The Lagrangian basis functions can be obtained by multiplying
        V and VInv. ---*/
  VInv.MatMatMult('R', V, lag);

  /*--- Check if the sum of the elements of the relevant rows of lag is 1. ---*/
  CheckRowSum(nIntTot, rDOFs.size(), 1.0, lag);
}

void CFEMStandardPyraBase::GradVandermondePyramid(const unsigned short          mPoly,
                                                  const vector<passivedouble>   &r,
                                                  const vector<passivedouble>   &s,
                                                  const vector<passivedouble>   &t,
                                                  ColMajorMatrix<passivedouble> &VDr,
                                                  ColMajorMatrix<passivedouble> &VDs,
                                                  ColMajorMatrix<passivedouble> &VDt) {

  /*--- For a pyramid the orthogonal basis for the reference element is
        obtained by a combination of Jacobi polynomials (of which the Legendre
        polynomials is a special case). This is the result of the
        orthonormalization of the monomial basis.
        Note that the sequence of the i, j and k loop must be identical to
        the evaluation of the Vandermonde matrix itself.  ---*/
  unsigned short ii = 0;
  for(unsigned short i=0; i<=mPoly; ++i) {
    for(unsigned short j=0; j<=mPoly; ++j) {
      unsigned short muij = max(i,j);
      for(unsigned short k=0; k<=(mPoly-muij); ++k, ++ii) {
        for(unsigned short l=0; l<r.size(); ++l) {

          /*--- Determine the coefficients a, b, c and d. ---*/
          passivedouble a, b;
          const passivedouble d = 1.0-t[l];
          if(fabs(d) < 1.e-8) a = b = 0.0;
          else {
            a = 2.0*r[l]/d;
            b = 2.0*s[l]/d;
          }

          const passivedouble c = t[l];

          /*--- Determine the value of the three 1D contributions to the 
                3D basis functions as well as the gradients of these
                contributions functions w.r.t. to their arguments. ---*/
          const passivedouble fa  = NormJacobi(i,0,         0,a);
          const passivedouble gb  = NormJacobi(j,0,         0,b);
          const passivedouble hc  = NormJacobi(k,2*(muij+1),0,c);
          const passivedouble dfa = GradNormJacobi(i,0,         0,a);
          const passivedouble dgb = GradNormJacobi(j,0,         0,b);
          const passivedouble dhc = GradNormJacobi(k,2*(muij+1),0,c);

          /*--- Computation of the powers of d that occur in the expressions for
                the gradients. Note the safeguard to avoid division by zero. This is
                allowed, because the implementation is such that when this clipping
                is active, it is multiplied by zero. ---*/
          const passivedouble tmpdmu   = muij > 0 ? pow(d, muij)   : 1.0;
          const passivedouble tmpdmum1 = muij > 1 ? pow(d, muij-1) : 1.0;

          /*--- Compute the derivatives of the basis function. ---*/
          VDr(l,ii) = 2.0*tmpdmum1* 2.0*dfa*gb*hc;
          VDs(l,ii) = 2.0*tmpdmum1* 2.0*fa*dgb*hc;
          VDt(l,ii) = 2.0*tmpdmum1*(a*dfa*gb*hc + b*fa*dgb*hc - muij*fa*gb*hc)
                    + 2.0*tmpdmu  * fa*gb*dhc;
        }
      }
    }
  }
}

void CFEMStandardPyraBase::HesVandermondePyramid(const unsigned short          mPoly,
                                                 const vector<passivedouble>   &r,
                                                 const vector<passivedouble>   &s,
                                                 const vector<passivedouble>   &t,
                                                 ColMajorMatrix<passivedouble> &VDr2,
                                                 ColMajorMatrix<passivedouble> &VDs2,
                                                 ColMajorMatrix<passivedouble> &VDt2,
                                                 ColMajorMatrix<passivedouble> &VDrs,
                                                 ColMajorMatrix<passivedouble> &VDrt,
                                                 ColMajorMatrix<passivedouble> &VDst) {

  /*--- For a pyramid the orthogonal basis for the reference element is
        obtained by a combination of Jacobi polynomials (of which the Legendre
        polynomials is a special case). This is the result of the
        orthonormalization of the monomial basis.
        Note that the sequence of the i, j and k loop must be identical to
        the evaluation of the Vandermonde matrix itself.  ---*/
  unsigned short ii = 0;
  for(unsigned short i=0; i<=mPoly; ++i) {
    for(unsigned short j=0; j<=mPoly; ++j) {
      unsigned short muij = max(i,j);
      for(unsigned short k=0; k<=(mPoly-muij); ++k, ++ii) {
        for(unsigned short l=0; l<r.size(); ++l) {

          /*--- Determine the coefficients a, b, c and d. ---*/
          const passivedouble d = 1.0-t[l];
          const passivedouble dInv = (fabs(d) < 1.e-8) ? 1.e+8 : 1.0/d; 
          const passivedouble a = 2.0*r[l]*dInv;
          const passivedouble b = 2.0*s[l]*dInv;
          const passivedouble c = t[l];

          /*--- Determine the value of the three 1D contributions to the 
                3D basis functions as well as the 1st and 2nd derivative of
                these contributions functions w.r.t. to their arguments. ---*/
          const passivedouble fa   = NormJacobi(i,0,         0,a);
          const passivedouble gb   = NormJacobi(j,0,         0,b);
          const passivedouble hc   = NormJacobi(k,2*(muij+1),0,c);
          const passivedouble dfa  = GradNormJacobi(i,0,         0,a);
          const passivedouble dgb  = GradNormJacobi(j,0,         0,b);
          const passivedouble dhc  = GradNormJacobi(k,2*(muij+1),0,c);
          const passivedouble d2fa = HesNormJacobi(i,0,         0,a);
          const passivedouble d2gb = HesNormJacobi(j,0,         0,b);
          const passivedouble d2hc = HesNormJacobi(k,2*(muij+1),0,c);

          /*--- Computation of the powers of d that occur in the expressions for
                the gradients. Note the safeguard to avoid division by zero. This is
                allowed, because the implementation is such that when this clipping
                is active, it is multiplied by zero. However, for the pyramid with
                its rational basis functions there are some terms in which a division
                by d cannot be avoided. This is handled via dInv. ---*/
          const passivedouble tmpdmu   = muij > 0 ? pow(d, muij)   : 1.0;
          const passivedouble tmpdmum1 = muij > 1 ? pow(d, muij-1) : 1.0;
          const passivedouble tmpdmum2 = muij > 2 ? pow(d, muij-2) : 1.0;

          /*--- Compute the 2nd derivative w.r.t. to r. ---*/
          VDr2(l,ii) = tmpdmum2*4.0*d2fa*gb*hc;

          /*--- Compute the 2nd derivative w.r.t. to s. ---*/
          VDs2(l,ii) = tmpdmum2*4.0*fa*d2gb*hc;

          /*--- Compute the 2nd derivative w.r.t. to r. ---*/
          VDt2(l,ii) = tmpdmum2*(muij*(muij-1)*fa*gb*hc
                     -           2.0*a*(muij-1)*dfa*gb*hc
                     -           2.0*b*(muij-1)*fa*dgb*hc
                     +           a*a*d2fa*gb*hc
                     +           b*b*fa*d2gb*hc)
                     + tmpdmum1*(2.0*a*dfa*gb*dhc
                     +           2.0*b*fa*dgb*dhc
                     -           2.0*muij*fa*gb*dhc)
                     + tmpdmu  *fa*gb*d2hc
                     + tmpdmum1*2.0*a*b*dfa*dgb*hc*dInv;

          /*--- Compute the cross derivative w.r.t. r and s. ---*/
          VDrs(l,ii) = tmpdmum1*4.0*dfa*dgb*hc*dInv;

          /*--- Compute the cross derivative w.r.t. r and t. ---*/
          VDrt(l,ii) = tmpdmum2*2.0*(a*d2fa*gb*hc - (muij-1)*dfa*gb*hc)
                     + tmpdmum1*2.0*(dfa*gb*dhc + b*dfa*dgb*hc*dInv);

          /*--- Compute the cross derivative w.r.t. s and t. ---*/
          VDst(l,ii) = tmpdmum2*2.0*(b*fa*d2gb*hc - (muij-1)*fa*dgb*hc)
                     + tmpdmum1*2.0*(fa*dgb*dhc + a*dfa*dgb*hc*dInv);

          /*--- Multiply all 2nd derivatives with 2 to obtain the
                correct expressions. ---*/
          VDr2(l,ii) *= 2.0;
          VDs2(l,ii) *= 2.0;
          VDt2(l,ii) *= 2.0;
          VDrs(l,ii) *= 2.0;
          VDrt(l,ii) *= 2.0;
          VDst(l,ii) *= 2.0;
        }
      }
    }
  }
}

void CFEMStandardPyraBase::VandermondePyramid(const unsigned short          mPoly,
                                              const vector<passivedouble>   &r,
                                              const vector<passivedouble>   &s,
                                              const vector<passivedouble>   &t,
                                              ColMajorMatrix<passivedouble> &V) {

  /*--- For a pyramid the orthogonal basis for the reference element is
        obtained by a combination of Jacobi polynomials (of which the Legendre
        polynomials is a special case). This is the result of the
        orthonormalization of the monomial basis. ---*/
  unsigned short ii = 0;
  for(unsigned short i=0; i<=mPoly; ++i) {
    for(unsigned short j=0; j<=mPoly; ++j) {
      unsigned short muij = max(i,j);
      for(unsigned short k=0; k<=(mPoly-muij); ++k, ++ii) {
        for(unsigned short l=0; l<r.size(); ++l) {

          /*--- Determine the coefficients a, b, c and d. ---*/
          passivedouble a, b;
          const passivedouble d = 1.0-t[l];
          if(fabs(d) < 1.e-8) a = b = 0.0;
          else {
            a = 2.0*r[l]/d;
            b = 2.0*s[l]/d;
          }

          const passivedouble c = t[l];

          /*--- Determine the value of the current basis function in this point. ---*/
          V(l,ii) = 2.0*pow(d,muij)*NormJacobi(i,0,0,a)*NormJacobi(j,0,0,b)
                  * NormJacobi(k,2*(muij+1),0,c);
        }
      }
    }
  }
}


void CFEMStandardPyraBase::LocationAllIntegrationPoints(vector<passivedouble> &rInt,
                                                        vector<passivedouble> &sInt,
                                                        vector<passivedouble> &tInt) {

  /*--- Determine the number of 1D integration points for GL and GJ. ---*/
  const unsigned short nGL = rLineIntGL.size();
  const unsigned short nGJ = rLineIntGJ.size();

  /*--- Determine the total number of integration points and determine the
        parametric coordinates of all integration points. ---*/
  const unsigned short nIntTot = nGL*nGL*nGJ;
  rInt.resize(nIntTot);
  sInt.resize(nIntTot);
  tInt.resize(nIntTot);

  unsigned short ii = 0;
  for(unsigned short k=0; k<nGJ; ++k) {
    const passivedouble zeta = rLineIntGJ[k];
    for(unsigned short j=0; j<nGL; ++j) {
      const passivedouble eta = rLineIntGL[j];
      for(unsigned short i=0; i<nGL; ++i, ++ii) {
        const passivedouble xi = rLineIntGL[i];

        rInt[ii] = 0.5*(1.0-zeta)*xi;
        sInt[ii] = 0.5*(1.0-zeta)*eta;
        tInt[ii] = zeta;
      }
    }
  }
}

void CFEMStandardPyraBase::LocationPyramidGridDOFsEquidistant(const unsigned short  mPoly,
                                                              vector<passivedouble> &rDOFs,
                                                              vector<passivedouble> &sDOFs,
                                                              vector<passivedouble> &tDOFs) {

  /*--- Determine the number of DOFs based on the polynomial degree. ---*/
  const unsigned short nD1D = mPoly + 1;
  const unsigned short nD   = nD1D*(nD1D+1)*(2*nD1D+1)/6;

  /*--- As the number of grid DOFs near the bottom is bigger
        than near the top of the pyramid, it is not possible
        to use a tensor product. The most convenient way is
        therefore to store all the grid DOFs explicitly.
        This is not a big issue, because the grid DOFs are
        only used in the pre- and post-processing steps.
        Allocate the memory for the parametric coordinates. ---*/
  rDOFs.resize(nD);
  sDOFs.resize(nD);
  tDOFs.resize(nD);

  /*--- Determine the spacing in t-direction, which is from the base quad
        to the top, and initialize the local counters mPoly and ii. ---*/
  const passivedouble dt = 2.0/mPoly;
  unsigned short nnPoly = mPoly, ii = 0;

  /*--- Outer loop in k-direction, which is from the base to top. ---*/
  for(unsigned int k=0; k<=mPoly; ++k, --nnPoly) {

    /*--- Determine the minimum and maximum value for r and s
          for this t-value. ---*/
    const passivedouble t     = -1.0 + k*dt;
    const passivedouble rsMin =  0.5*(t-1.0);
    const passivedouble rsMax = -rsMin;

    /*--- Determine the step size along the edges of the current quad.
          Take the exceptional situation nnPoly == 0 into account to avoid a
          division by zero.     ---*/
    const passivedouble dh = nnPoly ? (rsMax-rsMin)/nnPoly : 0.0;

    /*--- Loop over the vertices of the current quadrilateral. ---*/
    for(unsigned short j=0; j<=nnPoly; ++j) {
      const passivedouble s = rsMin + j*dh;

      for(unsigned short i=0; i<=nnPoly; ++i, ++ii) {
        const double r = rsMin + i*dh;
        rDOFs[ii] = r;
        sDOFs[ii] = s;
        tDOFs[ii] = t;
      }
    }
  }
}

void CFEMStandardPyraBase::LocationPyramidGridDOFsLGL(const unsigned short  mPoly,
                                                      vector<passivedouble> &rDOFs,
                                                      vector<passivedouble> &sDOFs,
                                                      vector<passivedouble> &tDOFs) {

  /*--- Determine the number of DOFs based on the polynomial degree. ---*/
  const unsigned short nD1D = mPoly + 1;
  const unsigned short nD   = nD1D*(nD1D+1)*(2*nD1D+1)/6;

  /*--- As the number of grid DOFs near the bottom is bigger
        than near the top of the pyramid, it is not possible
        to use a tensor product. The most convenient way is
        therefore to store all the grid DOFs explicitly.
        This is not a big issue, because the grid DOFs are
        only used in the pre- and post-processing steps.
        Allocate the memory for the parametric coordinates. ---*/
  rDOFs.resize(nD);
  sDOFs.resize(nD);
  tDOFs.resize(nD);

  /*--- Determine the location of the 1D LGL points along the standard line. ---*/
  vector<passivedouble> LGLPoints;
  Location1DGridDOFsLGL(mPoly, LGLPoints);

  /*--- Determine the location of the DOFs of the standard triangle. ---*/
  vector<passivedouble> xiTri, etaTri;
  LocationTriangleGridDOFsLGL(mPoly, xiTri, etaTri);

  /*--- Determine the location of the DOFs on the triangle 0-1-4
        of the pyramid. This is stored in a 2D vector for convenience.
        The mapping from the standard triangle to the triangle of
        the pyramid is given as: r =  0.5 + xi + 0.5*eta,
                                 s = -0.5      + 0.5*eta,
                                 t =                 eta,
        where (r,s,t) are the parametric coordinates of the pyramid
        and (xi, eta) are the parametric coordinates of the triangle.
        The (r,s,t) coordinates of the other faces can be obtained by
        swapping and negating the r and s coordinates and it is therefore
        not needed to compute those. ---*/
  vector<vector<passivedouble> > rDOFs_Face014, sDOFs_Face014, tDOFs_Face014;
  rDOFs_Face014.resize(mPoly+1);
  sDOFs_Face014.resize(mPoly+1);
  tDOFs_Face014.resize(mPoly+1);

  unsigned short ii = 0;
  for(unsigned short j=0; j<=mPoly; ++j) {
    for(unsigned short i=0; i<=(mPoly-j); ++i, ++ii) {
      rDOFs_Face014[j].push_back( 0.5 + xiTri[ii] + 0.5*etaTri[ii]);
      sDOFs_Face014[j].push_back(-0.5             + 0.5*etaTri[ii]);
      tDOFs_Face014[j].push_back(etaTri[ii]);
    }
  }

  /*--- Set the coordinates of the base quad. ---*/
  ii = 0;
  for(unsigned short j=0; j<=mPoly; ++j) {
    for(unsigned short i=0; i<=mPoly; ++i, ++ii) {
      rDOFs[ii] =  LGLPoints[i];
      sDOFs[ii] =  LGLPoints[j];
      tDOFs[ii] = -1.0;
    }
  }

  /*--- Loop over the interior planes of the pyramid. ---*/
  for(unsigned short k=1; k<mPoly; ++k) {

    /*--- Determine the approximte edge length of the curve along the triangular
          face for this k-value. Store the length values of the DOFs on this
          line, as they are converted later on to parametric values. ---*/
    vector<passivedouble> xi(rDOFs_Face014[k].size(), 0.0);
    for(unsigned short i=1; i<rDOFs_Face014[k].size(); ++i) {
      const passivedouble dr = rDOFs_Face014[k][i] - rDOFs_Face014[k][i-1];
      const passivedouble ds = sDOFs_Face014[k][i] - sDOFs_Face014[k][i-1];
      const passivedouble dt = tDOFs_Face014[k][i] - tDOFs_Face014[k][i-1];

      xi[i] = xi[i-1] + sqrt(dr*dr + ds*ds + dt*dt);
    }

    /*--- Convert the xi values to a scaled version between zero and one. ---*/
    for(unsigned short i=0; i<xi.size(); ++i) xi[i] /= xi.back();

    /*--- Store the coordinates of the corner points of the quad. As symmetry
          can be used, it is only necessary to store the coordinate of the
          lower left node. ---*/
    const passivedouble rLowerLeft = rDOFs_Face014[k][0];
    const passivedouble sLowerLeft = sDOFs_Face014[k][0];
    const passivedouble tLowerLeft = tDOFs_Face014[k][0];

    /*--- Carry out the transfinite interpolation to determine the coordinates
          of the current quad. First loop over the parametric v-direction. ---*/
    for(unsigned short j=0; j<rDOFs_Face014[k].size(); ++j) {

      /*--- Store the parametric v-coordinate and the (r,s,t) coordinates of
            the point, which corresponds to this v value and u = 0. The
            corresponding coordinate for u = 1 has the same s and t,
            but negative r. So no need to store it. ---*/
      const passivedouble v    = xi[j];
      const passivedouble omv  = 1.0 - v;
      const passivedouble r_c2 = sDOFs_Face014[k][j];
      const passivedouble s_c2 = rDOFs_Face014[k][j];
      const passivedouble t_c2 = tDOFs_Face014[k][j];

      /*--- Loop over the parametric u-direction. ---*/
      for(unsigned short i=0; i<rDOFs_Face014[k].size(); ++i, ++ii) {

        /*--- Store the parametric u-coordinate and the (r,s,t) coordinates of
              the point, which corresponds to this u value and v = 0. The
              corresponding coordinate for v = 1 has the same r and t,
              but negative s. So no need to store it. ---*/
        const passivedouble u    = xi[i];
        const passivedouble omu  = 1.0 - u;
        const passivedouble r_c1 = rDOFs_Face014[k][i];
        const passivedouble s_c1 = sDOFs_Face014[k][i];
        const passivedouble t_c1 = tDOFs_Face014[k][i];

        /*--- Interpolate the coordinate inside the pyramid. ---*/
        rDOFs[ii] = r_c1*(omv + v) + r_c2*(omu - u)
                  - rLowerLeft*(omu*omv - u*v - u*omv + omu*v);
        sDOFs[ii] = s_c1*(omv - v) + s_c2*(omu + u)
                  - sLowerLeft*(omu*omv - u*v + u*omv - omu*v);
        tDOFs[ii] = t_c1*(omv + v) + t_c2*(omu + u)
                  - tLowerLeft*(omu*omv + u*v + u*omv + omu*v);
      }
    }
  }

  /*--- Set the last coordinate to the top of the pyramid. ---*/
  rDOFs[ii] = 0.0;
  sDOFs[ii] = 0.0;
  tDOFs[ii] = 1.0;  
}

void CFEMStandardPyraBase::SubConnLinearElements(void) {

  /*--- The pyramid is split into several linear pyramids and tets.
        Set the VTK sub-types accordingly. ---*/
  VTK_SubType1 = PYRAMID;
  VTK_SubType2 = TETRAHEDRON;

  /*--- Initialize the number of DOFs for the current edges to the number of
        DOFs of the edges on the base of the pyramid. Also initialize the
        current k offset to zero.     ---*/
  unsigned short nDOFsCurrentEdges = nPoly + 1;
  unsigned short offCurrentK       = 0;

  /*--- Loop in the k-direction of the pyramid. ---*/
  for(unsigned short k=0; k<nPoly; ++k) {

    /*------------------------------------------------------------------------*/
    /*       Sub-pyramids in the same direction as the original pyramid.      */
    /*------------------------------------------------------------------------*/

    /*--- Determine the index of the first vertex of the quadrilateral of the
          next k value. ---*/
    unsigned short kk = offCurrentK + nDOFsCurrentEdges*nDOFsCurrentEdges;

    /*--- Loop in j-direction of the current quad. ---*/
    for(unsigned short j=0; j<(nDOFsCurrentEdges-1); ++j) {

      /*--- Index of the first vertex along the j-row of the current quad. ---*/
      const unsigned short jj = offCurrentK + j*nDOFsCurrentEdges;

      /*--- Loop in i-direction of the current quad. ---*/
      for(unsigned short i=0; i<(nDOFsCurrentEdges-1); ++i) {

        /*--- Determine the local indices of the corners of the quadrilateral
              of this subpyramid as well as the top of the subpyramid.
              Store the connectivity in subConn1ForPlotting.  ---*/
        const unsigned short n0 = jj + i;
        const unsigned short n1 = n0 + 1;
        const unsigned short n2 = n1 + nDOFsCurrentEdges;
        const unsigned short n3 = n0 + nDOFsCurrentEdges;
        const unsigned short n4 = kk + i;

        subConn1ForPlotting.push_back(n0);
        subConn1ForPlotting.push_back(n1);
        subConn1ForPlotting.push_back(n2);
        subConn1ForPlotting.push_back(n3);
        subConn1ForPlotting.push_back(n4);
      }

      /*--- Update kk for the next j-row. ---*/
      kk += nDOFsCurrentEdges - 1;
    }

    /*------------------------------------------------------------------------*/
    /*    Sub-pyramids in the opposite direction as the original pyramid.     */
    /*------------------------------------------------------------------------*/

    /*--- Reset the value of kk to the index of the first vertex of the
          quadrilateral of the next k-plane. ---*/
    kk = offCurrentK + nDOFsCurrentEdges*nDOFsCurrentEdges;

    /*--- Loop in j-direction of the current quad. Note that the starting index
          of this loop is 1. ---*/
    for(unsigned short j=1; j<(nDOFsCurrentEdges-1); ++j) {

      /*--- Index of the first vertex along the j-row of the current quad. ---*/
      const unsigned short jj = offCurrentK + j*nDOFsCurrentEdges;

      /*--- Loop in the i-direction of this quad. Again the starting index is 1. ---*/
      for(unsigned short i=1; i<(nDOFsCurrentEdges-1); ++i) {

        /*--- Determine the local indices of the corners of the quadrilateral
              of this subpyramid as well as the top of the subpyramid.  ---*/
        const unsigned short n0 = kk + i - 1;
        const unsigned short n1 = n0 + 1;
        const unsigned short n2 = n1 + nDOFsCurrentEdges-1;
        const unsigned short n3 = n0 + nDOFsCurrentEdges-1;
        const unsigned short n4 = jj + i;

        /*--- Store the connectivity of this subpyramid. Note that n1 and n3 are
              swapped, such that a positive volume is obtained according to the
              right hand rule. ---*/
        subConn1ForPlotting.push_back(n0);
        subConn1ForPlotting.push_back(n3);
        subConn1ForPlotting.push_back(n2);
        subConn1ForPlotting.push_back(n1);
        subConn1ForPlotting.push_back(n4);
      }

      /*--- Update kk for the next j-row. ---*/
      kk += nDOFsCurrentEdges - 1;
    }

    /*------------------------------------------------------------------------*/
    /*                   Sub-tetrahedra in the j-direction.                   */
    /*------------------------------------------------------------------------*/

    /*--- Reset the value of kk again. ---*/
    kk = offCurrentK + nDOFsCurrentEdges*nDOFsCurrentEdges;

    /*--- Loop in the i-direction of the current quad. Note that the starting
          index must be 1. ---*/
    for(unsigned short i=1; i<(nDOFsCurrentEdges-1); ++i) {

      /*--- Loop in the j-direction of the current quad. This loop starts at 0. ---*/
      for(unsigned short j=0; j<(nDOFsCurrentEdges-1); ++j) {

        /*--- Determine the local indices of the 4 corner points of this tet.
              Its connectivity is stored in subConn2ForPlotting, because in
              subConn1ForPlotting the pyramids are stored. ---*/
        const unsigned short n0 = kk + i-1 + j*(nDOFsCurrentEdges-1);
        const unsigned short n1 = n0 + 1;
        const unsigned short n2 = offCurrentK + i + j*nDOFsCurrentEdges;
        const unsigned short n3 = n2 + nDOFsCurrentEdges;

        subConn2ForPlotting.push_back(n0);
        subConn2ForPlotting.push_back(n1);
        subConn2ForPlotting.push_back(n2);
        subConn2ForPlotting.push_back(n3);
      }
    }

    /*------------------------------------------------------------------------*/
    /*                Sub-tetrahedra in the i-direction.                      */
    /*------------------------------------------------------------------------*/

    /*--- Loop in the j-direction of the current quad. Note that the starting
          index must be 1. ---*/
    for(unsigned short j=1; j<(nDOFsCurrentEdges-1); ++j) {

      /*--- Index of the first vertex along the j-row of the current quad. ---*/
      const unsigned short jj = offCurrentK + j*nDOFsCurrentEdges;

      /*--- Loop in the i-direction of the current quad. This loop starts at 0. ---*/
      for(unsigned short i=0; i<(nDOFsCurrentEdges-1); ++i) {

        /*--- Determine the local indices of the 4 corner points of this tet
              and store its connectivity in subConn2ForPlotting.   ---*/
        const unsigned short n0 = kk + i;
        const unsigned short n1 = jj + i;
        const unsigned short n2 = n0 + nDOFsCurrentEdges-1;
        const unsigned short n3 = n1 + 1;

        subConn2ForPlotting.push_back(n0);
        subConn2ForPlotting.push_back(n1);
        subConn2ForPlotting.push_back(n2);
        subConn2ForPlotting.push_back(n3);
      }

      /*--- Update kk for the next j-row. ---*/
      kk += nDOFsCurrentEdges - 1;
    }

    /*--- Update the value of offCurrentK with the amounts of DOFs present in the
          current quadrilateral plane and decrement nDOFsCurrentEdges, such that it
          contains the number of DOFs along an edge of the next quadrilateral plane. ---*/
    offCurrentK += nDOFsCurrentEdges*nDOFsCurrentEdges;
    --nDOFsCurrentEdges;
  }
}

void CFEMStandardPyraBase::LocalGridConnFaces(void) {

  /*--- Allocate the first index of gridConnFaces, which is equal to the number
        of faces of the pyramid, which is 5. Reserve memory for the second
        index afterwards. ---*/
  const unsigned short nDOFsQuad     = (nPoly+1)*(nPoly+1);
  const unsigned short nDOFsTriangle = (nPoly+1)*(nPoly+2)/2;
  gridConnFaces.resize(5);

  gridConnFaces[0].reserve(nDOFsQuad);
  gridConnFaces[1].reserve(nDOFsTriangle);
  gridConnFaces[2].reserve(nDOFsTriangle);
  gridConnFaces[3].reserve(nDOFsTriangle);
  gridConnFaces[4].reserve(nDOFsTriangle);

  /*--- Loop over all the nodes of the pyramid and pick the correct
        ones for the faces. ---*/
  unsigned short mPoly = nPoly;
  unsigned int ii = 0;
  for(unsigned short k=0; k<=nPoly; ++k, --mPoly) {
    for(unsigned short j=0; j<=mPoly; ++j) {
      for(unsigned short i=0; i<=mPoly; ++i, ++ii) {
        if(k == 0)     gridConnFaces[0].push_back(ii);
        if(j == 0)     gridConnFaces[1].push_back(ii);
        if(j == mPoly) gridConnFaces[2].push_back(ii);
        if(i == 0)     gridConnFaces[3].push_back(ii);
        if(i == mPoly) gridConnFaces[4].push_back(ii);
      }
    }
  }

  /*--- Make sure that the element is to the left of the faces. ---*/
  const unsigned short n0 = 0;
  const unsigned short n1 = nPoly;
  const unsigned short n2 = nDOFsQuad -1;
  const unsigned short n3 = n2 - nPoly;
  const unsigned short n4 = nDOFs -1;

  ChangeDirectionQuadConn(gridConnFaces[0], n0, n1, n2, n3);
  ChangeDirectionTriangleConn(gridConnFaces[1], n0, n4, n1);
  ChangeDirectionTriangleConn(gridConnFaces[2], n3, n2, n4);
  ChangeDirectionTriangleConn(gridConnFaces[3], n0, n3, n4);
  ChangeDirectionTriangleConn(gridConnFaces[4], n1, n4, n2);
}