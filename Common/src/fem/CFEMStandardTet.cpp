/*!
 * \file CFEMStandardTet.cpp
 * \brief Functions for the class CFEMStandardTet.
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

#include "../../include/fem/CFEMStandardTet.hpp"
#include "../../include/toolboxes/CGeneralSquareMatrixCM.hpp"

/*----------------------------------------------------------------------------------*/
/*               Public member functions of CFEMStandardTet.                        */
/*----------------------------------------------------------------------------------*/

CFEMStandardTet::CFEMStandardTet(const unsigned short val_nPoly,
                                 const unsigned short val_orderExact) {

  /*--- Store the command line arguments. ---*/
  VTK_Type   = TETRAHEDRON;
  nPoly      = val_nPoly;
  orderExact = val_orderExact;

  /*--- Determine the total number of DOFs and its padded version. ---*/
  nDOFs    = (nPoly+1)*(nPoly+2)*(nPoly+3)/6;
  nDOFsPad = ((nDOFs+baseVectorLen-1)/baseVectorLen)*baseVectorLen;

  /*--- Determine the parametric locations of the grid DOFs of the tetrahedron. ---*/
  LocationTetGridDOFsEquidistant(rTetDOFsEqui, sTetDOFsEqui, tTetDOFsEqui);
  LocationTetGridDOFsLGL();

  /*--- Determine the parametric location and weights of the
        integration rule of the tetrahedron. ---*/
  IntegrationPointsTetrahedron(rTetInt, sTetInt, tTetInt, wTetInt);

  /*--- Determine the total number of integration points
        and its padded version. ---*/
  nIntegration    = rTetInt.size();
  nIntegrationPad = ((nIntegration+baseVectorLen-1)/baseVectorLen)*baseVectorLen;

  /*--- Allocate the memory for the padded number of integration points and
        initialize the weights to zero. This is done such that the padded
        values are initialized appropriately. ---*/
  wIntegration.resize(nIntegrationPad);
  wIntegration.setConstant(0.0);

  /*--- Copy the values from wTetInt to wIntegration. ---*/
  for(unsigned short i=0; i<nIntegration; ++i)
    wIntegration(i) = wTetInt[i];
}

void CFEMStandardTet::LocationTetGridDOFsEquidistant(vector<passivedouble> &r,
                                                     vector<passivedouble> &s,
                                                     vector<passivedouble> &t) {

  /*--- For a tetrahedron it is not possible to apply
        a tensor product and therefore all DOFs are
        simply stored. Allocate the memory. ---*/
  r.resize(nDOFs);
  s.resize(nDOFs);
  t.resize(nDOFs);

  /*--- Determine the equidistant spacing along an edge. ---*/
  const passivedouble dh = 2.0/nPoly;

  /*--- Triple loop to compute the location of the grid DOFs. ---*/
  unsigned short ii = 0;
  for(unsigned short k=0; k<=nPoly; ++k) {
    const unsigned short uppBoundJ = nPoly - k;
    for(unsigned short j=0; j<=uppBoundJ; ++j) {
      const unsigned short uppBoundI = nPoly - k - j;
      for(unsigned short i=0; i<=uppBoundI; ++i, ++ii) {
        r[ii] = -1.0 + i*dh;
        s[ii] = -1.0 + j*dh;
        t[ii] = -1.0 + k*dh;
      }
    }
  }
}

void CFEMStandardTet::LocationTetGridDOFsLGL() {

  /*--- The code to determine the parametric coordinates of the DOFs of the
        tetrahedron is a translation of the Matlab code belonging to the book
        Nodal Discontinuous Galerkin Methods, Algorithms, Analysis and Applications,
        written by Jan S. Hesthaven and Tim Warburton. ---*/

  /*--- Local parameters. ---*/
  const passivedouble tol   = 1.e-10;
  const passivedouble sqrt3 = sqrt(3.0);
  const passivedouble sqrt6 = sqrt(6.0);
  const passivedouble alphaOpt[] = {0.0000, 0.0000, 0.0000, 0.0000, 0.1002, 1.1332,
                                    1.5608, 1.3413, 1.2577, 1.1603, 1.10153,
                                    0.6080, 0.4523, 0.8856, 0.8717, 0.9655};

  /*--- Create the equidistributed nodes. ---*/
  vector<passivedouble> r, s, t;
  LocationTetGridDOFsEquidistant(r, s, t);

  /*--- Create the barycentric coordinates. ---*/
  vector<passivedouble> L1(nDOFs), L2(nDOFs), L3(nDOFs), L4(nDOFs);

  for(unsigned short i=0; i<nDOFs; ++i) {
    L1[i] =  0.5*(1.0 + t[i]);
    L2[i] =  0.5*(1.0 + s[i]);
    L3[i] = -0.5*(1.0 + r[i] + s[i] + t[i]);
    L4[i] =  0.5*(1.0 + r[i]);
  }

  /*--- Set the coordinates of the vertices of the tetrahedron. ---*/
  const passivedouble v1[] = {-1.0, -1.0/sqrt3, -1.0/sqrt6};
  const passivedouble v2[] = { 1.0, -1.0/sqrt3, -1.0/sqrt6};
  const passivedouble v3[] = { 0.0,  2.0/sqrt3, -1.0/sqrt6};
  const passivedouble v4[] = { 0.0,  0.0,        3.0/sqrt6};

  /*--- Create for the four faces the orthogonal axis which
        are tangent to the face. ---*/
  passivedouble t1[4][3], t2[4][3];
  for(unsigned short i=0; i<3; ++i) {
    t1[0][i] = v2[i]-v1[i]; t2[0][i] = v3[i] - 0.5*(v1[i]+v2[i]);
    t1[1][i] = v2[i]-v1[i]; t2[1][i] = v4[i] - 0.5*(v1[i]+v2[i]);
    t1[2][i] = v3[i]-v2[i]; t2[2][i] = v4[i] - 0.5*(v2[i]+v3[i]);
    t1[3][i] = v3[i]-v1[i]; t2[3][i] = v4[i] - 0.5*(v1[i]+v3[i]);
  }

  /*--- Normalize the rows of t1 and t2. ---*/
  for(unsigned short i=0; i<4; ++i) {
    passivedouble norm;
    norm = sqrt(t1[i][0]*t1[i][0] + t1[i][1]*t1[i][1] + t1[i][2]*t1[i][2]);
    t1[i][0] /= norm; t1[i][1] /= norm; t1[i][2] /= norm;

    norm = sqrt(t2[i][0]*t2[i][0] + t2[i][1]*t2[i][1] + t2[i][2]*t2[i][2]);
    t2[i][0] /= norm; t2[i][1] /= norm; t2[i][2] /= norm;
  }

  /*--- Set the optimal value of alp for this case. ---*/
  const passivedouble alp = (nPoly < 16) ? alphaOpt[nPoly] : 1.0;

  /*--- Allocate the memory for the parametric coordinates. ---*/
  rTetDOFsLGL.resize(nDOFs);
  sTetDOFsLGL.resize(nDOFs);
  tTetDOFsLGL.resize(nDOFs);

  /*--- Initialize the parametric coordinates to the uniform distribution
        on the equilateral tetrahedron. ---*/
  for(unsigned short i=0; i<nDOFs; ++i) {
    rTetDOFsLGL[i] = L3[i]*v1[0] + L4[i]*v2[0] + L2[i]*v3[0] + L1[i]*v4[0];
    sTetDOFsLGL[i] = L3[i]*v1[1] + L4[i]*v2[1] + L2[i]*v3[1] + L1[i]*v4[1];
    tTetDOFsLGL[i] = L3[i]*v1[2] + L4[i]*v2[2] + L2[i]*v3[2] + L1[i]*v4[2];
  }

  /*--- Initialize the corrections, which are stored in r, s and t, to zero. ---*/
  r.assign(nDOFs, 0.0);
  s.assign(nDOFs, 0.0);
  t.assign(nDOFs, 0.0);

  /*--- Loop over the four faces of the tetrahedron. ---*/
  for(unsigned short face=0; face<4; ++face) {

    /*--- In order to use the same approach for each face, the barycentric
          coordinates are different for the faces. Store these La, Lb,
          Lc and Ld. ---*/
    vector<passivedouble> La, Lb, Lc, Ld;
    if(face == 0)      {La = L1; Lb = L2; Lc = L3; Ld = L4;}
    else if(face == 1) {La = L2; Lb = L1; Lc = L3; Ld = L4;}
    else if(face == 2) {La = L3; Lb = L1; Lc = L4; Ld = L2;}
    else               {La = L4; Lb = L1; Lc = L3; Ld = L2;}

    /*--- Compute the warping factors tangential to this face. ---*/
    vector<passivedouble> warp1, warp2;
    EvalShift(alp, Lb, Lc, Ld, warp1, warp2);

    /*--- Compute the volume blending. ---*/
    vector<passivedouble> blend(nDOFs);
    for(unsigned short i=0; i<nDOFs; ++i) {
      blend[i] = Lb[i]*Lc[i]*Ld[i];

      const passivedouble denom = (Lb[i]+0.5*La[i])*(Lc[i]+0.5*La[i])*(Lc[i]+0.5*La[i]);
      if(denom > tol)
        blend[i] *= (1.0 + alp*alp*La[i]*La[i])/denom;
    }

    /*--- Update the corrections. ---*/
    for(unsigned short i=0; i<nDOFs; ++i) {
      const passivedouble abv1 = blend[i]*warp1[i], abv2 = blend[i]*warp2[i];
      r[i] += abv1*t1[face][0] + abv2*t2[face][0];
      s[i] += abv1*t1[face][1] + abv2*t2[face][1];
      t[i] += abv1*t1[face][2] + abv2*t2[face][2];
    }

    /*--- Fix the face warping. ---*/
    for(unsigned short i=0; i<nDOFs; ++i) {
      unsigned short val = 0;
      if(Lb[i] > tol) ++val;
      if(Lc[i] > tol) ++val;
      if(Ld[i] > tol) ++val;
      if((La[i] < tol) && val < 3) {
        r[i] = warp1[i]*t1[face][0] + warp2[i]*t2[face][0];
        s[i] = warp1[i]*t1[face][1] + warp2[i]*t2[face][1];
        t[i] = warp1[i]*t1[face][2] + warp2[i]*t2[face][2];
      }
    }
  }

  /*--- Store the parametric coordinates of the equilateral tetrahedron in r,s,t. ---*/
  for(unsigned short i=0; i<nDOFs; ++i) {
    r[i] += rTetDOFsLGL[i];
    s[i] += sTetDOFsLGL[i];
    t[i] += tTetDOFsLGL[i];
  }

  /*--- Create the matrix needed to convert the coordinates to the standard tetrahedron. ---*/
  CGeneralSquareMatrixCM A(3);
  A(0,0) = v2[0]-v1[0]; A(1,0) = v2[1]-v1[1]; A(2,0) = v2[2]-v1[2];
  A(0,1) = v3[0]-v1[0]; A(1,1) = v3[1]-v1[1]; A(2,1) = v3[2]-v1[2];
  A(0,2) = v4[0]-v1[0]; A(1,2) = v4[1]-v1[1]; A(2,2) = v4[2]-v1[2];

  A.Invert();

  /*--- Convert the parametric coordinates from the equilateral tetrahedron to the
        parametric coordinates of the standard tetrahedron. ---*/
  for(unsigned short i=0; i<nDOFs; ++i) {
    const passivedouble rr = 2.0*r[i] + v1[0] - v2[0] - v3[0] - v4[0];
    const passivedouble ss = 2.0*s[i] + v1[1] - v2[1] - v3[1] - v4[1];
    const passivedouble tt = 2.0*t[i] + v1[2] - v2[2] - v3[2] - v4[2];

    rTetDOFsLGL[i] = A(0,0)*rr + A(0,1)*ss + A(0,2)*tt;
    sTetDOFsLGL[i] = A(1,0)*rr + A(1,1)*ss + A(1,2)*tt;
    tTetDOFsLGL[i] = A(2,0)*rr + A(2,1)*ss + A(2,2)*tt;
  }
}

void CFEMStandardTet::EvalShift(const passivedouble         alpha,
                                const vector<passivedouble> &L1,
                                const vector<passivedouble> &L2,
                                const vector<passivedouble> &L3,
                                vector<passivedouble>       &dx,
                                vector<passivedouble>       &dy) {

  /*--- The code for this function is a translation of the Matlab code belonging to
        the book Nodal Discontinuous Galerkin Methods, Algorithms, Analysis and
        Applications, written by Jan S. Hesthaven and Tim Warburton. ---*/

  /*--- Compute for each edge the nodal blending functions. ---*/
  vector<passivedouble> blend1(nDOFs), blend2(nDOFs), blend3(nDOFs);

  for(unsigned short i=0; i<nDOFs; ++i) {
    blend1[i] = L2[i]*L3[i];
    blend2[i] = L1[i]*L3[i];
    blend3[i] = L1[i]*L2[i];
  }

  /*--- Compute the warp factors for the 3 edges. ---*/
  vector<passivedouble> warp1, warp2, warp3;
  vector<passivedouble> tmp(nDOFs);

  for(unsigned short i=0; i<nDOFs; ++i) tmp[i] = L3[i] - L2[i];
  EvalWarp(tmp, warp1);

  for(unsigned short i=0; i<nDOFs; ++i) tmp[i] = L1[i] - L3[i];
  EvalWarp(tmp, warp2);

  for(unsigned short i=0; i<nDOFs; ++i) tmp[i] = L2[i] - L1[i];
  EvalWarp(tmp, warp3);

  /*--- Combine the blending and the warping. ---*/
  for(unsigned short i=0; i<nDOFs; ++i) {
    warp1[i] *= blend1[i]*(1.0 + alpha*alpha*L1[i]*L1[i]);
    warp2[i] *= blend2[i]*(1.0 + alpha*alpha*L2[i]*L2[i]);
    warp3[i] *= blend3[i]*(1.0 + alpha*alpha*L3[i]*L3[i]);
  }

  /*--- Evaluate the shift in the equilateral triangle, which is the
        final result to be stored in dx and dy. ---*/
  dx.resize(nDOFs);
  dy.resize(nDOFs);

  const passivedouble cos1 = cos(SU2_TYPE::GetValue(TWO3*PI_NUMBER));
  const passivedouble cos2 = cos(SU2_TYPE::GetValue(FOUR3*PI_NUMBER));
  const passivedouble sin1 = sin(SU2_TYPE::GetValue(TWO3*PI_NUMBER));
  const passivedouble sin2 = sin(SU2_TYPE::GetValue(FOUR3*PI_NUMBER));

  for(unsigned short i=0; i<nDOFs; ++i) {
    dx[i] = warp1[i] + cos1*warp2[i] + cos2*warp3[i];
    dy[i] =            sin1*warp2[i] + sin2*warp3[i];
  }
}

void CFEMStandardTet::EvalWarp(const vector<passivedouble> &xOut,
                               vector<passivedouble>       &warp) {

  /*--- The code for this function is a translation of the Matlab code belonging to
        the book Nodal Discontinuous Galerkin Methods, Algorithms, Analysis and
        Applications, written by Jan S. Hesthaven and Tim Warburton. ---*/

  /*--- Determine the 1D equidistributed and LGL nodes. Note that
        in this function the reversed definition is used. ---*/
  vector<passivedouble> xEq, xLGL;
  Location1DGridDOFsEquidistant(xEq);
  Location1DGridDOFsLGL(xLGL);

  reverse(xEq.begin(),  xEq.end());
  reverse(xLGL.begin(), xLGL.end());

  /*--- Initialize the warping vector to zero. ---*/
  warp.assign(nDOFs, 0.0);

  /*--- Loop over the internal points of the edge. ---*/
  for(unsigned short i=1; i<nPoly; ++i) {

    /*--- Initialize the vector d to the distance between the
          equidistributed and LGL node i. ---*/
    vector<passivedouble> d(nDOFs, xLGL[i]-xEq[i]);

    /*--- Loop over the other internal nodes of the edge to correct
          the entries of vector d. ---*/
    for(unsigned short j=1; j<nPoly; ++j) {
      if(j != i) {
        const passivedouble denom = 1.0/(xEq[i]-xEq[j]);
        for(unsigned short k=0; k<nDOFs; ++k)
          d[k] *= (xOut[k]-xEq[j])*denom; 
      }
    }

    /*--- Carry out the scaling relative to the end points of the edge. ---*/
    passivedouble denom = -1.0/(xEq[i]-xEq[0]);
    for(unsigned short k=0; k<nDOFs; ++k)
      d[k] *= denom;

    denom = 1.0/(xEq[i]-xEq[nPoly]);
    for(unsigned short k=0; k<nDOFs; ++k)
      d[k] *= denom;

    /*--- Add d to warp for the DOFs. Note that the factor 4 is included here instead
          of in EvalShift as in the Matlab code of Hesthaven and Warburton. ---*/
    for(unsigned short k=0; k<nDOFs; ++k)
      warp[k] += 4.0*d[k];
  }
}
