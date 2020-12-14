/*!
 * \file CFEMStandardPyra.cpp
 * \brief Functions for the class CFEMStandardPyra.
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

#include "../../include/fem/CFEMStandardPyra.hpp"
#include "../../include/fem/fem_gauss_jacobi_quadrature.hpp"

/*----------------------------------------------------------------------------------*/
/*              Public member functions of CFEMStandardPyra.                        */
/*----------------------------------------------------------------------------------*/

CFEMStandardPyra::CFEMStandardPyra(const unsigned short val_nPoly,
                                   const unsigned short val_orderExact) {

  /*--- Store the command line arguments. ---*/
  VTK_Type   = PYRAMID;
  nPoly      = val_nPoly;
  orderExact = val_orderExact;

  /*--- Determine the number of DOFs in 1D of the base quad and the total number of DOFs.
        Also determine the padded value of the latter. ---*/
  nDOFs1D  = nPoly + 1;
  nDOFs    = nDOFs1D*(nDOFs1D+1)*(2*nDOFs1D+1)/6;
  nDOFsPad = ((nDOFs+baseVectorLen-1)/baseVectorLen)*baseVectorLen;

  /*--- Determine the parametric locations of the grid DOFs of the pyramid. ---*/
  LocationPyramidGridDOFsEquidistant(rPyraDOFsEqui, sPyraDOFsEqui, tPyraDOFsEqui);
  LocationPyramidGridDOFsLGL(rPyraDOFsLGL, sPyraDOFsLGL, tPyraDOFsLGL);

  /*--- The 3D quadrature rule for a pyramid is obtained by transforming the
        standard pyramid into a standard hexahedron by means of the Duffy
        transformation. Hence the integration rule is a tensor product,
        albeit a special one. ---*/
  nInt1D          = orderExact/2 + 1;
  nIntegration    = nInt1D*nInt1D*nInt1D;
  nIntegrationPad = ((nIntegration+baseVectorLen-1)/baseVectorLen)*baseVectorLen;

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
        is 0.25*(1-t)^2. The Gauss-Jacobi rule only accounts for  (1-t)^2.---*/
  unsigned short ii = 0;
  for(unsigned short k=0; k<nInt1D; ++k)
    for(unsigned short j=0; j<nInt1D; ++j)
      for(unsigned short i=0; i<nInt1D; ++i, ++ii)
        wIntegration(ii) = 0.25*wLineIntGL[i]*wLineIntGL[j]*wLineIntGJ[k];
}

void CFEMStandardPyra::LocationAllIntegrationPoints(vector<passivedouble> &rInt,
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

void CFEMStandardPyra::LocationPyramidGridDOFsEquidistant(vector<passivedouble> &rDOFs,
                                                          vector<passivedouble> &sDOFs,
                                                          vector<passivedouble> &tDOFs) {

  /*--- Determine the number of DOFs based on the polynomial degree.
        This is done, because this function can be called for different
        values of nPoly. ---*/
  const unsigned short nD1D  = nPoly + 1;
  const unsigned short nD    = nD1D*(nD1D+1)*(2*nD1D+1)/6;

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
  const passivedouble dt = 2.0/nPoly;
  unsigned short mPoly = nPoly, ii = 0;

  /*--- Outer loop in k-direction, which is from the base to top. ---*/
  for(unsigned int k=0; k<=nPoly; ++k, --mPoly) {

    /*--- Determine the minimum and maximum value for r and s
          for this t-value. ---*/
    const passivedouble t     = -1.0 + k*dt;
    const passivedouble rsMin =  0.5*(t-1.0);
    const passivedouble rsMax = -rsMin;

    /*--- Determine the step size along the edges of the current quad.
          Take the exceptional situation mPoly == 0 into account to avoid a
          division by zero.     ---*/
    const passivedouble dh = mPoly ? (rsMax-rsMin)/mPoly : 0.0;

    /*--- Loop over the vertices of the current quadrilateral. ---*/
    for(unsigned short j=0; j<=mPoly; ++j) {
      const passivedouble s = rsMin + j*dh;

      for(unsigned short i=0; i<=mPoly; ++i, ++ii) {
        const double r = rsMin + i*dh;
        rDOFs[ii] = r;
        sDOFs[ii] = s;
        tDOFs[ii] = t;
      }
    }
  }
}

void CFEMStandardPyra::LocationPyramidGridDOFsLGL(vector<passivedouble> &rDOFs,
                                                  vector<passivedouble> &sDOFs,
                                                  vector<passivedouble> &tDOFs) {

  /*--- Determine the number of DOFs based on the polynomial degree.
        This is done, because this function can be called for different
        values of nPoly. ---*/
  const unsigned short nD1D  = nPoly + 1;
  const unsigned short nD    = nD1D*(nD1D+1)*(2*nD1D+1)/6;

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
  Location1DGridDOFsLGL(LGLPoints);

  /*--- Determine the location of the DOFs of the standard triangle. ---*/
  vector<passivedouble> xiTri, etaTri;
  LocationTriangleGridDOFsLGL(xiTri, etaTri);

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
  rDOFs_Face014.resize(nPoly+1);
  sDOFs_Face014.resize(nPoly+1);
  tDOFs_Face014.resize(nPoly+1);

  unsigned short ii = 0;
  for(unsigned short j=0; j<=nPoly; ++j) {
    for(unsigned short i=0; i<=(nPoly-j); ++i, ++ii) {
      rDOFs_Face014[j].push_back( 0.5 + xiTri[ii] + 0.5*etaTri[ii]);
      sDOFs_Face014[j].push_back(-0.5             + 0.5*etaTri[ii]);
      tDOFs_Face014[j].push_back(etaTri[ii]);
    }
  }

  /*--- Set the coordinates of the base quad. ---*/
  ii = 0;
  for(unsigned short j=0; j<=nPoly; ++j) {
    for(unsigned short i=0; i<=nPoly; ++i, ++ii) {
      rPyraDOFsLGL[ii] =  LGLPoints[i];
      sPyraDOFsLGL[ii] =  LGLPoints[j];
      tPyraDOFsLGL[ii] = -1.0;
    }
  }

  /*--- Loop over the interior planes of the pyramid. ---*/
  for(unsigned short k=1; k<nPoly; ++k) {

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
