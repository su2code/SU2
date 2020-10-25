/*!
 * \file CFEMStandardPyraGrid.cpp
 * \brief Functions for the class CFEMStandardPyraGrid.
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

#include "../../include/fem/CFEMStandardPyraGrid.hpp"
#include "../../include/toolboxes/CGeneralSquareMatrixCM.hpp"

/*----------------------------------------------------------------------------------*/
/*             Public member functions of CFEMStandardPyraGrid.                     */
/*----------------------------------------------------------------------------------*/

CFEMStandardPyraGrid::CFEMStandardPyraGrid(const unsigned short val_nPoly,
                                           const unsigned short val_orderExact)
  : CFEMStandardPyra(val_nPoly, val_orderExact) {

  /*--- Compute the values of the Lagrangian basis functions in the integration
        points for both the equidistant and LGL point distribution. ---*/
  LagBasisIntPointsPyra(rPyraDOFsEqui, sPyraDOFsEqui, tPyraDOFsEqui, lagBasisIntEqui);
  LagBasisIntPointsPyra(rPyraDOFsLGL,  sPyraDOFsLGL,  tPyraDOFsLGL,  lagBasisIntLGL);

  /*--- Compute the values of the derivatives of the Lagrangian basis functions in
        the integration points for both the equidistant and LGL point distribution. ---*/
  DerLagBasisIntPointsPyra(rPyraDOFsEqui, sPyraDOFsEqui, tPyraDOFsEqui, derLagBasisIntEqui);
  DerLagBasisIntPointsPyra(rPyraDOFsLGL,  sPyraDOFsLGL,  tPyraDOFsLGL,  derLagBasisIntLGL);

  /*--- Create the local grid connectivities of the faces of the volume element. ---*/
  LocalGridConnFaces();

  /*--- Determine the local subconnectivity of this standard element when split
        in several linear elements. Used for a.o. plotting and searcing. ---*/
  SubConnLinearElements();

  /*--- Set up the jitted gemm call, if supported. For this particular standard
        element the derivative of the coordinates are computed, which is 3. ---*/
  SetUpJittedGEMM(nIntegrationPad, 3, nDOFs);
}

void CFEMStandardPyraGrid::CoorIntPoints(const bool                LGLDistribution,
                                         ColMajorMatrix<su2double> &matCoorDOF,
                                         ColMajorMatrix<su2double> &matCoorInt) {

  /*--- Check for which point distribution the derivatives must be computed. ---*/
  if( LGLDistribution ) {

    /*--- LGL distribution. Call the function OwnGemm to compute the
          Cartesian coordinates in the integration points. ---*/
    OwnGemm(nIntegrationPad, 3, nDOFs, lagBasisIntLGL, matCoorDOF, matCoorInt, nullptr);
  }
  else {

    /*--- Equidistant distribution. Call the function OwnGemm to compute the
          Cartesian coordinates in the integration points. ---*/
    OwnGemm(nIntegrationPad, 3, nDOFs, lagBasisIntEqui, matCoorDOF, matCoorInt, nullptr);
  }
}

void CFEMStandardPyraGrid::DerivativesCoorIntPoints(const bool                         LGLDistribution,
                                                    ColMajorMatrix<su2double>          &matCoor,
                                                    vector<ColMajorMatrix<su2double> > &matDerCoor) {

  /*--- Check for which point distribution the derivatives must be computed. ---*/
  if( LGLDistribution ) {

    /*--- LGL distribution. Call the function OwnGemm 3 times to compute the derivatives
          of the Cartesian coordinates w.r.t. the three parametric coordinates. ---*/
    OwnGemm(nIntegrationPad, 3, nDOFs, derLagBasisIntLGL[0], matCoor, matDerCoor[0], nullptr);
    OwnGemm(nIntegrationPad, 3, nDOFs, derLagBasisIntLGL[1], matCoor, matDerCoor[1], nullptr);
    OwnGemm(nIntegrationPad, 3, nDOFs, derLagBasisIntLGL[2], matCoor, matDerCoor[2], nullptr);
  }
  else {

    /*--- Equidistant distribution. Call the function OwnGemm 3 times to compute the derivatives
          of the Cartesian coordinates w.r.t. the three parametric coordinates. ---*/
    OwnGemm(nIntegrationPad, 3, nDOFs, derLagBasisIntEqui[0], matCoor, matDerCoor[0], nullptr);
    OwnGemm(nIntegrationPad, 3, nDOFs, derLagBasisIntEqui[1], matCoor, matDerCoor[1], nullptr);
    OwnGemm(nIntegrationPad, 3, nDOFs, derLagBasisIntEqui[2], matCoor, matDerCoor[2], nullptr);
  }
}

passivedouble CFEMStandardPyraGrid::WorkEstimateVolume(CConfig *config) {

  /*--- TEMPORARY IMPLEMENTATION. ---*/
  return nIntegration + 0.1*nDOFs;
}

/*----------------------------------------------------------------------------------*/
/*             Private member functions of CFEMStandardPyraGrid.                    */
/*----------------------------------------------------------------------------------*/

void CFEMStandardPyraGrid::DerLagBasisIntPointsPyra(const vector<passivedouble>            &rDOFs,
                                                    const vector<passivedouble>            &sDOFs,
                                                    const vector<passivedouble>            &tDOFs,
                                                    vector<ColMajorMatrix<passivedouble> > &derLag) {

  /*--- Determine the parametric coordinates of all integration points
        of the pyramid. ---*/
  vector<passivedouble> rInt, sInt, tInt;
  LocationAllIntegrationPoints(rInt, sInt, tInt);

  /*--- Determine the padded number of the total number of integration points. ---*/
  const unsigned short nIntTot    = rInt.size();
  const unsigned short nIntTotPad = ((nIntTot+baseVectorLen-1)/baseVectorLen)*baseVectorLen;

  /*--- Determine the inverse of the Vandermonde matrix of the DOFs. ---*/
  CGeneralSquareMatrixCM VInv(rDOFs.size());
  VandermondePyramid(rDOFs, sDOFs, tDOFs, VInv.GetMat());
  VInv.Invert();

  /*--- Determine the gradient of the Vandermonde matrix of the integration points. Make
        sure to allocate the number of rows to nIntTotPad and initialize them to zero. ---*/
  ColMajorMatrix<passivedouble> VDr(nIntTotPad,rDOFs.size()),
                                VDs(nIntTotPad,rDOFs.size()),
                                VDt(nIntTotPad,rDOFs.size());
  VDr.setConstant(0.0);
  VDs.setConstant(0.0);
  VDt.setConstant(0.0);

  GradVandermondePyramid(rInt, sInt, tInt, VDr, VDs, VDt);

  /*--- The gradients of the Lagrangian basis functions can be obtained by
        multiplying VDr, VDs, VDt and VInv. ---*/
  derLag.resize(3);
  VInv.MatMatMult('R', VDr, derLag[0]);
  VInv.MatMatMult('R', VDs, derLag[1]);
  VInv.MatMatMult('R', VDt, derLag[2]);

  /*--- Check if the sum of the elements of the relevant rows of derLag is 0. ---*/
  for(unsigned short i=0; i<nIntTot; ++i) {
    passivedouble rowSumDr = 0.0, rowSumDs = 0.0, rowSumDt = 0.0;
    for(unsigned short j=0; j<rDOFs.size(); ++j) {
      rowSumDr += derLag[0](i,j);
      rowSumDs += derLag[1](i,j);
      rowSumDt += derLag[2](i,j);
    }

    assert(fabs(rowSumDr) < 1.e-6);
    assert(fabs(rowSumDs) < 1.e-6);
    assert(fabs(rowSumDt) < 1.e-6);
  }
}

void CFEMStandardPyraGrid::LagBasisIntPointsPyra(const vector<passivedouble>   &rDOFs,
                                                 const vector<passivedouble>   &sDOFs,
                                                 const vector<passivedouble>   &tDOFs,
                                                 ColMajorMatrix<passivedouble> &lag) {

  /*--- Determine the parametric coordinates of all integration points
        of the pyramid. ---*/
  vector<passivedouble> rInt, sInt, tInt;
  LocationAllIntegrationPoints(rInt, sInt, tInt);

  /*--- Determine the padded number of the total number of integration points. ---*/
  const unsigned short nIntTot    = rInt.size();
  const unsigned short nIntTotPad = ((nIntTot+baseVectorLen-1)/baseVectorLen)*baseVectorLen;

  /*--- Determine the inverse of the Vandermonde matrix of the DOFs. ---*/
  CGeneralSquareMatrixCM VInv(rDOFs.size());
  VandermondePyramid(rDOFs, sDOFs, tDOFs, VInv.GetMat());
  VInv.Invert();

  /*--- Determine the Vandermonde matrix of the integration points. Make sure to
        allocate the number of rows to nIntTotPad and initialize them to zero. ---*/ 
  ColMajorMatrix<passivedouble> V(nIntTotPad,rDOFs.size());
  V.setConstant(0.0);
  VandermondePyramid(rInt, sInt, tInt, V);

  /*--- The Lagrangian basis functions can be obtained by multiplying
        V and VInv. ---*/
  VInv.MatMatMult('R', V, lag);

  /*--- Check if the sum of the elements of the relevant rows of lag is 1. ---*/
  for(unsigned short i=0; i<nIntTot; ++i) {
    passivedouble rowSum = -1.0;
    for(unsigned short j=0; j<rDOFs.size(); ++j) rowSum += lag(i,j);
    assert(fabs(rowSum) < 1.e-6);
  }
}

void CFEMStandardPyraGrid::GradVandermondePyramid(const vector<passivedouble>   &r,
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
  for(unsigned short i=0; i<=nPoly; ++i) {
    for(unsigned short j=0; j<=nPoly; ++j) {
      unsigned short muij = max(i,j);
      const passivedouble scaleFact = pow(2,muij+1);
      for(unsigned short k=0; k<=(nPoly-muij); ++k, ++ii) {
        for(unsigned short l=0; l<r.size(); ++l) {

          /*--- Determine the coefficients a, b and c. ---*/
          passivedouble a, b;
          const passivedouble tmp = 0.5*(1.0-t[l]);
          if(fabs(tmp) < 1.e-8) a = b = 0.0;
          else {
            a = r[l]/tmp;
            b = s[l]/tmp;
          }

          const passivedouble c = t[l];

          /*--- Determine the value of the three 1D contributions to the 3D
                basis functions as well as the gradients of these basis
                functions w.r.t. to their arguments. ---*/
          const passivedouble fa  = NormJacobi(i,0,         0,a);
          const passivedouble gb  = NormJacobi(j,0,         0,b);
          const passivedouble hc  = NormJacobi(k,2*(muij+1),0,c);
          const passivedouble dfa = GradNormJacobi(i,0,         0,a);
          const passivedouble dgb = GradNormJacobi(j,0,         0,b);
          const passivedouble dhc = GradNormJacobi(k,2*(muij+1),0,c);

          /*--- Compute the derivative of the basis function w.r.t. r and s.
                As r is only present in the parameter a the derivative of
                the basis function w.r.t. a is multiplied by dadr. A similar
                argument holds for s, which is only present in the parameter b.
                Note that the implementation is such that all possible
                singularities are divided out of the expression.  ---*/
          VDr(l,ii) = dfa*gb*hc;
          VDs(l,ii) = fa*dgb*hc;
          if(muij > 1) {
            const passivedouble tmpt = pow(tmp, (muij-1));
            VDr(l,ii) *= tmpt;
            VDs(l,ii) *= tmpt;
          }

          /*--- Compute the derivative of the basis function w.r.t. t.
                As t is present in a, b and c, all parameters must be taken into
                account when the derivative is computed. Note that the
                implementation is such that all possible singularities are
                divided out of the expression. The first part is the derivative
                of the basis function w.r.t. c, which is equal to t. --*/
          VDt(l,ii) = dhc;
          if(muij > 0) VDt(l,ii) *= pow(tmp, muij);

          if(muij > 0) {
            passivedouble tmpt = 0.5*muij*hc;
            if(muij > 1) tmpt *= pow(tmp, (muij-1));
            VDt(l,ii) -= tmpt;
          }

          VDt(l,ii) *= fa*gb;

          /*--- Add the contribution from the derivative of the basis function
                w.r.t. a multiplied by dadt and the derivative w.r.t. b multiplied
                by dbdt.                      ---*/
          VDt(l,ii) += 0.5*a*VDr(l,ii) + 0.5*b*VDs(l,ii);

          /*--- Multiply the three derivatives with the scale factor to
                obtain the correct answers. See Vandermonde3D_Pyramid for
                the explanation. ---*/
          VDr(l,ii) *= scaleFact;
          VDs(l,ii) *= scaleFact;
          VDt(l,ii) *= scaleFact;
        }
      }
    }
  }
}

void CFEMStandardPyraGrid::VandermondePyramid(const vector<passivedouble>   &r,
                                              const vector<passivedouble>   &s,
                                              const vector<passivedouble>   &t,
                                              ColMajorMatrix<passivedouble> &V) {

  /*--- For a pyramid the orthogonal basis for the reference element is
        obtained by a combination of Jacobi polynomials (of which the Legendre
        polynomials is a special case). This is the result of the
        orthonormalization of the monomial basis. ---*/
  unsigned short ii = 0;
  for(unsigned short i=0; i<=nPoly; ++i) {
    for(unsigned short j=0; j<=nPoly; ++j) {
      unsigned short muij = max(i,j);
      const passivedouble scaleFact = pow(2,muij+1);
      for(unsigned short k=0; k<=(nPoly-muij); ++k, ++ii) {
        for(unsigned short l=0; l<r.size(); ++l) {

          /*--- Determine the coefficients a, b and c. ---*/
          passivedouble a, b;
          const passivedouble tmp = 0.5*(1.0-t[l]);
          if(fabs(tmp) < 1.e-8) a = b = 0.0;
          else {
            a = r[l]/tmp;
            b = s[l]/tmp;
          }

          const passivedouble c = t[l];

          /*--- Determine the value of the current basis function in this point.
                The multiplication with scaleFact is necessary, because this
                formulation is derived in literature for a t coordinate between
                0 and 1, while in this code t varies between -1 and 1. ---*/
          const passivedouble tmpt = muij ? pow(tmp,muij) : 1.0;

          V(l,ii) = scaleFact*tmpt*NormJacobi(i,0,0,a)*NormJacobi(j,0,0,b)
                  * NormJacobi(k,2*(muij+1),0,c);
        }
      }
    }
  }
}

void CFEMStandardPyraGrid::LocalGridConnFaces(void) {

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

void CFEMStandardPyraGrid::SubConnLinearElements(void) {

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
