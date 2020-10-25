/*!
 * \file CFEMStandardTetGrid.cpp
 * \brief Functions for the class CFEMStandardTetGrid.
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

#include "../../include/fem/CFEMStandardTetGrid.hpp"
#include "../../include/toolboxes/CGeneralSquareMatrixCM.hpp"

/*----------------------------------------------------------------------------------*/
/*             Public member functions of CFEMStandardTetGrid.                      */
/*----------------------------------------------------------------------------------*/

CFEMStandardTetGrid::CFEMStandardTetGrid(const unsigned short val_nPoly,
                                         const unsigned short val_orderExact)
  : CFEMStandardTet(val_nPoly, val_orderExact) {

  /*--- Compute the values of the Lagrangian basis functions in the integration
        points for both the equidistant and LGL point distribution. ---*/
  LagBasisIntPointsTet(rTetDOFsEqui, sTetDOFsEqui, tTetDOFsEqui, lagBasisIntEqui);
  LagBasisIntPointsTet(rTetDOFsLGL,  sTetDOFsLGL,  tTetDOFsLGL,  lagBasisIntLGL);

  /*--- Compute the values of the derivatives of the Lagrangian basis functions in
        the integration points for both the equidistant and LGL point distribution. ---*/
  DerLagBasisIntPointsTet(rTetDOFsEqui, sTetDOFsEqui, tTetDOFsEqui, derLagBasisIntEqui);
  DerLagBasisIntPointsTet(rTetDOFsLGL,  sTetDOFsLGL,  tTetDOFsLGL,  derLagBasisIntLGL);

  /*--- Create the local grid connectivities of the faces of the volume element. ---*/
  LocalGridConnFaces();

  /*--- Determine the local subconnectivity of this standard element when split
        in several linear elements. Used for a.o. plotting and searcing. ---*/
  SubConnLinearElements();

  /*--- Set up the jitted gemm call, if supported. For this particular standard
        element the derivative of the coordinates are computed, which is 3. ---*/
  SetUpJittedGEMM(nIntegrationPad, 3, nDOFs);
}

void CFEMStandardTetGrid::CoorIntPoints(const bool                LGLDistribution,
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

void CFEMStandardTetGrid::DerivativesCoorIntPoints(const bool                         LGLDistribution,
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

passivedouble CFEMStandardTetGrid::WorkEstimateVolume(CConfig *config) {

  /*--- TEMPORARY IMPLEMENTATION. ---*/
  return nIntegration + 0.1*nDOFs;
}

/*----------------------------------------------------------------------------------*/
/*             Private member functions of CFEMStandardTetGrid.                     */
/*----------------------------------------------------------------------------------*/

void CFEMStandardTetGrid::DerLagBasisIntPointsTet(const vector<passivedouble>            &rDOFs,
                                                  const vector<passivedouble>            &sDOFs,
                                                  const vector<passivedouble>            &tDOFs,
                                                  vector<ColMajorMatrix<passivedouble> > &derLag) {

  /*--- Determine the padded number of the total number of integration points. ---*/
  const unsigned short nIntTot    = rTetInt.size();
  const unsigned short nIntTotPad = ((nIntTot+baseVectorLen-1)/baseVectorLen)*baseVectorLen;

  /*--- Determine the inverse of the Vandermonde matrix of the DOFs. ---*/
  CGeneralSquareMatrixCM VInv(rDOFs.size());
  VandermondeTetrahedron(rDOFs, sDOFs, tDOFs, VInv.GetMat());
  VInv.Invert();

  /*--- Determine the gradient of the Vandermonde matrix of the integration points. Make
        sure to allocate the number of rows to nIntTotPad and initialize them to zero. ---*/
  ColMajorMatrix<passivedouble> VDr(nIntTotPad,rDOFs.size()),
                                VDs(nIntTotPad,rDOFs.size()),
                                VDt(nIntTotPad,rDOFs.size());
  VDr.setConstant(0.0);
  VDs.setConstant(0.0);
  VDt.setConstant(0.0);

  GradVandermondeTetrahedron(rTetInt, sTetInt, tTetInt, VDr, VDs, VDt);

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

void CFEMStandardTetGrid::LagBasisIntPointsTet(const vector<passivedouble>   &rDOFs,
                                               const vector<passivedouble>   &sDOFs,
                                               const vector<passivedouble>   &tDOFs,
                                               ColMajorMatrix<passivedouble> &lag) {

  /*--- Determine the padded number of the total number of integration points. ---*/
  const unsigned short nIntTot    = rTetInt.size();
  const unsigned short nIntTotPad = ((nIntTot+baseVectorLen-1)/baseVectorLen)*baseVectorLen;

  /*--- Determine the inverse of the Vandermonde matrix of the DOFs. ---*/
  CGeneralSquareMatrixCM VInv(rDOFs.size());
  VandermondeTetrahedron(rDOFs, sDOFs, tDOFs, VInv.GetMat());
  VInv.Invert();

  /*--- Determine the Vandermonde matrix of the integration points. Make sure to
        allocate the number of rows to nIntTotPad and initialize them to zero. ---*/ 
  ColMajorMatrix<passivedouble> V(nIntTotPad,rDOFs.size());
  V.setConstant(0.0);
  VandermondeTetrahedron(rTetInt, sTetInt, tTetInt, V);

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

void CFEMStandardTetGrid::GradVandermondeTetrahedron(const vector<passivedouble>   &r,
                                                     const vector<passivedouble>   &s,
                                                     const vector<passivedouble>   &t,
                                                     ColMajorMatrix<passivedouble> &VDr,
                                                     ColMajorMatrix<passivedouble> &VDs,
                                                     ColMajorMatrix<passivedouble> &VDt) {

  /*--- For a tetrahedron the orthogonal basis for the reference element is obtained by a
        combination of Jacobi polynomials (of which the Legendre polynomials is a special
        case). This is the result of the orthonormalization of the monomial basis.
        Note that the sequence of the i, j and k loop must be identical to
        the evaluation of the Vandermonde matrix itself. ---*/
  unsigned short ii = 0;
  for(unsigned short i=0; i<=nPoly; ++i) {
    for(unsigned short j=0; j<=(nPoly-i); ++j) {
      for(unsigned short k=0; k<=(nPoly-i-j); ++k, ++ii) {
        for(unsigned short l=0; l<r.size(); ++l) {

          /*--- Determine the coefficients a, b and c. ---*/
          passivedouble a, b;
          passivedouble tmp = s[l] + t[l];
          if(fabs(tmp) < 1.e-8) a = -1.0;
          else                  a = -1.0 - 2.0*(1.0+r[l])/tmp;

          tmp = 1.0 - t[l];
          if(fabs(tmp) < 1.e-8) b = -1.0;
          else                  b = -1.0 + 2.0*(1.0+s[l])/tmp;

          const passivedouble c = t[l];

          /*--- Determine the value of the three 1D contributions to the 3D basis functions as
                well as the gradients of these basis functions w.r.t. to their arguments. ---*/
          const passivedouble fa  = NormJacobi(i,0,    0,a);
          const passivedouble gb  = NormJacobi(j,2*i+1,0,b);
          const passivedouble hc  = NormJacobi(k,2*(i+j+1),0,c);
          const passivedouble dfa = GradNormJacobi(i,0,    0,a);
          const passivedouble dgb = GradNormJacobi(j,2*i+1,0,b);
          const passivedouble dhc = GradNormJacobi(k,2*(i+j+1),0,c);

          /*--- Compute the derivative of the basis function w.r.t. r. As r is only present in
                the parameter a the derivative of the basis function w.r.t. a is multiplied by
                dadr. Note that the implementation is such that all possible singularities are
                divided out of the expression. ---*/
          VDr(l,ii) = sqrt(8.0)*dfa*gb*hc;
          if(i > 0) {
            VDr(l,ii) *= 4.0;
            if(i > 1 ) VDr(l,ii) *= pow((1.0-b), (i-1));
          }

          if(i+j > 1) VDr(l,ii) *= pow((1.0-c), (i+j-1));

          /*--- Compute the derivative of the basis function w.r.t. s. As s is present in both
                the parameters a and b, both variables must be taken into account when the
                derivative is computed. Note that the implementation is such that all possible
                singularities are divided out of the expression. The first part is the derivative
                of the basis function w.r.t. b multiplied by dbds. This value is stored, because
                it is needed later on to compute the derivative w.r.t. t. ---*/
          VDs(l,ii) = dgb;
          if( i ) VDs(l,ii) *= pow((1.0-b), i);

          if(i > 0) {
            tmp = i*gb;
            if(i > 1) tmp *= pow((1.0-b), (i-1));
            VDs(l,ii) -= tmp;
          }

          if(i+j > 0) {
            VDs(l,ii) *= 2.0*sqrt(8.0)*fa*hc;
            if(i+j > 1) VDs(l,ii) *= pow((1.0-c), (i+j-1));
          }

          const passivedouble dPsidbXdbds = VDs(l,ii);

          /*--- Add the contribution from the derivative of the basis function
                w.r.t. a multiplied by dads. ---*/
          VDs(l,ii) += 0.5*(a+1.0)*VDr(l,ii);

          /*--- Compute the derivative of the basis function w.r.t. t. As t is present in a, b and c,
                all parameters must be taken into account when the derivative is computed. Note that
                the implementation is such that all possible singularities are divided out of the
                expression. The first part is the derivative of the basis function w.r.t. c,
                which is equal to t. ---*/
          VDt(l,ii) = dhc;
          if(i+j > 0) {
            VDt(l,ii) *= pow((1.0-c), (i+j));

            tmp = (i+j)*hc;
            if(i+j > 1) tmp *= pow((1.0-c), (i+j-1));

            VDt(l,ii) -= tmp;
          }

          VDt(l,ii) *= sqrt(8.0)*fa*gb;
          if( i ) VDt(l,ii) *= pow((1.0-b), i);

          /*--- Add the contribution from the derivative of the basis function w.r.t. a multiplied
                by dadt and the derivative w.r.t. b multiplied by dbdt. ---*/
          VDt(l,ii) += 0.5*(a+1.0)*VDr(l,ii) + 0.5*(b+1.0)*dPsidbXdbds;
        }
      }
    }
  }
}

void CFEMStandardTetGrid::VandermondeTetrahedron(const vector<passivedouble>   &r,
                                                 const vector<passivedouble>   &s,
                                                 const vector<passivedouble>   &t,
                                                 ColMajorMatrix<passivedouble> &V) {

  /*--- For a tetrahedron the orthogonal basis for the reference element is obtained by a
        combination of Jacobi polynomials (of which the Legendre polynomials is a special
        case). This is the result of the orthonormalization of the monomial basis. ---*/
  unsigned short ii = 0;
  for(unsigned short i=0; i<=nPoly; ++i) {
    for(unsigned short j=0; j<=(nPoly-i); ++j) {
      for(unsigned short k=0; k<=(nPoly-i-j); ++k, ++ii) {
        for(unsigned short l=0; l<r.size(); ++l) {

          /*--- Determine the coefficients a, b and c. ---*/
          passivedouble a, b;
          passivedouble tmp = s[l] + t[l];
          if(fabs(tmp) < 1.e-8) a = -1.0;
          else                  a = -1.0 - 2.0*(1.0+r[l])/tmp;

          tmp = 1.0 - t[l];
          if(fabs(tmp) < 1.e-8) b = -1.0;
          else                  b = -1.0 + 2.0*(1.0+s[l])/tmp;

          const passivedouble c = t[l];

          /*--- Determine the value of the current basis function in this point. ---*/
          passivedouble tmpb = 1.0;
          if( i )   tmpb = pow((1.0-b),i);

          passivedouble tmpc = 1.0;
          if( i+j ) tmpc = pow((1.0-c),i+j);

          V(l,ii) = sqrt(8.0)*tmpb*tmpc*NormJacobi(i,0,0,a)*NormJacobi(j,2*i+1,0,b)
                  * NormJacobi(k,2*(i+j+1),0,c);
        }
      }
    }
  }
}

void CFEMStandardTetGrid::LocalGridConnFaces(void) {

  /*--- Allocate the first index of gridConnFaces, which is equal to the number
        of faces of the tetrahedron, which is 4. Reserve memory for the second
        index afterwards. ---*/
  const unsigned short nDOFsTriangle = (nPoly+1)*(nPoly+2)/2;
  gridConnFaces.resize(4);

  gridConnFaces[0].reserve(nDOFsTriangle);
  gridConnFaces[1].reserve(nDOFsTriangle);
  gridConnFaces[2].reserve(nDOFsTriangle);
  gridConnFaces[3].reserve(nDOFsTriangle);

  /*--- Loop over all the nodes of the tetrahedron and pick the correct
        ones for the faces. ---*/
  unsigned int ii = 0;
  for(unsigned short k=0; k<=nPoly; ++k) {
    const unsigned short uppBoundJ = nPoly - k;
    for(unsigned short j=0; j<=uppBoundJ; ++j) {
      const unsigned short uppBoundI = nPoly - k - j;
      for(unsigned short i=0; i<=uppBoundI; ++i, ++ii) {
        if(k == 0)           gridConnFaces[0].push_back(ii);
        if(j == 0)           gridConnFaces[1].push_back(ii);
        if(i == 0)           gridConnFaces[2].push_back(ii);
        if((i+j+k) == nPoly) gridConnFaces[3].push_back(ii);
      }
    }
  }

  /*--- Make sure that the element is to the left of the faces. ---*/
  const unsigned short n0 = 0;
  const unsigned short n1 = nPoly;
  const unsigned short n2 = nDOFsTriangle -1;
  const unsigned short n3 = nDOFs -1;

  ChangeDirectionTriangleConn(gridConnFaces[0], n0, n1, n2);
  ChangeDirectionTriangleConn(gridConnFaces[1], n0, n3, n1);
  ChangeDirectionTriangleConn(gridConnFaces[2], n0, n2, n3);
  ChangeDirectionTriangleConn(gridConnFaces[3], n1, n3, n2);
}

void CFEMStandardTetGrid::SubConnLinearElements(void) {

  /*--- The tetrahedron is split into several linear tetrahedron.
        Set the VTK sub-types accordingly. ---*/
  VTK_SubType1 = TETRAHEDRON;
  VTK_SubType2 = NONE;

  /*--- Initialize the number of DOFs for the current edges to the number of
        DOFs of the edges present in the tetrahedron. Also initialize the
        current k offset to zero.    ---*/
  unsigned short nDOFsCurrentEdges = nPoly + 1;
  unsigned short offCurrentK       = 0;

  /*--- Loop in the k-direction of the tetrahedron, which is along the edge
        from the first vertex to the last vertex of the tet. ---*/
  for(unsigned short k=0; k<nPoly; ++k) {

    /*--- Determine the offset for the next k. ---*/
    const unsigned short offNextK = offCurrentK
                                  + nDOFsCurrentEdges*(nDOFsCurrentEdges+1)/2;

    /*------------------------------------------------------------------------*/
    /*      Step 1: The tetrahedron at the end of the current i-edge.         */
    /*------------------------------------------------------------------------*/

    unsigned short n0 = offCurrentK + nDOFsCurrentEdges - 2;
    unsigned short n1 = n0 + 1;
    unsigned short n2 = n0 + nDOFsCurrentEdges;
    unsigned short n3 = offNextK + nDOFsCurrentEdges - 2;

    subConn1ForPlotting.push_back(n0);
    subConn1ForPlotting.push_back(n1);
    subConn1ForPlotting.push_back(n2);
    subConn1ForPlotting.push_back(n3);

    /*------------------------------------------------------------------------*/
    /* Step 2: The prisms that run from the end of the j-edge to the base of  */
    /*         tet just created. These prisms are subdivided into tetrahedra. */
    /*------------------------------------------------------------------------*/

    for(unsigned short i=0; i<(nDOFsCurrentEdges-2); ++i) {

      /*--- Determine the lowest j-index on the current i-line that contributes
            to the subprism. Convert that index to the local vertex number in
            the tetrahedron. ---*/
      unsigned short j = nDOFsCurrentEdges-2 -i;
      n0 = j*nDOFsCurrentEdges + i - j*(j-1)/2 + offCurrentK;

      /*--- Increment the j index to obtain the n1 vertex of the subprism. ---*/
      ++j;
      n1 = j*nDOFsCurrentEdges + i - j*(j-1)/2 + offCurrentK;

      /*--- The n2 vertex of the prism is located on the next k-level.
            The i-index remains the same, but the j index must be adapted, because
            the number of DOFs on the edges on the next k-level is one less. ---*/
      j  = nDOFsCurrentEdges-2 -i;
      n2 = j*(nDOFsCurrentEdges-1) + i - j*(j-1)/2 + offNextK;

      /*--- The n3 vertex is part of the upper triangle of the prism. Hence the
            i-index must be incremented by one, stored in ii. The j-index must
            be computed accordingly and converted to the 1D numbering. ---*/
      unsigned short ii = i+1;
      j  = nDOFsCurrentEdges-2 -ii;
      n3 = j*nDOFsCurrentEdges + ii - j*(j-1)/2 + offCurrentK;

      /*--- Increment the j index to obtain the n4 vertex of the subprism. ---*/
      ++j;
      unsigned short n4 = j*nDOFsCurrentEdges + ii - j*(j-1)/2 + offCurrentK;

      /*--- The n5 vertex of the prism is located on the next k-level.
            The i-index remains the same, but the j index must be adapted, because
            the number of DOFs on the edges on the next k-level is one less. ---*/
      j  = nDOFsCurrentEdges-2 -ii;
      unsigned short n5 = j*(nDOFsCurrentEdges-1) + ii - j*(j-1)/2 + offNextK;

      /*--- Divide the subprism into 3 subtetrahedra. ---*/
      subConn1ForPlotting.push_back(n0);
      subConn1ForPlotting.push_back(n1);
      subConn1ForPlotting.push_back(n2);
      subConn1ForPlotting.push_back(n4);

      subConn1ForPlotting.push_back(n3);
      subConn1ForPlotting.push_back(n0);
      subConn1ForPlotting.push_back(n2);
      subConn1ForPlotting.push_back(n4);

      subConn1ForPlotting.push_back(n3);
      subConn1ForPlotting.push_back(n2);
      subConn1ForPlotting.push_back(n5);
      subConn1ForPlotting.push_back(n4);
    }

    /*------------------------------------------------------------------------*/
    /* Step 3: The remaining subelements for this k-level, which are prisms   */
    /*         and hexas. These are subdivided into tetrahedra again.         */
    /*------------------------------------------------------------------------*/

    for(unsigned short i=0; i<(nDOFsCurrentEdges-2); ++i) {

      /*--- Define the variables n4 to n7, because these indices will be used
            again for the prism treated after the hexahedra. ---*/
      unsigned short n4, n5, n6, n7;

      /*--- Initialize n3, n2, n7 and n6 to the quad on the line j = 0. ---*/
      n3 = offCurrentK + i; n2 = n3 + 1;
      n7 = offNextK    + i; n6 = n7 + 1;

      /*--- Loop in the j-direction for the hexahedra on this i-row. ---*/
      for(unsigned short j=0; j<(nDOFsCurrentEdges-3-i); ++j) {

        /*--- Set the values of n0, n1, n4 and n5 from the previous hex. ---*/
        n0 = n3; n1 = n2; n4 = n7; n5 = n6;

        /*--- Nodes n3 and n7 are the j-neighbors of n0 and n4 respectively.
              Convert these (i,j,k) indices to the 1D index again. ---*/
        n3 = (j+1)*nDOFsCurrentEdges + i - j*(j+1)/2 + offCurrentK;
        n7 = (j+1)*(nDOFsCurrentEdges-1) + i - j*(j+1)/2 + offNextK;

        /*--- Nodes n2 and n6 are the i-neighbors of n3 and n7 respectively.
              Just add an offset of 1 in the 1D numbering. ---*/
        n2 = n3 + 1;
        n6 = n7 + 1;

        /*--- Divide the hexahedron in 6 tetrahedra and add their connectivity
              to the vector to store the subtetrahedra. ---*/
        subConn1ForPlotting.push_back(n0);
        subConn1ForPlotting.push_back(n3);
        subConn1ForPlotting.push_back(n4);
        subConn1ForPlotting.push_back(n1);

        subConn1ForPlotting.push_back(n3);
        subConn1ForPlotting.push_back(n7);
        subConn1ForPlotting.push_back(n4);
        subConn1ForPlotting.push_back(n1);

        subConn1ForPlotting.push_back(n4);
        subConn1ForPlotting.push_back(n7);
        subConn1ForPlotting.push_back(n5);
        subConn1ForPlotting.push_back(n1);

        subConn1ForPlotting.push_back(n1);
        subConn1ForPlotting.push_back(n2);
        subConn1ForPlotting.push_back(n3);
        subConn1ForPlotting.push_back(n7);

        subConn1ForPlotting.push_back(n2);
        subConn1ForPlotting.push_back(n5);
        subConn1ForPlotting.push_back(n6);
        subConn1ForPlotting.push_back(n7);

        subConn1ForPlotting.push_back(n1);
        subConn1ForPlotting.push_back(n5);
        subConn1ForPlotting.push_back(n2);
        subConn1ForPlotting.push_back(n7);
      }

      /*--- The prism that is in between the last hexahedron treated above
            and one of the prisms treated in step 2.
            The ID's of four of the vertices of the prism have already been
            computed for the last hexahedron above. However, they must be
            put in the correct numbering of the prism, which is done here. ---*/
      n0 = n3;
      n1 = n2;
      n3 = n7;
      n4 = n6;

      /*--- Determine the j-index of the two remaining vertices and convert
            the (i,j,k) indices to the 1D index. ---*/
      unsigned short j = nDOFsCurrentEdges-2-i;

      n2 = j*nDOFsCurrentEdges + i - j*(j-1)/2 + offCurrentK;
      n5 = j*(nDOFsCurrentEdges-1) + i - j*(j-1)/2 + offNextK;

      /*--- Divide the prism in 3 tetrahedra and add their connectivity
            to the vector to store the subtetrahedra.     ---*/
      subConn1ForPlotting.push_back(n3);
      subConn1ForPlotting.push_back(n5);
      subConn1ForPlotting.push_back(n4);
      subConn1ForPlotting.push_back(n1);

      subConn1ForPlotting.push_back(n0);
      subConn1ForPlotting.push_back(n2);
      subConn1ForPlotting.push_back(n3);
      subConn1ForPlotting.push_back(n1);

      subConn1ForPlotting.push_back(n2);
      subConn1ForPlotting.push_back(n5);
      subConn1ForPlotting.push_back(n3);
      subConn1ForPlotting.push_back(n1);
    }

    /*--- Set offCurrentK to offNextK for the next k-value and decrement
          the value of nDOFsCurrentEdges, also for the next k-value. ---*/
    offCurrentK = offNextK;
    --nDOFsCurrentEdges;
  }
}
