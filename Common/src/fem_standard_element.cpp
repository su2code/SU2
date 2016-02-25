/*!
 * \file fem_standard_element.cpp
 * \brief Class for the FEM standard elements.
 * \author E. van der Weide
 * \version 4.1.0 "Cardinal"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
 *                      Dr. Thomas D. Economon (economon@stanford.edu).
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *
 * Copyright (C) 2012-2015 SU2, the open-source CFD code.
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

#include "../include/fem_standard_element.hpp"

/*----------------------------------------------------------------------------------*/
/*           Public member functions of FEMStandardElementClass.                    */
/*----------------------------------------------------------------------------------*/

FEMStandardElementClass::FEMStandardElementClass(unsigned short val_VTK_Type,
                                                 unsigned short val_nPoly,
                                                 bool           val_constJac,
                                                 CConfig        *config) {

  /*--- Copy the function arguments to the member variables. ---*/
  VTK_Type      = val_VTK_Type;
  nPoly         = val_nPoly;
  constJacobian = val_constJac;

  /*--- Determine the polynomial degree that must be integrated exactly by the
        integration rule and the corresponding number of integration points. ---*/
  if( constJacobian )
    orderExact = (unsigned short) ceil(nPoly*config->GetQuadrature_Factor_Straight());
  else
    orderExact = (unsigned short) ceil(nPoly*config->GetQuadrature_Factor_Curved());

  /*--- Determine the element type and compute the other member variables. ---*/
  switch( VTK_Type ) {
    case LINE:          DataStandardLine();          break;
    case TRIANGLE:      DataStandardTriangle();      break;
    case QUADRILATERAL: DataStandardQuadrilateral(); break;
    case TETRAHEDRON:   DataStandardTetrahedron();   break;
    case PYRAMID:       DataStandardPyramid();       break;
    case PRISM:         DataStandardPrism();         break;
    case HEXAHEDRON:    DataStandardHexahedron();    break;
  }

  /*--- To reduce the error due to round off in the Lagrangian basis functions,
        make sure that the row sum is 1. Also check if the difference is not
        too large to be solely caused by roundoff.   ---*/
  for(unsigned short j=0; j<nIntegration; ++j) {
    unsigned int jj = j*nDOFs;
    su2double val   = 0.0;
    for(unsigned short i=0; i<nDOFs; ++i)
      val += lagBasisIntegration[jj+i];

    if(fabs(val-1.0) > 1.e-6){
      cout << "In constructor FEMStandardElementClass::FEMStandardElementClass." << endl;
      cout << "Difference is too large to be caused by roundoff" << endl;
#ifndef HAVE_MPI
      exit(EXIT_FAILURE);
#else
      MPI_Abort(MPI_COMM_WORLD,1);
      MPI_Finalize();
#endif
    }

    val = 1.0/val;
    for(unsigned short i=0; i<nDOFs; ++i)
      lagBasisIntegration[jj+i] *= val;
  }

  /*--- Do a similar check for the derivatives of the Lagrangian basis functions.
        Only in this case there is no correction, because the sum must be zero
        in the integration points.                        --*/
  bool checkGradR = !drLagBasisIntegration.empty();
  bool checkGradS = !dsLagBasisIntegration.empty();
  bool checkGradT = !dtLagBasisIntegration.empty();

  for(unsigned short j=0; j<nIntegration; ++j) {
    unsigned int jj = j*nDOFs;
    su2double valR  = 0.0, valS = 0.0, valT = 0.0;
    for(unsigned short i=0; i<nDOFs; ++i) {
      if( checkGradR ) valR += drLagBasisIntegration[jj+i];
      if( checkGradS ) valS += dsLagBasisIntegration[jj+i];
      if( checkGradT ) valT += dtLagBasisIntegration[jj+i];
    }

    if(fabs(valR) > 1.e-6 || fabs(valS) > 1.e-6 || fabs(valT) > 1.e-6) {
      cout << "In constructor FEMStandardElementClass::FEMStandardElementClass." << endl;
      cout << "Difference is too large to be caused by roundoff" << endl;
#ifndef HAVE_MPI
      exit(EXIT_FAILURE);
#else
      MPI_Abort(MPI_COMM_WORLD,1);
      MPI_Finalize();
#endif
    }
  }
}

bool FEMStandardElementClass::SameStandardElement(unsigned short val_VTK_Type,
                                                  unsigned short val_nPoly,
                                                  bool           val_constJac) {
  if(val_VTK_Type != VTK_Type)      return false;
  if(val_nPoly    != nPoly)         return false;
  if(val_constJac != constJacobian) return false;

  return true;
}

/*----------------------------------------------------------------------------------*/
/*           Private member functions of FEMStandardElementClass.                   */
/*----------------------------------------------------------------------------------*/

void FEMStandardElementClass::Copy(const FEMStandardElementClass &other) {

  VTK_Type      = other.VTK_Type;
  nPoly         = other.nPoly;
  nDOFs         = other.nDOFs;
  nIntegration  = other.nIntegration;
  constJacobian = other.constJacobian;
  orderExact    = other.orderExact;

  rDOFs = other.rDOFs;
  sDOFs = other.sDOFs;
  tDOFs = other.tDOFs;

  rIntegration = other.rIntegration;
  sIntegration = other.sIntegration;
  tIntegration = other.tIntegration;
  wIntegration = other.wIntegration;

  lagBasisIntegration = other.lagBasisIntegration;

  drLagBasisIntegration = other.drLagBasisIntegration;
  dsLagBasisIntegration = other.dsLagBasisIntegration;
  dtLagBasisIntegration = other.dtLagBasisIntegration;

  connFace0 = other.connFace0;
  connFace1 = other.connFace1;
  connFace2 = other.connFace2;
  connFace3 = other.connFace3;
  connFace4 = other.connFace4;
  connFace5 = other.connFace5;

  subConn1ForPlotting = other.subConn1ForPlotting;
  subConn2ForPlotting = other.subConn2ForPlotting;
}

void FEMStandardElementClass::DataStandardLine() {

  /*--- Determine the location of the DOFs of the edge. ---*/
  nDOFs = nPoly + 1;
  rDOFs.resize(nDOFs);

  su2double dh = 2.0/nPoly;
  for(unsigned i=0; i<nDOFs; ++i)
    rDOFs[i] = -1.0 + i*dh;

  /*--- Allocate the memory for the integration points
        and weights and determine them.                ---*/
  nIntegration = orderExact/2 + 1;
  rIntegration.resize(nIntegration);
  wIntegration.resize(nIntegration);

  GaussLegendrePoints1D(rIntegration, wIntegration);

  /*--- Compute the inverse of the Vandermonde matrix in the DOFs and
        compute the Vandermonde matrix in the integration points. ---*/
  vector<su2double> VInv(nDOFs*nDOFs), V(nDOFs*nIntegration);

  Vandermonde1D(rDOFs, VInv);
  InverseMatrix(nDOFs, VInv);

  Vandermonde1D(rIntegration, V);

  /*--- Allocate the memory for lagBasisIntegration and determine its values.
        The Lagrange basis functions in the integration points are equal to
        the interpolation coefficients from the DOFs to the integration points
        and are obtained from the matrix product V*Vinv. Note that from a
        mathematical point of view the transpose of V*VInv is stored, because
        in this way the interpolation data for an integration point is
        contiguous in memory.                                              ---*/
  lagBasisIntegration.resize(nDOFs*nIntegration);
  MatMulTranspose(V, VInv, lagBasisIntegration);

  /*--- Compute the gradients of the 1D Vandermonde matrix in the integration
        points. The vector V can be used to store the data.     ---*/
  GradVandermonde1D(rIntegration, V);

  /*--- Allocate the memory to store the derivatives in r-direction of the
        Lagrange basis functions in the integration points and determine them.
        The derivatives of the Lagrange basis functions in the integration points
        are obtained from the matrix product V*Vinv. Note that from a
        mathematical point of view the transpose of V*VInv is stored, because
        in this way the gradient data for an integration point is
        contiguous in memory.                                              ---*/
  drLagBasisIntegration.resize(nDOFs*nIntegration);
  MatMulTranspose(V, VInv, drLagBasisIntegration);

  /*--- Determine the local connectivity of the two "faces" of the line element.
        For a line element the faces are just points.   ---*/
  connFace0.reserve(1); connFace0.push_back(0);
  connFace1.reserve(1); connFace1.push_back(nPoly);

  /*--- Determine the local subconnectivity of the line element used for plotting
        purposes. This is rather trivial, because the line element is subdivided
        into nPoly linear line elements.                    ---*/
  unsigned short nnPoly = max(nPoly,(unsigned short) 1);
  for(unsigned short i=0; i<nnPoly; ++i) {
    subConn1ForPlotting.push_back(i);
    subConn1ForPlotting.push_back(i+1);
  }
}

void FEMStandardElementClass::DataStandardTriangle(void) {

  /*--- Determine the location of the DOFs of the standard triangle. ---*/
  nDOFs = (nPoly+1)*(nPoly+2)/2;
  rDOFs.resize(nDOFs);
  sDOFs.resize(nDOFs);

  su2double dh = 2.0/nPoly;

  unsigned int ii = 0;
  for(unsigned short j=0; j<=nPoly; ++j) {
    su2double s = -1.0 + j*dh;
    unsigned short uppBoundI = nPoly - j;
    for(unsigned short i=0; i<=uppBoundI; ++i, ++ii) {
      su2double r = -1.0 + i*dh;
      rDOFs[ii]   = r;
      sDOFs[ii]   = s;
    }
  }

  /*--- Determine the integration points of the standard triangle. ---*/
  IntegrationPointsTriangle();

  /*--- Compute the inverse of the Vandermonde matrix in the DOFs and
        compute the Vandermonde matrix in the integration points. ---*/
  vector<su2double> VInv(nDOFs*nDOFs), V(nDOFs*nIntegration);

  Vandermonde2D_Triangle(rDOFs, sDOFs, VInv);
  InverseMatrix(nDOFs, VInv);

  Vandermonde2D_Triangle(rIntegration, sIntegration, V);

  /*--- Allocate the memory for lagBasisIntegration and determine its values.
        The Lagrange basis functions in the integration points are equal to
        the interpolation coefficients from the DOFs to the integration points
        and are obtained from the matrix product V*Vinv. Note that from a
        mathematical point of view the transpose of V*VInv is stored, because
        in this way the interpolation data for an integration point is
        contiguous in memory.                                              ---*/
  lagBasisIntegration.resize(nDOFs*nIntegration);
  MatMulTranspose(V, VInv, lagBasisIntegration);

  /*--- Compute the gradients of the 2D Vandermonde matrix in the integration points. ---*/
  vector<su2double> VDr(nDOFs*nIntegration), VDs(nDOFs*nIntegration);
  GradVandermonde2D_Triangle(rIntegration, sIntegration, VDr, VDs);

  /*--- Allocate the memory to store the derivatives in r- and s-direction of the
        Lagrange basis functions in the integration points and determine them.
        The derivatives of the Lagrange basis functions in the integration points
        are obtained from the matrix product VDr*Vinv and VDs*Vinv. Note that from
        a mathematical point of view the transpose of the result is stored, because
        in this way the gradient data for an integration point is
        contiguous in memory.                                              ---*/
  drLagBasisIntegration.resize(nDOFs*nIntegration);
  dsLagBasisIntegration.resize(nDOFs*nIntegration);

  MatMulTranspose(VDr, VInv, drLagBasisIntegration);
  MatMulTranspose(VDs, VInv, dsLagBasisIntegration);

  /*--- Determine the local connectivity of the three "faces" of the triangle.
        For a triangular element the faces are just lines. Make sure that the
        element is to the left of the face. ---*/
  connFace0.reserve(nPoly+1); connFace1.reserve(nPoly+1); connFace2.reserve(nPoly+1);

  for(signed short i=0; i<=nPoly; ++i) connFace0.push_back(i);
  for(signed short i=0; i<=nPoly; ++i) connFace1.push_back((i+1)*(nPoly+1) - i*(i+1)/2 -1);
  for(signed short i=nPoly; i>=0; --i) connFace2.push_back(i*(nPoly+1) - i*(i-1)/2);

  /*--- Determine the local subconnectivity of the triangular element used for
        plotting purposes. ---*/
  unsigned short jj = 0;

  /*--- Loop over subedges of the left boundary of the standard triangle. ---*/
  for(unsigned short j=0; j<nPoly; ++j) {

    /*--- Check if the "down" elements must be written. ---*/
    if( j ) {

      unsigned short kk = jj - (nPoly + 1 - j); // Offset of the relevant DOF on the previous row.

      for(unsigned short i=0; i<(nPoly-j); ++i) {
        unsigned short n0 = jj + i;
        unsigned short n1 = kk + i;
        unsigned short n2 = n0 + 1;

        subConn1ForPlotting.push_back(n0);
        subConn1ForPlotting.push_back(n1);
        subConn1ForPlotting.push_back(n2);
      }
    }

    /*--- The "upp" elements must always be written.
          Determine the offset of the DOF on the next row. ---*/
    unsigned short kk = jj + (nPoly + 1 - j);

    for(unsigned short i=0; i<(nPoly-j); ++i) {
      unsigned short n0 = jj + i;
      unsigned short n1 = n0 + 1;
      unsigned short n2 = kk + i;

      subConn1ForPlotting.push_back(n0);
      subConn1ForPlotting.push_back(n1);
      subConn1ForPlotting.push_back(n2);
    }

    /*--- Set jj to kk for the next edge. ---*/
    jj = kk;
  }
}

void FEMStandardElementClass::DataStandardQuadrilateral(void) {

  /*--- Determine the location of the DOFs of the standard quadrilateral. ---*/
  nDOFs = (nPoly+1)*(nPoly+1);
  rDOFs.resize(nDOFs);
  sDOFs.resize(nDOFs);

  su2double dh = 2.0/nPoly;

  unsigned int ii = 0;
  for(unsigned short j=0; j<=nPoly; ++j)
  {
    su2double s = -1.0 + j*dh;
    for(unsigned short i=0; i<=nPoly; ++i, ++ii)
    {
      su2double r = -1.0 + i*dh;
      rDOFs[ii]   = r;
      sDOFs[ii]   = s;
    }
  }

  /*--- The 2D quadrature rule is a tensor product of the 1D Gauss-Legendre
        quadrature rule. First determine the number of integration points in 1D,
        which is stored in M, and determine them.                    ---*/
  unsigned short M = orderExact/2 + 1;

  vector<su2double> GLPoints(M), GLWeights(M);
  GaussLegendrePoints1D(GLPoints, GLWeights);

  /*--- Allocate the memory for the integration points
        and weights and determine them.                ---*/
  nIntegration = M*M;
  rIntegration.resize(nIntegration);
  sIntegration.resize(nIntegration);
  wIntegration.resize(nIntegration);

  ii = 0;
  for(unsigned short j=0; j<M; ++j) {
    for(unsigned short i=0; i<M; ++i, ++ii) {
      rIntegration[ii] = GLPoints[i];
      sIntegration[ii] = GLPoints[j];
      wIntegration[ii] = GLWeights[i]*GLWeights[j];
    }
  }

  /*--- Compute the inverse of the Vandermonde matrix in the DOFs and
        compute the Vandermonde matrix in the integration points. ---*/
  vector<su2double> VInv(nDOFs*nDOFs), V(nDOFs*nIntegration);

  Vandermonde2D_Quadrilateral(rDOFs, sDOFs, VInv);
  InverseMatrix(nDOFs, VInv);

  Vandermonde2D_Quadrilateral(rIntegration, sIntegration, V);

  /*--- Allocate the memory for lagBasisIntegration and determine its values.
        The Lagrange basis functions in the integration points are equal to
        the interpolation coefficients from the DOFs to the integration points
        and are obtained from the matrix product V*Vinv. Note that from a
        mathematical point of view the transpose of V*VInv is stored, because
        in this way the interpolation data for an integration point is
        contiguous in memory.                                              ---*/
  lagBasisIntegration.resize(nDOFs*nIntegration);
  MatMulTranspose(V, VInv, lagBasisIntegration);

  /*--- Compute the gradients of the 2D Vandermonde matrix in the integration points. ---*/
  vector<su2double> VDr(nDOFs*nIntegration), VDs(nDOFs*nIntegration);
  GradVandermonde2D_Quadrilateral(rIntegration, sIntegration, VDr, VDs);

  /*--- Allocate the memory to store the derivatives in r- and s-direction of the
        Lagrange basis functions in the integration points and determine them.
        The derivatives of the Lagrange basis functions in the integration points
        are obtained from the matrix product VDr*Vinv and VDr*Vinv. Note that from
        a mathematical point of view the transpose of the result is stored, because
        in this way the gradient data for an integration point is
        contiguous in memory.                                              ---*/
  drLagBasisIntegration.resize(nDOFs*nIntegration);
  dsLagBasisIntegration.resize(nDOFs*nIntegration);

  MatMulTranspose(VDr, VInv, drLagBasisIntegration);
  MatMulTranspose(VDs, VInv, dsLagBasisIntegration);

  /*--- Determine the local connectivity of the four "faces" of the quad element.
        For a quad element the faces are just lines. Make sure that the element
        is to the left of the face. ---*/
  connFace0.reserve(nPoly+1); connFace1.reserve(nPoly+1);
  connFace2.reserve(nPoly+1); connFace3.reserve(nPoly+1);

  unsigned short n0 = 0, n1 = nPoly, n2 = nDOFs-1, n3 = nPoly*(nPoly+1);

  for(signed short i=n0; i<=n1; ++i)          connFace0.push_back(i);
  for(signed short i=n1; i<=n2; i+=(nPoly+1)) connFace1.push_back(i);
  for(signed short i=n2; i>=n3; --i)          connFace2.push_back(i);
  for(signed short i=n3; i>=n0; i-=(nPoly+1)) connFace3.push_back(i);

  /*--- Determine the local subconnectivity of the quadrilateral element used for
        plotting purposes. Note that the connectivity of the linear subelements
        obey the VTK connectivity rule of a quadrilateral, which is different
        from the connectivity for the high order quadrilateral. ---*/
  unsigned short nnPoly = max(nPoly,(unsigned short) 1);
  for(unsigned short j=0; j<nnPoly; ++j) {
    unsigned short jj = j*(nnPoly+1);
    for(unsigned short i=0; i<nnPoly; ++i) {
      n0 = jj + i;        subConn1ForPlotting.push_back(n0);
      n1 = n0 + 1;        subConn1ForPlotting.push_back(n1);
      n2 = n1 + nPoly+1;  subConn1ForPlotting.push_back(n2);
      n3 = n2 - 1;        subConn1ForPlotting.push_back(n3);
    }
  }
}

void FEMStandardElementClass::DataStandardTetrahedron(void) {

  /*--- Determine the location of the DOFs of the standard tetrahedron. ---*/
  nDOFs = (nPoly+1)*(nPoly+2)*(nPoly+3)/6;
  rDOFs.resize(nDOFs);
  sDOFs.resize(nDOFs);
  tDOFs.resize(nDOFs);

  su2double dh = 2.0/nPoly;

  unsigned int ii = 0;
  for(unsigned short k=0; k<=nPoly; ++k) {
    su2double t = -1.0 + k*dh;
    unsigned short uppBoundJ = nPoly - k;
    for(unsigned short j=0; j<=uppBoundJ; ++j) {
      su2double s = -1.0 + j*dh;
      unsigned short uppBoundI = nPoly - k - j;
      for(unsigned short i=0; i<=uppBoundI; ++i, ++ii) {
        su2double r = -1.0 + i*dh;
        rDOFs[ii] = r;
        sDOFs[ii] = s;
        tDOFs[ii] = t;
      }
    }
  }

  /*--- Determine the integration points of the standard tetrahedron. ---*/
  IntegrationPointsTetrahedron();

  /*--- Compute the inverse of the Vandermonde matrix in the DOFs and
        compute the Vandermonde matrix in the integration points. ---*/
  vector<su2double> VInv(nDOFs*nDOFs), V(nDOFs*nIntegration);

  Vandermonde3D_Tetrahedron(rDOFs, sDOFs, tDOFs, VInv);
  InverseMatrix(nDOFs, VInv);

  Vandermonde3D_Tetrahedron(rIntegration, sIntegration, tIntegration, V);

  /*--- Allocate the memory for lagBasisIntegration and determine its values.
        The Lagrange basis functions in the integration points are equal to
        the interpolation coefficients from the DOFs to the integration points
        and are obtained from the matrix product V*Vinv. Note that from a
        mathematical point of view the transpose of V*VInv is stored, because
        in this way the interpolation data for an integration point is
        contiguous in memory.                                              ---*/
  lagBasisIntegration.resize(nDOFs*nIntegration);
  MatMulTranspose(V, VInv, lagBasisIntegration);

  /*--- Compute the gradients of the 3D Vandermonde matrix in the integration points. ---*/
  vector<su2double> VDr(nDOFs*nIntegration), VDs(nDOFs*nIntegration), VDt(nDOFs*nIntegration);
  GradVandermonde3D_Tetrahedron(rIntegration, sIntegration, tIntegration, VDr, VDs, VDr);

  /*--- Allocate the memory to store the derivatives in r-, s- and t-direction of the
        Lagrange basis functions in the integration points and determine them.
        The derivatives of the Lagrange basis functions in the integration points
        are obtained from the matrix product VDr*Vinv, VDr*Vinv and VDt*Vinv. Note that
        from a mathematical point of view the transpose of the result is stored, because
        in this way the gradient data for an integration point is contiguous in memory. ---*/
  drLagBasisIntegration.resize(nDOFs*nIntegration);
  dsLagBasisIntegration.resize(nDOFs*nIntegration);
  dtLagBasisIntegration.resize(nDOFs*nIntegration);

  MatMulTranspose(VDr, VInv, drLagBasisIntegration);
  MatMulTranspose(VDs, VInv, dsLagBasisIntegration);
  MatMulTranspose(VDt, VInv, dtLagBasisIntegration);

  /*--- Determine the local connectivity of the four faces of the tetrahedron.
        For a tetrahedron the faces are triangles. ---*/
  unsigned short nDOFsTriangle = (nPoly+1)*(nPoly+2)/2;
  connFace0.reserve(nDOFsTriangle); connFace1.reserve(nDOFsTriangle);
  connFace2.reserve(nDOFsTriangle); connFace3.reserve(nDOFsTriangle);

  ii = 0;
  for(unsigned short k=0; k<=nPoly; ++k) {
    unsigned short uppBoundJ = nPoly - k;
    for(unsigned short j=0; j<=uppBoundJ; ++j) {
      unsigned short uppBoundI = nPoly - k - j;
      for(unsigned short i=0; i<=uppBoundI; ++i, ++ii) {
        if(k == 0)           connFace0.push_back(ii);
        if(j == 0)           connFace1.push_back(ii);
        if(i == 0)           connFace2.push_back(ii);
        if((i+j+k) == nPoly) connFace3.push_back(ii);
      }
    }
  }

  /*--- Make sure that the element is to the left of the faces. ---*/
  unsigned short n0 = 0;
  unsigned short n1 = nPoly;
  unsigned short n2 = nDOFsTriangle -1;
  unsigned short n3 = nDOFs -1;

  ChangeDirectionTriangleConn(connFace0, n0, n1, n2);
  ChangeDirectionTriangleConn(connFace1, n0, n3, n1);
  ChangeDirectionTriangleConn(connFace2, n0, n2, n3);
  ChangeDirectionTriangleConn(connFace3, n1, n3, n2);

  /*--- Determine the local subconnectivity of the tetrahedron used for
        plotting purposes. The high order tetrahedron is split in several
        linear subtetrahedra.           ---*/
  SubConnTetrahedron();
}

void FEMStandardElementClass::DataStandardPyramid(void) {

  /*--- Allocate the memory for the DOFs of the standard pyramid. ---*/
  unsigned short nDOFsEdge = nPoly+1;
  nDOFs = nDOFsEdge*(nDOFsEdge+1)*(2*nDOFsEdge+1)/6;
  rDOFs.resize(nDOFs);
  sDOFs.resize(nDOFs);
  tDOFs.resize(nDOFs);

  /*--- Determine the location of the DOFs of the standard pyramid.
        The outer loop is in the k-direction, which is from base to top. ---*/
  su2double dt         = 2.0/nPoly;
  unsigned short mPoly = nPoly, ii = 0;

  for(unsigned short k=0; k<=nPoly; ++k, --mPoly) {

    /*--- Determine the minimum and maximum value for r and s for this t-value. ---*/
    su2double t     = -1.0 + k*dt;
    su2double rsMin =  0.5*(t-1.0);
    su2double rsMax = -rsMin;

    /*--- Determine the step size along the edges of the current quad.
          Take the exceptional situation mPoly == 0 into account to avoid a
          division by zero.     ---*/
    su2double dh = mPoly ? (rsMax-rsMin)/mPoly : 0.0;

    /*--- Loop over the vertices of the current quadrilateral. ---*/
    for(unsigned short j=0; j<=mPoly; ++j) {
      su2double s = rsMin + j*dh;
      for(unsigned short i=0; i<=mPoly; ++i, ++ii) {
        su2double r = rsMin + i*dh;
        rDOFs[ii] = r;
        sDOFs[ii] = s;
        tDOFs[ii] = t;
      }
    }
  }

  /*--- Determine the integration points of the standard pyramid. ---*/
  IntegrationPointsPyramid();

  /*--- Compute the inverse of the Vandermonde matrix in the DOFs and
        compute the Vandermonde matrix in the integration points. ---*/
  vector<su2double> VInv(nDOFs*nDOFs), V(nDOFs*nIntegration);

  Vandermonde3D_Pyramid(rDOFs, sDOFs, tDOFs, VInv);
  InverseMatrix(nDOFs, VInv);

  Vandermonde3D_Pyramid(rIntegration, sIntegration, tIntegration, V);

  /*--- Allocate the memory for lagBasisIntegration and determine its values.
        The Lagrange basis functions in the integration points are equal to
        the interpolation coefficients from the DOFs to the integration points
        and are obtained from the matrix product V*Vinv. Note that from a
        mathematical point of view the transpose of V*VInv is stored, because
        in this way the interpolation data for an integration point is
        contiguous in memory.                                              ---*/
  lagBasisIntegration.resize(nDOFs*nIntegration);
  MatMulTranspose(V, VInv, lagBasisIntegration);

  /*--- Compute the gradients of the 3D Vandermonde matrix in the integration points. ---*/
  vector<su2double> VDr(nDOFs*nIntegration), VDs(nDOFs*nIntegration), VDt(nDOFs*nIntegration);
  GradVandermonde3D_Pyramid(rIntegration, sIntegration, tIntegration, VDr, VDs, VDr);

  /*--- Allocate the memory to store the derivatives in r-, s- and t-direction of the
        Lagrange basis functions in the integration points and determine them.
        The derivatives of the Lagrange basis functions in the integration points
        are obtained from the matrix product VDr*Vinv, VDr*Vinv and VDt*Vinv. Note that
        from a mathematical point of view the transpose of the result is stored, because
        in this way the gradient data for an integration point is contiguous in memory. ---*/
  drLagBasisIntegration.resize(nDOFs*nIntegration);
  dsLagBasisIntegration.resize(nDOFs*nIntegration);
  dtLagBasisIntegration.resize(nDOFs*nIntegration);

  MatMulTranspose(VDr, VInv, drLagBasisIntegration);
  MatMulTranspose(VDs, VInv, dsLagBasisIntegration);
  MatMulTranspose(VDt, VInv, dtLagBasisIntegration);

  /*--- Determine the local connectivity of the five faces of the pyramid.
        For a pyramid there are four triangular faces and one quadrilateral face. ---*/
  unsigned short nDOFsQuad     = (nPoly+1)*(nPoly+1);
  unsigned short nDOFsTriangle = (nPoly+1)*(nPoly+2)/2;

  connFace0.reserve(nDOFsQuad);
  connFace1.reserve(nDOFsTriangle);
  connFace2.reserve(nDOFsTriangle);
  connFace3.reserve(nDOFsTriangle);
  connFace4.reserve(nDOFsTriangle);

  mPoly = nPoly; ii = 0;
  for(unsigned short k=0; k<=nPoly; ++k, --mPoly) {
    for(unsigned short j=0; j<=mPoly; ++j) {
      for(unsigned short i=0; i<=mPoly; ++i, ++ii) {
        if(k == 0)     connFace0.push_back(ii);
        if(j == 0)     connFace1.push_back(ii);
        if(j == nPoly) connFace2.push_back(ii);
        if(i == 0)     connFace3.push_back(ii);
        if(i == nPoly) connFace4.push_back(ii);
      }
    }
  }

  /*--- Make sure that the element is to the left of the faces. ---*/
  unsigned short n0 = 0;
  unsigned short n1 = nPoly;
  unsigned short n2 = nDOFsQuad -1;
  unsigned short n3 = n2 - nPoly;
  unsigned short n4 = nDOFs -1;

  ChangeDirectionQuadConn(connFace0, n0, n1, n2, n3);
  ChangeDirectionTriangleConn(connFace1, n0, n4, n1);
  ChangeDirectionTriangleConn(connFace2, n3, n2, n4);
  ChangeDirectionTriangleConn(connFace3, n0, n3, n4);
  ChangeDirectionTriangleConn(connFace4, n1, n4, n2);

  /*--- Determine the local subconnectivity of the pyramid used for
        plotting purposes. The high order pyramid is split in several
        linear subpyramids and subtetrahedra, i.e. two element types. ---*/
  SubConnPyramid();
}

void FEMStandardElementClass::DataStandardPrism(void) {

  /*--- Allocate the memory for the DOFs of the standard prism
        and determine its locations. ---*/
  unsigned short nDOFsEdge = nPoly+1;
  nDOFs = nDOFsEdge*nDOFsEdge*(nDOFsEdge+1)/2;
  rDOFs.resize(nDOFs);
  sDOFs.resize(nDOFs);
  tDOFs.resize(nDOFs);

  su2double dh = 2.0/nPoly;

  unsigned short ii = 0;
  for(unsigned short k=0; k<=nPoly; ++k) {
    const su2double t = -1.0 + k*dh;

    for(unsigned short j=0; j<=nPoly; ++j) {
      su2double s = -1.0 + j*dh;
      unsigned short uppBoundI = nPoly - j;
      for(unsigned short i=0; i<=uppBoundI; ++i, ++ii) {
        su2double r = -1.0 + i*dh;
        rDOFs[ii] = r;
        sDOFs[ii] = s;
        tDOFs[ii] = t;
      }
    }
  }

  /*--- The 3D quadrature rule for a prism is a tensor product of the 1D Gauss-Legendre
        quadrature rule with the triangle quadrature rule. Determine the number of
        integration points in 1D, which is stored in M, and the actual integration
        1D integration points. ---*/
  unsigned short M = orderExact/2 + 1;

  vector<su2double> GLPoints(M), GLWeights(M);
  GaussLegendrePoints1D(GLPoints, GLWeights);

  /*--- Also determine the integration rule for a triangle. ---*/
  IntegrationPointsTriangle();

  unsigned short    nIntTriangle = nIntegration;
  vector<su2double> rTriangle    = rIntegration;
  vector<su2double> sTriangle    = sIntegration;
  vector<su2double> wTriangle    = wIntegration;

  /*--- Allocate the memory for the integration points and weights
        of the prism and determine them.                ---*/
  nIntegration = M*nIntTriangle;
  rIntegration.resize(nIntegration);
  sIntegration.resize(nIntegration);
  wIntegration.resize(nIntegration);

  ii = 0;
  for(unsigned short k=0; k<M; ++k) {
    for(unsigned short j=0; j<nIntTriangle; ++j, ++ii) {
      rIntegration[ii] = rTriangle[j];
      sIntegration[ii] = sTriangle[j];
      tIntegration[ii] = GLPoints[k];
      wIntegration[ii] = wTriangle[j]*GLWeights[k];
    }
  }

  /*--- Compute the inverse of the Vandermonde matrix in the DOFs and
        compute the Vandermonde matrix in the integration points. ---*/
  vector<su2double> VInv(nDOFs*nDOFs), V(nDOFs*nIntegration);

  Vandermonde3D_Prism(rDOFs, sDOFs, tDOFs, VInv);
  InverseMatrix(nDOFs, VInv);

  Vandermonde3D_Prism(rIntegration, sIntegration, tIntegration, V);

  /*--- Allocate the memory for lagBasisIntegration and determine its values.
        The Lagrange basis functions in the integration points are equal to
        the interpolation coefficients from the DOFs to the integration points
        and are obtained from the matrix product V*Vinv. Note that from a
        mathematical point of view the transpose of V*VInv is stored, because
        in this way the interpolation data for an integration point is
        contiguous in memory.                                              ---*/
  lagBasisIntegration.resize(nDOFs*nIntegration);
  MatMulTranspose(V, VInv, lagBasisIntegration);

  /*--- Compute the gradients of the 3D Vandermonde matrix in the integration points. ---*/
  vector<su2double> VDr(nDOFs*nIntegration), VDs(nDOFs*nIntegration), VDt(nDOFs*nIntegration);
  GradVandermonde3D_Prism(rIntegration, sIntegration, tIntegration, VDr, VDs, VDr);

  /*--- Allocate the memory to store the derivatives in r-, s- and t-direction of the
        Lagrange basis functions in the integration points and determine them.
        The derivatives of the Lagrange basis functions in the integration points
        are obtained from the matrix product VDr*Vinv, VDr*Vinv and VDt*Vinv. Note that
        from a mathematical point of view the transpose of the result is stored, because
        in this way the gradient data for an integration point is contiguous in memory. ---*/
  drLagBasisIntegration.resize(nDOFs*nIntegration);
  dsLagBasisIntegration.resize(nDOFs*nIntegration);
  dtLagBasisIntegration.resize(nDOFs*nIntegration);

  MatMulTranspose(VDr, VInv, drLagBasisIntegration);
  MatMulTranspose(VDs, VInv, dsLagBasisIntegration);
  MatMulTranspose(VDt, VInv, dtLagBasisIntegration);

  /*--- Determine the local connectivity of the five faces of the prism.
        For a prism there are two triangular faces and three quadrilateral faces. ---*/
  unsigned short nDOFsQuad     = (nPoly+1)*(nPoly+1);
  unsigned short nDOFsTriangle = (nPoly+1)*(nPoly+2)/2;

  connFace0.reserve(nDOFsTriangle);
  connFace1.reserve(nDOFsTriangle);
  connFace2.reserve(nDOFsQuad);
  connFace3.reserve(nDOFsQuad);
  connFace4.reserve(nDOFsQuad);

  ii = 0;
  for(unsigned short k=0; k<=nPoly; ++k) {
    for(unsigned short j=0; j<=nPoly; ++j) {
      unsigned short uppBoundI = nPoly - j;
      for(unsigned short i=0; i<=uppBoundI; ++i, ++ii) {
        if(k == 0)         connFace0.push_back(ii);
        if(k == nPoly)     connFace1.push_back(ii);
        if(j == 0)         connFace2.push_back(ii);
        if(i == 0)         connFace3.push_back(ii);
        if((i+j) == nPoly) connFace4.push_back(ii);
      }
    }
  }

  /*--- Make sure that the element is to the left of the faces. ---*/
  unsigned short n0 = 0;
  unsigned short n1 = nPoly;
  unsigned short n2 = nDOFsTriangle -1;
  unsigned short n3 = n0 + nDOFsTriangle*nPoly;
  unsigned short n4 = n1 + nDOFsTriangle*nPoly;
  unsigned short n5 = n2 + nDOFsTriangle*nPoly;

  ChangeDirectionTriangleConn(connFace0, n0, n1, n2);
  ChangeDirectionTriangleConn(connFace1, n3, n5, n4);
  ChangeDirectionQuadConn(connFace2, n0, n3, n4, n1);
  ChangeDirectionQuadConn(connFace3, n0, n2, n5, n3);
  ChangeDirectionQuadConn(connFace4, n1, n4, n5, n2);

  /*--- Determine the local subconnectivity of the prism used for
        plotting purposes. The high order prism is split in several
        linear subprisms.                      ---*/
  SubConnPrism();
}

void FEMStandardElementClass::DataStandardHexahedron(void) {

  /*--- Allocate the memory for the DOFs of the standard hexahedron
        and determine its locations. ---*/
  unsigned short nDOFsEdge = nPoly+1;
  nDOFs = nDOFsEdge*nDOFsEdge*nDOFsEdge;
  rDOFs.resize(nDOFs);
  sDOFs.resize(nDOFs);
  tDOFs.resize(nDOFs);

  su2double dh = 2.0/nPoly;

  unsigned short ii = 0;
  for(unsigned short k=0; k<=nPoly; ++k) {
    su2double t = -1.0 + k*dh;
    for(unsigned short j=0; j<=nPoly; ++j) {
      su2double s = -1.0 + j*dh;
      for(unsigned short i=0; i<=nPoly; ++i, ++ii) {
        su2double r = -1.0 + i*dh;
        rDOFs[ii] = r;
        sDOFs[ii] = s;
        tDOFs[ii] = t;
      }
    }
  }

  /*--- The 3D quadrature rule is a tensor product of the 1D Gauss-Legendre
        quadrature rule. Determine the number of integration points in 1D, which
        is stored in M, and the actual integration 1D integration points. ---*/
  unsigned short M = orderExact/2 + 1;

  vector<su2double> GLPoints(M), GLWeights(M);
  GaussLegendrePoints1D(GLPoints, GLWeights);

  /*--- Allocate the memory for the integration points and weights
        of the hexahedron and determine them.                ---*/
  nIntegration = M*M*M;
  rIntegration.resize(nIntegration);
  sIntegration.resize(nIntegration);
  wIntegration.resize(nIntegration);

  ii = 0;
  for(unsigned short k=0; k<M; ++k) {
    for(unsigned short j=0; j<M; ++j) {
      for(unsigned short i=0; i<M; ++i, ++ii) {
        rIntegration[ii] = GLPoints[i];
        sIntegration[ii] = GLPoints[j];
        tIntegration[ii] = GLPoints[k];
        wIntegration[ii] = GLWeights[i]*GLWeights[j]*GLWeights[k];
      }
    }
  }

  /*--- Compute the inverse of the Vandermonde matrix in the DOFs and
        compute the Vandermonde matrix in the integration points. ---*/
  vector<su2double> VInv(nDOFs*nDOFs), V(nDOFs*nIntegration);

  Vandermonde3D_Hexahedron(rDOFs, sDOFs, tDOFs, VInv);
  InverseMatrix(nDOFs, VInv);

  Vandermonde3D_Hexahedron(rIntegration, sIntegration, tIntegration, V);

  /*--- Allocate the memory for lagBasisIntegration and determine its values.
        The Lagrange basis functions in the integration points are equal to
        the interpolation coefficients from the DOFs to the integration points
        and are obtained from the matrix product V*Vinv. Note that from a
        mathematical point of view the transpose of V*VInv is stored, because
        in this way the interpolation data for an integration point is
        contiguous in memory.                                              ---*/
  lagBasisIntegration.resize(nDOFs*nIntegration);
  MatMulTranspose(V, VInv, lagBasisIntegration);

  /*--- Compute the gradients of the 3D Vandermonde matrix in the integration points. ---*/
  vector<su2double> VDr(nDOFs*nIntegration), VDs(nDOFs*nIntegration), VDt(nDOFs*nIntegration);
  GradVandermonde3D_Hexahedron(rIntegration, sIntegration, tIntegration, VDr, VDs, VDr);

  /*--- Allocate the memory to store the derivatives in r-, s- and t-direction of the
        Lagrange basis functions in the integration points and determine them.
        The derivatives of the Lagrange basis functions in the integration points
        are obtained from the matrix product VDr*Vinv, VDr*Vinv and VDt*Vinv. Note that
        from a mathematical point of view the transpose of the result is stored, because
        in this way the gradient data for an integration point is contiguous in memory. ---*/
  drLagBasisIntegration.resize(nDOFs*nIntegration);
  dsLagBasisIntegration.resize(nDOFs*nIntegration);
  dtLagBasisIntegration.resize(nDOFs*nIntegration);

  MatMulTranspose(VDr, VInv, drLagBasisIntegration);
  MatMulTranspose(VDs, VInv, dsLagBasisIntegration);
  MatMulTranspose(VDt, VInv, dtLagBasisIntegration);

  /*--- Determine the local connectivity of the six faces of the hexahedron.
        For a hexahedron the faces are all quadrilateral faces. ---*/
  unsigned short nDOFsQuad = (nPoly+1)*(nPoly+1);

  connFace0.reserve(nDOFsQuad);
  connFace1.reserve(nDOFsQuad);
  connFace2.reserve(nDOFsQuad);
  connFace3.reserve(nDOFsQuad);
  connFace4.reserve(nDOFsQuad);
  connFace5.reserve(nDOFsQuad);

  ii = 0;
  for(unsigned short k=0; k<=nPoly; ++k) {
    for(unsigned short j=0; j<=nPoly; ++j) {
      for(unsigned short i=0; i<=nPoly; ++i, ++ii) {
        if(k == 0)     connFace0.push_back(ii);
        if(k == nPoly) connFace1.push_back(ii);
        if(j == 0)     connFace2.push_back(ii);
        if(j == nPoly) connFace3.push_back(ii);
        if(i == 0)     connFace4.push_back(ii);
        if(i == nPoly) connFace5.push_back(ii);
      }
    }
  }

  /*--- Make sure that the element is to the left of the faces. ---*/
  unsigned short n0 = 0;
  unsigned short n1 = nPoly;
  unsigned short n2 = nDOFsQuad -1;
  unsigned short n3 = n2 - nPoly;
  unsigned short n4 = n0 + nDOFsQuad*nPoly;
  unsigned short n5 = n1 + nDOFsQuad*nPoly;
  unsigned short n6 = n2 + nDOFsQuad*nPoly;
  unsigned short n7 = n3 + nDOFsQuad*nPoly;

  ChangeDirectionQuadConn(connFace0, n0, n1, n2, n3);
  ChangeDirectionQuadConn(connFace1, n4, n7, n6, n5);
  ChangeDirectionQuadConn(connFace2, n0, n4, n5, n1);
  ChangeDirectionQuadConn(connFace3, n3, n2, n6, n7);
  ChangeDirectionQuadConn(connFace4, n0, n3, n7, n4);
  ChangeDirectionQuadConn(connFace5, n1, n5, n6, n2);

  /*--- Determine the local subconnectivity of the hexahedron used for
        plotting purposes. The high order hexahedron is split in several
        linear subhexahedra.                      ---*/
  SubConnHexahedron();
}

void FEMStandardElementClass::SubConnTetrahedron(void) {

  /*--- Initialize the number of DOFs for the current edges to the number of
        DOFs of the edges present in the tetrahedron. Also initialize the
        current k offset to zero.    ---*/
  unsigned short nDOFsCurrentEdges = nPoly + 1;
  unsigned short offCurrentK       = 0;

  /*--- Loop in the k-direction of the tetrahedron, which is along the edge
        from the first vertex to the last vertex of the tet.    ---*/
  for(unsigned short k=0; k<nPoly; ++k) {

    /*--- Determine the offset for the next k. ---*/
    unsigned short offNextK = offCurrentK
                            + nDOFsCurrentEdges*(nDOFsCurrentEdges+1)/2;

    /*--------------------------------------------------------------------------
       Step 1: The tetrahedron at the end of the current i-edge.
      ------------------------------------------------------------------------*/

    unsigned short n0 = offCurrentK + nDOFsCurrentEdges - 2;
    unsigned short n1 = n0 + 1;
    unsigned short n2 = n0 + nDOFsCurrentEdges;
    unsigned short n3 = offNextK + nDOFsCurrentEdges - 2;

    subConn1ForPlotting.push_back(n0);
    subConn1ForPlotting.push_back(n1);
    subConn1ForPlotting.push_back(n2);
    subConn1ForPlotting.push_back(n3);

    /*--------------------------------------------------------------------------
       Step 2: The prisms that run from the end of the j-edge to the base of tet
               just created. These prisms are subdivided into tetrahedra.
      ------------------------------------------------------------------------*/

    for(unsigned short i=0; i<(nDOFsCurrentEdges-2); ++i) {

      /*--- Determine the lowest j-index on the current i-line that contributes
            to the subprism. Convert that index to the local vertex number in
            the tetrahedron.     ---*/
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

    /*--------------------------------------------------------------------------
       Step 3: The remaining subelements for this k-level, which are prisms and
               hexas. These are subdivided into tetrahedra again.
      ------------------------------------------------------------------------*/

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
              to the vector to store the subtetrahedra.       ---*/
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
            the (i,j,k) indices to the 1D index.     ---*/
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

void FEMStandardElementClass::SubConnPyramid(void) {

  /*--- Initialize the number of DOFs for the current edges to the number of
        DOFs of the edges on the base of the pyramid. Also initialize the
        current k offset to zero.     ---*/
  unsigned short nDOFsCurrentEdges = nPoly + 1;
  unsigned short offCurrentK       = 0;

  /*--- Loop in the k-direction of the pyramid. ---*/
  for(unsigned short k=0; k<nPoly; ++k) {

    /*--------------------------------------------------------------------------
       Sub-pyramids in the same direction as the original pyramid.
      ------------------------------------------------------------------------*/

    /*--- Determine the index of the first vertex of the quadrilateral of the
          next k value.    ---*/
    unsigned short kk = offCurrentK + nDOFsCurrentEdges*nDOFsCurrentEdges;

    // Loop in j-direction of the current quad.
    for(unsigned short j=0; j<(nDOFsCurrentEdges-1); ++j) {

      /*--- Index of the first vertex along the j-row of the current quad. ---*/
      unsigned short jj = offCurrentK + j*nDOFsCurrentEdges;

      /*--- Loop in i-direction of the current quad. ---*/
      for(unsigned short i=0; i<(nDOFsCurrentEdges-1); ++i) {

        /*--- Determine the local indices of the corners of the quadrilateral
              of this subpyramid as well as the top of the subpyramid.
              Store the connectivity in subConn1ForPlotting.  ---*/
        unsigned short n0 = jj + i;
        unsigned short n1 = n0 + 1;
        unsigned short n2 = n1 + nDOFsCurrentEdges;
        unsigned short n3 = n0 + nDOFsCurrentEdges;
        unsigned short n4 = kk + i;

        subConn1ForPlotting.push_back(n0);
        subConn1ForPlotting.push_back(n1);
        subConn1ForPlotting.push_back(n2);
        subConn1ForPlotting.push_back(n3);
        subConn1ForPlotting.push_back(n4);
      }

      /*--- Update kk for the next j-row. ---*/
      kk += nDOFsCurrentEdges - 1;
    }

    /*--------------------------------------------------------------------------
       Sub-pyramids in the opposite direction as the original pyramid.
      ------------------------------------------------------------------------*/

    /*--- Reset the value of kk to the index of the first vertex of the
          quadrilateral of the next k-plane.       ---*/
    kk = offCurrentK + nDOFsCurrentEdges*nDOFsCurrentEdges;

    /*--- Loop in j-direction of the current quad. Note that the starting index
          of this loop is 1.                       ---*/
    for(unsigned short j=1; j<(nDOFsCurrentEdges-1); ++j) {

      /*--- Index of the first vertex along the j-row of the current quad. ---*/
      unsigned short jj = offCurrentK + j*nDOFsCurrentEdges;

      /*--- Loop in the i-direction of this quad. Again the starting index is 1. ---*/
      for(unsigned short i=1; i<(nDOFsCurrentEdges-1); ++i) {

        /*--- Determine the local indices of the corners of the quadrilateral
              of this subpyramid as well as the top of the subpyramid.  ---*/
        unsigned short n0 = kk + i - 1;
        unsigned short n1 = n0 + 1;
        unsigned short n2 = n1 + nDOFsCurrentEdges-1;
        unsigned short n3 = n0 + nDOFsCurrentEdges-1;
        unsigned short n4 = jj + i;

        /*--- Store the connectivity of this subpyramid. Note that n1 and n3 are
              swapped, such that a positive volume is obtained according to the
              right hand rule.                      ---*/
        subConn1ForPlotting.push_back(n0);
        subConn1ForPlotting.push_back(n3);
        subConn1ForPlotting.push_back(n2);
        subConn1ForPlotting.push_back(n1);
        subConn1ForPlotting.push_back(n4);
      }

      /*--- Update kk for the next j-row. ---*/
      kk += nDOFsCurrentEdges - 1;
    }

    /*--------------------------------------------------------------------------
       Sub-tetrahedra in the j-direction.
      ------------------------------------------------------------------------*/

    /*--- Reset the value of kk again. ---*/
    kk = offCurrentK + nDOFsCurrentEdges*nDOFsCurrentEdges;

    /*--- Loop in the i-direction of the current quad. Note that the starting
          index must be 1.                         ---*/
    for(unsigned short i=1; i<(nDOFsCurrentEdges-1); ++i) {

      /*--- Loop in the j-direction of the current quad. This loop starts at 0. ---*/
      for(unsigned short j=0; j<(nDOFsCurrentEdges-1); ++j) {

        /*--- Determine the local indices of the 4 corner points of this tet.
              Its connectivity is stored in subConn2ForPlotting, because in
              subConn1ForPlotting the pyramids are stored. ---*/
        unsigned short n0 = kk + i-1 + j*(nDOFsCurrentEdges-1);
        unsigned short n1 = n0 + 1;
        unsigned short n2 = offCurrentK + i + j*nDOFsCurrentEdges;
        unsigned short n3 = n2 + nDOFsCurrentEdges;

        subConn2ForPlotting.push_back(n0);
        subConn2ForPlotting.push_back(n1);
        subConn2ForPlotting.push_back(n2);
        subConn2ForPlotting.push_back(n3);
      }
    }

    /*--------------------------------------------------------------------------
       Sub-tetrahedra in the i-direction.
      ------------------------------------------------------------------------*/

    /*--- Loop in the j-direction of the current quad. Note that the starting
          index must be 1.                          ---*/
    for(unsigned short j=1; j<(nDOFsCurrentEdges-1); ++j) {

      /*--- Index of the first vertex along the j-row of the current quad. ---*/
      unsigned short jj = offCurrentK + j*nDOFsCurrentEdges;

      /*--- Loop in the i-direction of the current quad. This loop starts at 0. ---*/
      for(unsigned short i=0; i<(nDOFsCurrentEdges-1); ++i) {

        /*--- Determine the local indices of the 4 corner points of this tet
              and store its connectivity in subConn2ForPlotting.   ---*/
        unsigned short n0 = kk + i;
        unsigned short n1 = jj + i;
        unsigned short n2 = n0 + nDOFsCurrentEdges-1;
        unsigned short n3 = n1 + 1;

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

void FEMStandardElementClass::SubConnPrism(void) {

  /*--- Determine the number of DOFs for a triangle. This is the offset in
        k-direction, the structured direction of a prisms.    ---*/
  unsigned short nDOFTria = (nPoly+1)*(nPoly+2)/2;

  /*--- Loop in k-direction, which is the structured direction of the prism. ---*/
  for(unsigned short k=0; k<nPoly; ++k) {

    /*--- Initialize the counter jj to the ID of the first vertex for this
          k-value. jj contains the ID of the first vertex on the "0-2 edge"
          of the current bottom triangle.                 ---*/
    unsigned short jj = k*nDOFTria;

    /*--- Loop over subedges of the left boundary of the standard triangle. ---*/
    for(unsigned short j=0; j<nPoly; ++j) {

      /*--- Check if the "down" elements must be written. ---*/
      if( j ) {

        /*--- Determine the offset of the relevant DOF on the previous row. ---*/
        unsigned short kk = jj - (nPoly + 1 - j);

        /*--- Loop over the edges of this row. ---*/
        for(unsigned short i=0; i<(nPoly-j); ++i) {

          /*--- Determine the local connectivity of this subelement and add
                it to subConn1ForPlotting.           ---*/
          unsigned short n0 = jj + i;
          unsigned short n1 = kk + i;
          unsigned short n2 = n0 + 1;
          unsigned short n3 = n0 + nDOFTria;
          unsigned short n4 = n1 + nDOFTria;
          unsigned short n5 = n2 + nDOFTria;

          subConn1ForPlotting.push_back(n0);
          subConn1ForPlotting.push_back(n1);
          subConn1ForPlotting.push_back(n2);
          subConn1ForPlotting.push_back(n3);
          subConn1ForPlotting.push_back(n4);
          subConn1ForPlotting.push_back(n5);
        }
      }

      /*--- The "upp" elements must always be written.
            Determine the offset of the DOF on the next row. ---*/
      unsigned short kk = jj + (nPoly + 1 - j);

      /*--- Loop over the edges of this row. ---*/
      for(unsigned short i=0; i<(nPoly-j); ++i) {

        /*--- Determine the local connectivity of this subelement and add
              it to subConn1ForPlotting.           ---*/
        unsigned short n0 = jj + i;
        unsigned short n1 = n0 + 1;
        unsigned short n2 = kk + i;
        unsigned short n3 = n0 + nDOFTria;
        unsigned short n4 = n1 + nDOFTria;
        unsigned short n5 = n2 + nDOFTria;

        subConn1ForPlotting.push_back(n0);
        subConn1ForPlotting.push_back(n1);
        subConn1ForPlotting.push_back(n2);
        subConn1ForPlotting.push_back(n3);
        subConn1ForPlotting.push_back(n4);
        subConn1ForPlotting.push_back(n5);
      }

      /*--- Set jj to kk for the next edge. ---*/
      jj = kk;
    }
  }
}

void FEMStandardElementClass::SubConnHexahedron(void) {

  /*--- Determine the nodal offset in j- and k-direction. ---*/
  unsigned short jOff = nPoly+1;
  unsigned short kOff = jOff*jOff;

  /*--- Loop over the subelements in k-direction. ---*/
  for(unsigned short k=0; k<nPoly; ++k) {

    /*--- Abbreviate the offset in k-direction used in the connectivity. ---*/
    unsigned short kk = k*kOff;

    /*--- Loop over the subelements in j-direction. ---*/
    for(unsigned short j=0; j<nPoly; ++j) {

      /*--- Abbreviate the offset in j-direction used in the connectivity. ---*/
      unsigned short jj = j*jOff;

      /*--- Loop over the subelements in i-direction. ---*/
      for(unsigned short i=0; i<nPoly; ++i) {

        /*--- Determine the 8 vertices of this subhexahedron and store
              them in subConn1ForPlotting.           ---*/
        unsigned short n0 = kk + jj + i;
        unsigned short n1 = n0 + 1;
        unsigned short n2 = n1 + jOff;
        unsigned short n3 = n0 + jOff;
        unsigned short n4 = n0 + kOff;
        unsigned short n5 = n1 + kOff;
        unsigned short n6 = n2 + kOff;
        unsigned short n7 = n3 + kOff;

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

void FEMStandardElementClass::GaussLegendrePoints1D(vector<su2double> &GLPoints,
                                                    vector<su2double> &GLWeights) {

  /*--- Determine the number of integration points. Check if the number makes sense. ---*/
  unsigned short nIntPoints = GLPoints.size();
  if(nIntPoints < 1 || nIntPoints > 100) {
    cout << "Invalid number of Gauss Legendre integration points" << endl;
#ifndef HAVE_MPI
    exit(EXIT_FAILURE);
#else
    MPI_Abort(MPI_COMM_WORLD,1);
    MPI_Finalize();
#endif
  }

  /*--- The distribution of points is symmetric. Hence only half
        the number of integration points need to be computed.    ---*/
  unsigned short nn = nIntPoints/2;

  if(2*nn < nIntPoints) GLPoints[nn] = 0.0;

  /*--- The remaing points must be computed. These are the roots of P_n(x),
        P_n is the classis Legendre polynomial of order n.
        Loop over roots to be computed.                       ---*/
  unsigned short ii = nIntPoints -1;
  for(unsigned short i=0; i<nn; ++i, --ii) {

    /*--- Initial guess of this root and determine the Legendre
          Polynomials P_n and P_{n-1} and the value f = P_n.   ---*/
    su2double x = (1.0 - (nIntPoints-1)/(8.0*nIntPoints*nIntPoints*nIntPoints))
                * cos((4*i+3)*PI_NUMBER/(4.0*nIntPoints+2.0));

    su2double Pnm1, Pn;
    Legendre(x, nIntPoints, Pnm1, Pn);
    su2double f = Pn;

    /*--- Solve the root using Halley's method.
          Loop until machine precision has been reached. ---*/
    for(;;) {

      /*--- Determine the value of the first and second derivative of f. ---*/
      su2double df  = nIntPoints*(Pnm1 - x*Pn)/(1.0-x*x);
      su2double d2f = (2.0*x*df - nIntPoints*(nIntPoints+1)*Pn)/(1.0-x*x);

      /*--- Compute the new value of the root. ---*/
      x = x - 2.0*f*df/(2.0*df*df - f*d2f);

      /*--- Determine the new value of the Legendre polynomials and
            compute the new value of f. Store the old value.        ---*/
      su2double fOld = f;
      Legendre(x, nIntPoints, Pnm1, Pn);
      f = Pn;

      /*--- Convergence criterion. ---*/
      if(fabs(fOld) <= fabs(f)) break;
    }

    /*--- Store the symmetric equivalent as well. ---*/
    GLPoints[ii] =  x;
    GLPoints[i]  = -x;
  }

  /*--- Compute the integration weights of the points.
        Make sure the sum is exactly 2.               ---*/
  su2double f = 0.0;
  for(unsigned short i=0; i<nIntPoints; ++i) {
    su2double Pnm1, Pn;
    Legendre(GLPoints[i], nIntPoints, Pnm1, Pn);
    GLWeights[i] = 2.0*(1.0-GLPoints[i]*GLPoints[i])/(nIntPoints*nIntPoints*Pnm1*Pnm1);
    f           += GLWeights[i];
  }

  if(fabs(f-2.0) > 1.e-6) {
    cout << "Something wrong in computing the Gauss Legendre weights" << endl;
#ifndef HAVE_MPI
    exit(EXIT_FAILURE);
#else
    MPI_Abort(MPI_COMM_WORLD,1);
    MPI_Finalize();
#endif
  }

  f = 2.0/f;
  for(unsigned short i=0; i<nIntPoints; ++i)
    GLWeights[i] *= f;
}

void FEMStandardElementClass::Legendre(su2double      x,
                                       unsigned short n,
                                       su2double      &Pnm1,
                                       su2double      &Pn) {

  /*--- Initialization of the polynomials Pnm1 and Pn. ---*/
  Pnm1 = 1.0;
  Pn   = x;

  /*--- Recursive definition of Pn and Pnm1. ---*/
  for(unsigned i=2; i<=n; ++i) {
    su2double tmp = Pnm1;
    Pnm1          = Pn;
    Pn            = ((2*i-1)*x*Pn - (i-1)*tmp)/i;
  }
}

su2double FEMStandardElementClass::NormJacobi(unsigned short n,
                                              unsigned short alpha,
                                              unsigned short beta,
                                              su2double      x) {
  /*--- Some abbreviations. ---*/
  su2double ap1   = alpha + 1;
  su2double bp1   = beta  + 1;
  su2double apb   = alpha + beta;
  su2double apbp1 = apb + 1;
  su2double apbp2 = apb + 2;
  su2double apbp3 = apb + 3;
  su2double b2ma2 = beta*beta - alpha*alpha;

  /*--- Initialize the normalized polynomials. ---*/
  su2double Pnm1 = sqrt(pow(0.5,apbp1)*tgamma(apbp2)/(tgamma(ap1)*tgamma(bp1)));
  su2double Pn   = 0.5*Pnm1*(apbp2*x + alpha - beta)*sqrt(apbp3/(ap1*bp1));

  /*--- Take care of the special situation of n == 0. ---*/
  if(n == 0) Pn = Pnm1;
  else
  {
    /*--- The value of the normalized Legendre polynomial must be obtained via recursion. ---*/
    for(unsigned short i=2; i<=n; ++i)
    {
      /*--- Compute the coefficients a for i and i-1 and the coefficient bi. ---*/
      unsigned short j = i-1;
      su2double   tmp  = 2*j + apb;
      su2double   aim1 = 2.0*sqrt(j*(j+apb)*(j+alpha)*(j+beta)/((tmp-1.0)*(tmp+1.0)))
                       / tmp;

      su2double bi = b2ma2/(tmp*(tmp+2.0));

      tmp          = 2*i + apb;
      su2double ai = 2.0*sqrt(i*(i+apb)*(i+alpha)*(i+beta)/((tmp-1.0)*(tmp+1.0)))
                   / tmp;

      /*--- Compute the new value of Pn and make sure to store Pnm1 correctly. ---*/
      tmp  = Pnm1;
      Pnm1 = Pn;

      Pn = ((x-bi)*Pn - aim1*tmp)/ai;
    }
  }

  /*--- Return Pn. ---*/
  return Pn;
}

su2double FEMStandardElementClass::GradNormJacobi(unsigned short n,
                                                  unsigned short alpha,
                                                  unsigned short beta,
                                                  su2double      x) {

  /*--- Make a distinction for n == 0 and n > 0. For n == 0 the derivative is
        zero, because the polynomial itself is constant. ---*/
  su2double grad;
  if(n == 0) grad = 0.0;
  else
  {
    su2double tmp = n*(n+alpha+beta+1.0);
    grad          = sqrt(tmp)*NormJacobi(n-1, alpha+1, beta+1, x);
  }

  /*--- Return the gradient. ---*/
  return grad;
}

void FEMStandardElementClass::InverseMatrix(unsigned short    n,
                                            vector<su2double> &A) {

 /*--- Check the dimensions of A. ---*/
 if(A.size() != n*n) {
   cout << "Wrong size of the A matrix in InverseMatrix" << endl;
#ifndef HAVE_MPI
    exit(EXIT_FAILURE);
#else
    MPI_Abort(MPI_COMM_WORLD,1);
    MPI_Finalize();
#endif
  }

  /*--- Create a local matrix to carry out the actual inversion. ---*/
  vector<vector<su2double> > augmentedmatrix(n, vector<su2double>(2*n));

  /*--- Copy the data from A into the first part of augmentedmatrix. Note
        that A is stored in column major order, such that also Lapack
        routines can be used to invert the matrix.       ---*/
  unsigned int ii = 0;
  for(unsigned short j=0; j<n; ++j)
    for(unsigned short i=0; i<n; ++i, ++ii)
      augmentedmatrix[i][j] = A[ii];

  /*--- Augmenting with identity matrix of similar dimensions ---*/
  for(unsigned short j=0; j<n; ++j)
    for(unsigned short i=0; i<n; ++i)
      augmentedmatrix[i][j+n] = i == j ? 1 : 0;

  /*--- Outer loop of the Gauss-Jordan elimination. ---*/
  for(unsigned short j=0; j<n; ++j) {

    /*--- Find the pivot in the current column. ---*/
    unsigned short jj = j;
    su2double  valMax = fabs(augmentedmatrix[j][j]);
    for(unsigned short i=j+1; i<n; ++i) {
      su2double val = fabs(augmentedmatrix[i][j]);
      if(val > valMax){
        jj = i;
        valMax = val;
      }
    }

    /* Swap the rows j and jj, if needed. */
    if(jj > j) {
      for(unsigned short k=j; k<2*n; ++k) {
        su2double valTmp       = augmentedmatrix[j][k];
        augmentedmatrix[j][k]  = augmentedmatrix[jj][k];
        augmentedmatrix[jj][k] = valTmp;
      }
    }

    /*--- Performing row operations to form required identity matrix out
          of the input matrix.              ---*/
    for(unsigned i=0; i<n; ++i) {
      if(i != j) {
        valMax = augmentedmatrix[i][j]/augmentedmatrix[j][j];
        for(unsigned short k=j; k<2*n; ++k)
          augmentedmatrix[i][k] -= valMax*augmentedmatrix[j][k];
      }
    }

    valMax = 1.0/augmentedmatrix[j][j];
    for(unsigned short k=j; k<2*n; ++k)
      augmentedmatrix[j][k] *= valMax;
  }

  /*--- Store the inverse in A. Again column major order is used. ---*/
  ii = 0;
  for(unsigned short j=0; j<n; ++j)
    for(unsigned short i=0; i<n; ++i, ++ii)
      A[ii] = augmentedmatrix[i][j+n];
}

void FEMStandardElementClass::Vandermonde1D(vector<su2double> &r,
                                            vector<su2double> &V) {

  /*--- Determine the number or rows of the Vandermonde matrix and check
        if the dimension of V is correct.     ---*/
  unsigned short nRows = r.size();
  if(V.size() != nRows*nDOFs) {
    cout << "Wrong size of the V matrix in Vandermonde1D" << endl;
#ifndef HAVE_MPI
    exit(EXIT_FAILURE);
#else
    MPI_Abort(MPI_COMM_WORLD,1);
    MPI_Finalize();
#endif
  }

  /*--- Compute the Vandermonde matrix. ---*/
  unsigned int ii = 0;
  for(unsigned short i=0; i<nDOFs; ++i) {
    for(unsigned short k=0; k<nRows; ++k, ++ii) {
      V[ii] = NormJacobi(i, 0, 0, r[k]);
    }
  }
}

void FEMStandardElementClass::GradVandermonde1D(vector<su2double> &r,
                                                vector<su2double> &VDr) {

  /*--- Determine the number or rows of the gradient of the Vandermonde matrix
        and check if the dimension of VDr is correct.     ---*/
  unsigned short nRows = r.size();
  if(VDr.size() != nRows*nDOFs) {
    cout << "Wrong size of the VDr matrix in GradVandermonde1D" << endl;
#ifndef HAVE_MPI
    exit(EXIT_FAILURE);
#else
    MPI_Abort(MPI_COMM_WORLD,1);
    MPI_Finalize();
#endif
  }

  /*--- Compute the gradient of the Vandermonde matrix. ---*/
  unsigned int ii = 0;
  for(unsigned short i=0; i<nDOFs; ++i) {
    for(unsigned short k=0; k<nRows; ++k, ++ii) {
      VDr[ii] = GradNormJacobi(i, 0, 0, r[k]);
    }
  }
}

void FEMStandardElementClass::Vandermonde2D_Triangle(vector<su2double> &r,
                                                     vector<su2double> &s,
                                                     vector<su2double> &V) {

  /*--- Determine the number or rows of the Vandermonde matrix and check
        if the dimension of V is correct.     ---*/
  unsigned short nRows = r.size();
  if(V.size() != nRows*nDOFs) {
    cout << "Wrong size of the V matrix in Vandermonde2D_Triangle" << endl;
#ifndef HAVE_MPI
    exit(EXIT_FAILURE);
#else
    MPI_Abort(MPI_COMM_WORLD,1);
    MPI_Finalize();
#endif
  }

  /*--- For a triangle the orthogonal basis for the reference element is obtained
        by a combination of a Jacobi polynomial and a Legendre polynomial. This
        is the result of the orthonormalization of the monomial basis. ---*/
  unsigned int ii = 0;
  for(unsigned short i=0; i<=nPoly; ++i) {
    for(unsigned short j=0; j<=(nPoly-i); ++j) {
      for(unsigned short k=0; k<nRows; ++k, ++ii) {

        /*--- Determine the coefficients a and b. ---*/
        su2double a;
        if(fabs(s[k]-1.0) < 1.e-8) a = -1.0;
        else a = 2.0*(1.0+r[k])/(1.0-s[k]) - 1.0;

        su2double b = s[k];

        /*--- Determine the value of the current basis function in this point. ---*/
        su2double tmp = pow((1.0-b),i);
        V[ii] = sqrt(2.0)*tmp*NormJacobi(i,0,0,a)*NormJacobi(j,2*i+1,0,b);
      }
    }
  }
}

void FEMStandardElementClass::GradVandermonde2D_Triangle(vector<su2double> &r,
                                                         vector<su2double> &s,
                                                         vector<su2double> &VDr,
                                                         vector<su2double> &VDs) {

  /*--- Determine the number or rows of the gradient of the Vandermonde matrix
        and check if the dimensions of VDr and VDs are correct.     ---*/
  unsigned short nRows = r.size();
  if(VDr.size() != nRows*nDOFs || VDs.size() != nRows*nDOFs) {
    cout << "Wrong size of the VDr and/or VDs matrices in GradVandermonde2D_Triangle" << endl;
#ifndef HAVE_MPI
    exit(EXIT_FAILURE);
#else
    MPI_Abort(MPI_COMM_WORLD,1);
    MPI_Finalize();
#endif
  }

  /*--- For a triangle the orthogonal basis for the reference element is obtained
        by a combination of a Jacobi polynomial and a Legendre polynomial. This
        is the result of the orthonormalization of the monomial basis. ---*/
  unsigned int ii = 0;
  for(unsigned short i=0; i<=nPoly; ++i) {
    for(unsigned short j=0; j<=(nPoly-i); ++j) {
      for(unsigned k=0; k<nRows; ++k, ++ii) {

        /*--- Determine the coefficients a and b. ---*/
        su2double a;
        if(fabs(s[k]-1.0) < 1.e-8) a = -1.0;
        else a = 2.0*(1.0+r[k])/(1.0-s[k]) - 1.0;

        su2double b = s[k];

        /*--- Determine the value of the two 1D contributions to the 2D
              basis functions as well as the gradients of these basis
              functions w.r.t. to their arguments. ---*/
        su2double fa  = NormJacobi(i,0,    0,a);
        su2double gb  = NormJacobi(j,2*i+1,0,b);
        su2double dfa = GradNormJacobi(i,0,    0,a);
        su2double dgb = GradNormJacobi(j,2*i+1,0,b);

        /*--- Determine the gradients of the basis functions w.r.t. the
              coordinates r and s. The product rule must be used in order
              to change the derivative of a to the derivative of r and s. ---*/
        VDr[ii] = sqrt(2.0)*dfa*gb;
        VDs[ii] = VDr[ii];
        if(i > 0)
        {
          su2double tmp = pow((1.0-b), (i-1));
          VDr[ii]       = 2.0*tmp*VDr[ii];
          VDs[ii]       = (a+1.0)*tmp*VDs[ii] - i*tmp*sqrt(2.0)*fa*gb;
        }

        su2double tmp = pow((1.0-b), i);
        VDs[ii] += sqrt(2.0)*fa*dgb*tmp;
      }
    }
  }
}

void FEMStandardElementClass::Vandermonde2D_Quadrilateral(vector<su2double> &r,
                                                          vector<su2double> &s,
                                                          vector<su2double> &V) {

  /*--- Determine the number or rows of the Vandermonde matrix and check
        if the dimension of V is correct.     ---*/
  unsigned short nRows = r.size();
  if(V.size() != nRows*nDOFs) {
    cout << "Wrong size of the V matrix in Vandermonde2D_Quadrilateral" << endl;
#ifndef HAVE_MPI
    exit(EXIT_FAILURE);
#else
    MPI_Abort(MPI_COMM_WORLD,1);
    MPI_Finalize();
#endif
  }

  /*--- For a quadrilateral the basis functions are the product of the 1D
        basis functions, which are the normalized Legendre polynomials.
        The Legendre polynomials are implemented via Jacobi polynomials. ---*/
  unsigned int ii = 0;
  for(unsigned short i=0; i<=nPoly; ++i) {
    for(unsigned short j=0; j<=nPoly; ++j) {
      for(unsigned short k=0; k<nRows; ++k, ++ii) {
        V[ii] = NormJacobi(i,0,0,r[k])*NormJacobi(j,0,0,s[k]);
      }
    }
  }
}

void FEMStandardElementClass::GradVandermonde2D_Quadrilateral(vector<su2double> &r,
                                                              vector<su2double> &s,
                                                              vector<su2double> &VDr,
                                                              vector<su2double> &VDs) {

  /*--- Determine the number or rows of the gradient of the Vandermonde matrix
        and check if the dimensions of VDr and VDs are correct.     ---*/
  unsigned short nRows = r.size();
  if(VDr.size() != nRows*nDOFs || VDs.size() != nRows*nDOFs) {
    cout << "Wrong size of the VDr and/or VDs matrices in GradVandermonde2D_Quadrilateral" << endl;
#ifndef HAVE_MPI
    exit(EXIT_FAILURE);
#else
    MPI_Abort(MPI_COMM_WORLD,1);
    MPI_Finalize();
#endif
  }

  /*--- For a quadrilateral the basis functions are the product of the 1D
        basis functions, which are the normalized Legendre polynomials.
        The Legendre polynomials are implemented via Jacobi polynomials.
        Hence the derivatives in r- and s-direction can be computed easily. ---*/
  unsigned int ii = 0;
  for(unsigned short i=0; i<=nPoly; ++i) {
    for(unsigned short j=0; j<=nPoly; ++j) {
      for(unsigned short k=0; k<nRows; ++k, ++ii) {
        VDr[ii] = GradNormJacobi(i,0,0,r[k])*NormJacobi(j,0,0,s[k]);
        VDs[ii] = GradNormJacobi(j,0,0,s[k])*NormJacobi(i,0,0,r[k]);
      }
    }
  }
}

void FEMStandardElementClass::Vandermonde3D_Tetrahedron(vector<su2double> &r,
                                                        vector<su2double> &s,
                                                        vector<su2double> &t,
                                                        vector<su2double> &V) {

  /*--- Determine the number or rows of the Vandermonde matrix and check
        if the dimension of V is correct.     ---*/
  unsigned short nRows = r.size();
  if(V.size() != nRows*nDOFs) {
    cout << "Wrong size of the V matrix in Vandermonde3D_Tetrahedron" << endl;
#ifndef HAVE_MPI
    exit(EXIT_FAILURE);
#else
    MPI_Abort(MPI_COMM_WORLD,1);
    MPI_Finalize();
#endif
  }

  /*--- For a tetrahedron the orthogonal basis for the reference element is obtained by a
        combination of Jacobi polynomials (of which the Legendre polynomials is a special
        case). This is the result of the orthonormalization of the monomial basis. ---*/
  unsigned int ii = 0;
  for(unsigned short i=0; i<=nPoly; ++i) {
    for(unsigned short j=0; j<=(nPoly-i); ++j) {
      for(unsigned short k=0; k<=(nPoly-i-j); ++k) {
        for(unsigned short l=0; l<nRows; ++l, ++ii) {

          /*--- Determine the coefficients a, b and c. ---*/
          su2double a, b;
          su2double tmp = s[l] + t[l];
          if(fabs(tmp) < 1.e-8) a = -1.0;
          else                  a = -1.0 - 2.0*(1.0+r[l])/tmp;

          tmp = 1.0 - t[l];
          if(fabs(tmp) < 1.e-8) b = -1.0;
          else                  b = -1.0 + 2.0*(1.0+s[l])/tmp;

          su2double c = t[l];

          /*--- Determine the value of the current basis function in this point. ---*/
          su2double tmpb = pow((1.0-b),i);
          su2double tmpc = pow((1.0-c),i+j);
          V[ii] = sqrt(8.0)*tmpb*tmpc*NormJacobi(i,0,0,a)*NormJacobi(j,2*i+1,0,b)
                * NormJacobi(k,2*(i+j+1),0,c);
        }
      }
    }
  }
}

void FEMStandardElementClass::GradVandermonde3D_Tetrahedron(vector<su2double> &r,
                                                            vector<su2double> &s,
                                                            vector<su2double> &t,
                                                            vector<su2double> &VDr,
                                                            vector<su2double> &VDs,
                                                            vector<su2double> &VDt) {

  /*--- Determine the number or rows of the gradient of the Vandermonde matrix
        and check if the dimensions of VDr, VDs and VDt are correct.     ---*/
  unsigned short nRows = r.size();
  if(VDr.size() != nRows*nDOFs || VDs.size() != nRows*nDOFs || VDt.size() != nRows*nDOFs) {
    cout << "Wrong size of the VDr, VDs and VDt matrices in GradVandermonde3D_Tetrahedron" << endl;
#ifndef HAVE_MPI
    exit(EXIT_FAILURE);
#else
    MPI_Abort(MPI_COMM_WORLD,1);
    MPI_Finalize();
#endif
  }

  /*--- For a tetrahedron the orthogonal basis for the reference element is obtained by a
        combination of Jacobi polynomials (of which the Legendre polynomials is a special
        case). This is the result of the orthonormalization of the monomial basis.
        Note that the sequence of the i, j and k loop must be identical to
        the evaluation of the Vandermonde matrix itself.                ---*/
  unsigned int ii = 0;
  for(unsigned short i=0; i<=nPoly; ++i) {
    for(unsigned short j=0; j<=(nPoly-i); ++j) {
      for(unsigned short k=0; k<=(nPoly-i-j); ++k) {
        for(unsigned short l=0; l<nRows; ++l, ++ii) {

          /*--- Determine the coefficients a, b and c. */
          su2double a, b;
          su2double tmp = s[l] + t[l];
          if(fabs(tmp) < 1.e-8) a = -1.0;
          else                  a = -1.0 - 2.0*(1.0+r[l])/tmp;

          tmp = 1.0 - t[l];
          if(fabs(tmp) < 1.e-8) b = -1.0;
          else                  b = -1.0 + 2.0*(1.0+s[l])/tmp;

          su2double c = t[l];

          /*--- Determine the value of the three 1D contributions to the 3D basis functions as
                well as the gradients of these basis functions w.r.t. to their arguments. ---*/
          su2double fa  = NormJacobi(i,0,    0,a);
          su2double gb  = NormJacobi(j,2*i+1,0,b);
          su2double hc  = NormJacobi(k,2*(i+j+1),0,c);
          su2double dfa = GradNormJacobi(i,0,    0,a);
          su2double dgb = GradNormJacobi(j,2*i+1,0,b);
          su2double dhc = GradNormJacobi(k,2*(i+j+1),0,c);

          /*--- Compute the derivative of the basis function w.r.t. r. As r is only present in
                the parameter a the derivative of the basis function w.r.t. a is multiplied by
                dadr. Note that the implementation is such that all possible singularities are
                divided out of the expression.                                  ---*/
          VDr[ii] = sqrt(8.0)*dfa*gb*hc;
          if(i   > 0) VDr[ii] *= 4.0*pow((1.0-b), (i-1));
          if(i+j > 0) VDr[ii] *=     pow((1.0-c), (i+j-1));

          /*--- Compute the derivative of the basis function w.r.t. s. As s is present in both
                the parameters a and b, both variables must be taken into account when the
                derivative is computed. Note that the implementation is such that all possible
                singularities are divided out of the expression. The first part is the derivative
                of the basis function w.r.t. b multiplied by dbds. This value is stored, because
                it is needed later on to compute the derivative w.r.t. t.       ---*/
          VDs[ii] = dgb*pow((1.0-b), i);
          if(i   > 0) VDs[ii] -= i*gb*pow((1.0-b), (i-1));
          if(i+j > 0) VDs[ii] *= 2.0*sqrt(8.0)*fa*hc*pow((1.0-c), (i+j-1));

          su2double dPsidbXdbds = VDs[ii];

          /*--- Add the contribution from the derivative of the basis function
                w.r.t. a multiplied by dads.           ---*/
          VDs[ii] += 0.5*(a+1.0)*VDr[ii];

          /*--- Compute the derivative of the basis function w.r.t. t. As t is present in a, b and c,
                all parameters must be taken into account when the derivative is computed. Note that
                the implementation is such that all possible singularities are divided out of the
                expression. The first part is the derivative of the basis function w.r.t. c,
                which is equal to t.                                     ---*/
          VDt[ii] = dhc*pow((1.0-c), (i+j));
          if(i+j > 0) VDt[ii] -= (i+j)*hc*pow((1.0-c), (i+j-1));
          VDt[ii] *= sqrt(8.0)*fa*gb*pow((1.0-b), i);

          /*--- Add the contribution from the derivative of the basis function w.r.t. a multiplied
                by dadt and the derivative w.r.t. b multiplied by dbdt.           ---*/
          VDt[ii] += 0.5*(a+1.0)*VDr[ii] + 0.5*(b+1.0)*dPsidbXdbds;
        }
      }
    }
  }
}

void FEMStandardElementClass::Vandermonde3D_Pyramid(vector<su2double> &r,
                                                    vector<su2double> &s,
                                                    vector<su2double> &t,
                                                    vector<su2double> &V) {

  /*--- Determine the number or rows of the Vandermonde matrix and check
        if the dimension of V is correct.     ---*/
  unsigned short nRows = r.size();
  if(V.size() != nRows*nDOFs) {
    cout << "Wrong size of the V matrix in Vandermonde3D_Pyramid" << endl;
#ifndef HAVE_MPI
    exit(EXIT_FAILURE);
#else
    MPI_Abort(MPI_COMM_WORLD,1);
    MPI_Finalize();
#endif
  }

  /*--- For a pyramid the orthogonal basis for the reference element is
        obtained by a combination of Jacobi polynomials (of which the Legendre
        polynomials is a special case). This is the result of the
        orthonormalization of the monomial basis. ---*/
  unsigned int ii = 0;
  for(unsigned short i=0; i<=nPoly; ++i) {
    for(unsigned short j=0; j<=nPoly; ++j) {
      unsigned short muij = max(i,j);
      for(unsigned short k=0; k<=(nPoly-muij); ++k) {
        for(unsigned short l=0; l<nRows; ++l, ++ii) {

          /*--- Determine the coefficients a, b and c. ---*/
          su2double a, b;
          su2double tmp = 0.5*(1.0-t[l]);
          if(fabs(tmp) < 1.e-8) a = b = 0.0;
          else {
            a = r[l]/tmp;
            b = s[l]/tmp;
          }

          su2double c = t[l];

          /*--- Determine the value of the current basis function in this point. ---*/
          su2double tmpt = pow(tmp,muij);
          V[ii] = tmpt*NormJacobi(i,0,0,a)*NormJacobi(j,0,0,b)
                * NormJacobi(k,2*(muij+1),0,c);
        }
      }
    }
  }
}

void FEMStandardElementClass::GradVandermonde3D_Pyramid(vector<su2double> &r,
                                                        vector<su2double> &s,
                                                        vector<su2double> &t,
                                                        vector<su2double> &VDr,
                                                        vector<su2double> &VDs,
                                                        vector<su2double> &VDt) {

  /*--- Determine the number or rows of the gradient of the Vandermonde matrix
        and check if the dimensions of VDr, VDs and VDt are correct.     ---*/
  unsigned short nRows = r.size();
  if(VDr.size() != nRows*nDOFs || VDs.size() != nRows*nDOFs || VDt.size() != nRows*nDOFs) {
    cout << "Wrong size of the VDr, VDs and VDt matrices in GradVandermonde3D_Pyramid" << endl;
#ifndef HAVE_MPI
    exit(EXIT_FAILURE);
#else
    MPI_Abort(MPI_COMM_WORLD,1);
    MPI_Finalize();
#endif
  }

  /*--- For a pyramid the orthogonal basis for the reference element is
        obtained by a combination of Jacobi polynomials (of which the Legendre
        polynomials is a special case). This is the result of the
        orthonormalization of the monomial basis.
        Note that the sequence of the i, j and k loop must be identical to
        the evaluation of the Vandermonde matrix itself.  ---*/
  unsigned int ii = 0;
  for(unsigned short i=0; i<=nPoly; ++i) {
    for(unsigned short j=0; j<=nPoly; ++j) {
      unsigned short muij = max(i,j);
      for(unsigned short k=0; k<=(nPoly-muij); ++k) {
        for(unsigned short l=0; l<nRows; ++l, ++ii) {

          /*--- Determine the coefficients a, b and c. ---*/
          su2double a, b;
          su2double tmp = 0.5*(1.0-t[l]);
          if(fabs(tmp) < 1.e-8) a = b = 0.0;
          else {
            a = r[l]/tmp;
            b = s[l]/tmp;
          }

          su2double c = t[l];

          /*--- Determine the value of the three 1D contributions to the 3D
                basis functions as well as the gradients of these basis
                functions w.r.t. to their arguments. ---*/
          su2double fa  = NormJacobi(i,0,         0,a);
          su2double gb  = NormJacobi(j,0,         0,b);
          su2double hc  = NormJacobi(k,2*(muij+1),0,c);
          su2double dfa = GradNormJacobi(i,0,         0,a);
          su2double dgb = GradNormJacobi(j,0,         0,b);
          su2double dhc = GradNormJacobi(k,2*(muij+1),0,c);

          /*--- Compute the derivative of the basis function w.r.t. r and s.
                As r is only present in the parameter a the derivative of
                the basis function w.r.t. a is multiplied by dadr. A similar
                argument holds for s, which is only present in the parameter b.
                Note that the implementation is such that all possible
                singularities are divided out of the expression.  ---*/
          VDr[ii] = dfa*gb*hc;
          VDs[ii] = fa*dgb*hc;
          if(muij > 0)
          {
            su2double tmpt = pow(tmp, (muij-1));
            VDr[ii] *= tmpt;
            VDs[ii] *= tmpt;
          }

          /*--- Compute the derivative of the basis function w.r.t. t.
                As t is present in a, b and c, all parameters must be taken into
                account when the derivative is computed. Note that the
                implementation is such that all possible singularities are
                divided out of the expression.
                The first part is the derivative of the basis function w.r.t. c,
                which is equal to t.       --*/
          VDt[ii] = dhc*pow(tmp, muij);
          if(muij > 0) VDt[ii] -= 0.5*muij*hc*pow(tmp, (muij-1));
          VDt[ii] *= fa*gb;

          /*--- Add the contribution from the derivative of the basis function
                w.r.t. a multiplied by dadt and the derivative w.r.t. b multiplied
                by dbdt.                      ---*/
          VDt[ii] += 0.5*a*VDr[ii] + 0.5*b*VDs[ii];
        }
      }
    }
  }
}

void FEMStandardElementClass::Vandermonde3D_Prism(vector<su2double> &r,
                                                  vector<su2double> &s,
                                                  vector<su2double> &t,
                                                  vector<su2double> &V) {

  /*--- Determine the number or rows of the Vandermonde matrix and check
        if the dimension of V is correct.     ---*/
  unsigned short nRows = r.size();
  if(V.size() != nRows*nDOFs) {
    cout << "Wrong size of the V matrix in Vandermonde3D_Prism" << endl;
#ifndef HAVE_MPI
    exit(EXIT_FAILURE);
#else
    MPI_Abort(MPI_COMM_WORLD,1);
    MPI_Finalize();
#endif
  }

  /*--- For a prism the orthogonal basis for the reference element is a tensor
        product of the 1D basis functions in the structured direction of the prism
        and the basis functions of a triangle. For that triangle the orthogonal
        basis is obtained by a combination of a Jacobi polynomial and a Legendre
        polynomial. This is the result of the orthonormalization of the
        monomial basis.                   ---*/
  unsigned int ii = 0;
  for(unsigned short i=0; i<=nPoly; ++i) {
    for(unsigned short j=0; j<=(nPoly-i); ++j) {
      for(unsigned short k=0; k<=nPoly; ++k) {
        for(unsigned short l=0; l<nRows; ++l, ++ii) {

          /*--- Determine the coefficients a and b. ---*/
          su2double a;
          if(fabs(s[l]-1.0) < 1.e-8) a = -1.0;
          else a = 2.0*(1.0+r[l])/(1.0-s[l]) - 1.0;

          su2double b = s[l];

          /*--- Determine the value of the current basis function in this point. ---*/
          su2double tmp = pow((1.0-b),i);
          V[ii] = sqrt(2.0)*tmp*NormJacobi(i,0,0,a)*NormJacobi(j,2*i+1,0,b)
                * NormJacobi(k,0,0,t[l]);
        }
      }
    }
  }
}

void FEMStandardElementClass::GradVandermonde3D_Prism(vector<su2double> &r,
                                                      vector<su2double> &s,
                                                      vector<su2double> &t,
                                                      vector<su2double> &VDr,
                                                      vector<su2double> &VDs,
                                                      vector<su2double> &VDt) {

  /*--- Determine the number or rows of the gradient of the Vandermonde matrix
        and check if the dimensions of VDr, VDs and VDt are correct.     ---*/
  unsigned short nRows = r.size();
  if(VDr.size() != nRows*nDOFs || VDs.size() != nRows*nDOFs || VDt.size() != nRows*nDOFs) {
    cout << "Wrong size of the VDr, VDs and VDt matrices in GradVandermonde3D_Prism" << endl;
#ifndef HAVE_MPI
    exit(EXIT_FAILURE);
#else
    MPI_Abort(MPI_COMM_WORLD,1);
    MPI_Finalize();
#endif
  }

  /*--- For a prism the orthogonal basis for the reference element is a tensor
        product of the 1D basis functions in the structured direction of the prism
        and the basis functions of a triangle. Hence the derivative matrices also
        follows this tensor product rule.
        Note that the sequence of the i, j and k loop must be identical to
        the evaluation of the Vandermonde matrix itself.          ---*/
  unsigned int ii = 0;
  for(unsigned short i=0; i<=nPoly; ++i) {
    for(unsigned short j=0; j<=(nPoly-i); ++j) {
      for(unsigned short k=0; k<=nPoly; ++k) {
        for(unsigned short l=0; l<nRows; ++l, ++ii) {

          /*--- Determine the coefficients a and b. ---*/
          su2double a;
          if(fabs(s[l]-1.0) < 1.e-8) a = -1.0;
          else a = 2.0*(1.0+r[l])/(1.0-s[l]) - 1.0;

          su2double b = s[l];

          /*--- Determine the value of the two 1D contributions to the 2D
                basis functions of the triangle as well as the gradients of
                these basis functions w.r.t. to its argument. ---*/
          su2double fa  = NormJacobi(i,0,    0,a);
          su2double gb  = NormJacobi(j,2*i+1,0,b);
          su2double dfa = GradNormJacobi(i,0,    0,a);
          su2double dgb = GradNormJacobi(j,2*i+1,0,b);

          /*--- Determine the gradients of the basis functions w.r.t. the
                coordinates r and s. The product rule must be used in order
                to change the derivative of a to the derivative of r and s. ---*/
          VDr[ii] = sqrt(2.0)*dfa*gb;
          VDs[ii] = VDr[ii];
          if(i > 0)
          {
            su2double tmp = pow((1.0-b), (i-1));
            VDr[ii]       = 2.0*tmp*VDr[ii];
            VDs[ii]       = (a+1.0)*tmp*VDs[ii] - i*tmp*sqrt(2.0)*fa*gb;
          }

          su2double tmp = pow((1.0-b), i);
          VDs[ii] += sqrt(2.0)*fa*dgb*tmp;

          /*--- Multiply VDr and VDs with the contribution from the structured
                direction of the prism.                 ---*/
          VDr[ii] *= NormJacobi(k,0,0,t[l]);
          VDs[ii] *= NormJacobi(k,0,0,t[l]);

          /*--- Compute the derivative of the basis function in the t-direction,
                which is the structured direction.             ---*/
          VDt[ii] = sqrt(2.0)*tmp*fa*gb*GradNormJacobi(k,0,0,t[l]);
        }
      }
    }
  }
}

void FEMStandardElementClass::Vandermonde3D_Hexahedron(vector<su2double> &r,
                                                       vector<su2double> &s,
                                                       vector<su2double> &t,
                                                       vector<su2double> &V) {

  /*--- Determine the number or rows of the Vandermonde matrix and check
        if the dimension of V is correct.     ---*/
  unsigned short nRows = r.size();
  if(V.size() != nRows*nDOFs) {
    cout << "Wrong size of the V matrix in Vandermonde3D_Hexahedron" << endl;
#ifndef HAVE_MPI
    exit(EXIT_FAILURE);
#else
    MPI_Abort(MPI_COMM_WORLD,1);
    MPI_Finalize();
#endif
  }

  /*--- For a hexahedron the basis functions are the tensor product of the 1D
        basis functions, which are the normalized Legendre polynomials. Note
        that the Legendre polynomials are a special kind of Jacobi polynomials. ---*/
  unsigned int ii = 0;
  for(unsigned short i=0; i<=nPoly; ++i) {
    for(unsigned short j=0; j<=nPoly; ++j) {
      for(unsigned short k=0; k<=nPoly; ++k) {
        for(unsigned short l=0; l<nRows; ++l, ++ii) {
          V[ii] = NormJacobi(i,0,0,r[l])*NormJacobi(j,0,0,s[l])
                * NormJacobi(k,0,0,t[l]);
        }
      }
    }
  }
}

void FEMStandardElementClass::GradVandermonde3D_Hexahedron(vector<su2double> &r,
                                                           vector<su2double> &s,
                                                           vector<su2double> &t,
                                                           vector<su2double> &VDr,
                                                           vector<su2double> &VDs,
                                                           vector<su2double> &VDt) {

  /*--- Determine the number or rows of the gradient of the Vandermonde matrix
        and check if the dimensions of VDr, VDs and VDt are correct.     ---*/
  unsigned short nRows = r.size();
  if(VDr.size() != nRows*nDOFs || VDs.size() != nRows*nDOFs || VDt.size() != nRows*nDOFs) {
    cout << "Wrong size of the VDr, VDs and VDt matrices in GradVandermonde3D_Hexahedron" << endl;
#ifndef HAVE_MPI
    exit(EXIT_FAILURE);
#else
    MPI_Abort(MPI_COMM_WORLD,1);
    MPI_Finalize();
#endif
  }

  /*--- For a hexahedron the basis functions are the tensor product of the 1D
        basis functions, which are the normalized Legendre polynomials.
        The derivatives are therefore easy to compute.
        Note that the Legendre polynomials are a special kind of Jacobi polynomials.
        Also note that the sequence of the i, j and k loop must be identical to
        the evaluation of the Vandermonde matrix itself.      ---*/
  unsigned int ii = 0;
  for(unsigned short i=0; i<=nPoly; ++i) {
    for(unsigned short j=0; j<=nPoly; ++j) {
      for(unsigned short k=0; k<=nPoly; ++k) {
        for(unsigned short l=0; l<nRows; ++l, ++ii) {
          VDr[ii] = NormJacobi(j,0,0,s[l])*NormJacobi(k,0,0,t[l])
                  * GradNormJacobi(i,0,0,r[l]);
          VDs[ii] = NormJacobi(i,0,0,r[l])*NormJacobi(k,0,0,t[l])
                  * GradNormJacobi(j,0,0,s[l]);
          VDt[ii] = NormJacobi(i,0,0,r[l])*NormJacobi(j,0,0,s[l])
                  * GradNormJacobi(k,0,0,t[l]);
        }
      }
    }
  }
}

void FEMStandardElementClass::MatMulTranspose(vector<su2double> &A,
                                              vector<su2double> &B,
                                              vector<su2double> &C) {

  /*--- Check if the dimensions of the matrices correspond to the
        assumptions made in this function.                    ---*/
  unsigned int dimA = nDOFs*nIntegration;
  unsigned int dimB = nDOFs*nDOFs;

  if(A.size() != dimA || B.size() != dimB || C.size() != dimA) {
    cout << "Unexpected size of the matrices in FEMStandardElementClass::MatMulTranspose" << endl;
#ifndef HAVE_MPI
    exit(EXIT_FAILURE);
#else
    MPI_Abort(MPI_COMM_WORLD,1);
    MPI_Finalize();
#endif
  }

  /*--- Carry out the actual matrix matrix multiplication and
        store the transpose of the result.                    ---*/
  for(unsigned short j=0; j<nDOFs; ++j) {
    for(unsigned short i=0; i<nIntegration; ++i) {
      unsigned int ii = i*nDOFs + j;
      C[ii] = 0.0;

      for(unsigned short k=0; k<nDOFs; ++k) {
        unsigned int indA = k*nIntegration + i;
        unsigned int indB = j*nDOFs + k;

        C[ii] += A[indA]*B[indB];
      }
    }
  }
}

void FEMStandardElementClass::ChangeDirectionQuadConn(std::vector<unsigned short> &connQuad,
                                                      unsigned short              vert0,
                                                      unsigned short              vert1,
                                                      unsigned short              vert2,
                                                      unsigned short              vert3) {

  /*--- Determine the indices of the 4 corner vertices of the quad. ---*/
  unsigned short ind0 = 0;
  unsigned short ind1 = nPoly;
  unsigned short ind2 = nDOFs - 1;
  unsigned short ind3 = ind2 - nPoly;

  /*--- There exists a linear mapping from the indices of the numbering used in the
        connectivity of this face to the indices of the target numbering. This
        mapping is of the form ii = a + b*i + c*j and jj = d + e*i + f*j, where
        ii,jj are the indices of the target numbering and i,j the indices of the
        numbering used for this face. The values of the coefficients a,b,c,d,e,f
        depend on how the corner points coincide with each other. This is
        determined below. The bool verticesDontMatch is there to check if vertices
        do not match. This should not happen, but it is checked for security. ---*/
  signed short a, b, c, d, e, f;
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
  if( verticesDontMatch ) {
    cout << "In function FEMStandardElementClass::ChangeDirectionQuadConn." << endl;
    cout << "Corner vertices do not match. This should not happen." << endl;
#ifndef HAVE_MPI
    exit(EXIT_FAILURE);
#else
    MPI_Abort(MPI_COMM_WORLD,1);
    MPI_Finalize();
#endif
  }

  /*--- Copy the connectivity, such that things works out correctly when carrying
        out the renumbering.      ---*/
  vector<unsigned short> connQuadOr = connQuad;

  /*--- Loop over the vertices of the original face to copy the connectivity data. ---*/
  unsigned short ind = 0;
  for(unsigned short j=0; j<=nPoly; ++j) {
    for(unsigned short i=0; i<=nPoly; ++i, ++ind) {

      /*--- Determine the ii and jj indices of the target, convert it to
            a 1D index and shore the modified index in connQuad. ---*/
      unsigned short ii = a + i*b + j*c;
      unsigned short jj = d + i*e + j*f;

      unsigned short iind = jj*(nPoly+1) + ii;

      connQuad[iind] = connQuadOr[ind];
    }
  }
}

void FEMStandardElementClass::ChangeDirectionTriangleConn(std::vector<unsigned short> &connTriangle,
                                                          unsigned short              vert0,
                                                          unsigned short              vert1,
                                                          unsigned short              vert2) {

  /*--- Determine the indices of the 3 corner vertices of the triangle. ---*/
  unsigned short ind0 = 0;
  unsigned short ind1 = nPoly;
  unsigned short ind2 = nDOFs-1;

  /*--- There exists a linear mapping from the indices of the numbering used in the
        connectivity of this face to the indices of the target numbering. This
        mapping is of the form ii = a + b*i + c*j and jj = d + e*i + f*j, where
        ii,jj are the indices of the target numbering and i,j the indices of the
        numbering used for this face. The values of the coefficients a,b,c,d,e,f
        depend on how the corner points coincide with each other. This is
        determined below. The bool verticesDontMatch is there to check if vertices
        do not match. This should not happen, but it is checked for security. ---*/
  signed short a, b, c, d, e, f;
  bool verticesDontMatch = false;

  if(vert0 == connTriangle[ind0]) {

    /*--- Vert0 coincides with the first vertex of the face connectivity.
          Check the situation for the neighboring vertices.  ---*/
    if(vert1 == connTriangle[ind1]) {

      /*--- The vertex numbering is the same for both faces. ---*/
      if(vert2 != connTriangle[ind2]) verticesDontMatch = true;

      a = 0; b = 1; c = 0; d = 0; e = 0; f = 1;
    }
    else if(vert1 == connTriangle[ind2]) {

      /*--- The ii-index corresponds to the j-index and the
            jj-index corresponds to i-index.    ---*/
      if(vert2 != connTriangle[ind1]) verticesDontMatch = true;

      a = 0; b = 0; c = 1; d = 0; e = 1; f = 0;
    }
    else {
      verticesDontMatch = true;
    }
  }
  else if(vert0 == connTriangle[ind1]) {

    /*--- Vert0 coincides with the second vertex of the face connectivity.
          Check the situation for the neighboring vertices.    ---*/
    if(vert1 == connTriangle[ind0]) {

      /*--- The ii-index corresponds to a combination of the i and j index
            and the jj-index corresponds to the j-index.   ---*/
      if(vert2 != connTriangle[ind2]) verticesDontMatch = true;

      a = nPoly; b = -1; c = -1; d = 0; e = 0; f = 1;
    }
    else if(vert1 == connTriangle[ind2]) {

      /*--- The jj-index corresponds to a combination of the i and j index
            and the ii-index corresponds to the j-index.  ---*/
      if(vert2 != connTriangle[ind0]) verticesDontMatch = true;

      a = 0; b = 0; c = 1; d = nPoly; e = -1; f = -1;
    }
    else {
      verticesDontMatch = true;
    }
  }
  else if(vert0 == connTriangle[ind2]) {

    /*--- Vert0 coincides with the third vertex of the face connectivity.
          Check the situation for the neighboring vertices.    ---*/
    if(vert1 == connTriangle[ind0]) {

      /*--- The ii-index corresponds to a combination of the i and j index
            and the jj-index corresponds with the i-index.  ---*/
      if(vert2 != connTriangle[ind1]) verticesDontMatch = true;

      a = nPoly; b = -1; c = -1; d = 0; e = 1; f = 0;
    }
    else if(vert1 == connTriangle[ind1]) {

      /*--- The jj-index corresponds to a combination of the i and j index
            and the ii-index corresponds with the i-index.   ---*/
      if(vert2 != connTriangle[ind0]) verticesDontMatch = true;

      a = 0; b = 1; c = 0; d = nPoly; e = -1; f = -1;
    }
    else {
      verticesDontMatch = true;
    }
  }
  else {
    verticesDontMatch = true;
  }

  /*--- If non-matching vertices have been found, terminate with an error message. ---*/
  if( verticesDontMatch ) {
    cout << "In function FEMStandardElementClass::ChangeDirectionTriangleConn." << endl;
    cout << "Corner vertices do not match. This should not happen." << endl;
#ifndef HAVE_MPI
    exit(EXIT_FAILURE);
#else
    MPI_Abort(MPI_COMM_WORLD,1);
    MPI_Finalize();
#endif
  }

  /*--- Copy the connectivity, such that things works out correctly when carrying
        out the renumbering.  ---*/
  vector<unsigned short> connTriangleOr = connTriangle;

  /*--- Loop over the vertices of the original face to copy the connectivity data. ---*/
  unsigned short ind = 0;
  for(unsigned short j=0; j<=nPoly; ++j) {
    for(unsigned short i=0; i<=(nPoly-j); ++i, ++ind) {

      /*--- Determine the ii and jj indices of the target, convert it to
            a 1D index and shore the modified index in connTriangle. ---*/
      unsigned short ii = a + i*b + j*c;
      unsigned short jj = d + i*e + j*f;

      unsigned short iind = jj*(nPoly+1) + ii - jj*(jj-1)/2;

      connTriangle[iind] = connTriangleOr[ind];
    }
  }
}

void FEMStandardElementClass::IntegrationPointsTriangle(void) {

  /*--- Set the number of integration points, depending on the order of
        polynomials that must be integrated exactly. ---*/
  switch( orderExact ) {
    case  1: nIntegration =   1; break;
    case  2: nIntegration =   3; break;
    case  3: nIntegration =   6; break;
    case  4: nIntegration =   6; break;
    case  5: nIntegration =   7; break;
    case  6: nIntegration =  12; break;
    case  7: nIntegration =  15; break;
    case  8: nIntegration =  16; break;
    case  9: nIntegration =  19; break;
    case 10: nIntegration =  25; break;
    case 11: nIntegration =  28; break;
    case 12: nIntegration =  36; break;
    case 13: nIntegration =  40; break;
    case 14: nIntegration =  46; break;
    case 15: nIntegration =  54; break;
    case 16: nIntegration =  58; break;
    case 17: nIntegration =  66; break;
    case 18: nIntegration =  73; break;
    case 19: nIntegration =  82; break;
    case 20: nIntegration =  85; break;
    case 21: nIntegration =  93; break;
    case 22: nIntegration = 100; break;
    case 23: nIntegration = 106; break;
    case 24: nIntegration = 118; break;
    case 25: nIntegration = 126; break;
    case 26: nIntegration = 138; break;
    case 27: nIntegration = 145; break;
    case 28: nIntegration = 225; break;
    default:
      cout << "FEMStandardElementClass::IntegrationPointsTriangle: Polynomial order not supported" << endl;
#ifndef HAVE_MPI
    exit(EXIT_FAILURE);
#else
    MPI_Abort(MPI_COMM_WORLD,1);
    MPI_Finalize();
#endif
  }

  /*--- Allocate the memory for the integration points and their weights. ---*/
  rIntegration.resize(nIntegration);
  sIntegration.resize(nIntegration);
  wIntegration.resize(nIntegration);

  /*--- Set the pointers to the data arrays of the variables just allocated, such
        that the names are shorter. This is useful for the code below. ---*/
  su2double *r = rIntegration.data();
  su2double *s = sIntegration.data();
  su2double *w = wIntegration.data();
 
  /*--- Set the data for the integration points, depending on the order.
        These integration rules come from the Matlab codes corresponding to
        the book "Nodal Discontinuous Methods: Algorithms, Analysis, and Applications"
        by Jan S. Hesthaven and Tim Warburton. ---*/
  switch( orderExact ) {
    case  1: {
      r[0] = -0.333333333333333; s[0] = -0.333333333333333; w[0] = 2.000000000000000;
    
      break;
    }
    case  2: {
      r[0] = -0.666666666666667; s[0] = -0.666666666666667; w[0] = 0.666666666666667;
      r[1] =  0.333333333333333; s[1] = -0.666666666666667; w[1] = 0.666666666666667;
      r[2] = -0.666666666666667; s[2] =  0.333333333333333; w[2] = 0.666666666666667;
    
      break;
    }
    case  3: {
      r[0] = -0.816847572980458; s[0] = -0.816847572980458; w[0] = 0.219903487310644;
      r[1] =  0.633695145960917; s[1] = -0.816847572980459; w[1] = 0.219903487310644;
      r[2] = -0.816847572980459; s[2] =  0.633695145960917; w[2] = 0.219903487310644;
      r[3] = -0.108103018168070; s[3] = -0.108103018168070; w[3] = 0.446763179356023;
      r[4] = -0.783793963663860; s[4] = -0.108103018168070; w[4] = 0.446763179356023;
      r[5] = -0.108103018168070; s[5] = -0.783793963663860; w[5] = 0.446763179356023;
    
      break;
    }
    case  4: {
      r[0] = -0.816847572980458; s[0] = -0.816847572980458; w[0] = 0.219903487310644;
      r[1] =  0.633695145960917; s[1] = -0.816847572980459; w[1] = 0.219903487310644;
      r[2] = -0.816847572980459; s[2] =  0.633695145960917; w[2] = 0.219903487310644;
      r[3] = -0.108103018168070; s[3] = -0.108103018168070; w[3] = 0.446763179356023;
      r[4] = -0.783793963663860; s[4] = -0.108103018168070; w[4] = 0.446763179356023;
      r[5] = -0.108103018168070; s[5] = -0.783793963663860; w[5] = 0.446763179356023;
    
      break;
    }
    case  5: {
      r[0] = -0.333333333333333; s[0] = -0.333333333333333; w[0] = 0.450000000000000;
      r[1] = -0.059715871789770; s[1] = -0.059715871789770; w[1] = 0.264788305577012;
      r[2] = -0.880568256420460; s[2] = -0.059715871789770; w[2] = 0.264788305577012;
      r[3] = -0.059715871789770; s[3] = -0.880568256420460; w[3] = 0.264788305577012;
      r[4] = -0.797426985353087; s[4] = -0.797426985353087; w[4] = 0.251878361089654;
      r[5] =  0.594853970706175; s[5] = -0.797426985353087; w[5] = 0.251878361089654;
      r[6] = -0.797426985353087; s[6] =  0.594853970706175; w[6] = 0.251878361089654;
    
      break;
    }
    case  6: {
      r[0]  = -0.501426509658179; s[0]  = -0.501426509658179; w[0]  = 0.233572551452759;
      r[1]  =  0.002853019316358; s[1]  = -0.501426509658179; w[1]  = 0.233572551452759;
      r[2]  = -0.501426509658179; s[2]  =  0.002853019316358; w[2]  = 0.233572551452759;
      r[3]  = -0.873821971016996; s[3]  = -0.873821971016996; w[3]  = 0.101689812740414;
      r[4]  =  0.747643942033991; s[4]  = -0.873821971016996; w[4]  = 0.101689812740414;
      r[5]  = -0.873821971016996; s[5]  =  0.747643942033991; w[5]  = 0.101689812740414;
      r[6]  = -0.379295097932431; s[6]  = -0.893709900310366; w[6]  = 0.165702151236747;
      r[7]  = -0.893709900310366; s[7]  = -0.379295097932431; w[7]  = 0.165702151236747;
      r[8]  =  0.273004998242797; s[8]  = -0.893709900310366; w[8]  = 0.165702151236747;
      r[9]  = -0.893709900310366; s[9]  =  0.273004998242797; w[9]  = 0.165702151236747;
      r[10] =  0.273004998242797; s[10] = -0.379295097932431; w[10] = 0.165702151236747;
      r[11] = -0.379295097932431; s[11] =  0.273004998242797; w[11] = 0.165702151236747;
    
      break;
    }
    case  7: {
      r[0]  = -0.158765070723907; s[0]  = -0.158765070723907; w[0]  = 0.278402959350494;
      r[1]  = -0.682469858552186; s[1]  = -0.158765070723907; w[1]  = 0.278402959350494;
      r[2]  = -0.158765070723907; s[2]  = -0.682469858552186; w[2]  = 0.278402959350494;
      r[3]  = -0.901937632421608; s[3]  = -0.901937632421608; w[3]  = 0.063653882657145;
      r[4]  =  0.803875264843215; s[4]  = -0.901937632421608; w[4]  = 0.063653882657145;
      r[5]  = -0.901937632421608; s[5]  =  0.803875264843215; w[5]  = 0.063653882657145;
      r[6]  = -0.696696002590528; s[6]  = -0.696696002590528; w[6]  = 0.168642646759703;
      r[7]  =  0.393392005181056; s[7]  = -0.696696002590528; w[7]  = 0.168642646759703;
      r[8]  = -0.696696002590528; s[8]  =  0.393392005181056; w[8]  = 0.168642646759703;
      r[9]  =  0.348275575962283; s[9]  = -0.375762055856278; w[9]  = 0.077983588949662;
      r[10] = -0.375762055856278; s[10] =  0.348275575962283; w[10] = 0.077983588949662;
      r[11] = -0.972513520106005; s[11] = -0.375762055856278; w[11] = 0.077983588949662;
      r[12] = -0.375762055856278; s[12] = -0.972513520106005; w[12] = 0.077983588949662;
      r[13] = -0.972513520106005; s[13] =  0.348275575962283; w[13] = 0.077983588949662;
      r[14] =  0.348275575962283; s[14] = -0.972513520106005; w[14] = 0.077983588949662;
    
      break;
    }
    case  8: {
      r[0]  = -0.333333333333333; s[0]  = -0.333333333333333; w[0]  = 0.288631215355574;
      r[1]  = -0.081414823414554; s[1]  = -0.081414823414554; w[1]  = 0.190183268534569;
      r[2]  = -0.837170353170893; s[2]  = -0.081414823414554; w[2]  = 0.190183268534569;
      r[3]  = -0.081414823414554; s[3]  = -0.837170353170893; w[3]  = 0.190183268534569;
      r[4]  = -0.658861384496480; s[4]  = -0.658861384496480; w[4]  = 0.206434741069437;
      r[5]  =  0.317722768992959; s[5]  = -0.658861384496480; w[5]  = 0.206434741069437;
      r[6]  = -0.658861384496480; s[6]  =  0.317722768992959; w[6]  = 0.206434741069437;
      r[7]  = -0.898905543365938; s[7]  = -0.898905543365938; w[7]  = 0.064916995246396;
      r[8]  =  0.797811086731876; s[8]  = -0.898905543365938; w[8]  = 0.064916995246396;
      r[9]  = -0.898905543365938; s[9]  =  0.797811086731876; w[9]  = 0.064916995246396;
      r[10] = -0.473774340730724; s[10] =  0.456984785910809; w[10] = 0.054460628348870;
      r[11] =  0.456984785910809; s[11] = -0.473774340730724; w[11] = 0.054460628348870;
      r[12] = -0.983210445180085; s[12] =  0.456984785910809; w[12] = 0.054460628348870;
      r[13] =  0.456984785910809; s[13] = -0.983210445180085; w[13] = 0.054460628348870;
      r[14] = -0.983210445180085; s[14] = -0.473774340730724; w[14] = 0.054460628348870;
      r[15] = -0.473774340730724; s[15] = -0.983210445180085; w[15] = 0.054460628348870;
    
      break;
    }
    case  9: {
      r[0]  = -0.333333333333333; s[0]  = -0.333333333333333; w[0]  = 0.194271592565598;
      r[1]  = -0.020634961602525; s[1]  = -0.020634961602525; w[1]  = 0.062669400454278;
      r[2]  = -0.958730076794950; s[2]  = -0.020634961602525; w[2]  = 0.062669400454278;
      r[3]  = -0.020634961602525; s[3]  = -0.958730076794950; w[3]  = 0.062669400454278;
      r[4]  = -0.125820817014127; s[4]  = -0.125820817014127; w[4]  = 0.155655082009549;
      r[5]  = -0.748358365971746; s[5]  = -0.125820817014127; w[5]  = 0.155655082009549;
      r[6]  = -0.125820817014127; s[6]  = -0.748358365971746; w[6]  = 0.155655082009549;
      r[7]  = -0.623592928761935; s[7]  = -0.623592928761935; w[7]  = 0.159295477854420;
      r[8]  =  0.247185857523869; s[8]  = -0.623592928761935; w[8]  = 0.159295477854420;
      r[9]  = -0.623592928761935; s[9]  =  0.247185857523869; w[9]  = 0.159295477854420;
      r[10] = -0.910540973211095; s[10] = -0.910540973211095; w[10] = 0.051155351317396;
      r[11] =  0.821081946422189; s[11] = -0.910540973211095; w[11] = 0.051155351317396;
      r[12] = -0.910540973211095; s[12] =  0.821081946422189; w[12] = 0.051155351317396;
      r[13] =  0.482397197568996; s[13] = -0.556074021678469; w[13] = 0.086567078754579;
      r[14] = -0.556074021678469; s[14] =  0.482397197568996; w[14] = 0.086567078754579;
      r[15] = -0.926323175890527; s[15] = -0.556074021678469; w[15] = 0.086567078754579;
      r[16] = -0.556074021678469; s[16] = -0.926323175890527; w[16] = 0.086567078754579;
      r[17] = -0.926323175890527; s[17] =  0.482397197568996; w[17] = 0.086567078754579;
      r[18] =  0.482397197568996; s[18] = -0.926323175890527; w[18] = 0.086567078754579;
    
      break;
    }
    case 10: {
      r[0]  = -0.333333333333333; s[0]  = -0.333333333333333; w[0]  = 0.167046799610393;
      r[1]  = -0.004269134091050; s[1]  = -0.004269134091050; w[1]  = 0.014459701184113;
      r[2]  = -0.991461731817899; s[2]  = -0.004269134091050; w[2]  = 0.014459701184113;
      r[3]  = -0.004269134091050; s[3]  = -0.991461731817899; w[3]  = 0.014459701184113;
      r[4]  = -0.143975100541888; s[4]  = -0.143975100541888; w[4]  = 0.148984355841961;
      r[5]  = -0.712049798916225; s[5]  = -0.143975100541888; w[5]  = 0.148984355841961;
      r[6]  = -0.143975100541888; s[6]  = -0.712049798916225; w[6]  = 0.148984355841961;
      r[7]  = -0.630487174513551; s[7]  = -0.630487174513551; w[7]  = 0.157292946806217;
      r[8]  =  0.260974349027102; s[8]  = -0.630487174513551; w[8]  = 0.157292946806217;
      r[9]  = -0.630487174513551; s[9]  =  0.260974349027102; w[9]  = 0.157292946806217;
      r[10] = -0.959037562856645; s[10] = -0.959037562856645; w[10] = 0.013856646174215;
      r[11] =  0.918075125713290; s[11] = -0.959037562856645; w[11] = 0.013856646174215;
      r[12] = -0.959037562856645; s[12] =  0.918075125713290; w[12] = 0.013856646174215;
      r[13] = -0.726852847487933; s[13] =  0.656846867693389; w[13] = 0.059036640669559;
      r[14] =  0.656846867693389; s[14] = -0.726852847487933; w[14] = 0.059036640669559;
      r[15] = -0.929994020205456; s[15] =  0.656846867693389; w[15] = 0.059036640669559;
      r[16] =  0.656846867693389; s[16] = -0.929994020205456; w[16] = 0.059036640669559;
      r[17] = -0.929994020205456; s[17] = -0.726852847487933; w[17] = 0.059036640669559;
      r[18] = -0.726852847487933; s[18] = -0.929994020205456; w[18] = 0.059036640669559;
      r[19] = -0.334512798822723; s[19] =  0.259414658305837; w[19] = 0.079158734392122;
      r[20] =  0.259414658305837; s[20] = -0.334512798822723; w[20] = 0.079158734392122;
      r[21] = -0.924901859483115; s[21] =  0.259414658305837; w[21] = 0.079158734392122;
      r[22] =  0.259414658305837; s[22] = -0.924901859483115; w[22] = 0.079158734392122;
      r[23] = -0.924901859483115; s[23] = -0.334512798822723; w[23] = 0.079158734392122;
      r[24] = -0.334512798822723; s[24] = -0.924901859483115; w[24] = 0.079158734392122;
    
      break;
    }
    case 11: {
      r[0]  = -0.333333333333333; s[0]  = -0.333333333333333; w[0]  = 0.171542269659139;
      r[1]  = -0.008227986519701; s[1]  = -0.008227986519701; w[1]  = 0.033232309648099;
      r[2]  = -0.983544026960598; s[2]  = -0.008227986519701; w[2]  = 0.033232309648099;
      r[3]  = -0.008227986519701; s[3]  = -0.983544026960598; w[3]  = 0.033232309648099;
      r[4]  = -0.123062767323470; s[4]  = -0.123062767323470; w[4]  = 0.134650934561697;
      r[5]  = -0.753874465353061; s[5]  = -0.123062767323470; w[5]  = 0.134650934561697;
      r[6]  = -0.123062767323470; s[6]  = -0.753874465353061; w[6]  = 0.134650934561697;
      r[7]  = -0.579583807314522; s[7]  = -0.579583807314522; w[7]  = 0.141045284213931;
      r[8]  =  0.159167614629044; s[8]  = -0.579583807314522; w[8]  = 0.141045284213931;
      r[9]  = -0.579583807314522; s[9]  =  0.159167614629044; w[9]  = 0.141045284213931;
      r[10] = -0.794793191450311; s[10] = -0.794793191450311; w[10] = 0.077255352409353;
      r[11] =  0.589586382900622; s[11] = -0.794793191450311; w[11] = 0.077255352409353;
      r[12] = -0.794793191450311; s[12] =  0.589586382900622; w[12] = 0.077255352409353;
      r[13] = -0.943045809035533; s[13] = -0.943045809035533; w[13] = 0.020851693311810;
      r[14] =  0.886091618071066; s[14] = -0.943045809035533; w[14] = 0.020851693311810;
      r[15] = -0.943045809035533; s[15] =  0.886091618071066; w[15] = 0.020851693311810;
      r[16] = -0.701412405190965; s[16] =  0.686811031831200; w[16] = 0.020556454130106;
      r[17] =  0.686811031831200; s[17] = -0.701412405190965; w[17] = 0.020556454130106;
      r[18] = -0.985398626640235; s[18] =  0.686811031831200; w[18] = 0.020556454130106;
      r[19] =  0.686811031831200; s[19] = -0.985398626640235; w[19] = 0.020556454130106;
      r[20] = -0.985398626640235; s[20] = -0.701412405190966; w[20] = 0.020556454130106;
      r[21] = -0.701412405190966; s[21] = -0.985398626640235; w[21] = 0.020556454130106;
      r[22] = -0.420944470627256; s[22] =  0.328932876191794; w[22] = 0.080668713854259;
      r[23] =  0.328932876191794; s[23] = -0.420944470627256; w[23] = 0.080668713854259;
      r[24] = -0.907988405564538; s[24] =  0.328932876191794; w[24] = 0.080668713854259;
      r[25] =  0.328932876191794; s[25] = -0.907988405564538; w[25] = 0.080668713854259;
      r[26] = -0.907988405564538; s[26] = -0.420944470627256; w[26] = 0.080668713854259;
      r[27] = -0.420944470627256; s[27] = -0.907988405564538; w[27] = 0.080668713854259;
    
      break;
    }
    case 12: {
      r[0]  = -0.115127288420454; s[0]  = -0.115127288420454; w[0]  = 0.088694050812375;
      r[1]  = -0.769745423159092; s[1]  = -0.115127288420454; w[1]  = 0.088694050812375;
      r[2]  = -0.115127288420454; s[2]  = -0.769745423159092; w[2]  = 0.088694050812375;
      r[3]  = -0.241326933740970; s[3]  = -0.241326933740970; w[3]  = 0.083486638970495;
      r[4]  = -0.517346132518060; s[4]  = -0.241326933740970; w[4]  = 0.083486638970495;
      r[5]  = -0.241326933740970; s[5]  = -0.517346132518060; w[5]  = 0.083486638970495;
      r[6]  = -0.773435586772427; s[6]  = -0.773435586772427; w[6]  = 0.061949370456958;
      r[7]  =  0.546871173544855; s[7]  = -0.773435586772427; w[7]  = 0.061949370456958;
      r[8]  = -0.773435586772427; s[8]  =  0.546871173544855; w[8]  = 0.061949370456958;
      r[9]  = -0.550360151853422; s[9]  = -0.550360151853422; w[9]  = 0.096926744607805;
      r[10] =  0.100720303706844; s[10] = -0.550360151853422; w[10] = 0.096926744607805;
      r[11] = -0.550360151853422; s[11] =  0.100720303706844; w[11] = 0.096926744607805;
      r[12] = -0.949850839686919; s[12] = -0.949850839686919; w[12] = 0.016284715813394;
      r[13] =  0.899701679373837; s[13] = -0.949850839686919; w[13] = 0.016284715813394;
      r[14] = -0.949850839686919; s[14] =  0.899701679373837; w[14] = 0.016284715813394;
      r[15] = -0.023247583079632; s[15] = -0.023247583079632; w[15] = 0.047730162209543;
      r[16] = -0.953504833840737; s[16] = -0.023247583079632; w[16] = 0.047730162209543;
      r[17] = -0.023247583079632; s[17] = -0.953504833840737; w[17] = 0.047730162209543;
      r[18] =  0.701114686742190; s[18] = -0.956341185792092; w[18] = 0.030740478172077;
      r[19] = -0.956341185792092; s[19] =  0.701114686742190; w[19] = 0.030740478172077;
      r[20] = -0.744773500950097; s[20] = -0.956341185792092; w[20] = 0.030740478172077;
      r[21] = -0.956341185792092; s[21] = -0.744773500950097; w[21] = 0.030740478172077;
      r[22] = -0.744773500950097; s[22] =  0.701114686742189; w[22] = 0.030740478172077;
      r[23] =  0.701114686742189; s[23] = -0.744773500950097; w[23] = 0.030740478172077;
      r[24] =  0.280974603820889; s[24] = -0.820754518409649; w[24] = 0.073364295183933;
      r[25] = -0.820754518409649; s[25] =  0.280974603820889; w[25] = 0.073364295183933;
      r[26] = -0.460220085411240; s[26] = -0.820754518409649; w[26] = 0.073364295183933;
      r[27] = -0.820754518409649; s[27] = -0.460220085411240; w[27] = 0.073364295183933;
      r[28] = -0.460220085411240; s[28] =  0.280974603820889; w[28] = 0.073364295183933;
      r[29] =  0.280974603820889; s[29] = -0.460220085411240; w[29] = 0.073364295183933;
      r[30] =  0.380092839755442; s[30] = -0.412599349012509; w[30] = 0.031692718542038;
      r[31] = -0.412599349012509; s[31] =  0.380092839755442; w[31] = 0.031692718542038;
      r[32] = -0.967493490742933; s[32] = -0.412599349012509; w[32] = 0.031692718542038;
      r[33] = -0.412599349012509; s[33] = -0.967493490742933; w[33] = 0.031692718542038;
      r[34] = -0.967493490742933; s[34] =  0.380092839755442; w[34] = 0.031692718542038;
      r[35] =  0.380092839755442; s[35] = -0.967493490742933; w[35] = 0.031692718542038;
    
      break;
    }
    case 13: {
      r[0]  = -0.333333333333333; s[0]  = -0.333333333333333; w[0]  = 0.097649880256618;
      r[1]  = -0.178369512906405; s[1]  = -0.178369512906405; w[1]  = 0.090830392435257;
      r[2]  = -0.643260974187191; s[2]  = -0.178369512906405; w[2]  = 0.090830392435257;
      r[3]  = -0.178369512906405; s[3]  = -0.643260974187191; w[3]  = 0.090830392435257;
      r[4]  = -0.542160866418819; s[4]  = -0.542160866418819; w[4]  = 0.093688172298398;
      r[5]  =  0.084321732837639; s[5]  = -0.542160866418819; w[5]  = 0.093688172298398;
      r[6]  = -0.542160866418819; s[6]  =  0.084321732837639; w[6]  = 0.093688172298398;
      r[7]  = -0.771083702251175; s[7]  = -0.771083702251175; w[7]  = 0.062636625766914;
      r[8]  =  0.542167404502350; s[8]  = -0.771083702251175; w[8]  = 0.062636625766914;
      r[9]  = -0.771083702251175; s[9]  =  0.542167404502350; w[9]  = 0.062636625766914;
      r[10] = -0.951222762153923; s[10] = -0.951222762153923; w[10] = 0.015384123058757;
      r[11] =  0.902445524307847; s[11] = -0.951222762153923; w[11] = 0.015384123058757;
      r[12] = -0.951222762153923; s[12] =  0.902445524307847; w[12] = 0.015384123058757;
      r[13] = -0.071124044893491; s[13] = -0.071124044893491; w[13] = 0.063336411850675;
      r[14] = -0.857751910213018; s[14] = -0.071124044893491; w[14] = 0.063336411850675;
      r[15] = -0.071124044893491; s[15] = -0.857751910213018; w[15] = 0.063336411850675;
      r[16] =  0.089497598920545; s[16] = -0.972567007080596; w[16] = 0.017363112895214;
      r[17] = -0.972567007080596; s[17] =  0.089497598920545; w[17] = 0.017363112895214;
      r[18] = -0.116930591839948; s[18] = -0.972567007080596; w[18] = 0.017363112895214;
      r[19] = -0.972567007080596; s[19] = -0.116930591839948; w[19] = 0.017363112895214;
      r[20] = -0.116930591839948; s[20] =  0.089497598920545; w[20] = 0.017363112895214;
      r[21] =  0.089497598920545; s[21] = -0.116930591839948; w[21] = 0.017363112895214;
      r[22] =  0.708303576893352; s[22] = -0.752827405974001; w[22] = 0.030203519313196;
      r[23] = -0.752827405974001; s[23] =  0.708303576893352; w[23] = 0.030203519313196;
      r[24] = -0.955476170919351; s[24] = -0.752827405974001; w[24] = 0.030203519313196;
      r[25] = -0.752827405974001; s[25] = -0.955476170919351; w[25] = 0.030203519313196;
      r[26] = -0.955476170919351; s[26] =  0.708303576893352; w[26] = 0.030203519313196;
      r[27] =  0.708303576893352; s[27] = -0.955476170919351; w[27] = 0.030203519313196;
      r[28] =  0.400693795223418; s[28] = -0.438009873048077; w[28] = 0.032954120635442;
      r[29] = -0.438009873048077; s[29] =  0.400693795223418; w[29] = 0.032954120635442;
      r[30] = -0.962683922175342; s[30] = -0.438009873048076; w[30] = 0.032954120635442;
      r[31] = -0.438009873048076; s[31] = -0.962683922175342; w[31] = 0.032954120635442;
      r[32] = -0.962683922175342; s[32] =  0.400693795223419; w[32] = 0.032954120635442;
      r[33] =  0.400693795223419; s[33] = -0.962683922175342; w[33] = 0.032954120635442;
      r[34] = -0.461147708177324; s[34] = -0.809703416305633; w[34] = 0.073599737741711;
      r[35] = -0.809703416305633; s[35] = -0.461147708177324; w[35] = 0.073599737741711;
      r[36] =  0.270851124482956; s[36] = -0.809703416305633; w[36] = 0.073599737741711;
      r[37] = -0.809703416305633; s[37] =  0.270851124482956; w[37] = 0.073599737741711;
      r[38] =  0.270851124482956; s[38] = -0.461147708177324; w[38] = 0.073599737741711;
      r[39] = -0.461147708177324; s[39] =  0.270851124482956; w[39] = 0.073599737741711;
    
      break;
    }
    case 14: {
      r[0]  = -0.333333333333333; s[0]  = -0.333333333333333; w[0]  = 0.086016161940337;
      r[1]  = -0.031437623162377; s[1]  = -0.031437623162377; w[1]  = 0.020861042581505;
      r[2]  = -0.937124753675245; s[2]  = -0.031437623162377; w[2]  = 0.020861042581505;
      r[3]  = -0.031437623162377; s[3]  = -0.937124753675245; w[3]  = 0.020861042581505;
      r[4]  = -0.190179067155837; s[4]  = -0.190179067155837; w[4]  = 0.084906577433761;
      r[5]  = -0.619641865688325; s[5]  = -0.190179067155837; w[5]  = 0.084906577433761;
      r[6]  = -0.190179067155837; s[6]  = -0.619641865688325; w[6]  = 0.084906577433761;
      r[7]  = -0.544945017724160; s[7]  = -0.544945017724160; w[7]  = 0.094412955694385;
      r[8]  =  0.089890035448320; s[8]  = -0.544945017724160; w[8]  = 0.094412955694385;
      r[9]  = -0.544945017724160; s[9]  =  0.089890035448320; w[9]  = 0.094412955694385;
      r[10] = -0.829084434528895; s[10] = -0.829084434528895; w[10] = 0.036931420530108;
      r[11] =  0.658168869057790; s[11] = -0.829084434528895; w[11] = 0.036931420530108;
      r[12] = -0.829084434528895; s[12] =  0.658168869057790; w[12] = 0.036931420530108;
      r[13] = -0.961731843229187; s[13] = -0.961731843229187; w[13] = 0.009566859223433;
      r[14] =  0.923463686458374; s[14] = -0.961731843229187; w[14] = 0.009566859223433;
      r[15] = -0.961731843229187; s[15] =  0.923463686458374; w[15] = 0.009566859223433;
      r[16] = -0.228055909294146; s[16] =  0.200257489743243; w[16] = 0.023873332102599;
      r[17] =  0.200257489743243; s[17] = -0.228055909294146; w[17] = 0.023873332102599;
      r[18] = -0.972201580449098; s[18] =  0.200257489743243; w[18] = 0.023873332102599;
      r[19] =  0.200257489743243; s[19] = -0.972201580449098; w[19] = 0.023873332102599;
      r[20] = -0.972201580449097; s[20] = -0.228055909294146; w[20] = 0.023873332102599;
      r[21] = -0.228055909294146; s[21] = -0.972201580449097; w[21] = 0.023873332102599;
      r[22] = -0.801926305237733; s[22] =  0.769518568627872; w[22] = 0.018006541968931;
      r[23] =  0.769518568627872; s[23] = -0.801926305237733; w[23] = 0.018006541968931;
      r[24] = -0.967592263390139; s[24] =  0.769518568627872; w[24] = 0.018006541968931;
      r[25] =  0.769518568627872; s[25] = -0.967592263390139; w[25] = 0.018006541968931;
      r[26] = -0.967592263390139; s[26] = -0.801926305237733; w[26] = 0.018006541968931;
      r[27] = -0.801926305237733; s[27] = -0.967592263390139; w[27] = 0.018006541968931;
      r[28] = -0.544893784145096; s[28] =  0.504355701557747; w[28] = 0.030117913357456;
      r[29] =  0.504355701557747; s[29] = -0.544893784145096; w[29] = 0.030117913357456;
      r[30] = -0.959461917412651; s[30] =  0.504355701557747; w[30] = 0.030117913357456;
      r[31] =  0.504355701557747; s[31] = -0.959461917412651; w[31] = 0.030117913357456;
      r[32] = -0.959461917412651; s[32] = -0.544893784145095; w[32] = 0.030117913357456;
      r[33] = -0.544893784145095; s[33] = -0.959461917412651; w[33] = 0.030117913357456;
      r[34] = -0.288429073961769; s[34] =  0.109399611316741; w[34] = 0.067035335144205;
      r[35] =  0.109399611316741; s[35] = -0.288429073961769; w[35] = 0.067035335144205;
      r[36] = -0.820970537354972; s[36] =  0.109399611316741; w[36] = 0.067035335144205;
      r[37] =  0.109399611316741; s[37] = -0.820970537354972; w[37] = 0.067035335144205;
      r[38] = -0.820970537354972; s[38] = -0.288429073961769; w[38] = 0.067035335144205;
      r[39] = -0.288429073961769; s[39] = -0.820970537354972; w[39] = 0.067035335144205;
      r[40] = -0.593507885280605; s[40] =  0.397658756992448; w[40] = 0.056624756038490;
      r[41] =  0.397658756992448; s[41] = -0.593507885280605; w[41] = 0.056624756038490;
      r[42] = -0.804150871711842; s[42] =  0.397658756992447; w[42] = 0.056624756038490;
      r[43] =  0.397658756992447; s[43] = -0.804150871711842; w[43] = 0.056624756038490;
      r[44] = -0.804150871711842; s[44] = -0.593507885280605; w[44] = 0.056624756038490;
      r[45] = -0.593507885280605; s[45] = -0.804150871711842; w[45] = 0.056624756038490;
    
      break;
    }
    case 15: {
      r[0]  = -0.083438407261750; s[0]  = -0.083438407261750; w[0]  = 0.065323637697611;
      r[1]  = -0.833123185476500; s[1]  = -0.083438407261750; w[1]  = 0.065323637697611;
      r[2]  = -0.083438407261750; s[2]  = -0.833123185476500; w[2]  = 0.065323637697611;
      r[3]  = -0.192779070841739; s[3]  = -0.192779070841739; w[3]  = 0.054825636062729;
      r[4]  = -0.614441858316522; s[4]  = -0.192779070841739; w[4]  = 0.054825636062729;
      r[5]  = -0.192779070841739; s[5]  = -0.614441858316522; w[5]  = 0.054825636062729;
      r[6]  = -0.413605664173949; s[6]  = -0.413605664173949; w[6]  = 0.053020073197407;
      r[7]  = -0.172788671652101; s[7]  = -0.413605664173949; w[7]  = 0.053020073197407;
      r[8]  = -0.413605664173949; s[8]  = -0.172788671652101; w[8]  = 0.053020073197407;
      r[9]  = -0.707064426114454; s[9]  = -0.707064426114454; w[9]  = 0.058431924272972;
      r[10] =  0.414128852228908; s[10] = -0.707064426114454; w[10] = 0.058431924272972;
      r[11] = -0.707064426114454; s[11] =  0.414128852228908; w[11] = 0.058431924272972;
      r[12] = -0.887274264668793; s[12] = -0.887274264668793; w[12] = 0.021169216132488;
      r[13] =  0.774548529337586; s[13] = -0.887274264668793; w[13] = 0.021169216132488;
      r[14] = -0.887274264668793; s[14] =  0.774548529337586; w[14] = 0.021169216132488;
      r[15] = -0.966849746283259; s[15] = -0.966849746283259; w[15] = 0.007229286128184;
      r[16] =  0.933699492566519; s[16] = -0.966849746283259; w[16] = 0.007229286128184;
      r[17] = -0.966849746283259; s[17] =  0.933699492566519; w[17] = 0.007229286128184;
      r[18] = -0.520930891690411; s[18] =  0.501106485071962; w[18] = 0.017055496203419;
      r[19] =  0.501106485071962; s[19] = -0.520930891690411; w[19] = 0.017055496203419;
      r[20] = -0.980175593381550; s[20] =  0.501106485071962; w[20] = 0.017055496203419;
      r[21] =  0.501106485071962; s[21] = -0.980175593381550; w[21] = 0.017055496203419;
      r[22] = -0.980175593381550; s[22] = -0.520930891690411; w[22] = 0.017055496203419;
      r[23] = -0.520930891690411; s[23] = -0.980175593381550; w[23] = 0.017055496203419;
      r[24] = -0.190242385363320; s[24] =  0.158634844102864; w[24] = 0.027832353033384;
      r[25] =  0.158634844102864; s[25] = -0.190242385363320; w[25] = 0.027832353033384;
      r[26] = -0.968392458739544; s[26] =  0.158634844102864; w[26] = 0.027832353033384;
      r[27] =  0.158634844102864; s[27] = -0.968392458739544; w[27] = 0.027832353033384;
      r[28] = -0.968392458739544; s[28] = -0.190242385363320; w[28] = 0.027832353033384;
      r[29] = -0.190242385363320; s[29] = -0.968392458739544; w[29] = 0.027832353033384;
      r[30] = -0.809995773773910; s[30] =  0.799708556139969; w[30] = 0.008583865881470;
      r[31] =  0.799708556139969; s[31] = -0.809995773773910; w[31] = 0.008583865881470;
      r[32] = -0.989712782366059; s[32] =  0.799708556139969; w[32] = 0.008583865881470;
      r[33] =  0.799708556139969; s[33] = -0.989712782366059; w[33] = 0.008583865881470;
      r[34] = -0.989712782366059; s[34] = -0.809995773773910; w[34] = 0.008583865881470;
      r[35] = -0.809995773773910; s[35] = -0.989712782366059; w[35] = 0.008583865881470;
      r[36] = -0.700493785355452; s[36] =  0.602649133849474; w[36] = 0.032470658563550;
      r[37] =  0.602649133849474; s[37] = -0.700493785355452; w[37] = 0.032470658563550;
      r[38] = -0.902155348494022; s[38] =  0.602649133849474; w[38] = 0.032470658563550;
      r[39] =  0.602649133849474; s[39] = -0.902155348494022; w[39] = 0.032470658563550;
      r[40] = -0.902155348494023; s[40] = -0.700493785355452; w[40] = 0.032470658563550;
      r[41] = -0.700493785355452; s[41] = -0.902155348494023; w[41] = 0.032470658563550;
      r[42] = -0.426160775117330; s[42] =  0.288623277852292; w[42] = 0.051214681842525;
      r[43] =  0.288623277852292; s[43] = -0.426160775117330; w[43] = 0.051214681842525;
      r[44] = -0.862462502734962; s[44] =  0.288623277852292; w[44] = 0.051214681842525;
      r[45] =  0.288623277852292; s[45] = -0.862462502734962; w[45] = 0.051214681842525;
      r[46] = -0.862462502734962; s[46] = -0.426160775117330; w[46] = 0.051214681842525;
      r[47] = -0.426160775117330; s[47] = -0.862462502734962; w[47] = 0.051214681842525;
      r[48] = -0.436328663801831; s[48] =  0.099519827552432; w[48] = 0.066176391063291;
      r[49] =  0.099519827552432; s[49] = -0.436328663801831; w[49] = 0.066176391063291;
      r[50] = -0.663191163750602; s[50] =  0.099519827552432; w[50] = 0.066176391063291;
      r[51] =  0.099519827552432; s[51] = -0.663191163750602; w[51] = 0.066176391063291;
      r[52] = -0.663191163750602; s[52] = -0.436328663801831; w[52] = 0.066176391063291;
      r[53] = -0.436328663801831; s[53] = -0.663191163750602; w[53] = 0.066176391063291;
    
      break;
    }
    case 16: {
      r[0]  = -0.333333333333333; s[0]  = -0.333333333333333; w[0]  = 0.092421202321849;
      r[1]  = -0.015796436959213; s[1]  = -0.015796436959213; w[1]  = 0.027991652734724;
      r[2]  = -0.968407126081574; s[2]  = -0.015796436959213; w[2]  = 0.027991652734724;
      r[3]  = -0.015796436959213; s[3]  = -0.968407126081574; w[3]  = 0.027991652734724;
      r[4]  = -0.087376299046785; s[4]  = -0.087376299046785; w[4]  = 0.036406763596424;
      r[5]  = -0.825247401906431; s[5]  = -0.087376299046785; w[5]  = 0.036406763596424;
      r[6]  = -0.087376299046785; s[6]  = -0.825247401906431; w[6]  = 0.036406763596424;
      r[7]  = -0.640498098534045; s[7]  = -0.640498098534045; w[7]  = 0.062729423327251;
      r[8]  =  0.280996197068091; s[8]  = -0.640498098534045; w[8]  = 0.062729423327251;
      r[9]  = -0.640498098534045; s[9]  =  0.280996197068091; w[9]  = 0.062729423327251;
      r[10] = -0.828211591878549; s[10] = -0.828211591878549; w[10] = 0.032217170192673;
      r[11] =  0.656423183757099; s[11] = -0.828211591878549; w[11] = 0.032217170192673;
      r[12] = -0.828211591878549; s[12] =  0.656423183757099; w[12] = 0.032217170192673;
      r[13] = -0.977882767295650; s[13] = -0.977882767295650; w[13] = 0.003676792575648;
      r[14] =  0.955765534591299; s[14] = -0.977882767295650; w[14] = 0.003676792575648;
      r[15] = -0.977882767295650; s[15] =  0.955765534591299; w[15] = 0.003676792575648;
      r[16] = -0.648594849699253; s[16] =  0.619534954455185; w[16] = 0.019075686313307;
      r[17] =  0.619534954455185; s[17] = -0.648594849699253; w[17] = 0.019075686313307;
      r[18] = -0.970940104755933; s[18] =  0.619534954455185; w[18] = 0.019075686313307;
      r[19] =  0.619534954455185; s[19] = -0.970940104755933; w[19] = 0.019075686313307;
      r[20] = -0.970940104755933; s[20] = -0.648594849699253; w[20] = 0.019075686313307;
      r[21] = -0.648594849699253; s[21] = -0.970940104755933; w[21] = 0.019075686313307;
      r[22] = -0.355089058161032; s[22] =  0.324703448463838; w[22] = 0.025111336723323;
      r[23] =  0.324703448463838; s[23] = -0.355089058161032; w[23] = 0.025111336723323;
      r[24] = -0.969614390302806; s[24] =  0.324703448463838; w[24] = 0.025111336723323;
      r[25] =  0.324703448463838; s[25] = -0.969614390302806; w[25] = 0.025111336723323;
      r[26] = -0.969614390302806; s[26] = -0.355089058161032; w[26] = 0.025111336723323;
      r[27] = -0.355089058161032; s[27] = -0.969614390302806; w[27] = 0.025111336723323;
      r[28] = -0.863082725186344; s[28] =  0.829810165174004; w[28] = 0.014134510743040;
      r[29] =  0.829810165174004; s[29] = -0.863082725186344; w[29] = 0.014134510743040;
      r[30] = -0.966727439987660; s[30] =  0.829810165174004; w[30] = 0.014134510743040;
      r[31] =  0.829810165174004; s[31] = -0.966727439987660; w[31] = 0.014134510743040;
      r[32] = -0.966727439987660; s[32] = -0.863082725186344; w[32] = 0.014134510743040;
      r[33] = -0.863082725186344; s[33] = -0.966727439987660; w[33] = 0.014134510743040;
      r[34] = -0.650110471151428; s[34] =  0.511910525190604; w[34] = 0.025255543672496;
      r[35] =  0.511910525190604; s[35] = -0.650110471151428; w[35] = 0.025255543672496;
      r[36] = -0.861800054039176; s[36] =  0.511910525190604; w[36] = 0.025255543672496;
      r[37] =  0.511910525190604; s[37] = -0.861800054039176; w[37] = 0.025255543672496;
      r[38] = -0.861800054039176; s[38] = -0.650110471151428; w[38] = 0.025255543672496;
      r[39] = -0.650110471151428; s[39] = -0.861800054039176; w[39] = 0.025255543672496;
      r[40] = -0.282624366842035; s[40] =  0.131065359625217; w[40] = 0.038432852520849;
      r[41] =  0.131065359625217; s[41] = -0.282624366842035; w[41] = 0.038432852520849;
      r[42] = -0.848440992783182; s[42] =  0.131065359625217; w[42] = 0.038432852520849;
      r[43] =  0.131065359625217; s[43] = -0.848440992783182; w[43] = 0.038432852520849;
      r[44] = -0.848440992783182; s[44] = -0.282624366842035; w[44] = 0.038432852520849;
      r[45] = -0.282624366842035; s[45] = -0.848440992783182; w[45] = 0.038432852520849;
      r[46] = -0.500121711290305; s[46] =  0.333566061453817; w[46] = 0.032993291211606;
      r[47] =  0.333566061453817; s[47] = -0.500121711290305; w[47] = 0.032993291211606;
      r[48] = -0.833444350163512; s[48] =  0.333566061453817; w[48] = 0.032993291211606;
      r[49] =  0.333566061453817; s[49] = -0.833444350163512; w[49] = 0.032993291211606;
      r[50] = -0.833444350163512; s[50] = -0.500121711290305; w[50] = 0.032993291211606;
      r[51] = -0.500121711290305; s[51] = -0.833444350163512; w[51] = 0.032993291211606;
      r[52] = -0.353308953197778; s[52] = -0.027538833736648; w[52] = 0.081415677215044;
      r[53] = -0.027538833736648; s[53] = -0.353308953197778; w[53] = 0.081415677215044;
      r[54] = -0.619152213065573; s[54] = -0.027538833736648; w[54] = 0.081415677215044;
      r[55] = -0.027538833736648; s[55] = -0.619152213065573; w[55] = 0.081415677215044;
      r[56] = -0.619152213065573; s[56] = -0.353308953197778; w[56] = 0.081415677215044;
      r[57] = -0.353308953197778; s[57] = -0.619152213065573; w[57] = 0.081415677215044;
    
      break;
    }
    case 17: {
      r[0]  = -0.013565183960096; s[0]  = -0.013565183960096; w[0]  = 0.022448424077156;
      r[1]  = -0.972869632079808; s[1]  = -0.013565183960096; w[1]  = 0.022448424077156;
      r[2]  = -0.013565183960096; s[2]  = -0.972869632079808; w[2]  = 0.022448424077156;
      r[3]  = -0.918012666387151; s[3]  = -0.918012666387151; w[3]  = 0.011377816766740;
      r[4]  =  0.836025332774302; s[4]  = -0.918012666387151; w[4]  = 0.011377816766740;
      r[5]  = -0.918012666387151; s[5]  =  0.836025332774302; w[5]  = 0.011377816766740;
      r[6]  = -0.469265478050798; s[6]  = -0.469265478050798; w[6]  = 0.056148078797670;
      r[7]  = -0.061469043898404; s[7]  = -0.469265478050798; w[7]  = 0.056148078797670;
      r[8]  = -0.469265478050798; s[8]  = -0.061469043898404; w[8]  = 0.056148078797670;
      r[9]  = -0.278161974512876; s[9]  = -0.278161974512876; w[9]  = 0.031804753050686;
      r[10] = -0.443676050974249; s[10] = -0.278161974512876; w[10] = 0.031804753050686;
      r[11] = -0.278161974512876; s[11] = -0.443676050974249; w[11] = 0.031804753050686;
      r[12] = -0.743831619059579; s[12] = -0.743831619059579; w[12] = 0.042944604159838;
      r[13] =  0.487663238119159; s[13] = -0.743831619059579; w[13] = 0.042944604159838;
      r[14] = -0.743831619059579; s[14] =  0.487663238119159; w[14] = 0.042944604159838;
      r[15] = -0.975086972700342; s[15] = -0.975086972700342; w[15] = 0.004138336267131;
      r[16] =  0.950173945400684; s[16] = -0.975086972700342; w[16] = 0.004138336267131;
      r[17] = -0.975086972700342; s[17] =  0.950173945400684; w[17] = 0.004138336267131;
      r[18] = -0.331485915788080; s[18] =  0.305875449252880; w[18] = 0.020094111670483;
      r[19] =  0.305875449252880; s[19] = -0.331485915788080; w[19] = 0.020094111670483;
      r[20] = -0.974389533464800; s[20] =  0.305875449252880; w[20] = 0.020094111670483;
      r[21] =  0.305875449252880; s[21] = -0.974389533464800; w[21] = 0.020094111670483;
      r[22] = -0.974389533464800; s[22] = -0.331485915788080; w[22] = 0.020094111670483;
      r[23] = -0.331485915788080; s[23] = -0.974389533464800; w[23] = 0.020094111670483;
      r[24] = -0.617368143234940; s[24] =  0.598534426856625; w[24] = 0.013023443943004;
      r[25] =  0.598534426856625; s[25] = -0.617368143234940; w[25] = 0.013023443943004;
      r[26] = -0.981166283621686; s[26] =  0.598534426856625; w[26] = 0.013023443943004;
      r[27] =  0.598534426856625; s[27] = -0.981166283621686; w[27] = 0.013023443943004;
      r[28] = -0.981166283621686; s[28] = -0.617368143234940; w[28] = 0.013023443943004;
      r[29] = -0.617368143234940; s[29] = -0.981166283621686; w[29] = 0.013023443943004;
      r[30] = -0.846737941679184; s[30] =  0.837351785097968; w[30] = 0.005856954369545;
      r[31] =  0.837351785097968; s[31] = -0.846737941679184; w[31] = 0.005856954369545;
      r[32] = -0.990613843418784; s[32] =  0.837351785097968; w[32] = 0.005856954369545;
      r[33] =  0.837351785097968; s[33] = -0.990613843418784; w[33] = 0.005856954369545;
      r[34] = -0.990613843418784; s[34] = -0.846737941679184; w[34] = 0.005856954369545;
      r[35] = -0.846737941679184; s[35] = -0.990613843418784; w[35] = 0.005856954369545;
      r[36] = -0.530582382227094; s[36] =  0.412831673770520; w[36] = 0.037746172643499;
      r[37] =  0.412831673770520; s[37] = -0.530582382227094; w[37] = 0.037746172643499;
      r[38] = -0.882249291543426; s[38] =  0.412831673770520; w[38] = 0.037746172643499;
      r[39] =  0.412831673770520; s[39] = -0.882249291543426; w[39] = 0.037746172643499;
      r[40] = -0.882249291543426; s[40] = -0.530582382227094; w[40] = 0.037746172643499;
      r[41] = -0.530582382227094; s[41] = -0.882249291543426; w[41] = 0.037746172643499;
      r[42] = -0.230167008731505; s[42] =  0.092655855430620; w[42] = 0.048362752988591;
      r[43] =  0.092655855430620; s[43] = -0.230167008731505; w[43] = 0.048362752988591;
      r[44] = -0.862488846699115; s[44] =  0.092655855430620; w[44] = 0.048362752988591;
      r[45] =  0.092655855430620; s[45] = -0.862488846699115; w[45] = 0.048362752988591;
      r[46] = -0.862488846699115; s[46] = -0.230167008731505; w[46] = 0.048362752988591;
      r[47] = -0.230167008731505; s[47] = -0.862488846699115; w[47] = 0.048362752988591;
      r[48] = -0.769841956364772; s[48] =  0.681416462088999; w[48] = 0.023467126169388;
      r[49] =  0.681416462088999; s[49] = -0.769841956364772; w[49] = 0.023467126169388;
      r[50] = -0.911574505724228; s[50] =  0.681416462088999; w[50] = 0.023467126169388;
      r[51] =  0.681416462088999; s[51] = -0.911574505724228; w[51] = 0.023467126169388;
      r[52] = -0.911574505724228; s[52] = -0.769841956364772; w[52] = 0.023467126169388;
      r[53] = -0.769841956364772; s[53] = -0.911574505724228; w[53] = 0.023467126169388;
      r[54] = -0.502901760016855; s[54] =  0.206946937904245; w[54] = 0.054703428836887;
      r[55] =  0.206946937904245; s[55] = -0.502901760016855; w[55] = 0.054703428836887;
      r[56] = -0.704045177887389; s[56] =  0.206946937904245; w[56] = 0.054703428836887;
      r[57] =  0.206946937904245; s[57] = -0.704045177887389; w[57] = 0.054703428836887;
      r[58] = -0.704045177887389; s[58] = -0.502901760016855; w[58] = 0.054703428836887;
      r[59] = -0.502901760016855; s[59] = -0.704045177887389; w[59] = 0.054703428836887;
      r[60] = -0.258115377755090; s[60] = -0.067105414124881; w[60] = 0.045648336152327;
      r[61] = -0.067105414124881; s[61] = -0.258115377755090; w[61] = 0.045648336152327;
      r[62] = -0.674779208120030; s[62] = -0.067105414124880; w[62] = 0.045648336152327;
      r[63] = -0.067105414124880; s[63] = -0.674779208120030; w[63] = 0.045648336152327;
      r[64] = -0.674779208120030; s[64] = -0.258115377755090; w[64] = 0.045648336152327;
      r[65] = -0.258115377755090; s[65] = -0.674779208120030; w[65] = 0.045648336152327;
    
      break;
    }
    case 18: {
      r[0]  = -0.333333333333333; s[0]  = -0.333333333333333; w[0]  = 0.044365924378076;
      r[1]  = -0.012131216528029; s[1]  = -0.012131216528029; w[1]  = 0.019118866558615;
      r[2]  = -0.975737566943942; s[2]  = -0.012131216528029; w[2]  = 0.019118866558615;
      r[3]  = -0.012131216528029; s[3]  = -0.975737566943942; w[3]  = 0.019118866558615;
      r[4]  = -0.230136518749796; s[4]  = -0.230136518749796; w[4]  = 0.049294282264207;
      r[5]  = -0.539726962500408; s[5]  = -0.230136518749796; w[5]  = 0.049294282264207;
      r[6]  = -0.230136518749796; s[6]  = -0.539726962500408; w[6]  = 0.049294282264207;
      r[7]  = -0.499402807850562; s[7]  = -0.499402807850562; w[7]  = 0.060654125981259;
      r[8]  = -0.001194384298875; s[8]  = -0.499402807850563; w[8]  = 0.060654125981259;
      r[9]  = -0.499402807850563; s[9]  = -0.001194384298875; w[9]  = 0.060654125981259;
      r[10] = -0.710045241598598; s[10] = -0.710045241598598; w[10] = 0.031366231844497;
      r[11] =  0.420090483197196; s[11] = -0.710045241598598; w[11] = 0.031366231844497;
      r[12] = -0.710045241598598; s[12] =  0.420090483197196; w[12] = 0.031366231844497;
      r[13] = -0.903019171569391; s[13] = -0.903019171569391; w[13] = 0.016173035699059;
      r[14] =  0.806038343138783; s[14] = -0.903019171569391; w[14] = 0.016173035699059;
      r[15] = -0.903019171569391; s[15] =  0.806038343138783; w[15] = 0.016173035699059;
      r[16] = -0.972866921676664; s[16] = -0.972866921676664; w[16] = 0.004750076743416;
      r[17] =  0.945733843353328; s[17] = -0.972866921676664; w[17] = 0.004750076743416;
      r[18] = -0.972866921676664; s[18] =  0.945733843353328; w[18] = 0.004750076743416;
      r[19] = -0.859848143426920; s[19] =  0.852682278520227; w[19] = 0.004329998364046;
      r[20] =  0.852682278520227; s[20] = -0.859848143426920; w[20] = 0.004329998364046;
      r[21] = -0.992834135093306; s[21] =  0.852682278520227; w[21] = 0.004329998364046;
      r[22] =  0.852682278520227; s[22] = -0.992834135093306; w[22] = 0.004329998364046;
      r[23] = -0.992834135093306; s[23] = -0.859848143426920; w[23] = 0.004329998364046;
      r[24] = -0.859848143426920; s[24] = -0.992834135093306; w[24] = 0.004329998364046;
      r[25] = -0.592208428296456; s[25] =  0.573110619333459; w[25] = 0.012300136592151;
      r[26] =  0.573110619333459; s[26] = -0.592208428296456; w[26] = 0.012300136592151;
      r[27] = -0.980902191037003; s[27] =  0.573110619333459; w[27] = 0.012300136592151;
      r[28] =  0.573110619333459; s[28] = -0.980902191037003; w[28] = 0.012300136592151;
      r[29] = -0.980902191037003; s[29] = -0.592208428296456; w[29] = 0.012300136592151;
      r[30] = -0.592208428296456; s[30] = -0.980902191037003; w[30] = 0.012300136592151;
      r[31] = -0.316660705097116; s[31] =  0.294011442575493; w[31] = 0.017146550496094;
      r[32] =  0.294011442575493; s[32] = -0.316660705097116; w[32] = 0.017146550496094;
      r[33] = -0.977350737478378; s[33] =  0.294011442575493; w[33] = 0.017146550496094;
      r[34] =  0.294011442575493; s[34] = -0.977350737478378; w[34] = 0.017146550496094;
      r[35] = -0.977350737478378; s[35] = -0.316660705097116; w[35] = 0.017146550496094;
      r[36] = -0.316660705097116; s[36] = -0.977350737478378; w[36] = 0.017146550496094;
      r[37] = -0.762708275393600; s[37] =  0.720805143081058; w[37] = 0.010992784367686;
      r[38] =  0.720805143081058; s[38] = -0.762708275393600; w[38] = 0.010992784367686;
      r[39] = -0.958096867687458; s[39] =  0.720805143081058; w[39] = 0.010992784367686;
      r[40] =  0.720805143081058; s[40] = -0.958096867687458; w[40] = 0.010992784367686;
      r[41] = -0.958096867687458; s[41] = -0.762708275393600; w[41] = 0.010992784367686;
      r[42] = -0.762708275393600; s[42] = -0.958096867687458; w[42] = 0.010992784367686;
      r[43] = -0.499464753852078; s[43] =  0.389673957377540; w[43] = 0.033988423627451;
      r[44] =  0.389673957377540; s[44] = -0.499464753852078; w[44] = 0.033988423627451;
      r[45] = -0.890209203525462; s[45] =  0.389673957377540; w[45] = 0.033988423627451;
      r[46] =  0.389673957377540; s[46] = -0.890209203525462; w[46] = 0.033988423627451;
      r[47] = -0.890209203525462; s[47] = -0.499464753852078; w[47] = 0.033988423627451;
      r[48] = -0.499464753852078; s[48] = -0.890209203525462; w[48] = 0.033988423627451;
      r[49] = -0.730193961187213; s[49] =  0.595225312221749; w[49] = 0.025768775250048;
      r[50] =  0.595225312221749; s[50] = -0.730193961187213; w[50] = 0.025768775250048;
      r[51] = -0.865031351034536; s[51] =  0.595225312221749; w[51] = 0.025768775250048;
      r[52] =  0.595225312221749; s[52] = -0.865031351034536; w[52] = 0.025768775250048;
      r[53] = -0.865031351034535; s[53] = -0.730193961187213; w[53] = 0.025768775250048;
      r[54] = -0.730193961187213; s[54] = -0.865031351034535; w[54] = 0.025768775250048;
      r[55] = -0.213662181077459; s[55] =  0.092159471784540; w[55] = 0.040023698716305;
      r[56] =  0.092159471784540; s[56] = -0.213662181077459; w[56] = 0.040023698716305;
      r[57] = -0.878497290707081; s[57] =  0.092159471784540; w[57] = 0.040023698716305;
      r[58] =  0.092159471784540; s[58] = -0.878497290707081; w[58] = 0.040023698716305;
      r[59] = -0.878497290707081; s[59] = -0.213662181077459; w[59] = 0.040023698716305;
      r[60] = -0.213662181077459; s[60] = -0.878497290707081; w[60] = 0.040023698716305;
      r[61] = -0.268906542375276; s[61] = -0.010964677486621; w[61] = 0.045277248097616;
      r[62] = -0.010964677486621; s[62] = -0.268906542375276; w[62] = 0.045277248097616;
      r[63] = -0.720128780138103; s[63] = -0.010964677486621; w[63] = 0.045277248097616;
      r[64] = -0.010964677486621; s[64] = -0.720128780138103; w[64] = 0.045277248097616;
      r[65] = -0.720128780138103; s[65] = -0.268906542375276; w[65] = 0.045277248097616;
      r[66] = -0.268906542375276; s[66] = -0.720128780138103; w[66] = 0.045277248097616;
      r[67] = -0.507963230765286; s[67] =  0.233062121320884; w[67] = 0.045433087546729;
      r[68] =  0.233062121320884; s[68] = -0.507963230765286; w[68] = 0.045433087546729;
      r[69] = -0.725098890555598; s[69] =  0.233062121320884; w[69] = 0.045433087546729;
      r[70] =  0.233062121320884; s[70] = -0.725098890555598; w[70] = 0.045433087546729;
      r[71] = -0.725098890555598; s[71] = -0.507963230765286; w[71] = 0.045433087546729;
      r[72] = -0.507963230765286; s[72] = -0.725098890555598; w[72] = 0.045433087546729;
    
      break;
    }
    case 19: {
      r[0]  = -0.333333333333333; s[0]  = -0.333333333333333; w[0]  = 0.041424760048010;
      r[1]  = -0.121193890022463; s[1]  = -0.121193890022463; w[1]  = 0.031784855337689;
      r[2]  = -0.757612219955075; s[2]  = -0.121193890022463; w[2]  = 0.031784855337689;
      r[3]  = -0.121193890022463; s[3]  = -0.757612219955075; w[3]  = 0.031784855337689;
      r[4]  = -0.227921589064499; s[4]  = -0.227921589064499; w[4]  = 0.048406599909930;
      r[5]  = -0.544156821871001; s[5]  = -0.227921589064499; w[5]  = 0.048406599909930;
      r[6]  = -0.227921589064499; s[6]  = -0.544156821871001; w[6]  = 0.048406599909930;
      r[7]  = -0.480832787187258; s[7]  = -0.480832787187258; w[7]  = 0.048177661858475;
      r[8]  = -0.038334425625484; s[8]  = -0.480832787187258; w[8]  = 0.048177661858475;
      r[9]  = -0.480832787187258; s[9]  = -0.038334425625484; w[9]  = 0.048177661858475;
      r[10] = -0.607587000462791; s[10] = -0.607587000462791; w[10] = 0.042691162490028;
      r[11] =  0.215174000925582; s[11] = -0.607587000462791; w[11] = 0.042691162490028;
      r[12] = -0.607587000462791; s[12] =  0.215174000925582; w[12] = 0.042691162490028;
      r[13] = -0.728032884013154; s[13] = -0.728032884013154; w[13] = 0.031889013978028;
      r[14] =  0.456065768026307; s[14] = -0.728032884013154; w[14] = 0.031889013978028;
      r[15] = -0.728032884013154; s[15] =  0.456065768026307; w[15] = 0.031889013978028;
      r[16] = -0.884999963315497; s[16] = -0.884999963315497; w[16] = 0.015529272923518;
      r[17] =  0.769999926630994; s[17] = -0.884999963315497; w[17] = 0.015529272923518;
      r[18] = -0.884999963315497; s[18] =  0.769999926630994; w[18] = 0.015529272923518;
      r[19] = -0.973927705240279; s[19] = -0.973927705240279; w[19] = 0.004442852871415;
      r[20] =  0.947855410480559; s[20] = -0.973927705240279; w[20] = 0.004442852871415;
      r[21] = -0.973927705240279; s[21] =  0.947855410480559; w[21] = 0.004442852871415;
      r[22] = -0.419340397282846; s[22] =  0.418549313572482; w[22] = 0.004599409024561;
      r[23] =  0.418549313572482; s[23] = -0.419340397282846; w[23] = 0.004599409024561;
      r[24] = -0.999208916289636; s[24] =  0.418549313572482; w[24] = 0.004599409024561;
      r[25] =  0.418549313572482; s[25] = -0.999208916289636; w[25] = 0.004599409024561;
      r[26] = -0.999208916289636; s[26] = -0.419340397282846; w[26] = 0.004599409024561;
      r[27] = -0.419340397282846; s[27] = -0.999208916289636; w[27] = 0.004599409024561;
      r[28] = -0.154285565801825; s[28] =  0.131285016869205; w[28] = 0.016608327001253;
      r[29] =  0.131285016869205; s[29] = -0.154285565801825; w[29] = 0.016608327001253;
      r[30] = -0.976999451067380; s[30] =  0.131285016869205; w[30] = 0.016608327001253;
      r[31] =  0.131285016869205; s[31] = -0.976999451067380; w[31] = 0.016608327001253;
      r[32] = -0.976999451067380; s[32] = -0.154285565801825; w[32] = 0.016608327001253;
      r[33] = -0.154285565801825; s[33] = -0.976999451067380; w[33] = 0.016608327001253;
      r[34] = -0.864328042937615; s[34] =  0.842262778125134; w[34] = 0.008508350664254;
      r[35] =  0.842262778125134; s[35] = -0.864328042937615; w[35] = 0.008508350664254;
      r[36] = -0.977934735187519; s[36] =  0.842262778125134; w[36] = 0.008508350664254;
      r[37] =  0.842262778125134; s[37] = -0.977934735187519; w[37] = 0.008508350664254;
      r[38] = -0.977934735187519; s[38] = -0.864328042937615; w[38] = 0.008508350664254;
      r[39] = -0.864328042937615; s[39] = -0.977934735187519; w[39] = 0.008508350664254;
      r[40] = -0.676980679779850; s[40] =  0.654351228666662; w[40] = 0.012752760516072;
      r[41] =  0.654351228666662; s[41] = -0.676980679779850; w[41] = 0.012752760516072;
      r[42] = -0.977370548886812; s[42] =  0.654351228666662; w[42] = 0.012752760516072;
      r[43] =  0.654351228666662; s[43] = -0.977370548886812; w[43] = 0.012752760516072;
      r[44] = -0.977370548886812; s[44] = -0.676980679779850; w[44] = 0.012752760516072;
      r[45] = -0.676980679779850; s[45] = -0.977370548886812; w[45] = 0.012752760516072;
      r[46] = -0.451038034994882; s[46] =  0.391933493150074; w[46] = 0.024864731137205;
      r[47] =  0.391933493150074; s[47] = -0.451038034994882; w[47] = 0.024864731137205;
      r[48] = -0.940895458155192; s[48] =  0.391933493150074; w[48] = 0.024864731137205;
      r[49] =  0.391933493150074; s[49] = -0.940895458155192; w[49] = 0.024864731137205;
      r[50] = -0.940895458155192; s[50] = -0.451038034994882; w[50] = 0.024864731137205;
      r[51] = -0.451038034994882; s[51] = -0.940895458155192; w[51] = 0.024864731137205;
      r[52] = -0.755079003035426; s[52] =  0.643775758022122; w[52] = 0.012779469878381;
      r[53] =  0.643775758022122; s[53] = -0.755079003035426; w[53] = 0.012779469878381;
      r[54] = -0.888696754986697; s[54] =  0.643775758022122; w[54] = 0.012779469878381;
      r[55] =  0.643775758022122; s[55] = -0.888696754986697; w[55] = 0.012779469878381;
      r[56] = -0.888696754986697; s[56] = -0.755079003035426; w[56] = 0.012779469878381;
      r[57] = -0.755079003035426; s[57] = -0.888696754986697; w[57] = 0.012779469878381;
      r[58] = -0.668500025831376; s[58] =  0.544622467980326; w[58] = 0.017801502950615;
      r[59] =  0.544622467980326; s[59] = -0.668500025831376; w[59] = 0.017801502950615;
      r[60] = -0.876122442148950; s[60] =  0.544622467980326; w[60] = 0.017801502950615;
      r[61] =  0.544622467980326; s[61] = -0.876122442148950; w[61] = 0.017801502950615;
      r[62] = -0.876122442148950; s[62] = -0.668500025831376; w[62] = 0.017801502950615;
      r[63] = -0.668500025831376; s[63] = -0.876122442148950; w[63] = 0.017801502950615;
      r[64] = -0.205954158284095; s[64] =  0.090524206162520; w[64] = 0.035738233266717;
      r[65] =  0.090524206162520; s[65] = -0.205954158284095; w[65] = 0.035738233266717;
      r[66] = -0.884570047878425; s[66] =  0.090524206162520; w[66] = 0.035738233266717;
      r[67] =  0.090524206162520; s[67] = -0.884570047878425; w[67] = 0.035738233266717;
      r[68] = -0.884570047878425; s[68] = -0.205954158284095; w[68] = 0.035738233266717;
      r[69] = -0.205954158284095; s[69] = -0.884570047878425; w[69] = 0.035738233266717;
      r[70] = -0.488476686471250; s[70] =  0.295757669634562; w[70] = 0.038344199404453;
      r[71] =  0.295757669634562; s[71] = -0.488476686471250; w[71] = 0.038344199404453;
      r[72] = -0.807280983163311; s[72] =  0.295757669634562; w[72] = 0.038344199404453;
      r[73] =  0.295757669634562; s[73] = -0.807280983163311; w[73] = 0.038344199404453;
      r[74] = -0.807280983163311; s[74] = -0.488476686471250; w[74] = 0.038344199404453;
      r[75] = -0.488476686471250; s[75] = -0.807280983163311; w[75] = 0.038344199404453;
      r[76] = -0.342507806805715; s[76] =  0.042214391480752; w[76] = 0.042971513130613;
      r[77] =  0.042214391480752; s[77] = -0.342507806805715; w[77] = 0.042971513130613;
      r[78] = -0.699706584675037; s[78] =  0.042214391480751; w[78] = 0.042971513130613;
      r[79] =  0.042214391480751; s[79] = -0.699706584675037; w[79] = 0.042971513130613;
      r[80] = -0.699706584675037; s[80] = -0.342507806805715; w[80] = 0.042971513130613;
      r[81] = -0.342507806805715; s[81] = -0.699706584675037; w[81] = 0.042971513130613;
    
      break;
    }
    case 20: {
      r[0]  = -0.333333333333333; s[0]  = -0.333333333333333; w[0]  = 0.055220853995399;
      r[1]  = -0.001500649324429; s[1]  = -0.001500649324429; w[1]  = 0.003558059094653;
      r[2]  = -0.996998701351142; s[2]  = -0.001500649324429; w[2]  = 0.003558059094653;
      r[3]  = -0.001500649324429; s[3]  = -0.996998701351142; w[3]  = 0.003558059094653;
      r[4]  = -0.094139751938951; s[4]  = -0.094139751938951; w[4]  = 0.040224796227922;
      r[5]  = -0.811720496122098; s[5]  = -0.094139751938951; w[5]  = 0.040224796227922;
      r[6]  = -0.094139751938951; s[6]  = -0.811720496122098; w[6]  = 0.040224796227922;
      r[7]  = -0.204472124089526; s[7]  = -0.204472124089526; w[7]  = 0.053635694518663;
      r[8]  = -0.591055751820947; s[8]  = -0.204472124089526; w[8]  = 0.053635694518663;
      r[9]  = -0.204472124089526; s[9]  = -0.591055751820947; w[9]  = 0.053635694518663;
      r[10] = -0.470999594934425; s[10] = -0.470999594934425; w[10] = 0.049046267603004;
      r[11] = -0.058000810131149; s[11] = -0.470999594934425; w[11] = 0.049046267603004;
      r[12] = -0.470999594934425; s[12] = -0.058000810131149; w[12] = 0.049046267603004;
      r[13] = -0.577962071815846; s[13] = -0.577962071815846; w[13] = 0.032789156821391;
      r[14] =  0.155924143631693; s[14] = -0.577962071815846; w[14] = 0.032789156821391;
      r[15] = -0.577962071815846; s[15] =  0.155924143631693; w[15] = 0.032789156821391;
      r[16] = -0.784528785657457; s[16] = -0.784528785657457; w[16] = 0.029591814797299;
      r[17] =  0.569057571314915; s[17] = -0.784528785657457; w[17] = 0.029591814797299;
      r[18] = -0.784528785657457; s[18] =  0.569057571314915; w[18] = 0.029591814797299;
      r[19] = -0.921861824324395; s[19] = -0.921861824324395; w[19] = 0.009158564555409;
      r[20] =  0.843723648648789; s[20] = -0.921861824324395; w[20] = 0.009158564555409;
      r[21] = -0.921861824324395; s[21] =  0.843723648648789; w[21] = 0.009158564555409;
      r[22] = -0.977651240541341; s[22] = -0.977651240541341; w[22] = 0.003303653031152;
      r[23] =  0.955302481082681; s[23] = -0.977651240541341; w[23] = 0.003303653031152;
      r[24] = -0.977651240541341; s[24] =  0.955302481082681; w[24] = 0.003303653031152;
      r[25] = -0.872900668183296; s[25] =  0.862201431808621; w[25] = 0.004698341817151;
      r[26] =  0.862201431808621; s[26] = -0.872900668183296; w[26] = 0.004698341817151;
      r[27] = -0.989300763625325; s[27] =  0.862201431808621; w[27] = 0.004698341817151;
      r[28] =  0.862201431808621; s[28] = -0.989300763625325; w[28] = 0.004698341817151;
      r[29] = -0.989300763625325; s[29] = -0.872900668183295; w[29] = 0.004698341817151;
      r[30] = -0.872900668183295; s[30] = -0.989300763625325; w[30] = 0.004698341817151;
      r[31] = -0.685786162118586; s[31] =  0.669876527986188; w[31] = 0.008931851508364;
      r[32] =  0.669876527986188; s[32] = -0.685786162118586; w[32] = 0.008931851508364;
      r[33] = -0.984090365867602; s[33] =  0.669876527986188; w[33] = 0.008931851508364;
      r[34] =  0.669876527986188; s[34] = -0.984090365867602; w[34] = 0.008931851508364;
      r[35] = -0.984090365867602; s[35] = -0.685786162118586; w[35] = 0.008931851508364;
      r[36] = -0.685786162118586; s[36] = -0.984090365867602; w[36] = 0.008931851508364;
      r[37] = -0.208715771271252; s[37] =  0.187870974708724; w[37] = 0.012199133615816;
      r[38] =  0.187870974708724; s[38] = -0.208715771271252; w[38] = 0.012199133615816;
      r[39] = -0.979155203437472; s[39] =  0.187870974708724; w[39] = 0.012199133615816;
      r[40] =  0.187870974708724; s[40] = -0.979155203437472; w[40] = 0.012199133615816;
      r[41] = -0.979155203437472; s[41] = -0.208715771271252; w[41] = 0.012199133615816;
      r[42] = -0.208715771271252; s[42] = -0.979155203437472; w[42] = 0.012199133615816;
      r[43] = -0.453664858574179; s[43] =  0.431736028981932; w[43] = 0.013782162654376;
      r[44] =  0.431736028981932; s[44] = -0.453664858574179; w[44] = 0.013782162654376;
      r[45] = -0.978071170407753; s[45] =  0.431736028981932; w[45] = 0.013782162654376;
      r[46] =  0.431736028981932; s[46] = -0.978071170407753; w[46] = 0.013782162654376;
      r[47] = -0.978071170407753; s[47] = -0.453664858574179; w[47] = 0.013782162654376;
      r[48] = -0.453664858574179; s[48] = -0.978071170407753; w[48] = 0.013782162654376;
      r[49] = -0.796429235029966; s[49] =  0.719295810859041; w[49] = 0.015994950144956;
      r[50] =  0.719295810859041; s[50] = -0.796429235029966; w[50] = 0.015994950144956;
      r[51] = -0.922866575829075; s[51] =  0.719295810859041; w[51] = 0.015994950144956;
      r[52] =  0.719295810859041; s[52] = -0.922866575829075; w[52] = 0.015994950144956;
      r[53] = -0.922866575829075; s[53] = -0.796429235029966; w[53] = 0.015994950144956;
      r[54] = -0.796429235029966; s[54] = -0.922866575829075; w[54] = 0.015994950144956;
      r[55] = -0.106682901647172; s[55] =  0.035521886012736; w[55] = 0.014772268570672;
      r[56] =  0.035521886012736; s[56] = -0.106682901647172; w[56] = 0.014772268570672;
      r[57] = -0.928838984365563; s[57] =  0.035521886012736; w[57] = 0.014772268570672;
      r[58] =  0.035521886012736; s[58] = -0.928838984365563; w[58] = 0.014772268570672;
      r[59] = -0.928838984365564; s[59] = -0.106682901647172; w[59] = 0.014772268570672;
      r[60] = -0.106682901647172; s[60] = -0.928838984365564; w[60] = 0.014772268570672;
      r[61] = -0.601978411700994; s[61] =  0.502636778975466; w[61] = 0.025598663757297;
      r[62] =  0.502636778975466; s[62] = -0.601978411700994; w[62] = 0.025598663757297;
      r[63] = -0.900658367274472; s[63] =  0.502636778975466; w[63] = 0.025598663757297;
      r[64] =  0.502636778975466; s[64] = -0.900658367274472; w[64] = 0.025598663757297;
      r[65] = -0.900658367274472; s[65] = -0.601978411700994; w[65] = 0.025598663757297;
      r[66] = -0.601978411700994; s[66] = -0.900658367274472; w[66] = 0.025598663757297;
      r[67] = -0.351477632615435; s[67] =  0.234438182446771; w[67] = 0.034516142351393;
      r[68] =  0.234438182446771; s[68] = -0.351477632615435; w[68] = 0.034516142351393;
      r[69] = -0.882960549831337; s[69] =  0.234438182446771; w[69] = 0.034516142351393;
      r[70] =  0.234438182446771; s[70] = -0.882960549831337; w[70] = 0.034516142351393;
      r[71] = -0.882960549831337; s[71] = -0.351477632615435; w[71] = 0.034516142351393;
      r[72] = -0.351477632615435; s[72] = -0.882960549831337; w[72] = 0.034516142351393;
      r[73] = -0.582937273579734; s[73] =  0.339941699570946; w[73] = 0.037345891805871;
      r[74] =  0.339941699570946; s[74] = -0.582937273579734; w[74] = 0.037345891805871;
      r[75] = -0.757004425991211; s[75] =  0.339941699570946; w[75] = 0.037345891805871;
      r[76] =  0.339941699570946; s[76] = -0.757004425991211; w[76] = 0.037345891805871;
      r[77] = -0.757004425991211; s[77] = -0.582937273579734; w[77] = 0.037345891805871;
      r[78] = -0.582937273579734; s[78] = -0.757004425991211; w[78] = 0.037345891805871;
      r[79] = -0.353658866927485; s[79] =  0.072237177039608; w[79] = 0.045636448116791;
      r[80] =  0.072237177039608; s[80] = -0.353658866927485; w[80] = 0.045636448116791;
      r[81] = -0.718578310112123; s[81] =  0.072237177039608; w[81] = 0.045636448116791;
      r[82] =  0.072237177039608; s[82] = -0.718578310112123; w[82] = 0.045636448116791;
      r[83] = -0.718578310112123; s[83] = -0.353658866927485; w[83] = 0.045636448116791;
      r[84] = -0.353658866927485; s[84] = -0.718578310112123; w[84] = 0.045636448116791;
    
      break;
    }
    case 21: {
      r[0]  = -0.006427416686680; s[0]  = -0.006427416686680; w[0]  = 0.009411977990382;
      r[1]  = -0.987145166626641; s[1]  = -0.006427416686680; w[1]  = 0.009411977990382;
      r[2]  = -0.006427416686680; s[2]  = -0.987145166626641; w[2]  = 0.009411977990382;
      r[3]  = -0.037371238381685; s[3]  = -0.037371238381685; w[3]  = 0.023551229711283;
      r[4]  = -0.925257523236629; s[4]  = -0.037371238381685; w[4]  = 0.023551229711283;
      r[5]  = -0.037371238381685; s[5]  = -0.925257523236629; w[5]  = 0.023551229711283;
      r[6]  = -0.101972641308225; s[6]  = -0.101972641308225; w[6]  = 0.041353961937512;
      r[7]  = -0.796054717383551; s[7]  = -0.101972641308225; w[7]  = 0.041353961937512;
      r[8]  = -0.101972641308225; s[8]  = -0.796054717383551; w[8]  = 0.041353961937512;
      r[9]  = -0.405469125577621; s[9]  = -0.405469125577621; w[9]  = 0.045947520102707;
      r[10] = -0.189061748844758; s[10] = -0.405469125577621; w[10] = 0.045947520102707;
      r[11] = -0.405469125577621; s[11] = -0.189061748844758; w[11] = 0.045947520102707;
      r[12] = -0.556282620757932; s[12] = -0.556282620757932; w[12] = 0.047724273523632;
      r[13] =  0.112565241515865; s[13] = -0.556282620757932; w[13] = 0.047724273523632;
      r[14] = -0.556282620757932; s[14] =  0.112565241515865; w[14] = 0.047724273523632;
      r[15] = -0.787490889149774; s[15] = -0.787490889149774; w[15] = 0.019906824503099;
      r[16] =  0.574981778299548; s[16] = -0.787490889149774; w[16] = 0.019906824503099;
      r[17] = -0.787490889149774; s[17] =  0.574981778299548; w[17] = 0.019906824503099;
      r[18] = -0.892713038713217; s[18] = -0.892713038713217; w[18] = 0.014287994054933;
      r[19] =  0.785426077426434; s[19] = -0.892713038713217; w[19] = 0.014287994054933;
      r[20] = -0.892713038713217; s[20] =  0.785426077426434; w[20] = 0.014287994054933;
      r[21] = -0.975095902008750; s[21] =  0.956581360464913; w[21] = 0.001538784565604;
      r[22] =  0.956581360464913; s[22] = -0.975095902008750; w[22] = 0.001538784565604;
      r[23] = -0.981485458456163; s[23] =  0.956581360464913; w[23] = 0.001538784565604;
      r[24] =  0.956581360464913; s[24] = -0.981485458456163; w[24] = 0.001538784565604;
      r[25] = -0.981485458456163; s[25] = -0.975095902008750; w[25] = 0.001538784565604;
      r[26] = -0.975095902008750; s[26] = -0.981485458456163; w[26] = 0.001538784565604;
      r[27] = -0.734190363492767; s[27] =  0.718759150923754; w[27] = 0.007257720730941;
      r[28] =  0.718759150923754; s[28] = -0.734190363492767; w[28] = 0.007257720730941;
      r[29] = -0.984568787430988; s[29] =  0.718759150923754; w[29] = 0.007257720730941;
      r[30] =  0.718759150923754; s[30] = -0.984568787430988; w[30] = 0.007257720730941;
      r[31] = -0.984568787430988; s[31] = -0.734190363492767; w[31] = 0.007257720730941;
      r[32] = -0.734190363492767; s[32] = -0.984568787430988; w[32] = 0.007257720730941;
      r[33] = -0.279909645448207; s[33] =  0.260066127858194; w[33] = 0.013303360825119;
      r[34] =  0.260066127858194; s[34] = -0.279909645448207; w[34] = 0.013303360825119;
      r[35] = -0.980156482409987; s[35] =  0.260066127858194; w[35] = 0.013303360825119;
      r[36] =  0.260066127858194; s[36] = -0.980156482409987; w[36] = 0.013303360825119;
      r[37] = -0.980156482409987; s[37] = -0.279909645448207; w[37] = 0.013303360825119;
      r[38] = -0.279909645448207; s[38] = -0.980156482409987; w[38] = 0.013303360825119;
      r[39] = -0.527626383829046; s[39] =  0.508213994286434; w[39] = 0.011413497063052;
      r[40] =  0.508213994286434; s[40] = -0.527626383829046; w[40] = 0.011413497063052;
      r[41] = -0.980587610457388; s[41] =  0.508213994286434; w[41] = 0.011413497063052;
      r[42] =  0.508213994286434; s[42] = -0.980587610457388; w[42] = 0.011413497063052;
      r[43] = -0.980587610457388; s[43] = -0.527626383829046; w[43] = 0.011413497063052;
      r[44] = -0.527626383829046; s[44] = -0.980587610457388; w[44] = 0.011413497063052;
      r[45] = -0.887945885434766; s[45] =  0.867094084690673; w[45] = 0.006552521291655;
      r[46] =  0.867094084690673; s[46] = -0.887945885434766; w[46] = 0.006552521291655;
      r[47] = -0.979148199255907; s[47] =  0.867094084690673; w[47] = 0.006552521291655;
      r[48] =  0.867094084690673; s[48] = -0.979148199255907; w[48] = 0.006552521291655;
      r[49] = -0.979148199255907; s[49] = -0.887945885434766; w[49] = 0.006552521291655;
      r[50] = -0.887945885434766; s[50] = -0.979148199255907; w[50] = 0.006552521291655;
      r[51] = -0.738154400021390; s[51] =  0.655149643679198; w[51] = 0.017432381135646;
      r[52] =  0.655149643679198; s[52] = -0.738154400021390; w[52] = 0.017432381135646;
      r[53] = -0.916995243657807; s[53] =  0.655149643679198; w[53] = 0.017432381135646;
      r[54] =  0.655149643679198; s[54] = -0.916995243657807; w[54] = 0.017432381135646;
      r[55] = -0.916995243657807; s[55] = -0.738154400021390; w[55] = 0.017432381135646;
      r[56] = -0.738154400021390; s[56] = -0.916995243657807; w[56] = 0.017432381135646;
      r[57] = -0.296840355384531; s[57] =  0.194393499586184; w[57] = 0.027402609499462;
      r[58] =  0.194393499586184; s[58] = -0.296840355384531; w[58] = 0.027402609499462;
      r[59] = -0.897553144201652; s[59] =  0.194393499586184; w[59] = 0.027402609499462;
      r[60] =  0.194393499586184; s[60] = -0.897553144201652; w[60] = 0.027402609499462;
      r[61] = -0.897553144201652; s[61] = -0.296840355384531; w[61] = 0.027402609499462;
      r[62] = -0.296840355384531; s[62] = -0.897553144201652; w[62] = 0.027402609499462;
      r[63] = -0.531632641825188; s[63] =  0.431035654635683; w[63] = 0.025014173853794;
      r[64] =  0.431035654635683; s[64] = -0.531632641825188; w[64] = 0.025014173853794;
      r[65] = -0.899403012810496; s[65] =  0.431035654635683; w[65] = 0.025014173853794;
      r[66] =  0.431035654635683; s[66] = -0.899403012810496; w[66] = 0.025014173853794;
      r[67] = -0.899403012810496; s[67] = -0.531632641825188; w[67] = 0.025014173853794;
      r[68] = -0.531632641825188; s[68] = -0.899403012810496; w[68] = 0.025014173853794;
      r[69] = -0.682130525363537; s[69] =  0.500223407469218; w[69] = 0.008690112508894;
      r[70] =  0.500223407469218; s[70] = -0.682130525363537; w[70] = 0.008690112508894;
      r[71] = -0.818092882105681; s[71] =  0.500223407469218; w[71] = 0.008690112508894;
      r[72] =  0.500223407469218; s[72] = -0.818092882105681; w[72] = 0.008690112508894;
      r[73] = -0.818092882105681; s[73] = -0.682130525363537; w[73] = 0.008690112508894;
      r[74] = -0.682130525363537; s[74] = -0.818092882105681; w[74] = 0.008690112508894;
      r[75] = -0.376270450471865; s[75] =  0.138207573987544; w[75] = 0.037386369243416;
      r[76] =  0.138207573987544; s[76] = -0.376270450471865; w[76] = 0.037386369243416;
      r[77] = -0.761937123515678; s[77] =  0.138207573987544; w[77] = 0.037386369243416;
      r[78] =  0.138207573987544; s[78] = -0.761937123515678; w[78] = 0.037386369243416;
      r[79] = -0.761937123515678; s[79] = -0.376270450471865; w[79] = 0.037386369243416;
      r[80] = -0.376270450471865; s[80] = -0.761937123515678; w[80] = 0.037386369243416;
      r[81] = -0.590719081005038; s[81] =  0.339375331943743; w[81] = 0.032509611084921;
      r[82] =  0.339375331943743; s[82] = -0.590719081005038; w[82] = 0.032509611084921;
      r[83] = -0.748656250938705; s[83] =  0.339375331943743; w[83] = 0.032509611084921;
      r[84] =  0.339375331943743; s[84] = -0.748656250938705; w[84] = 0.032509611084921;
      r[85] = -0.748656250938705; s[85] = -0.590719081005038; w[85] = 0.032509611084921;
      r[86] = -0.590719081005038; s[86] = -0.748656250938705; w[86] = 0.032509611084921;
      r[87] = -0.305527033262063; s[87] = -0.081590222620926; w[87] = 0.043740300619054;
      r[88] = -0.081590222620926; s[88] = -0.305527033262063; w[88] = 0.043740300619054;
      r[89] = -0.612882744117011; s[89] = -0.081590222620926; w[89] = 0.043740300619054;
      r[90] = -0.081590222620926; s[90] = -0.612882744117011; w[90] = 0.043740300619054;
      r[91] = -0.612882744117011; s[91] = -0.305527033262063; w[91] = 0.043740300619054;
      r[92] = -0.305527033262063; s[92] = -0.612882744117011; w[92] = 0.043740300619054;
    
      break;
    }
    case 22: {
      r[0]   = -0.333333333333333; s[0]   = -0.333333333333333; w[0]   = 0.052134891986775;
      r[1]   = -0.005980679319437; s[1]   = -0.005980679319437; w[1]   = 0.002127559342619;
      r[2]   = -0.988038641361126; s[2]   = -0.005980679319437; w[2]   = 0.002127559342619;
      r[3]   = -0.005980679319437; s[3]   = -0.988038641361126; w[3]   = 0.002127559342619;
      r[4]   = -0.112289615572167; s[4]   = -0.112289615572167; w[4]   = 0.040243429730281;
      r[5]   = -0.775420768855665; s[5]   = -0.112289615572167; w[5]   = 0.040243429730281;
      r[6]   = -0.112289615572167; s[6]   = -0.775420768855665; w[6]   = 0.040243429730281;
      r[7]   = -0.209082673923590; s[7]   = -0.209082673923590; w[7]   = 0.048340550896214;
      r[8]   = -0.581834652152820; s[8]   = -0.209082673923590; w[8]   = 0.048340550896214;
      r[9]   = -0.209082673923590; s[9]   = -0.581834652152820; w[9]   = 0.048340550896214;
      r[10]  = -0.467712894350161; s[10]  = -0.467712894350161; w[10]  = 0.047378470721283;
      r[11]  = -0.064574211299679; s[11]  = -0.467712894350161; w[11]  = 0.047378470721283;
      r[12]  = -0.467712894350161; s[12]  = -0.064574211299679; w[12]  = 0.047378470721283;
      r[13]  = -0.616276252022817; s[13]  = -0.616276252022817; w[13]  = 0.040231109011798;
      r[14]  =  0.232552504045635; s[14]  = -0.616276252022817; w[14]  = 0.040231109011798;
      r[15]  = -0.616276252022817; s[15]  =  0.232552504045635; w[15]  = 0.040231109011798;
      r[16]  = -0.755246756631224; s[16]  = -0.755246756631224; w[16]  = 0.027899919845537;
      r[17]  =  0.510493513262448; s[17]  = -0.755246756631224; w[17]  = 0.027899919845537;
      r[18]  = -0.755246756631224; s[18]  =  0.510493513262448; w[18]  = 0.027899919845537;
      r[19]  = -0.893562558892192; s[19]  = -0.893562558892192; w[19]  = 0.012350949744993;
      r[20]  =  0.787125117784384; s[20]  = -0.893562558892192; w[20]  = 0.012350949744993;
      r[21]  = -0.893562558892192; s[21]  =  0.787125117784384; w[21]  = 0.012350949744993;
      r[22]  = -0.943105536827199; s[22]  = -0.943105536827199; w[22]  = 0.003039705316805;
      r[23]  =  0.886211073654399; s[23]  = -0.943105536827199; w[23]  = 0.003039705316805;
      r[24]  = -0.943105536827199; s[24]  =  0.886211073654399; w[24]  = 0.003039705316805;
      r[25]  = -0.991786296535875; s[25]  = -0.991786296535875; w[25]  = 0.000691671630889;
      r[26]  =  0.983572593071750; s[26]  = -0.991786296535875; w[26]  = 0.000691671630889;
      r[27]  = -0.991786296535875; s[27]  =  0.983572593071750; w[27]  = 0.000691671630889;
      r[28]  = -0.413667191986728; s[28]  =  0.401491377065593; w[28]  = 0.008022198405286;
      r[29]  =  0.401491377065593; s[29]  = -0.413667191986728; w[29]  = 0.008022198405286;
      r[30]  = -0.987824185078864; s[30]  =  0.401491377065593; w[30]  = 0.008022198405286;
      r[31]  =  0.401491377065593; s[31]  = -0.987824185078864; w[31]  = 0.008022198405286;
      r[32]  = -0.987824185078865; s[32]  = -0.413667191986728; w[32]  = 0.008022198405286;
      r[33]  = -0.413667191986728; s[33]  = -0.987824185078865; w[33]  = 0.008022198405286;
      r[34]  = -0.932052869826067; s[34]  =  0.918490219898898; w[34]  = 0.003090974933637;
      r[35]  =  0.918490219898898; s[35]  = -0.932052869826067; w[35]  = 0.003090974933637;
      r[36]  = -0.986437350072830; s[36]  =  0.918490219898898; w[36]  = 0.003090974933637;
      r[37]  =  0.918490219898898; s[37]  = -0.986437350072830; w[37]  = 0.003090974933637;
      r[38]  = -0.986437350072830; s[38]  = -0.932052869826067; w[38]  = 0.003090974933637;
      r[39]  = -0.932052869826067; s[39]  = -0.986437350072830; w[39]  = 0.003090974933637;
      r[40]  = -0.640599464469824; s[40]  =  0.623868058857117; w[40]  = 0.008631214056178;
      r[41]  =  0.623868058857117; s[41]  = -0.640599464469824; w[41]  = 0.008631214056178;
      r[42]  = -0.983268594387293; s[42]  =  0.623868058857116; w[42]  = 0.008631214056178;
      r[43]  =  0.623868058857116; s[43]  = -0.983268594387293; w[43]  = 0.008631214056178;
      r[44]  = -0.983268594387293; s[44]  = -0.640599464469824; w[44]  = 0.008631214056178;
      r[45]  = -0.640599464469824; s[45]  = -0.983268594387293; w[45]  = 0.008631214056178;
      r[46]  = -0.158411570812104; s[46]  =  0.140890429939385; w[46]  = 0.010897457944070;
      r[47]  =  0.140890429939385; s[47]  = -0.158411570812104; w[47]  = 0.010897457944070;
      r[48]  = -0.982478859127282; s[48]  =  0.140890429939385; w[48]  = 0.010897457944070;
      r[49]  =  0.140890429939385; s[49]  = -0.982478859127282; w[49]  = 0.010897457944070;
      r[50]  = -0.982478859127282; s[50]  = -0.158411570812104; w[50]  = 0.010897457944070;
      r[51]  = -0.158411570812104; s[51]  = -0.982478859127282; w[51]  = 0.010897457944070;
      r[52]  = -0.814806216564922; s[52]  =  0.793459188393794; w[52]  = 0.008001551263950;
      r[53]  =  0.793459188393794; s[53]  = -0.814806216564922; w[53]  = 0.008001551263950;
      r[54]  = -0.978652971828873; s[54]  =  0.793459188393795; w[54]  = 0.008001551263950;
      r[55]  =  0.793459188393795; s[55]  = -0.978652971828873; w[55]  = 0.008001551263950;
      r[56]  = -0.978652971828873; s[56]  = -0.814806216564922; w[56]  = 0.008001551263950;
      r[57]  = -0.814806216564922; s[57]  = -0.978652971828873; w[57]  = 0.008001551263950;
      r[58]  = -0.363960589331557; s[58]  =  0.295112627690632; w[58]  = 0.017745434082141;
      r[59]  =  0.295112627690632; s[59]  = -0.363960589331557; w[59]  = 0.017745434082141;
      r[60]  = -0.931152038359075; s[60]  =  0.295112627690632; w[60]  = 0.017745434082141;
      r[61]  =  0.295112627690632; s[61]  = -0.931152038359075; w[61]  = 0.017745434082141;
      r[62]  = -0.931152038359075; s[62]  = -0.363960589331557; w[62]  = 0.017745434082141;
      r[63]  = -0.363960589331557; s[63]  = -0.931152038359075; w[63]  = 0.017745434082141;
      r[64]  = -0.570323384314285; s[64]  =  0.491002374048080; w[64]  = 0.017950672107247;
      r[65]  =  0.491002374048080; s[65]  = -0.570323384314285; w[65]  = 0.017950672107247;
      r[66]  = -0.920678989733795; s[66]  =  0.491002374048080; w[66]  = 0.017950672107247;
      r[67]  =  0.491002374048080; s[67]  = -0.920678989733795; w[67]  = 0.017950672107247;
      r[68]  = -0.920678989733795; s[68]  = -0.570323384314285; w[68]  = 0.017950672107247;
      r[69]  = -0.570323384314285; s[69]  = -0.920678989733795; w[69]  = 0.017950672107247;
      r[70]  = -0.134636707603169; s[70]  =  0.047157943371727; w[70]  = 0.019228676191702;
      r[71]  =  0.047157943371727; s[71]  = -0.134636707603169; w[71]  = 0.019228676191702;
      r[72]  = -0.912521235768557; s[72]  =  0.047157943371727; w[72]  = 0.019228676191702;
      r[73]  =  0.047157943371727; s[73]  = -0.912521235768557; w[73]  = 0.019228676191702;
      r[74]  = -0.912521235768557; s[74]  = -0.134636707603169; w[74]  = 0.019228676191702;
      r[75]  = -0.134636707603169; s[75]  = -0.912521235768557; w[75]  = 0.019228676191702;
      r[76]  = -0.752804662376914; s[76]  =  0.649201355096884; w[76]  = 0.019100455942623;
      r[77]  =  0.649201355096884; s[77]  = -0.752804662376914; w[77]  = 0.019100455942623;
      r[78]  = -0.896396692719970; s[78]  =  0.649201355096884; w[78]  = 0.019100455942623;
      r[79]  =  0.649201355096884; s[79]  = -0.896396692719970; w[79]  = 0.019100455942623;
      r[80]  = -0.896396692719970; s[80]  = -0.752804662376914; w[80]  = 0.019100455942623;
      r[81]  = -0.752804662376914; s[81]  = -0.896396692719970; w[81]  = 0.019100455942623;
      r[82]  = -0.328638645186542; s[82]  =  0.161176552903625; w[82]  = 0.028083739150636;
      r[83]  =  0.161176552903625; s[83]  = -0.328638645186542; w[83]  = 0.028083739150636;
      r[84]  = -0.832537907717083; s[84]  =  0.161176552903625; w[84]  = 0.028083739150636;
      r[85]  =  0.161176552903625; s[85]  = -0.832537907717083; w[85]  = 0.028083739150636;
      r[86]  = -0.832537907717083; s[86]  = -0.328638645186542; w[86]  = 0.028083739150636;
      r[87]  = -0.328638645186542; s[87]  = -0.832537907717083; w[87]  = 0.028083739150636;
      r[88]  = -0.559939069093864; s[88]  =  0.361861076343249; w[88]  = 0.031874915632607;
      r[89]  =  0.361861076343249; s[89]  = -0.559939069093864; w[89]  = 0.031874915632607;
      r[90]  = -0.801922007249385; s[90]  =  0.361861076343249; w[90]  = 0.031874915632607;
      r[91]  =  0.361861076343249; s[91]  = -0.801922007249385; w[91]  = 0.031874915632607;
      r[92]  = -0.801922007249385; s[92]  = -0.559939069093864; w[92]  = 0.031874915632607;
      r[93]  = -0.559939069093864; s[93]  = -0.801922007249385; w[93]  = 0.031874915632607;
      r[94]  = -0.388196135129182; s[94]  =  0.068617773977901; w[94]  = 0.040865211838584;
      r[95]  =  0.068617773977901; s[95]  = -0.388196135129182; w[95]  = 0.040865211838584;
      r[96]  = -0.680421638848718; s[96]  =  0.068617773977901; w[96]  = 0.040865211838584;
      r[97]  =  0.068617773977901; s[97]  = -0.680421638848718; w[97]  = 0.040865211838584;
      r[98]  = -0.680421638848718; s[98]  = -0.388196135129182; w[98]  = 0.040865211838584;
      r[99]  = -0.388196135129182; s[99]  = -0.680421638848718; w[99]  = 0.040865211838584;
    
      break;
    }
    case 23: {
      r[0]   = -0.333333333333333; s[0]   = -0.333333333333333; w[0]   = 0.049993627092845;
      r[1]   = -0.023538910955543; s[1]   = -0.023538910955543; w[1]   = 0.008227775544650;
      r[2]   = -0.952922178088914; s[2]   = -0.023538910955543; w[2]   = 0.008227775544650;
      r[3]   = -0.023538910955543; s[3]   = -0.952922178088914; w[3]   = 0.008227775544650;
      r[4]   = -0.112255282121933; s[4]   = -0.112255282121933; w[4]   = 0.037662477311797;
      r[5]   = -0.775489435756134; s[5]   = -0.112255282121933; w[5]   = 0.037662477311797;
      r[6]   = -0.112255282121933; s[6]   = -0.775489435756134; w[6]   = 0.037662477311797;
      r[7]   = -0.212088442838205; s[7]   = -0.212088442838205; w[7]   = 0.046916667243031;
      r[8]   = -0.575823114323591; s[8]   = -0.212088442838205; w[8]   = 0.046916667243031;
      r[9]   = -0.212088442838205; s[9]   = -0.575823114323591; w[9]   = 0.046916667243031;
      r[10]  = -0.467032671065063; s[10]  = -0.467032671065063; w[10]  = 0.047191879210718;
      r[11]  = -0.065934657869875; s[11]  = -0.467032671065063; w[11]  = 0.047191879210718;
      r[12]  = -0.467032671065063; s[12]  = -0.065934657869875; w[12]  = 0.047191879210718;
      r[13]  = -0.602458675644835; s[13]  = -0.602458675644835; w[13]  = 0.039689081561875;
      r[14]  =  0.204917351289670; s[14]  = -0.602458675644835; w[14]  = 0.039689081561875;
      r[15]  = -0.602458675644835; s[15]  =  0.204917351289670; w[15]  = 0.039689081561875;
      r[16]  = -0.729340080155585; s[16]  = -0.729340080155585; w[16]  = 0.029223692430473;
      r[17]  =  0.458680160311169; s[17]  = -0.729340080155585; w[17]  = 0.029223692430473;
      r[18]  = -0.729340080155585; s[18]  =  0.458680160311169; w[18]  = 0.029223692430473;
      r[19]  = -0.835032835950571; s[19]  = -0.835032835950571; w[19]  = 0.018566057479642;
      r[20]  =  0.670065671901142; s[20]  = -0.835032835950571; w[20]  = 0.018566057479642;
      r[21]  = -0.835032835950571; s[21]  =  0.670065671901142; w[21]  = 0.018566057479642;
      r[22]  = -0.916650365183127; s[22]  = -0.916650365183127; w[22]  = 0.008784960019827;
      r[23]  =  0.833300730366254; s[23]  = -0.916650365183127; w[23]  = 0.008784960019827;
      r[24]  = -0.916650365183127; s[24]  =  0.833300730366254; w[24]  = 0.008784960019827;
      r[25]  = -0.982250743407210; s[25]  = -0.982250743407210; w[25]  = 0.002061552741126;
      r[26]  =  0.964501486814420; s[26]  = -0.982250743407210; w[26]  = 0.002061552741126;
      r[27]  = -0.982250743407210; s[27]  =  0.964501486814420; w[27]  = 0.002061552741126;
      r[28]  = -0.105601712649592; s[28]  =  0.098292547426640; w[28]  = 0.004448377162851;
      r[29]  =  0.098292547426640; s[29]  = -0.105601712649592; w[29]  = 0.004448377162851;
      r[30]  = -0.992690834777047; s[30]  =  0.098292547426640; w[30]  = 0.004448377162851;
      r[31]  =  0.098292547426640; s[31]  = -0.992690834777047; w[31]  = 0.004448377162851;
      r[32]  = -0.992690834777047; s[32]  = -0.105601712649592; w[32]  = 0.004448377162851;
      r[33]  = -0.105601712649592; s[33]  = -0.992690834777047; w[33]  = 0.004448377162851;
      r[34]  = -0.775227228315076; s[34]  =  0.773362241907340; w[34]  = 0.002318973002522;
      r[35]  =  0.773362241907340; s[35]  = -0.775227228315076; w[35]  = 0.002318973002522;
      r[36]  = -0.998135013592264; s[36]  =  0.773362241907340; w[36]  = 0.002318973002522;
      r[37]  =  0.773362241907340; s[37]  = -0.998135013592264; w[37]  = 0.002318973002522;
      r[38]  = -0.998135013592264; s[38]  = -0.775227228315076; w[38]  = 0.002318973002522;
      r[39]  = -0.775227228315076; s[39]  = -0.998135013592264; w[39]  = 0.002318973002522;
      r[40]  = -0.399056075238285; s[40]  =  0.392811698367078; w[40]  = 0.005127488569570;
      r[41]  =  0.392811698367078; s[41]  = -0.399056075238285; w[41]  = 0.005127488569570;
      r[42]  = -0.993755623128794; s[42]  =  0.392811698367078; w[42]  = 0.005127488569570;
      r[43]  =  0.392811698367078; s[43]  = -0.993755623128794; w[43]  = 0.005127488569570;
      r[44]  = -0.993755623128794; s[44]  = -0.399056075238285; w[44]  = 0.005127488569570;
      r[45]  = -0.399056075238285; s[45]  = -0.993755623128794; w[45]  = 0.005127488569570;
      r[46]  = -0.907633572417643; s[46]  =  0.891831436837026; w[46]  = 0.004167853762943;
      r[47]  =  0.891831436837026; s[47]  = -0.907633572417643; w[47]  = 0.004167853762943;
      r[48]  = -0.984197864419384; s[48]  =  0.891831436837026; w[48]  = 0.004167853762943;
      r[49]  =  0.891831436837026; s[49]  = -0.984197864419384; w[49]  = 0.004167853762943;
      r[50]  = -0.984197864419383; s[50]  = -0.907633572417642; w[50]  = 0.004167853762943;
      r[51]  = -0.907633572417642; s[51]  = -0.984197864419383; w[51]  = 0.004167853762943;
      r[52]  = -0.610008084404048; s[52]  =  0.588812673873118; w[52]  = 0.010245094978581;
      r[53]  =  0.588812673873118; s[53]  = -0.610008084404048; w[53]  = 0.010245094978581;
      r[54]  = -0.978804589469070; s[54]  =  0.588812673873118; w[54]  = 0.010245094978581;
      r[55]  =  0.588812673873118; s[55]  = -0.978804589469070; w[55]  = 0.010245094978581;
      r[56]  = -0.978804589469070; s[56]  = -0.610008084404048; w[56]  = 0.010245094978581;
      r[57]  = -0.610008084404048; s[57]  = -0.978804589469070; w[57]  = 0.010245094978581;
      r[58]  = -0.243171741354056; s[58]  =  0.207960140279272; w[58]  = 0.012823590639281;
      r[59]  =  0.207960140279272; s[59]  = -0.243171741354056; w[59]  = 0.012823590639281;
      r[60]  = -0.964788398925216; s[60]  =  0.207960140279272; w[60]  = 0.012823590639281;
      r[61]  =  0.207960140279272; s[61]  = -0.964788398925216; w[61]  = 0.012823590639281;
      r[62]  = -0.964788398925216; s[62]  = -0.243171741354055; w[62]  = 0.012823590639281;
      r[63]  = -0.243171741354055; s[63]  = -0.964788398925216; w[63]  = 0.012823590639281;
      r[64]  = -0.789506197850739; s[64]  =  0.739925802118596; w[64]  = 0.011635516794817;
      r[65]  =  0.739925802118596; s[65]  = -0.789506197850739; w[65]  = 0.011635516794817;
      r[66]  = -0.950419604267857; s[66]  =  0.739925802118596; w[66]  = 0.011635516794817;
      r[67]  =  0.739925802118596; s[67]  = -0.950419604267857; w[67]  = 0.011635516794817;
      r[68]  = -0.950419604267857; s[68]  = -0.789506197850739; w[68]  = 0.011635516794817;
      r[69]  = -0.789506197850739; s[69]  = -0.950419604267857; w[69]  = 0.011635516794817;
      r[70]  = -0.455884762353173; s[70]  =  0.385117794610837; w[70]  = 0.019844273900077;
      r[71]  =  0.385117794610837; s[71]  = -0.455884762353173; w[71]  = 0.019844273900077;
      r[72]  = -0.929233032257664; s[72]  =  0.385117794610837; w[72]  = 0.019844273900077;
      r[73]  =  0.385117794610837; s[73]  = -0.929233032257664; w[73]  = 0.019844273900077;
      r[74]  = -0.929233032257664; s[74]  = -0.455884762353173; w[74]  = 0.019844273900077;
      r[75]  = -0.455884762353173; s[75]  = -0.929233032257664; w[75]  = 0.019844273900077;
      r[76]  = -0.112482048478797; s[76]  =  0.019740291196956; w[76]  = 0.011754407449188;
      r[77]  =  0.019740291196956; s[77]  = -0.112482048478797; w[77]  = 0.011754407449188;
      r[78]  = -0.907258242718159; s[78]  =  0.019740291196956; w[78]  = 0.011754407449188;
      r[79]  =  0.019740291196956; s[79]  = -0.907258242718159; w[79]  = 0.011754407449188;
      r[80]  = -0.907258242718159; s[80]  = -0.112482048478797; w[80]  = 0.011754407449188;
      r[81]  = -0.112482048478797; s[81]  = -0.907258242718159; w[81]  = 0.011754407449188;
      r[82]  = -0.664205324633765; s[82]  =  0.550530931038120; w[82]  = 0.021600208906600;
      r[83]  =  0.550530931038120; s[83]  = -0.664205324633765; w[83]  = 0.021600208906600;
      r[84]  = -0.886325606404355; s[84]  =  0.550530931038120; w[84]  = 0.021600208906600;
      r[85]  =  0.550530931038120; s[85]  = -0.886325606404355; w[85]  = 0.021600208906600;
      r[86]  = -0.886325606404354; s[86]  = -0.664205324633765; w[86]  = 0.021600208906600;
      r[87]  = -0.664205324633765; s[87]  = -0.886325606404354; w[87]  = 0.021600208906600;
      r[88]  = -0.295594853518405; s[88]  =  0.152149871080588; w[88]  = 0.028495796363599;
      r[89]  =  0.152149871080588; s[89]  = -0.295594853518405; w[89]  = 0.028495796363599;
      r[90]  = -0.856555017562183; s[90]  =  0.152149871080588; w[90]  = 0.028495796363599;
      r[91]  =  0.152149871080588; s[91]  = -0.856555017562183; w[91]  = 0.028495796363599;
      r[92]  = -0.856555017562183; s[92]  = -0.295594853518405; w[92]  = 0.028495796363599;
      r[93]  = -0.295594853518405; s[93]  = -0.856555017562183; w[93]  = 0.028495796363599;
      r[94]  = -0.525935282564036; s[94]  =  0.319991164262401; w[94]  = 0.031988381018369;
      r[95]  =  0.319991164262401; s[95]  = -0.525935282564036; w[95]  = 0.031988381018369;
      r[96]  = -0.794055881698366; s[96]  =  0.319991164262401; w[96]  = 0.031988381018369;
      r[97]  =  0.319991164262401; s[97]  = -0.794055881698366; w[97]  = 0.031988381018369;
      r[98]  = -0.794055881698366; s[98]  = -0.525935282564036; w[98]  = 0.031988381018369;
      r[99]  = -0.525935282564036; s[99]  = -0.794055881698366; w[99]  = 0.031988381018369;
      r[100] = -0.368673747852155; s[100] =  0.056982821883363; w[100] = 0.041389027831224;
      r[101] =  0.056982821883363; s[101] = -0.368673747852155; w[101] = 0.041389027831224;
      r[102] = -0.688309074031208; s[102] =  0.056982821883363; w[102] = 0.041389027831224;
      r[103] =  0.056982821883363; s[103] = -0.688309074031208; w[103] = 0.041389027831224;
      r[104] = -0.688309074031208; s[104] = -0.368673747852155; w[104] = 0.041389027831224;
      r[105] = -0.368673747852155; s[105] = -0.688309074031208; w[105] = 0.041389027831224;
    
      break;
    }
    case 24: {
      r[0]   = -0.333333333333333; s[0]   = -0.333333333333333; w[0]   = 0.031330242175782;
      r[1]   = -0.031382760283421; s[1]   = -0.031382760283421; w[1]   = 0.014717653847587;
      r[2]   = -0.937234479433159; s[2]   = -0.031382760283421; w[2]   = 0.014717653847587;
      r[3]   = -0.031382760283421; s[3]   = -0.937234479433159; w[3]   = 0.014717653847587;
      r[4]   = -0.120096233850662; s[4]   = -0.120096233850662; w[4]   = 0.030963348398970;
      r[5]   = -0.759807532298676; s[5]   = -0.120096233850662; w[5]   = 0.030963348398970;
      r[6]   = -0.120096233850662; s[6]   = -0.759807532298676; w[6]   = 0.030963348398970;
      r[7]   = -0.217343938317536; s[7]   = -0.217343938317536; w[7]   = 0.037986940843648;
      r[8]   = -0.565312123364928; s[8]   = -0.217343938317536; w[8]   = 0.037986940843648;
      r[9]   = -0.217343938317536; s[9]   = -0.565312123364928; w[9]   = 0.037986940843648;
      r[10]  = -0.429818284439074; s[10]  = -0.429818284439074; w[10]  = 0.033466254011049;
      r[11]  = -0.140363431121852; s[11]  = -0.429818284439074; w[11]  = 0.033466254011049;
      r[12]  = -0.429818284439074; s[12]  = -0.140363431121852; w[12]  = 0.033466254011049;
      r[13]  = -0.522392801883300; s[13]  = -0.522392801883300; w[13]  = 0.029203686754054;
      r[14]  =  0.044785603766601; s[14]  = -0.522392801883300; w[14]  = 0.029203686754054;
      r[15]  = -0.522392801883300; s[15]  =  0.044785603766601; w[15]  = 0.029203686754054;
      r[16]  = -0.704274510649686; s[16]  = -0.704274510649686; w[16]  = 0.022147482605586;
      r[17]  =  0.408549021299373; s[17]  = -0.704274510649686; w[17]  = 0.022147482605586;
      r[18]  = -0.704274510649686; s[18]  =  0.408549021299373; w[18]  = 0.022147482605586;
      r[19]  = -0.843105586121546; s[19]  = -0.843105586121546; w[19]  = 0.015827753677335;
      r[20]  =  0.686211172243092; s[20]  = -0.843105586121546; w[20]  = 0.015827753677335;
      r[21]  = -0.843105586121546; s[21]  =  0.686211172243092; w[21]  = 0.015827753677335;
      r[22]  = -0.936285318267427; s[22]  = -0.936285318267427; w[22]  = 0.005897827695084;
      r[23]  =  0.872570636534854; s[23]  = -0.936285318267427; w[23]  = 0.005897827695084;
      r[24]  = -0.936285318267427; s[24]  =  0.872570636534854; w[24]  = 0.005897827695084;
      r[25]  = -0.983008712604130; s[25]  = -0.983008712604130; w[25]  = 0.001905131297570;
      r[26]  =  0.966017425208259; s[26]  = -0.983008712604130; w[26]  = 0.001905131297570;
      r[27]  = -0.983008712604130; s[27]  =  0.966017425208259; w[27]  = 0.001905131297570;
      r[28]  = -0.133936561710478; s[28]  =  0.120327965916980; w[28]  = 0.008627443386199;
      r[29]  =  0.120327965916980; s[29]  = -0.133936561710478; w[29]  = 0.008627443386199;
      r[30]  = -0.986391404206501; s[30]  =  0.120327965916980; w[30]  = 0.008627443386199;
      r[31]  =  0.120327965916980; s[31]  = -0.986391404206501; w[31]  = 0.008627443386199;
      r[32]  = -0.986391404206501; s[32]  = -0.133936561710478; w[32]  = 0.008627443386199;
      r[33]  = -0.133936561710478; s[33]  = -0.986391404206501; w[33]  = 0.008627443386199;
      r[34]  = -0.589559941056661; s[34]  =  0.574991625876557; w[34]  = 0.007292016855628;
      r[35]  =  0.574991625876557; s[35]  = -0.589559941056661; w[35]  = 0.007292016855628;
      r[36]  = -0.985431684819897; s[36]  =  0.574991625876557; w[36]  = 0.007292016855628;
      r[37]  =  0.574991625876557; s[37]  = -0.985431684819897; w[37]  = 0.007292016855628;
      r[38]  = -0.985431684819897; s[38]  = -0.589559941056661; w[38]  = 0.007292016855628;
      r[39]  = -0.589559941056661; s[39]  = -0.985431684819897; w[39]  = 0.007292016855628;
      r[40]  = -0.375079684130172; s[40]  =  0.361371712031196; w[40]  = 0.008117703702749;
      r[41]  =  0.361371712031196; s[41]  = -0.375079684130172; w[41]  = 0.008117703702749;
      r[42]  = -0.986292027901024; s[42]  =  0.361371712031196; w[42]  = 0.008117703702749;
      r[43]  =  0.361371712031196; s[43]  = -0.986292027901024; w[43]  = 0.008117703702749;
      r[44]  = -0.986292027901024; s[44]  = -0.375079684130172; w[44]  = 0.008117703702749;
      r[45]  = -0.375079684130172; s[45]  = -0.986292027901024; w[45]  = 0.008117703702749;
      r[46]  = -0.906832105591067; s[46]  =  0.897244066241012; w[46]  = 0.002858359515171;
      r[47]  =  0.897244066241012; s[47]  = -0.906832105591067; w[47]  = 0.002858359515171;
      r[48]  = -0.990411960649945; s[48]  =  0.897244066241012; w[48]  = 0.002858359515171;
      r[49]  =  0.897244066241012; s[49]  = -0.990411960649945; w[49]  = 0.002858359515171;
      r[50]  = -0.990411960649945; s[50]  = -0.906832105591067; w[50]  = 0.002858359515171;
      r[51]  = -0.906832105591067; s[51]  = -0.990411960649945; w[51]  = 0.002858359515171;
      r[52]  = -0.768281564078976; s[52]  =  0.757477531703546; w[52]  = 0.004596788633153;
      r[53]  =  0.757477531703546; s[53]  = -0.768281564078976; w[53]  = 0.004596788633153;
      r[54]  = -0.989195967624570; s[54]  =  0.757477531703546; w[54]  = 0.004596788633153;
      r[55]  =  0.757477531703546; s[55]  = -0.989195967624570; w[55]  = 0.004596788633153;
      r[56]  = -0.989195967624570; s[56]  = -0.768281564078975; w[56]  = 0.004596788633153;
      r[57]  = -0.768281564078975; s[57]  = -0.989195967624570; w[57]  = 0.004596788633153;
      r[58]  = -0.498678726055310; s[58]  =  0.424585822926150; w[58]  = 0.017097477649229;
      r[59]  =  0.424585822926150; s[59]  = -0.498678726055310; w[59]  = 0.017097477649229;
      r[60]  = -0.925907096870839; s[60]  =  0.424585822926150; w[60]  = 0.017097477649229;
      r[61]  =  0.424585822926150; s[61]  = -0.925907096870839; w[61]  = 0.017097477649229;
      r[62]  = -0.925907096870839; s[62]  = -0.498678726055310; w[62]  = 0.017097477649229;
      r[63]  = -0.498678726055310; s[63]  = -0.925907096870839; w[63]  = 0.017097477649229;
      r[64]  = -0.279607664483490; s[64]  =  0.208396882316209; w[64]  = 0.018575372163249;
      r[65]  =  0.208396882316209; s[65]  = -0.279607664483490; w[65]  = 0.018575372163249;
      r[66]  = -0.928789217832719; s[66]  =  0.208396882316209; w[66]  = 0.018575372163249;
      r[67]  =  0.208396882316209; s[67]  = -0.928789217832719; w[67]  = 0.018575372163249;
      r[68]  = -0.928789217832719; s[68]  = -0.279607664483490; w[68]  = 0.018575372163249;
      r[69]  = -0.279607664483490; s[69]  = -0.928789217832719; w[69]  = 0.018575372163249;
      r[70]  = -0.837737372642508; s[70]  =  0.780926628976898; w[70]  = 0.008711688567270;
      r[71]  =  0.780926628976898; s[71]  = -0.837737372642508; w[71]  = 0.008711688567270;
      r[72]  = -0.943189256334389; s[72]  =  0.780926628976898; w[72]  = 0.008711688567270;
      r[73]  =  0.780926628976898; s[73]  = -0.943189256334389; w[73]  = 0.008711688567270;
      r[74]  = -0.943189256334389; s[74]  = -0.837737372642508; w[74]  = 0.008711688567270;
      r[75]  = -0.837737372642508; s[75]  = -0.943189256334389; w[75]  = 0.008711688567270;
      r[76]  = -0.693711793201882; s[76]  =  0.625024627714433; w[76]  = 0.013292883787775;
      r[77]  =  0.625024627714433; s[77]  = -0.693711793201882; w[77]  = 0.013292883787775;
      r[78]  = -0.931312834512551; s[78]  =  0.625024627714433; w[78]  = 0.013292883787775;
      r[79]  =  0.625024627714433; s[79]  = -0.931312834512551; w[79]  = 0.013292883787775;
      r[80]  = -0.931312834512551; s[80]  = -0.693711793201882; w[80]  = 0.013292883787775;
      r[81]  = -0.693711793201882; s[81]  = -0.931312834512551; w[81]  = 0.013292883787775;
      r[82]  =  0.004700993829951; s[82]  = -0.135567929734853; w[82]  = 0.012281001510740;
      r[83]  = -0.135567929734853; s[83]  =  0.004700993829951; w[83]  = 0.012281001510740;
      r[84]  = -0.869133064095098; s[84]  = -0.135567929734853; w[84]  = 0.012281001510740;
      r[85]  = -0.135567929734853; s[85]  = -0.869133064095098; w[85]  = 0.012281001510740;
      r[86]  = -0.869133064095098; s[86]  =  0.004700993829951; w[86]  = 0.012281001510740;
      r[87]  =  0.004700993829951; s[87]  = -0.869133064095098; w[87]  = 0.012281001510740;
      r[88]  = -0.521976982650353; s[88]  =  0.349105045944945; w[88]  = 0.021820885555420;
      r[89]  =  0.349105045944945; s[89]  = -0.521976982650353; w[89]  = 0.021820885555420;
      r[90]  = -0.827128063294591; s[90]  =  0.349105045944945; w[90]  = 0.021820885555420;
      r[91]  =  0.349105045944945; s[91]  = -0.827128063294591; w[91]  = 0.021820885555420;
      r[92]  = -0.827128063294591; s[92]  = -0.521976982650353; w[92]  = 0.021820885555420;
      r[93]  = -0.521976982650353; s[93]  = -0.827128063294591; w[93]  = 0.021820885555420;
      r[94]  = -0.326968721136217; s[94]  =  0.144732602043697; w[94]  = 0.027371090249668;
      r[95]  =  0.144732602043697; s[95]  = -0.326968721136217; w[95]  = 0.027371090249668;
      r[96]  = -0.817763880907480; s[96]  =  0.144732602043697; w[96]  = 0.027371090249668;
      r[97]  =  0.144732602043697; s[97]  = -0.817763880907480; w[97]  = 0.027371090249668;
      r[98]  = -0.817763880907480; s[98]  = -0.326968721136217; w[98]  = 0.027371090249668;
      r[99]  = -0.326968721136217; s[99]  = -0.817763880907480; w[99]  = 0.027371090249668;
      r[100] = -0.696116803866709; s[100] =  0.531011336134871; w[100] = 0.018737543889147;
      r[101] =  0.531011336134871; s[101] = -0.696116803866709; w[101] = 0.018737543889147;
      r[102] = -0.834894532268162; s[102] =  0.531011336134871; w[102] = 0.018737543889147;
      r[103] =  0.531011336134871; s[103] = -0.834894532268162; w[103] = 0.018737543889147;
      r[104] = -0.834894532268162; s[104] = -0.696116803866709; w[104] = 0.018737543889147;
      r[105] = -0.696116803866709; s[105] = -0.834894532268162; w[105] = 0.018737543889147;
      r[106] = -0.536404917748559; s[106] =  0.226734558438523; w[106] = 0.028700757241388;
      r[107] =  0.226734558438523; s[107] = -0.536404917748559; w[107] = 0.028700757241388;
      r[108] = -0.690329640689964; s[108] =  0.226734558438523; w[108] = 0.028700757241388;
      r[109] =  0.226734558438523; s[109] = -0.690329640689964; w[109] = 0.028700757241388;
      r[110] = -0.690329640689964; s[110] = -0.536404917748559; w[110] = 0.028700757241388;
      r[111] = -0.536404917748559; s[111] = -0.690329640689964; w[111] = 0.028700757241388;
      r[112] = -0.344489563789357; s[112] =  0.009282127094255; w[112] = 0.033972574031807;
      r[113] =  0.009282127094255; s[113] = -0.344489563789357; w[113] = 0.033972574031807;
      r[114] = -0.664792563304898; s[114] =  0.009282127094255; w[114] = 0.033972574031807;
      r[115] =  0.009282127094255; s[115] = -0.664792563304898; w[115] = 0.033972574031807;
      r[116] = -0.664792563304898; s[116] = -0.344489563789357; w[116] = 0.033972574031807;
      r[117] = -0.344489563789357; s[117] = -0.664792563304898; w[117] = 0.033972574031807;
    
      break;
    }
    case 25: {
      r[0]   = -0.027946483073174; s[0]   = -0.027946483073174; w[0]   = 0.016011163760041;
      r[1]   = -0.944107033853652; s[1]   = -0.027946483073174; w[1]   = 0.016011163760041;
      r[2]   = -0.027946483073174; s[2]   = -0.944107033853652; w[2]   = 0.016011163760041;
      r[3]   = -0.131178601327651; s[3]   = -0.131178601327651; w[3]   = 0.031894153664781;
      r[4]   = -0.737642797344697; s[4]   = -0.131178601327651; w[4]   = 0.031894153664781;
      r[5]   = -0.131178601327651; s[5]   = -0.737642797344697; w[5]   = 0.031894153664781;
      r[6]   = -0.220221729512072; s[6]   = -0.220221729512072; w[6]   = 0.026218282461591;
      r[7]   = -0.559556540975855; s[7]   = -0.220221729512072; w[7]   = 0.026218282461591;
      r[8]   = -0.220221729512072; s[8]   = -0.559556540975855; w[8]   = 0.026218282461591;
      r[9]   = -0.403113531960391; s[9]   = -0.403113531960391; w[9]   = 0.039166001931271;
      r[10]  = -0.193772936079218; s[10]  = -0.403113531960391; w[10]  = 0.039166001931271;
      r[11]  = -0.403113531960391; s[11]  = -0.193772936079218; w[11]  = 0.039166001931271;
      r[12]  = -0.531911655325256; s[12]  = -0.531911655325256; w[12]  = 0.032941770883075;
      r[13]  =  0.063823310650513; s[13]  = -0.531911655325256; w[13]  = 0.032941770883075;
      r[14]  = -0.531911655325256; s[14]  =  0.063823310650513; w[14]  = 0.032941770883075;
      r[15]  = -0.697063330781965; s[15]  = -0.697063330781965; w[15]  = 0.017094558148184;
      r[16]  =  0.394126661563930; s[16]  = -0.697063330781965; w[16]  = 0.017094558148184;
      r[17]  = -0.697063330781965; s[17]  =  0.394126661563930; w[17]  = 0.017094558148184;
      r[18]  = -0.774532212908013; s[18]  = -0.774532212908013; w[18]  = 0.016323771714453;
      r[19]  =  0.549064425816026; s[19]  = -0.774532212908013; w[19]  = 0.016323771714453;
      r[20]  = -0.774532212908013; s[20]  =  0.549064425816026; w[20]  = 0.016323771714453;
      r[21]  = -0.844568615816947; s[21]  = -0.844568615816947; w[21]  = 0.012242293079968;
      r[22]  =  0.689137231633895; s[22]  = -0.844568615816948; w[22]  = 0.012242293079968;
      r[23]  = -0.844568615816948; s[23]  =  0.689137231633895; w[23]  = 0.012242293079968;
      r[24]  = -0.930213812771406; s[24]  = -0.930213812771406; w[24]  = 0.005816996529873;
      r[25]  =  0.860427625542812; s[25]  = -0.930213812771406; w[25]  = 0.005816996529873;
      r[26]  = -0.930213812771406; s[26]  =  0.860427625542812; w[26]  = 0.005816996529873;
      r[27]  = -0.985483630758135; s[27]  = -0.985483630758135; w[27]  = 0.001384550491324;
      r[28]  =  0.970967261516271; s[28]  = -0.985483630758135; w[28]  = 0.001384550491324;
      r[29]  = -0.985483630758135; s[29]  =  0.970967261516271; w[29]  = 0.001384550491324;
      r[30]  = -0.545571095693272; s[30]  =  0.542986390284387; w[30]  = 0.002496578398555;
      r[31]  =  0.542986390284387; s[31]  = -0.545571095693272; w[31]  = 0.002496578398555;
      r[32]  = -0.997415294591116; s[32]  =  0.542986390284387; w[32]  = 0.002496578398555;
      r[33]  =  0.542986390284387; s[33]  = -0.997415294591116; w[33]  = 0.002496578398555;
      r[34]  = -0.997415294591116; s[34]  = -0.545571095693272; w[34]  = 0.002496578398555;
      r[35]  = -0.545571095693272; s[35]  = -0.997415294591116; w[35]  = 0.002496578398555;
      r[36]  = -0.129978890292857; s[36]  =  0.119179487748624; w[36]  = 0.006809505817606;
      r[37]  =  0.119179487748624; s[37]  = -0.129978890292857; w[37]  = 0.006809505817606;
      r[38]  = -0.989200597455768; s[38]  =  0.119179487748624; w[38]  = 0.006809505817606;
      r[39]  =  0.119179487748624; s[39]  = -0.989200597455768; w[39]  = 0.006809505817606;
      r[40]  = -0.989200597455768; s[40]  = -0.129978890292857; w[40]  = 0.006809505817606;
      r[41]  = -0.129978890292857; s[41]  = -0.989200597455768; w[41]  = 0.006809505817606;
      r[42]  = -0.359380801455591; s[42]  =  0.346612795387641; w[42]  = 0.006719308652128;
      r[43]  =  0.346612795387641; s[43]  = -0.359380801455591; w[43]  = 0.006719308652128;
      r[44]  = -0.987231993932050; s[44]  =  0.346612795387641; w[44]  = 0.006719308652128;
      r[45]  =  0.346612795387641; s[45]  = -0.987231993932050; w[45]  = 0.006719308652128;
      r[46]  = -0.987231993932050; s[46]  = -0.359380801455591; w[46]  = 0.006719308652128;
      r[47]  = -0.359380801455591; s[47]  = -0.987231993932050; w[47]  = 0.006719308652128;
      r[48]  = -0.816499355439990; s[48]  =  0.806442932436004; w[48]  = 0.003432313078994;
      r[49]  =  0.806442932436004; s[49]  = -0.816499355439990; w[49]  = 0.003432313078994;
      r[50]  = -0.989943576996014; s[50]  =  0.806442932436004; w[50]  = 0.003432313078994;
      r[51]  =  0.806442932436004; s[51]  = -0.989943576996014; w[51]  = 0.003432313078994;
      r[52]  = -0.989943576996014; s[52]  = -0.816499355439990; w[52]  = 0.003432313078994;
      r[53]  = -0.816499355439990; s[53]  = -0.989943576996014; w[53]  = 0.003432313078994;
      r[54]  = -0.923978328282551; s[54]  =  0.910324811038988; w[54]  = 0.002961712633431;
      r[55]  =  0.910324811038988; s[55]  = -0.923978328282551; w[55]  = 0.002961712633431;
      r[56]  = -0.986346482756436; s[56]  =  0.910324811038988; w[56]  = 0.002961712633431;
      r[57]  =  0.910324811038988; s[57]  = -0.986346482756436; w[57]  = 0.002961712633431;
      r[58]  = -0.986346482756436; s[58]  = -0.923978328282551; w[58]  = 0.002961712633431;
      r[59]  = -0.923978328282551; s[59]  = -0.986346482756436; w[59]  = 0.002961712633431;
      r[60]  = -0.685149563029377; s[60]  =  0.665117163749518; w[60]  = 0.007022625221457;
      r[61]  =  0.665117163749518; s[61]  = -0.685149563029377; w[61]  = 0.007022625221457;
      r[62]  = -0.979967600720141; s[62]  =  0.665117163749517; w[62]  = 0.007022625221457;
      r[63]  =  0.665117163749517; s[63]  = -0.979967600720141; w[63]  = 0.007022625221457;
      r[64]  = -0.979967600720141; s[64]  = -0.685149563029376; w[64]  = 0.007022625221457;
      r[65]  = -0.685149563029376; s[65]  = -0.979967600720141; w[65]  = 0.007022625221457;
      r[66]  = -0.520220680442934; s[66]  =  0.468705054096154; w[66]  = 0.014787100299413;
      r[67]  =  0.468705054096154; s[67]  = -0.520220680442934; w[67]  = 0.014787100299413;
      r[68]  = -0.948484373653220; s[68]  =  0.468705054096154; w[68]  = 0.014787100299413;
      r[69]  =  0.468705054096154; s[69]  = -0.948484373653220; w[69]  = 0.014787100299413;
      r[70]  = -0.948484373653220; s[70]  = -0.520220680442934; w[70]  = 0.014787100299413;
      r[71]  = -0.520220680442934; s[71]  = -0.948484373653220; w[71]  = 0.014787100299413;
      r[72]  = -0.276113763747879; s[72]  =  0.215657967508047; w[72]  = 0.015966174954753;
      r[73]  =  0.215657967508047; s[73]  = -0.276113763747879; w[73]  = 0.015966174954753;
      r[74]  = -0.939544203760168; s[74]  =  0.215657967508047; w[74]  = 0.015966174954753;
      r[75]  =  0.215657967508047; s[75]  = -0.939544203760168; w[75]  = 0.015966174954753;
      r[76]  = -0.939544203760168; s[76]  = -0.276113763747879; w[76]  = 0.015966174954753;
      r[77]  = -0.276113763747879; s[77]  = -0.939544203760168; w[77]  = 0.015966174954753;
      r[78]  = -0.832896078090343; s[78]  =  0.771886097876019; w[78]  = 0.008711925226316;
      r[79]  =  0.771886097876019; s[79]  = -0.832896078090343; w[79]  = 0.008711925226316;
      r[80]  = -0.938990019785676; s[80]  =  0.771886097876019; w[80]  = 0.008711925226316;
      r[81]  =  0.771886097876019; s[81]  = -0.938990019785676; w[81]  = 0.008711925226316;
      r[82]  = -0.938990019785676; s[82]  = -0.832896078090343; w[82]  = 0.008711925226316;
      r[83]  = -0.832896078090343; s[83]  = -0.938990019785676; w[83]  = 0.008711925226316;
      r[84]  = -0.703113558535164; s[84]  =  0.611200463810025; w[84]  = 0.014730113402836;
      r[85]  =  0.611200463810025; s[85]  = -0.703113558535164; w[85]  = 0.014730113402836;
      r[86]  = -0.908086905274861; s[86]  =  0.611200463810025; w[86]  = 0.014730113402836;
      r[87]  =  0.611200463810025; s[87]  = -0.908086905274861; w[87]  = 0.014730113402836;
      r[88]  = -0.908086905274861; s[88]  = -0.703113558535164; w[88]  = 0.014730113402836;
      r[89]  = -0.703113558535164; s[89]  = -0.908086905274861; w[89]  = 0.014730113402836;
      r[90]  = -0.432520582544930; s[90]  =  0.297634981464375; w[90]  = 0.021927145692839;
      r[91]  =  0.297634981464375; s[91]  = -0.432520582544930; w[91]  = 0.021927145692839;
      r[92]  = -0.865114398919445; s[92]  =  0.297634981464375; w[92]  = 0.021927145692839;
      r[93]  =  0.297634981464375; s[93]  = -0.865114398919445; w[93]  = 0.021927145692839;
      r[94]  = -0.865114398919445; s[94]  = -0.432520582544930; w[94]  = 0.021927145692839;
      r[95]  = -0.432520582544930; s[95]  = -0.865114398919445; w[95]  = 0.021927145692839;
      r[96]  = -0.186201249762425; s[96]  =  0.046111066930604; w[96]  = 0.023499923487082;
      r[97]  =  0.046111066930604; s[97]  = -0.186201249762425; w[97]  = 0.023499923487082;
      r[98]  = -0.859909817168179; s[98]  =  0.046111066930604; w[98]  = 0.023499923487082;
      r[99]  =  0.046111066930604; s[99]  = -0.859909817168179; w[99]  = 0.023499923487082;
      r[100] = -0.859909817168179; s[100] = -0.186201249762425; w[100] = 0.023499923487082;
      r[101] = -0.186201249762425; s[101] = -0.859909817168179; w[101] = 0.023499923487082;
      r[102] = -0.611772025950215; s[102] =  0.443948976669982; w[102] = 0.020031201427597;
      r[103] =  0.443948976669982; s[103] = -0.611772025950215; w[103] = 0.020031201427597;
      r[104] = -0.832176950719767; s[104] =  0.443948976669982; w[104] = 0.020031201427597;
      r[105] =  0.443948976669982; s[105] = -0.832176950719767; w[105] = 0.020031201427597;
      r[106] = -0.832176950719767; s[106] = -0.611772025950215; w[106] = 0.020031201427597;
      r[107] = -0.611772025950215; s[107] = -0.832176950719767; w[107] = 0.020031201427597;
      r[108] = -0.351731305998594; s[108] =  0.110980234644288; w[108] = 0.026619281575257;
      r[109] =  0.110980234644288; s[109] = -0.351731305998594; w[109] = 0.026619281575257;
      r[110] = -0.759248928645695; s[110] =  0.110980234644288; w[110] = 0.026619281575257;
      r[111] =  0.110980234644288; s[111] = -0.759248928645695; w[111] = 0.026619281575257;
      r[112] = -0.759248928645695; s[112] = -0.351731305998594; w[112] = 0.026619281575257;
      r[113] = -0.351731305998594; s[113] = -0.759248928645695; w[113] = 0.026619281575257;
      r[114] = -0.541445032888038; s[114] =  0.245311234573305; w[114] = 0.028308893010452;
      r[115] =  0.245311234573305; s[115] = -0.541445032888038; w[115] = 0.028308893010452;
      r[116] = -0.703866201685267; s[116] =  0.245311234573305; w[116] = 0.028308893010452;
      r[117] =  0.245311234573305; s[117] = -0.703866201685267; w[117] = 0.028308893010452;
      r[118] = -0.703866201685267; s[118] = -0.541445032888038; w[118] = 0.028308893010452;
      r[119] = -0.541445032888038; s[119] = -0.703866201685267; w[119] = 0.028308893010452;
      r[120] = -0.348763754808032; s[120] = -0.034779976926618; w[120] = 0.029762759122336;
      r[121] = -0.034779976926618; s[121] = -0.348763754808032; w[121] = 0.029762759122336;
      r[122] = -0.616456268265350; s[122] = -0.034779976926618; w[122] = 0.029762759122336;
      r[123] = -0.034779976926618; s[123] = -0.616456268265350; w[123] = 0.029762759122336;
      r[124] = -0.616456268265350; s[124] = -0.348763754808032; w[124] = 0.029762759122336;
      r[125] = -0.348763754808032; s[125] = -0.616456268265350; w[125] = 0.029762759122336;
    
      break;
    }
    case 26: {
      r[0]   = -0.027285103318006; s[0]   = -0.027285103318006; w[0]   = 0.005361616431413;
      r[1]   = -0.945429793363989; s[1]   = -0.027285103318006; w[1]   = 0.005361616431413;
      r[2]   = -0.027285103318006; s[2]   = -0.945429793363989; w[2]   = 0.005361616431413;
      r[3]   = -0.068629018183452; s[3]   = -0.068629018183452; w[3]   = 0.012732191427031;
      r[4]   = -0.862741963633096; s[4]   = -0.068629018183452; w[4]   = 0.012732191427031;
      r[5]   = -0.068629018183452; s[5]   = -0.862741963633096; w[5]   = 0.012732191427031;
      r[6]   = -0.136174624247252; s[6]   = -0.136174624247252; w[6]   = 0.029707048808635;
      r[7]   = -0.727650751505495; s[7]   = -0.136174624247252; w[7]   = 0.029707048808635;
      r[8]   = -0.136174624247252; s[8]   = -0.727650751505495; w[8]   = 0.029707048808635;
      r[9]   = -0.214440393597162; s[9]   = -0.214440393597162; w[9]   = 0.026709487358142;
      r[10]  = -0.571119212805675; s[10]  = -0.214440393597162; w[10]  = 0.026709487358142;
      r[11]  = -0.214440393597162; s[11]  = -0.571119212805675; w[11]  = 0.026709487358142;
      r[12]  = -0.392927303733924; s[12]  = -0.392927303733924; w[12]  = 0.031691197689589;
      r[13]  = -0.214145392532152; s[13]  = -0.392927303733924; w[13]  = 0.031691197689589;
      r[14]  = -0.392927303733924; s[14]  = -0.214145392532152; w[14]  = 0.031691197689589;
      r[15]  = -0.541283541373341; s[15]  = -0.541283541373341; w[15]  = 0.028281596577704;
      r[16]  =  0.082567082746682; s[16]  = -0.541283541373341; w[16]  = 0.028281596577704;
      r[17]  = -0.541283541373341; s[17]  =  0.082567082746682; w[17]  = 0.028281596577704;
      r[18]  = -0.710171670270578; s[18]  = -0.710171670270578; w[18]  = 0.022460561269283;
      r[19]  =  0.420343340541156; s[19]  = -0.710171670270578; w[19]  = 0.022460561269283;
      r[20]  = -0.710171670270578; s[20]  =  0.420343340541156; w[20]  = 0.022460561269283;
      r[21]  = -0.832251071789951; s[21]  = -0.832251071789951; w[21]  = 0.013160382318625;
      r[22]  =  0.664502143579902; s[22]  = -0.832251071789951; w[22]  = 0.013160382318625;
      r[23]  = -0.832251071789951; s[23]  =  0.664502143579902; w[23]  = 0.013160382318625;
      r[24]  = -0.931928355163542; s[24]  = -0.931928355163542; w[24]  = 0.005828465225510;
      r[25]  =  0.863856710327084; s[25]  = -0.931928355163542; w[25]  = 0.005828465225510;
      r[26]  = -0.931928355163542; s[26]  =  0.863856710327084; w[26]  = 0.005828465225510;
      r[27]  = -0.987078078630829; s[27]  = -0.987078078630829; w[27]  = 0.001098762787251;
      r[28]  =  0.974156157261657; s[28]  = -0.987078078630829; w[28]  = 0.001098762787251;
      r[29]  = -0.987078078630829; s[29]  =  0.974156157261657; w[29]  = 0.001098762787251;
      r[30]  = -0.932113855389598; s[30]  =  0.919072286067742; w[30]  = 0.002557622023150;
      r[31]  =  0.919072286067742; s[31]  = -0.932113855389598; w[31]  = 0.002557622023150;
      r[32]  = -0.986958430678145; s[32]  =  0.919072286067742; w[32]  = 0.002557622023150;
      r[33]  =  0.919072286067742; s[33]  = -0.986958430678145; w[33]  = 0.002557622023150;
      r[34]  = -0.986958430678145; s[34]  = -0.932113855389597; w[34]  = 0.002557622023150;
      r[35]  = -0.932113855389597; s[35]  = -0.986958430678145; w[35]  = 0.002557622023150;
      r[36]  = -0.331832115907669; s[36]  =  0.318412220533634; w[36]  = 0.007088900226242;
      r[37]  =  0.318412220533634; s[37]  = -0.331832115907669; w[37]  = 0.007088900226242;
      r[38]  = -0.986580104625965; s[38]  =  0.318412220533634; w[38]  = 0.007088900226242;
      r[39]  =  0.318412220533634; s[39]  = -0.986580104625965; w[39]  = 0.007088900226242;
      r[40]  = -0.986580104625965; s[40]  = -0.331832115907669; w[40]  = 0.007088900226242;
      r[41]  = -0.331832115907669; s[41]  = -0.986580104625965; w[41]  = 0.007088900226242;
      r[42]  = -0.117284648816543; s[42]  =  0.105635324925259; w[42]  = 0.006588000446609;
      r[43]  =  0.105635324925259; s[43]  = -0.117284648816543; w[43]  = 0.006588000446609;
      r[44]  = -0.988350676108717; s[44]  =  0.105635324925259; w[44]  = 0.006588000446609;
      r[45]  =  0.105635324925259; s[45]  = -0.988350676108717; w[45]  = 0.006588000446609;
      r[46]  = -0.988350676108717; s[46]  = -0.117284648816543; w[46]  = 0.006588000446609;
      r[47]  = -0.117284648816543; s[47]  = -0.988350676108717; w[47]  = 0.006588000446609;
      r[48]  = -0.528179583523521; s[48]  =  0.515363521190658; w[48]  = 0.006069708653929;
      r[49]  =  0.515363521190658; s[49]  = -0.528179583523521; w[49]  = 0.006069708653929;
      r[50]  = -0.987183937667137; s[50]  =  0.515363521190658; w[50]  = 0.006069708653929;
      r[51]  =  0.515363521190658; s[51]  = -0.987183937667137; w[51]  = 0.006069708653929;
      r[52]  = -0.987183937667137; s[52]  = -0.528179583523521; w[52]  = 0.006069708653929;
      r[53]  = -0.528179583523521; s[53]  = -0.987183937667137; w[53]  = 0.006069708653929;
      r[54]  = -0.698183069548291; s[54]  =  0.685666877859863; w[54]  = 0.004942128706706;
      r[55]  =  0.685666877859863; s[55]  = -0.698183069548291; w[55]  = 0.004942128706706;
      r[56]  = -0.987483808311572; s[56]  =  0.685666877859863; w[56]  = 0.004942128706706;
      r[57]  =  0.685666877859863; s[57]  = -0.987483808311572; w[57]  = 0.004942128706706;
      r[58]  = -0.987483808311572; s[58]  = -0.698183069548291; w[58]  = 0.004942128706706;
      r[59]  = -0.698183069548291; s[59]  = -0.987483808311572; w[59]  = 0.004942128706706;
      r[60]  = -0.834663613179069; s[60]  =  0.821127985248516; w[60]  = 0.004076340070443;
      r[61]  =  0.821127985248516; s[61]  = -0.834663613179069; w[61]  = 0.004076340070443;
      r[62]  = -0.986464372069447; s[62]  =  0.821127985248516; w[62]  = 0.004076340070443;
      r[63]  =  0.821127985248516; s[63]  = -0.986464372069447; w[63]  = 0.004076340070443;
      r[64]  = -0.986464372069447; s[64]  = -0.834663613179069; w[64]  = 0.004076340070443;
      r[65]  = -0.834663613179069; s[65]  = -0.986464372069447; w[65]  = 0.004076340070443;
      r[66]  = -0.700166221086294; s[66]  =  0.634916644768372; w[66]  = 0.011120280239045;
      r[67]  =  0.634916644768372; s[67]  = -0.700166221086294; w[67]  = 0.011120280239045;
      r[68]  = -0.934750423682077; s[68]  =  0.634916644768372; w[68]  = 0.011120280239045;
      r[69]  =  0.634916644768372; s[69]  = -0.934750423682077; w[69]  = 0.011120280239045;
      r[70]  = -0.934750423682078; s[70]  = -0.700166221086294; w[70]  = 0.011120280239045;
      r[71]  = -0.700166221086294; s[71]  = -0.934750423682078; w[71]  = 0.011120280239045;
      r[72]  = -0.160245113591115; s[72]  =  0.096801406370527; w[72]  = 0.012755358640694;
      r[73]  =  0.096801406370527; s[73]  = -0.160245113591115; w[73]  = 0.012755358640694;
      r[74]  = -0.936556292779412; s[74]  =  0.096801406370527; w[74]  = 0.012755358640694;
      r[75]  =  0.096801406370527; s[75]  = -0.936556292779412; w[75]  = 0.012755358640694;
      r[76]  = -0.936556292779412; s[76]  = -0.160245113591115; w[76]  = 0.012755358640694;
      r[77]  = -0.160245113591115; s[77]  = -0.936556292779412; w[77]  = 0.012755358640694;
      r[78]  = -0.348814552390868; s[78]  =  0.279444440434072; w[78]  = 0.014549157045387;
      r[79]  =  0.279444440434072; s[79]  = -0.348814552390868; w[79]  = 0.014549157045387;
      r[80]  = -0.930629888043204; s[80]  =  0.279444440434072; w[80]  = 0.014549157045387;
      r[81]  =  0.279444440434072; s[81]  = -0.930629888043204; w[81]  = 0.014549157045387;
      r[82]  = -0.930629888043204; s[82]  = -0.348814552390868; w[82]  = 0.014549157045387;
      r[83]  = -0.348814552390868; s[83]  = -0.930629888043204; w[83]  = 0.014549157045387;
      r[84]  = -0.835149343614232; s[84]  =  0.764834325997890; w[84]  = 0.009121880355576;
      r[85]  =  0.764834325997890; s[85]  = -0.835149343614232; w[85]  = 0.009121880355576;
      r[86]  = -0.929684982383659; s[86]  =  0.764834325997890; w[86]  = 0.009121880355576;
      r[87]  =  0.764834325997890; s[87]  = -0.929684982383659; w[87]  = 0.009121880355576;
      r[88]  = -0.929684982383659; s[88]  = -0.835149343614232; w[88]  = 0.009121880355576;
      r[89]  = -0.835149343614232; s[89]  = -0.929684982383659; w[89]  = 0.009121880355576;
      r[90]  = -0.534046892746188; s[90]  =  0.466663082505331; w[90]  = 0.013582808914065;
      r[91]  =  0.466663082505331; s[91]  = -0.534046892746188; w[91]  = 0.013582808914065;
      r[92]  = -0.932616189759143; s[92]  =  0.466663082505331; w[92]  = 0.013582808914065;
      r[93]  =  0.466663082505331; s[93]  = -0.932616189759143; w[93]  = 0.013582808914065;
      r[94]  = -0.932616189759143; s[94]  = -0.534046892746188; w[94]  = 0.013582808914065;
      r[95]  = -0.534046892746188; s[95]  = -0.932616189759143; w[95]  = 0.013582808914065;
      r[96]  = -0.704373154566392; s[96]  =  0.546726033293271; w[96]  = 0.015585248954434;
      r[97]  =  0.546726033293271; s[97]  = -0.704373154566392; w[97]  = 0.015585248954434;
      r[98]  = -0.842352878726879; s[98]  =  0.546726033293271; w[98]  = 0.015585248954434;
      r[99]  =  0.546726033293271; s[99]  = -0.842352878726879; w[99]  = 0.015585248954434;
      r[100] = -0.842352878726879; s[100] = -0.704373154566392; w[100] = 0.015585248954434;
      r[101] = -0.704373154566392; s[101] = -0.842352878726879; w[101] = 0.015585248954434;
      r[102] = -0.220224863280788; s[102] =  0.060870104244893; w[102] = 0.018960151369388;
      r[103] =  0.060870104244893; s[103] = -0.220224863280788; w[103] = 0.018960151369388;
      r[104] = -0.840645240964104; s[104] =  0.060870104244893; w[104] = 0.018960151369388;
      r[105] =  0.060870104244893; s[105] = -0.840645240964104; w[105] = 0.018960151369388;
      r[106] = -0.840645240964104; s[106] = -0.220224863280788; w[106] = 0.018960151369388;
      r[107] = -0.220224863280788; s[107] = -0.840645240964104; w[107] = 0.018960151369388;
      r[108] = -0.394989254221112; s[108] =  0.230609434358343; w[108] = 0.019206200353033;
      r[109] =  0.230609434358343; s[109] = -0.394989254221112; w[109] = 0.019206200353033;
      r[110] = -0.835620180137232; s[110] =  0.230609434358343; w[110] = 0.019206200353033;
      r[111] =  0.230609434358343; s[111] = -0.835620180137232; w[111] = 0.019206200353033;
      r[112] = -0.835620180137232; s[112] = -0.394989254221112; w[112] = 0.019206200353033;
      r[113] = -0.394989254221112; s[113] = -0.835620180137232; w[113] = 0.019206200353033;
      r[114] = -0.557563751773578; s[114] =  0.393274056332418; w[114] = 0.018121839059318;
      r[115] =  0.393274056332418; s[115] = -0.557563751773578; w[115] = 0.018121839059318;
      r[116] = -0.835710304558840; s[116] =  0.393274056332418; w[116] = 0.018121839059318;
      r[117] =  0.393274056332418; s[117] = -0.835710304558840; w[117] = 0.018121839059318;
      r[118] = -0.835710304558840; s[118] = -0.557563751773578; w[118] = 0.018121839059318;
      r[119] = -0.557563751773578; s[119] = -0.835710304558840; w[119] = 0.018121839059318;
      r[120] = -0.350522450436603; s[120] =  0.060630141977535; w[120] = 0.028403705951863;
      r[121] =  0.060630141977535; s[121] = -0.350522450436603; w[121] = 0.028403705951863;
      r[122] = -0.710107691540932; s[122] =  0.060630141977535; w[122] = 0.028403705951863;
      r[123] =  0.060630141977535; s[123] = -0.710107691540932; w[123] = 0.028403705951863;
      r[124] = -0.710107691540932; s[124] = -0.350522450436603; w[124] = 0.028403705951863;
      r[125] = -0.350522450436603; s[125] = -0.710107691540932; w[125] = 0.028403705951863;
      r[126] = -0.544246685468649; s[126] =  0.246315777233803; w[126] = 0.026986162760930;
      r[127] =  0.246315777233803; s[127] = -0.544246685468649; w[127] = 0.026986162760930;
      r[128] = -0.702069091765154; s[128] =  0.246315777233803; w[128] = 0.026986162760930;
      r[129] =  0.246315777233803; s[129] = -0.702069091765154; w[129] = 0.026986162760930;
      r[130] = -0.702069091765154; s[130] = -0.544246685468649; w[130] = 0.026986162760930;
      r[131] = -0.544246685468649; s[131] = -0.702069091765154; w[131] = 0.026986162760930;
      r[132] = -0.373996265063702; s[132] = -0.064072130507071; w[132] = 0.025102184615930;
      r[133] = -0.064072130507071; s[133] = -0.373996265063702; w[133] = 0.025102184615930;
      r[134] = -0.561931604429227; s[134] = -0.064072130507071; w[134] = 0.025102184615930;
      r[135] = -0.064072130507071; s[135] = -0.561931604429227; w[135] = 0.025102184615930;
      r[136] = -0.561931604429227; s[136] = -0.373996265063702; w[136] = 0.025102184615930;
      r[137] = -0.373996265063702; s[137] = -0.561931604429227; w[137] = 0.025102184615930;
    
      break;
    }
    case 27: {
      r[0]   = -0.333333333333333; s[0]   = -0.333333333333333; w[0]   = 0.029215421142604;
      r[1]   = -0.029244446609137; s[1]   = -0.029244446609137; w[1]   = 0.013615273717612;
      r[2]   = -0.941511106781726; s[2]   = -0.029244446609137; w[2]   = 0.013615273717612;
      r[3]   = -0.029244446609137; s[3]   = -0.941511106781726; w[3]   = 0.013615273717612;
      r[4]   = -0.132843721814340; s[4]   = -0.132843721814340; w[4]   = 0.030951102870660;
      r[5]   = -0.734312556371321; s[5]   = -0.132843721814340; w[5]   = 0.030951102870660;
      r[6]   = -0.132843721814340; s[6]   = -0.734312556371321; w[6]   = 0.030951102870660;
      r[7]   = -0.224940014502969; s[7]   = -0.224940014502969; w[7]   = 0.031733186724762;
      r[8]   = -0.550119970994063; s[8]   = -0.224940014502969; w[8]   = 0.031733186724762;
      r[9]   = -0.224940014502969; s[9]   = -0.550119970994063; w[9]   = 0.031733186724762;
      r[10]  = -0.427786783949924; s[10]  = -0.427786783949924; w[10]  = 0.031559855928489;
      r[11]  = -0.144426432100151; s[11]  = -0.427786783949924; w[11]  = 0.031559855928489;
      r[12]  = -0.427786783949924; s[12]  = -0.144426432100151; w[12]  = 0.031559855928489;
      r[13]  = -0.532776206162091; s[13]  = -0.532776206162091; w[13]  = 0.027863700370354;
      r[14]  =  0.065552412324182; s[14]  = -0.532776206162091; w[14]  = 0.027863700370354;
      r[15]  = -0.532776206162091; s[15]  =  0.065552412324182; w[15]  = 0.027863700370354;
      r[16]  = -0.690368113175391; s[16]  = -0.690368113175391; w[16]  = 0.019750305737275;
      r[17]  =  0.380736226350782; s[17]  = -0.690368113175391; w[17]  = 0.019750305737275;
      r[18]  = -0.690368113175391; s[18]  =  0.380736226350782; w[18]  = 0.019750305737275;
      r[19]  = -0.781499070170348; s[19]  = -0.781499070170348; w[19]  = 0.015635498945203;
      r[20]  =  0.562998140340695; s[20]  = -0.781499070170348; w[20]  = 0.015635498945203;
      r[21]  = -0.781499070170348; s[21]  =  0.562998140340695; w[21]  = 0.015635498945203;
      r[22]  = -0.862704525220276; s[22]  = -0.862704525220276; w[22]  = 0.010432865446986;
      r[23]  =  0.725409050440551; s[23]  = -0.862704525220276; w[23]  = 0.010432865446986;
      r[24]  = -0.862704525220276; s[24]  =  0.725409050440551; w[24]  = 0.010432865446986;
      r[25]  = -0.932961091535336; s[25]  = -0.932961091535336; w[25]  = 0.005852724472600;
      r[26]  =  0.865922183070671; s[26]  = -0.932961091535336; w[26]  = 0.005852724472600;
      r[27]  = -0.932961091535336; s[27]  =  0.865922183070671; w[27]  = 0.005852724472600;
      r[28]  = -0.986811381964774; s[28]  = -0.986811381964774; w[28]  = 0.001134879007165;
      r[29]  =  0.973622763929549; s[29]  = -0.986811381964774; w[29]  = 0.001134879007165;
      r[30]  = -0.986811381964774; s[30]  =  0.973622763929549; w[30]  = 0.001134879007165;
      r[31]  = -0.717055503863097; s[31]  =  0.716419063026126; w[31]  = 0.001375795459783;
      r[32]  =  0.716419063026126; s[32]  = -0.717055503863097; w[32]  = 0.001375795459783;
      r[33]  = -0.999363559163030; s[33]  =  0.716419063026126; w[33]  = 0.001375795459783;
      r[34]  =  0.716419063026126; s[34]  = -0.999363559163030; w[34]  = 0.001375795459783;
      r[35]  = -0.999363559163030; s[35]  = -0.717055503863097; w[35]  = 0.001375795459783;
      r[36]  = -0.717055503863097; s[36]  = -0.999363559163030; w[36]  = 0.001375795459783;
      r[37]  = -0.838822221132805; s[37]  =  0.827544031058903; w[37]  = 0.003088879415937;
      r[38]  =  0.827544031058903; s[38]  = -0.838822221132805; w[38]  = 0.003088879415937;
      r[39]  = -0.988721809926098; s[39]  =  0.827544031058903; w[39]  = 0.003088879415937;
      r[40]  =  0.827544031058903; s[40]  = -0.988721809926098; w[40]  = 0.003088879415937;
      r[41]  = -0.988721809926098; s[41]  = -0.838822221132805; w[41]  = 0.003088879415937;
      r[42]  = -0.838822221132805; s[42]  = -0.988721809926098; w[42]  = 0.003088879415937;
      r[43]  = -0.537630806497654; s[43]  =  0.526128943172145; w[43]  = 0.005625149324346;
      r[44]  =  0.526128943172145; s[44]  = -0.537630806497654; w[44]  = 0.005625149324346;
      r[45]  = -0.988498136674491; s[45]  =  0.526128943172145; w[45]  = 0.005625149324346;
      r[46]  =  0.526128943172145; s[46]  = -0.988498136674491; w[46]  = 0.005625149324346;
      r[47]  = -0.988498136674491; s[47]  = -0.537630806497654; w[47]  = 0.005625149324346;
      r[48]  = -0.537630806497654; s[48]  = -0.988498136674491; w[48]  = 0.005625149324346;
      r[49]  = -0.113842988120976; s[49]  =  0.103273580134309; w[49]  = 0.005913105372899;
      r[50]  =  0.103273580134309; s[50]  = -0.113842988120976; w[50]  = 0.005913105372899;
      r[51]  = -0.989430592013333; s[51]  =  0.103273580134310; w[51]  = 0.005913105372899;
      r[52]  =  0.103273580134310; s[52]  = -0.989430592013333; w[52]  = 0.005913105372899;
      r[53]  = -0.989430592013333; s[53]  = -0.113842988120976; w[53]  = 0.005913105372899;
      r[54]  = -0.113842988120976; s[54]  = -0.989430592013333; w[54]  = 0.005913105372899;
      r[55]  = -0.330849908839651; s[55]  =  0.320515152715889; w[55]  = 0.005565783438844;
      r[56]  =  0.320515152715889; s[56]  = -0.330849908839651; w[56]  = 0.005565783438844;
      r[57]  = -0.989665243876238; s[57]  =  0.320515152715889; w[57]  = 0.005565783438844;
      r[58]  =  0.320515152715889; s[58]  = -0.989665243876238; w[58]  = 0.005565783438844;
      r[59]  = -0.989665243876238; s[59]  = -0.330849908839651; w[59]  = 0.005565783438844;
      r[60]  = -0.330849908839651; s[60]  = -0.989665243876238; w[60]  = 0.005565783438844;
      r[61]  = -0.931869179489789; s[61]  =  0.919222453867072; w[61]  = 0.002449540456172;
      r[62]  =  0.919222453867072; s[62]  = -0.931869179489789; w[62]  = 0.002449540456172;
      r[63]  = -0.987353274377283; s[63]  =  0.919222453867072; w[63]  = 0.002449540456172;
      r[64]  =  0.919222453867072; s[64]  = -0.987353274377283; w[64]  = 0.002449540456172;
      r[65]  = -0.987353274377283; s[65]  = -0.931869179489789; w[65]  = 0.002449540456172;
      r[66]  = -0.931869179489789; s[66]  = -0.987353274377283; w[66]  = 0.002449540456172;
      r[67]  = -0.695123977280935; s[67]  =  0.663705924287710; w[67]  = 0.007906643758942;
      r[68]  =  0.663705924287710; s[68]  = -0.695123977280935; w[68]  = 0.007906643758942;
      r[69]  = -0.968581947006775; s[69]  =  0.663705924287710; w[69]  = 0.007906643758942;
      r[70]  =  0.663705924287710; s[70]  = -0.968581947006775; w[70]  = 0.007906643758942;
      r[71]  = -0.968581947006775; s[71]  = -0.695123977280935; w[71]  = 0.007906643758942;
      r[72]  = -0.695123977280935; s[72]  = -0.968581947006775; w[72]  = 0.007906643758942;
      r[73]  = -0.831869809255826; s[73]  =  0.777736654895860; w[73]  = 0.007221969465327;
      r[74]  =  0.777736654895860; s[74]  = -0.831869809255826; w[74]  = 0.007221969465327;
      r[75]  = -0.945866845640035; s[75]  =  0.777736654895860; w[75]  = 0.007221969465327;
      r[76]  =  0.777736654895860; s[76]  = -0.945866845640035; w[76]  = 0.007221969465327;
      r[77]  = -0.945866845640035; s[77]  = -0.831869809255826; w[77]  = 0.007221969465327;
      r[78]  = -0.831869809255826; s[78]  = -0.945866845640035; w[78]  = 0.007221969465327;
      r[79]  = -0.423211420319719; s[79]  =  0.374995160760160; w[79]  = 0.008707930088036;
      r[80]  =  0.374995160760160; s[80]  = -0.423211420319719; w[80]  = 0.008707930088036;
      r[81]  = -0.951783740440441; s[81]  =  0.374995160760160; w[81]  = 0.008707930088036;
      r[82]  =  0.374995160760160; s[82]  = -0.951783740440441; w[82]  = 0.008707930088036;
      r[83]  = -0.951783740440441; s[83]  = -0.423211420319719; w[83]  = 0.008707930088036;
      r[84]  = -0.423211420319719; s[84]  = -0.951783740440441; w[84]  = 0.008707930088036;
      r[85]  = -0.223962554104902; s[85]  =  0.171509638849087; w[85]  = 0.011256464740666;
      r[86]  =  0.171509638849087; s[86]  = -0.223962554104902; w[86]  = 0.011256464740666;
      r[87]  = -0.947547084744185; s[87]  =  0.171509638849087; w[87]  = 0.011256464740666;
      r[88]  =  0.171509638849087; s[88]  = -0.947547084744185; w[88]  = 0.011256464740666;
      r[89]  = -0.947547084744185; s[89]  = -0.223962554104902; w[89]  = 0.011256464740666;
      r[90]  = -0.223962554104902; s[90]  = -0.947547084744185; w[90]  = 0.011256464740666;
      r[91]  = -0.566926820194262; s[91]  =  0.490923426439871; w[91]  = 0.013269437724243;
      r[92]  =  0.490923426439871; s[92]  = -0.566926820194262; w[92]  = 0.013269437724243;
      r[93]  = -0.923996606245609; s[93]  =  0.490923426439871; w[93]  = 0.013269437724243;
      r[94]  =  0.490923426439871; s[94]  = -0.923996606245609; w[94]  = 0.013269437724243;
      r[95]  = -0.923996606245609; s[95]  = -0.566926820194262; w[95]  = 0.013269437724243;
      r[96]  = -0.566926820194262; s[96]  = -0.923996606245609; w[96]  = 0.013269437724243;
      r[97]  = -0.350787336349322; s[97]  =  0.245677528944182; w[97]  = 0.013046146795576;
      r[98]  =  0.245677528944182; s[98]  = -0.350787336349322; w[98]  = 0.013046146795576;
      r[99]  = -0.894890192594860; s[99]  =  0.245677528944182; w[99]  = 0.013046146795576;
      r[100] =  0.245677528944182; s[100] = -0.894890192594860; w[100] = 0.013046146795576;
      r[101] = -0.894890192594860; s[101] = -0.350787336349322; w[101] = 0.013046146795576;
      r[102] = -0.350787336349322; s[102] = -0.894890192594860; w[102] = 0.013046146795576;
      r[103] = -0.729613610936203; s[103] =  0.622255871516285; w[103] = 0.012947418329903;
      r[104] =  0.622255871516285; s[104] = -0.729613610936203; w[104] = 0.012947418329903;
      r[105] = -0.892642260580081; s[105] =  0.622255871516285; w[105] = 0.012947418329903;
      r[106] =  0.622255871516285; s[106] = -0.892642260580081; w[106] = 0.012947418329903;
      r[107] = -0.892642260580081; s[107] = -0.729613610936203; w[107] = 0.012947418329903;
      r[108] = -0.729613610936203; s[108] = -0.892642260580081; w[108] = 0.012947418329903;
      r[109] = -0.176209680235251; s[109] =  0.035410146457445; w[109] = 0.021313811627966;
      r[110] =  0.035410146457445; s[110] = -0.176209680235251; w[110] = 0.021313811627966;
      r[111] = -0.859200466222194; s[111] =  0.035410146457445; w[111] = 0.021313811627966;
      r[112] =  0.035410146457445; s[112] = -0.859200466222194; w[112] = 0.021313811627966;
      r[113] = -0.859200466222194; s[113] = -0.176209680235251; w[113] = 0.021313811627966;
      r[114] = -0.176209680235251; s[114] = -0.859200466222194; w[114] = 0.021313811627966;
      r[115] = -0.469171656648510; s[115] =  0.304686011244002; w[115] = 0.017176891528936;
      r[116] =  0.304686011244002; s[116] = -0.469171656648510; w[116] = 0.017176891528936;
      r[117] = -0.835514354595492; s[117] =  0.304686011244002; w[117] = 0.017176891528936;
      r[118] =  0.304686011244002; s[118] = -0.835514354595492; w[118] = 0.017176891528936;
      r[119] = -0.835514354595492; s[119] = -0.469171656648510; w[119] = 0.017176891528936;
      r[120] = -0.469171656648510; s[120] = -0.835514354595492; w[120] = 0.017176891528936;
      r[121] = -0.623736093891471; s[121] =  0.443479545725970; w[121] = 0.017763399981969;
      r[122] =  0.443479545725970; s[122] = -0.623736093891471; w[122] = 0.017763399981969;
      r[123] = -0.819743451834499; s[123] =  0.443479545725970; w[123] = 0.017763399981969;
      r[124] =  0.443479545725970; s[124] = -0.819743451834499; w[124] = 0.017763399981969;
      r[125] = -0.819743451834499; s[125] = -0.623736093891471; w[125] = 0.017763399981969;
      r[126] = -0.623736093891471; s[126] = -0.819743451834499; w[126] = 0.017763399981969;
      r[127] = -0.336632062394500; s[127] =  0.104847567508036; w[127] = 0.023828398561090;
      r[128] =  0.104847567508036; s[128] = -0.336632062394500; w[128] = 0.023828398561090;
      r[129] = -0.768215505113536; s[129] =  0.104847567508036; w[129] = 0.023828398561090;
      r[130] =  0.104847567508036; s[130] = -0.768215505113536; w[130] = 0.023828398561090;
      r[131] = -0.768215505113536; s[131] = -0.336632062394500; w[131] = 0.023828398561090;
      r[132] = -0.336632062394500; s[132] = -0.768215505113536; w[132] = 0.023828398561090;
      r[133] = -0.530404689974552; s[133] =  0.224195004863778; w[133] = 0.026764788789016;
      r[134] =  0.224195004863778; s[134] = -0.530404689974552; w[134] = 0.026764788789016;
      r[135] = -0.693790314889227; s[135] =  0.224195004863778; w[135] = 0.026764788789016;
      r[136] =  0.224195004863778; s[136] = -0.693790314889227; w[136] = 0.026764788789016;
      r[137] = -0.693790314889227; s[137] = -0.530404689974552; w[137] = 0.026764788789016;
      r[138] = -0.530404689974552; s[138] = -0.693790314889227; w[138] = 0.026764788789016;
      r[139] = -0.350831684137242; s[139] = -0.019508133770036; w[139] = 0.028977845006028;
      r[140] = -0.019508133770036; s[140] = -0.350831684137242; w[140] = 0.028977845006028;
      r[141] = -0.629660182092722; s[141] = -0.019508133770036; w[141] = 0.028977845006028;
      r[142] = -0.019508133770036; s[142] = -0.629660182092722; w[142] = 0.028977845006028;
      r[143] = -0.629660182092722; s[143] = -0.350831684137242; w[143] = 0.028977845006028;
      r[144] = -0.350831684137242; s[144] = -0.629660182092722; w[144] = 0.028977845006028;
    
      break;
    }
    case 28: {
      r[0]   = -0.988064607832230; s[0]   =  0.987992518020485; w[0]   = 0.000005678109445;
      r[1]   = -0.937649986705898; s[1]   =  0.937273392400706; w[1]   = 0.000067869690329;
      r[2]   = -0.849117911767581; s[2]   =  0.848206583410427; w[2]   = 0.000250117104238;
      r[3]   = -0.726072255922453; s[3]   =  0.724417731360170; w[3]   = 0.000591434111178;
      r[4]   = -0.573547944561595; s[4]   =  0.570972172608539; w[4]   = 0.001096877668021;
      r[5]   = -0.397788705468703; s[5]   =  0.394151347077563; w[5]   = 0.001734258212849;
      r[6]   = -0.205989917758162; s[6]   =  0.201194093997435; w[6]   = 0.002437321168702;
      r[7]   = -0.006003740989757; s[7]   =  0.000000000000000; w[7]   = 0.003114968848543;
      r[8]   =  0.193982435778648; s[8]   = -0.201194093997435; w[8]   = 0.003665090319212;
      r[9]   =  0.385781223489189; s[9]   = -0.394151347077563; w[9]   = 0.003990796070869;
      r[10]  =  0.561540462582081; s[10]  = -0.570972172608539; w[10]  = 0.004016439454972;
      r[11]  =  0.714064773942939; s[11]  = -0.724417731360170; w[11]  = 0.003700816722647;
      r[12]  =  0.837110429788067; s[12]  = -0.848206583410427; w[12]  = 0.003045376334904;
      r[13]  =  0.925642504726384; s[13]  = -0.937273392400706; w[13]  = 0.002096114396383;
      r[14]  =  0.976057125852715; s[14]  = -0.987992518020485; w[14]  = 0.000940083783827;
      r[15]  = -0.988369112325678; s[15]  =  0.987992518020485; w[15]  = 0.000012991999963;
      r[16]  = -0.939240706051164; s[16]  =  0.937273392400706; w[16]  = 0.000155291655212;
      r[17]  = -0.852967326449712; s[17]  =  0.848206583410427; w[17]  = 0.000572289322752;
      r[18]  = -0.733060901773317; s[18]  =  0.724417731360170; w[18]  = 0.001353251821662;
      r[19]  = -0.584427902697520; s[19]  =  0.570972172608539; w[19]  = 0.002509749901695;
      r[20]  = -0.413152762435777; s[20]  =  0.394151347077563; w[20]  = 0.003968131092565;
      r[21]  = -0.226247286304345; s[21]  =  0.201194093997435; w[21]  = 0.005576799256556;
      r[22]  = -0.031363303799647; s[22]  =  0.000000000000000; w[22]  = 0.007127315095696;
      r[23]  =  0.163520678705051; s[23]  = -0.201194093997435; w[23]  = 0.008386040063107;
      r[24]  =  0.350426154836483; s[24]  = -0.394151347077563; w[24]  = 0.009131282674963;
      r[25]  =  0.521701295098226; s[25]  = -0.570972172608539; w[25]  = 0.009189956930633;
      r[26]  =  0.670334294174023; s[26]  = -0.724417731360170; w[26]  = 0.008467785129237;
      r[27]  =  0.790240718850418; s[27]  = -0.848206583410427; w[27]  = 0.006968081473429;
      r[28]  =  0.876514098451870; s[28]  = -0.937273392400706; w[28]  = 0.004796088983887;
      r[29]  =  0.925642504726384; s[29]  = -0.987992518020485; w[29]  = 0.002150992086750;
      r[30]  = -0.988903846377640; s[30]  =  0.987992518020485; w[30]  = 0.000019785289042;
      r[31]  = -0.942034135439991; s[31]  =  0.937273392400706; w[31]  = 0.000236490940043;
      r[32]  = -0.859727204070395; s[32]  =  0.848206583410427; w[32]  = 0.000871529379538;
      r[33]  = -0.745333518414343; s[33]  =  0.724417731360170; w[33]  = 0.002060843481789;
      r[34]  = -0.603533972474414; s[34]  =  0.570972172608539; w[34]  = 0.003822054138805;
      r[35]  = -0.440133265559207; s[35]  =  0.394151347077563; w[35]  = 0.006042997294437;
      r[36]  = -0.261820832829464; s[36]  =  0.201194093997435; w[36]  = 0.008492809847469;
      r[37]  = -0.075896708294786; s[37]  =  0.000000000000000; w[37]  = 0.010854063244177;
      r[38]  =  0.110027416239891; s[38]  = -0.201194093997435; w[38]  = 0.012770953436328;
      r[39]  =  0.288339848969634; s[39]  = -0.394151347077563; w[39]  = 0.013905870348620;
      r[40]  =  0.451740555884842; s[40]  = -0.570972172608539; w[40]  = 0.013995224344240;
      r[41]  =  0.593540101824770; s[41]  = -0.724417731360170; w[41]  = 0.012895441564853;
      r[42]  =  0.707933787480822; s[42]  = -0.848206583410427; w[42]  = 0.010611569151594;
      r[43]  =  0.790240718850418; s[43]  = -0.937273392400706; w[43]  = 0.007303879856139;
      r[44]  =  0.837110429788067; s[44]  = -0.987992518020485; w[44]  = 0.003275708150100;
      r[45]  = -0.989647042582769; s[45]  =  0.987992518020485; w[45]  = 0.000025769562269;
      r[46]  = -0.945916562813853; s[46]  =  0.937273392400706; w[46]  = 0.000308020165519;
      r[47]  = -0.869122370464600; s[47]  =  0.848206583410427; w[47]  = 0.001135132803359;
      r[48]  = -0.762390524754508; s[48]  =  0.724417731360170; w[48]  = 0.002684167732828;
      r[49]  = -0.630088403599617; s[49]  =  0.570972172608539; w[49]  = 0.004978075474027;
      r[50]  = -0.477631920189938; s[50]  =  0.394151347077563; w[50]  = 0.007870766747028;
      r[51]  = -0.311262465886975; s[51]  =  0.201194093997435; w[51]  = 0.011061551425453;
      r[52]  = -0.137791134319915; s[52]  =  0.000000000000000; w[52]  = 0.014136991279319;
      r[53]  =  0.035680197247145; s[53]  = -0.201194093997435; w[53]  = 0.016633665503545;
      r[54]  =  0.202049651550108; s[54]  = -0.394151347077563; w[54]  = 0.018111850228555;
      r[55]  =  0.354506134959787; s[55]  = -0.570972172608539; w[55]  = 0.018228230300094;
      r[56]  =  0.486808256114678; s[56]  = -0.724417731360170; w[56]  = 0.016795806403939;
      r[57]  =  0.593540101824770; s[57]  = -0.848206583410427; w[57]  = 0.013821152243283;
      r[58]  =  0.670334294174023; s[58]  = -0.937273392400706; w[58]  = 0.009513016785380;
      r[59]  =  0.714064773942939; s[59]  = -0.987992518020485; w[59]  = 0.004266481271556;
      r[60]  = -0.990568289973542; s[60]  =  0.987992518020485; w[60]  = 0.000030699031605;
      r[61]  = -0.950729122489687; s[61]  =  0.937273392400706; w[61]  = 0.000366941459749;
      r[62]  = -0.880768383276303; s[62]  =  0.848206583410427; w[62]  = 0.001352272787635;
      r[63]  = -0.783533962351248; s[63]  =  0.724417731360170; w[63]  = 0.003197623195991;
      r[64]  = -0.663004610946658; s[64]  =  0.570972172608539; w[64]  = 0.005930333418609;
      r[65]  = -0.524114312723242; s[65]  =  0.394151347077563; w[65]  = 0.009376368701824;
      r[66]  = -0.372549075177309; s[66]  =  0.201194093997435; w[66]  = 0.013177519790991;
      r[67]  = -0.214513913695731; s[67]  =  0.000000000000000; w[67]  = 0.016841261700383;
      r[68]  = -0.056478752214152; s[68]  = -0.201194093997435; w[68]  = 0.019815525683434;
      r[69]  =  0.095086485331781; s[69]  = -0.394151347077563; w[69]  = 0.021576472924861;
      r[70]  =  0.233976783555196; s[70]  = -0.570972172608539; w[70]  = 0.021715115384404;
      r[71]  =  0.354506134959787; s[71]  = -0.724417731360170; w[71]  = 0.020008682578130;
      r[72]  =  0.451740555884842; s[72]  = -0.848206583410427; w[72]  = 0.016465005695410;
      r[73]  =  0.521701295098226; s[73]  = -0.937273392400706; w[73]  = 0.011332765372579;
      r[74]  =  0.561540462582081; s[74]  = -0.987992518020485; w[74]  = 0.005082618091388;
      r[75]  = -0.991629876411625; s[75]  =  0.987992518020485; w[75]  = 0.000034371743072;
      r[76]  = -0.956274807758919; s[76]  =  0.937273392400706; w[76]  = 0.000410840893588;
      r[77]  = -0.894188501892071; s[77]  =  0.848206583410427; w[77]  = 0.001514053388316;
      r[78]  = -0.807898304472545; s[78]  =  0.724417731360170; w[78]  = 0.003580174265665;
      r[79]  = -0.700935138254217; s[79]  =  0.570972172608539; w[79]  = 0.006639815197344;
      r[80]  = -0.577677642201529; s[80]  =  0.394151347077563; w[80]  = 0.010498120595869;
      r[81]  = -0.443171835046505; s[81]  =  0.201194093997435; w[81]  = 0.014754026459451;
      r[82]  = -0.302924326461218; s[82]  =  0.000000000000000; w[82]  = 0.018856084049129;
      r[83]  = -0.162676817875932; s[83]  = -0.201194093997435; w[83]  = 0.022186177283618;
      r[84]  = -0.028171010720908; s[84]  = -0.394151347077563; w[84]  = 0.024157797330925;
      r[85]  =  0.095086485331781; s[85]  = -0.570972172608539; w[85]  = 0.024313026429341;
      r[86]  =  0.202049651550108; s[86]  = -0.724417731360170; w[86]  = 0.022402442709917;
      r[87]  =  0.288339848969634; s[87]  = -0.848206583410427; w[87]  = 0.018434814254740;
      r[88]  =  0.350426154836483; s[88]  = -0.937273392400706; w[88]  = 0.012688572873940;
      r[89]  =  0.385781223489189; s[89]  = -0.987992518020485; w[89]  = 0.005690682540646;
      r[90]  = -0.992788341781213; s[90]  =  0.987992518020485; w[90]  = 0.000036637297986;
      r[91]  = -0.962326584707617; s[91]  =  0.937273392400706; w[91]  = 0.000437920771488;
      r[92]  = -0.908833322242456; s[92]  =  0.848206583410427; w[92]  = 0.001613849639200;
      r[93]  = -0.834486103249711; s[93]  =  0.724417731360170; w[93]  = 0.003816155355883;
      r[94]  = -0.742327153788413; s[94]  =  0.570972172608539; w[94]  = 0.007077467309460;
      r[95]  = -0.636129088126634; s[95]  =  0.394151347077563; w[95]  = 0.011190086338209;
      r[96]  = -0.520239531729724; s[96]  =  0.201194093997435; w[96]  = 0.015726512989614;
      r[97]  = -0.399402953001283; s[97]  =  0.000000000000000; w[97]  = 0.020098950720122;
      r[98]  = -0.278566374272841; s[98]  = -0.201194093997435; w[98]  = 0.023648541379510;
      r[99]  = -0.162676817875932; s[99]  = -0.394151347077563; w[99]  = 0.025750117404859;
      r[100] = -0.056478752214152; s[100] = -0.570972172608539; w[100] = 0.025915578164965;
      r[101] =  0.035680197247145; s[101] = -0.724417731360170; w[101] = 0.023879061573116;
      r[102] =  0.110027416239891; s[102] = -0.848206583410427; w[102] = 0.019649913644597;
      r[103] =  0.163520678705051; s[103] = -0.937273392400706; w[103] = 0.013524918548175;
      r[104] =  0.193982435778648; s[104] = -0.987992518020485; w[104] = 0.006065774189928;
      r[105] = -0.993996259010243; s[105] =  0.987992518020485; w[105] = 0.000037402932316;
      r[106] = -0.968636696200353; s[106] =  0.937273392400706; w[106] = 0.000447072297244;
      r[107] = -0.924103291705214; s[107] =  0.848206583410427; w[107] = 0.001647575343713;
      r[108] = -0.862208865680085; s[108] =  0.724417731360170; w[108] = 0.003895904128496;
      r[109] = -0.785486086304269; s[109] =  0.570972172608539; w[109] = 0.007225369917846;
      r[110] = -0.697075673538782; s[110] =  0.394151347077563; w[110] = 0.011423933120557;
      r[111] = -0.600597046998717; s[111] =  0.201194093997435; w[111] = 0.016055160539688;
      r[112] = -0.500000000000000; s[112] =  0.000000000000000; w[112] = 0.020518972050826;
      r[113] = -0.399402953001283; s[113] = -0.201194093997435; w[113] = 0.024142740900556;
      r[114] = -0.302924326461218; s[114] = -0.394151347077563; w[114] = 0.026288234977700;
      r[115] = -0.214513913695731; s[115] = -0.570972172608539; w[115] = 0.026457153482920;
      r[116] = -0.137791134319915; s[116] = -0.724417731360170; w[116] = 0.024378078430143;
      r[117] = -0.075896708294786; s[117] = -0.848206583410427; w[117] = 0.020060551144640;
      r[118] = -0.031363303799647; s[118] = -0.937273392400706; w[118] = 0.013807557894147;
      r[119] = -0.006003740989757; s[119] = -0.987992518020485; w[119] = 0.006192534764769;
      r[120] = -0.995204176239272; s[120] =  0.987992518020485; w[120] = 0.000036637297986;
      r[121] = -0.974946807693089; s[121] =  0.937273392400706; w[121] = 0.000437920771488;
      r[122] = -0.939373261167971; s[122] =  0.848206583410427; w[122] = 0.001613849639200;
      r[123] = -0.889931628110459; s[123] =  0.724417731360170; w[123] = 0.003816155355883;
      r[124] = -0.828645018820126; s[124] =  0.570972172608539; w[124] = 0.007077467309460;
      r[125] = -0.758022258950930; s[125] =  0.394151347077563; w[125] = 0.011190086338209;
      r[126] = -0.680954562267710; s[126] =  0.201194093997435; w[126] = 0.015726512989614;
      r[127] = -0.600597046998717; s[127] =  0.000000000000000; w[127] = 0.020098950720122;
      r[128] = -0.520239531729724; s[128] = -0.201194093997435; w[128] = 0.023648541379510;
      r[129] = -0.443171835046505; s[129] = -0.394151347077563; w[129] = 0.025750117404859;
      r[130] = -0.372549075177309; s[130] = -0.570972172608539; w[130] = 0.025915578164965;
      r[131] = -0.311262465886975; s[131] = -0.724417731360170; w[131] = 0.023879061573116;
      r[132] = -0.261820832829464; s[132] = -0.848206583410427; w[132] = 0.019649913644597;
      r[133] = -0.226247286304345; s[133] = -0.937273392400706; w[133] = 0.013524918548175;
      r[134] = -0.205989917758162; s[134] = -0.987992518020485; w[134] = 0.006065774189928;
      r[135] = -0.996362641608860; s[135] =  0.987992518020485; w[135] = 0.000034371743072;
      r[136] = -0.980998584641787; s[136] =  0.937273392400706; w[136] = 0.000410840893588;
      r[137] = -0.954018081518357; s[137] =  0.848206583410427; w[137] = 0.001514053388316;
      r[138] = -0.916519426887625; s[138] =  0.724417731360170; w[138] = 0.003580174265665;
      r[139] = -0.870037034354322; s[139] =  0.570972172608539; w[139] = 0.006639815197344;
      r[140] = -0.816473704876034; s[140] =  0.394151347077563; w[140] = 0.010498120595869;
      r[141] = -0.758022258950930; s[141] =  0.201194093997435; w[141] = 0.014754026459451;
      r[142] = -0.697075673538782; s[142] =  0.000000000000000; w[142] = 0.018856084049129;
      r[143] = -0.636129088126634; s[143] = -0.201194093997435; w[143] = 0.022186177283618;
      r[144] = -0.577677642201529; s[144] = -0.394151347077563; w[144] = 0.024157797330925;
      r[145] = -0.524114312723242; s[145] = -0.570972172608539; w[145] = 0.024313026429341;
      r[146] = -0.477631920189938; s[146] = -0.724417731360170; w[146] = 0.022402442709917;
      r[147] = -0.440133265559207; s[147] = -0.848206583410427; w[147] = 0.018434814254740;
      r[148] = -0.413152762435777; s[148] = -0.937273392400706; w[148] = 0.012688572873940;
      r[149] = -0.397788705468703; s[149] = -0.987992518020485; w[149] = 0.005690682540646;
      r[150] = -0.997424228046943; s[150] =  0.987992518020485; w[150] = 0.000030699031605;
      r[151] = -0.986544269911019; s[151] =  0.937273392400706; w[151] = 0.000366941459749;
      r[152] = -0.967438200134124; s[152] =  0.848206583410427; w[152] = 0.001352272787635;
      r[153] = -0.940883769008922; s[153] =  0.724417731360170; w[153] = 0.003197623195991;
      r[154] = -0.907967561661881; s[154] =  0.570972172608539; w[154] = 0.005930333418609;
      r[155] = -0.870037034354322; s[155] =  0.394151347077563; w[155] = 0.009376368701824;
      r[156] = -0.828645018820126; s[156] =  0.201194093997435; w[156] = 0.013177519790991;
      r[157] = -0.785486086304269; s[157] =  0.000000000000000; w[157] = 0.016841261700383;
      r[158] = -0.742327153788413; s[158] = -0.201194093997435; w[158] = 0.019815525683434;
      r[159] = -0.700935138254217; s[159] = -0.394151347077563; w[159] = 0.021576472924861;
      r[160] = -0.663004610946658; s[160] = -0.570972172608539; w[160] = 0.021715115384404;
      r[161] = -0.630088403599617; s[161] = -0.724417731360170; w[161] = 0.020008682578130;
      r[162] = -0.603533972474414; s[162] = -0.848206583410427; w[162] = 0.016465005695410;
      r[163] = -0.584427902697520; s[163] = -0.937273392400706; w[163] = 0.011332765372579;
      r[164] = -0.573547944561595; s[164] = -0.987992518020485; w[164] = 0.005082618091388;
      r[165] = -0.998345475437717; s[165] =  0.987992518020485; w[165] = 0.000025769562269;
      r[166] = -0.991356829586853; s[166] =  0.937273392400706; w[166] = 0.000308020165519;
      r[167] = -0.979084212945827; s[167] =  0.848206583410427; w[167] = 0.001135132803359;
      r[168] = -0.962027206605662; s[168] =  0.724417731360170; w[168] = 0.002684167732828;
      r[169] = -0.940883769008922; s[169] =  0.570972172608539; w[169] = 0.004978075474027;
      r[170] = -0.916519426887625; s[170] =  0.394151347077563; w[170] = 0.007870766747028;
      r[171] = -0.889931628110459; s[171] =  0.201194093997435; w[171] = 0.011061551425453;
      r[172] = -0.862208865680085; s[172] =  0.000000000000000; w[172] = 0.014136991279319;
      r[173] = -0.834486103249711; s[173] = -0.201194093997435; w[173] = 0.016633665503545;
      r[174] = -0.807898304472545; s[174] = -0.394151347077563; w[174] = 0.018111850228555;
      r[175] = -0.783533962351248; s[175] = -0.570972172608539; w[175] = 0.018228230300094;
      r[176] = -0.762390524754508; s[176] = -0.724417731360170; w[176] = 0.016795806403939;
      r[177] = -0.745333518414343; s[177] = -0.848206583410427; w[177] = 0.013821152243283;
      r[178] = -0.733060901773317; s[178] = -0.937273392400706; w[178] = 0.009513016785380;
      r[179] = -0.726072255922453; s[179] = -0.987992518020485; w[179] = 0.004266481271556;
      r[180] = -0.999088671642846; s[180] =  0.987992518020485; w[180] = 0.000019785289042;
      r[181] = -0.995239256960715; s[181] =  0.937273392400706; w[181] = 0.000236490940043;
      r[182] = -0.988479379340032; s[182] =  0.848206583410427; w[182] = 0.000871529379538;
      r[183] = -0.979084212945827; s[183] =  0.724417731360170; w[183] = 0.002060843481789;
      r[184] = -0.967438200134124; s[184] =  0.570972172608539; w[184] = 0.003822054138805;
      r[185] = -0.954018081518357; s[185] =  0.394151347077563; w[185] = 0.006042997294437;
      r[186] = -0.939373261167971; s[186] =  0.201194093997435; w[186] = 0.008492809847469;
      r[187] = -0.924103291705214; s[187] =  0.000000000000000; w[187] = 0.010854063244177;
      r[188] = -0.908833322242456; s[188] = -0.201194093997435; w[188] = 0.012770953436328;
      r[189] = -0.894188501892071; s[189] = -0.394151347077563; w[189] = 0.013905870348620;
      r[190] = -0.880768383276303; s[190] = -0.570972172608539; w[190] = 0.013995224344240;
      r[191] = -0.869122370464600; s[191] = -0.724417731360170; w[191] = 0.012895441564853;
      r[192] = -0.859727204070395; s[192] = -0.848206583410427; w[192] = 0.010611569151594;
      r[193] = -0.852967326449712; s[193] = -0.937273392400706; w[193] = 0.007303879856139;
      r[194] = -0.849117911767581; s[194] = -0.987992518020485; w[194] = 0.003275708150100;
      r[195] = -0.999623405694808; s[195] =  0.987992518020485; w[195] = 0.000012991999963;
      r[196] = -0.998032686349542; s[196] =  0.937273392400706; w[196] = 0.000155291655212;
      r[197] = -0.995239256960715; s[197] =  0.848206583410427; w[197] = 0.000572289322752;
      r[198] = -0.991356829586853; s[198] =  0.724417731360170; w[198] = 0.001353251821662;
      r[199] = -0.986544269911019; s[199] =  0.570972172608539; w[199] = 0.002509749901695;
      r[200] = -0.980998584641787; s[200] =  0.394151347077563; w[200] = 0.003968131092565;
      r[201] = -0.974946807693089; s[201] =  0.201194093997435; w[201] = 0.005576799256556;
      r[202] = -0.968636696200353; s[202] =  0.000000000000000; w[202] = 0.007127315095696;
      r[203] = -0.962326584707617; s[203] = -0.201194093997435; w[203] = 0.008386040063107;
      r[204] = -0.956274807758919; s[204] = -0.394151347077563; w[204] = 0.009131282674963;
      r[205] = -0.950729122489687; s[205] = -0.570972172608539; w[205] = 0.009189956930633;
      r[206] = -0.945916562813853; s[206] = -0.724417731360170; w[206] = 0.008467785129237;
      r[207] = -0.942034135439991; s[207] = -0.848206583410427; w[207] = 0.006968081473429;
      r[208] = -0.939240706051164; s[208] = -0.937273392400706; w[208] = 0.004796088983887;
      r[209] = -0.937649986705898; s[209] = -0.987992518020485; w[209] = 0.002150992086750;
      r[210] = -0.999927910188256; s[210] =  0.987992518020485; w[210] = 0.000005678109445;
      r[211] = -0.999623405694808; s[211] =  0.937273392400706; w[211] = 0.000067869690329;
      r[212] = -0.999088671642846; s[212] =  0.848206583410427; w[212] = 0.000250117104238;
      r[213] = -0.998345475437717; s[213] =  0.724417731360170; w[213] = 0.000591434111178;
      r[214] = -0.997424228046943; s[214] =  0.570972172608539; w[214] = 0.001096877668021;
      r[215] = -0.996362641608860; s[215] =  0.394151347077563; w[215] = 0.001734258212849;
      r[216] = -0.995204176239272; s[216] =  0.201194093997435; w[216] = 0.002437321168702;
      r[217] = -0.993996259010243; s[217] =  0.000000000000000; w[217] = 0.003114968848543;
      r[218] = -0.992788341781213; s[218] = -0.201194093997435; w[218] = 0.003665090319212;
      r[219] = -0.991629876411625; s[219] = -0.394151347077563; w[219] = 0.003990796070869;
      r[220] = -0.990568289973542; s[220] = -0.570972172608539; w[220] = 0.004016439454972;
      r[221] = -0.989647042582769; s[221] = -0.724417731360170; w[221] = 0.003700816722647;
      r[222] = -0.988903846377640; s[222] = -0.848206583410427; w[222] = 0.003045376334904;
      r[223] = -0.988369112325678; s[223] = -0.937273392400706; w[223] = 0.002096114396383;
      r[224] = -0.988064607832230; s[224] = -0.987992518020485; w[224] = 0.000940083783827;
    
      break;
    }
  }
}

void FEMStandardElementClass::IntegrationPointsTetrahedron(void) {

  /*--- Set the number of integration points, depending on the order of
        polynomials that must be integrated exactly. ---*/
  switch( orderExact ) {
    case  0: nIntegration =   1; break;
    case  1: nIntegration =   1; break;
    case  2: nIntegration =   4; break;
    case  3: nIntegration =   8; break;
    case  4: nIntegration =  14; break;
    case  5: nIntegration =  14; break;
    case  6: nIntegration =  24; break;
    case  7: nIntegration =  35; break;
    case  8: nIntegration =  46; break;
    case  9: nIntegration =  59; break;
    case 10: nIntegration =  81; break;
    case 11: nIntegration = 105; break;
    case 12: nIntegration = 132; break;
    default:
      cout << "FEMStandardElementClass::IntegrationPointsTetrahedron: Polynomial order not supported" << endl;
#ifndef HAVE_MPI
    exit(EXIT_FAILURE);
#else
    MPI_Abort(MPI_COMM_WORLD,1);
    MPI_Finalize();
#endif
  }

  /*--- Allocate the memory for the integration points and their weights. ---*/
  rIntegration.resize(nIntegration);
  sIntegration.resize(nIntegration);
  tIntegration.resize(nIntegration);
  wIntegration.resize(nIntegration);

  /*--- Set the pointers to the data arrays of the variables just allocated, such
        that the names are shorter. This is useful for the code below. ---*/
  su2double *r = rIntegration.data();
  su2double *s = sIntegration.data();
  su2double *t = tIntegration.data();
  su2double *w = wIntegration.data();

  /*--- Set the data for the integration points, depending on the order.
        These integration rules are obtained with the open source program
        Polyquad, developed by Freddie Witherden.           ---*/
  switch( orderExact )
  {
    case  0:
    case  1: {
      r[0] = -0.500000000000000; s[0] = -0.500000000000000; t[0] = -0.500000000000000; w[0] = 1.333333333333333;
    
      break;
    }
    case  2: {
      r[0] =-0.723606797749979; s[0] = -0.723606797749979; t[0] = 0.1708203932499369; w[0] = 0.3333333333333333;
      r[1] =-0.723606797749979; s[1] = 0.1708203932499369; t[1] = -0.723606797749979; w[1] = 0.3333333333333333;
      r[2] = 0.170820393249937; s[2] = -0.723606797749979; t[2] = -0.723606797749979; w[2] = 0.3333333333333333;
      r[3] =-0.723606797749979; s[3] = -0.723606797749979; t[3] = -0.723606797749979; w[3] = 0.3333333333333333;

      break;
    }
    case  3: {
      r[0] =-0.7860119704589262; s[0] =-0.7860119704589262; t[0] = 0.3580359113767786; w[0] = 0.1487451997767295;
      r[1] =-0.7860119704589262; s[1] = 0.3580359113767786; t[1] =-0.7860119704589262; w[1] = 0.1487451997767295;
      r[2] = 0.3580359113767786; s[2] =-0.7860119704589262; t[2] =-0.7860119704589262; w[2] = 0.1487451997767295;
      r[3] =-0.7860119704589262; s[3] =-0.7860119704589262; t[3] =-0.7860119704589262; w[3] = 0.1487451997767295;
      r[4] =-0.3438828566749746; s[4] =-0.3438828566749746; t[4] =-0.9683514299750762; w[4] = 0.1845881335566038;
      r[5] =-0.3438828566749746; s[5] =-0.9683514299750762; t[5] =-0.3438828566749746; w[5] = 0.1845881335566038;
      r[6] =-0.9683514299750762; s[6] =-0.3438828566749746; t[6] =-0.3438828566749746; w[6] = 0.1845881335566038;
      r[7] =-0.3438828566749746; s[7] =-0.3438828566749746; t[7] =-0.3438828566749746; w[7] = 0.1845881335566038;

      break;
    }
    case  4:
    case  5: {
      r[0]   = -0.8145294993782175; s[0]  = -0.8145294993782175; t[0]  =  0.4435884981346525; w[0]  = 0.09799072415514934;
      r[1]   = -0.8145294993782175; s[1]  =  0.4435884981346525; t[1]  = -0.8145294993782175; w[1]  = 0.09799072415514934;
      r[2]   =  0.4435884981346525; s[2]  = -0.8145294993782175; t[2]  = -0.8145294993782175; w[2]  = 0.09799072415514934;
      r[3]   = -0.8145294993782175; s[3]  = -0.8145294993782175; t[3]  = -0.8145294993782175; w[3]  = 0.09799072415514934;
      r[4]   = -0.3782281614733985; s[4]  = -0.3782281614733985; t[4]  = -0.8653155155798045; w[4]  =  0.1502505676240211;
      r[5]   = -0.3782281614733985; s[5]  = -0.8653155155798045; t[5]  = -0.3782281614733985; w[5]  =  0.1502505676240211;
      r[6]   = -0.8653155155798045; s[6]  = -0.3782281614733985; t[6]  = -0.3782281614733985; w[6]  =  0.1502505676240211;
      r[7]   = -0.3782281614733985; s[7]  = -0.3782281614733985; t[7]  = -0.3782281614733985; w[7]  =  0.1502505676240211;
      r[8]   =  -0.908992591748701; s[8]  =-0.09100740825129899; t[8]  =-0.09100740825129899; w[8]  =  0.0567280277027752;
      r[9]   =-0.09100740825129902; s[9]  =  -0.908992591748701; t[9]  =-0.09100740825129902; w[9]  =  0.0567280277027752;
      r[10]  =  -0.908992591748701; s[10] =  -0.908992591748701; t[10] =-0.09100740825129905; w[10] =  0.0567280277027752;
      r[11]  =  -0.908992591748701; s[11] =-0.09100740825129905; t[11] =  -0.908992591748701; w[11] =  0.0567280277027752;
      r[12]  =-0.09100740825129899; s[12] =  -0.908992591748701; t[12] =  -0.908992591748701; w[12] =  0.0567280277027752;
      r[13]  =-0.09100740825129902; s[13] =-0.09100740825129902; t[13] =  -0.908992591748701; w[13] =  0.0567280277027752;

      break;
    }
    case  6: {
      r[0]  =-0.5707942574816954; s[0]  =-0.5707942574816954; t[0]  =  -0.2876172275549138; w[0]  =  0.05323033367755636;
      r[1]  =-0.5707942574816954; s[1]  =-0.2876172275549138; t[1]  =  -0.5707942574816954; w[1]  =  0.05323033367755636;
      r[2]  =-0.2876172275549138; s[2]  =-0.5707942574816954; t[2]  =  -0.5707942574816954; w[2]  =  0.05323033367755636;
      r[3]  =-0.5707942574816954; s[3]  =-0.5707942574816954; t[3]  =  -0.5707942574816954; w[3]  =  0.05323033367755636;
      r[4]  =-0.3553242197154491; s[4]  =-0.3553242197154491; t[4]  =  -0.9340273408536528; w[4]  =  0.07380957539153983;
      r[5]  =-0.3553242197154491; s[5]  =-0.9340273408536528; t[5]  =  -0.3553242197154491; w[5]  =  0.07380957539153983;
      r[6]  =-0.9340273408536528; s[6]  =-0.3553242197154491; t[6]  =  -0.3553242197154491; w[6]  =  0.07380957539153983;
      r[7]  =-0.3553242197154491; s[7]  =-0.3553242197154491; t[7]  =  -0.3553242197154491; w[7]  =  0.07380957539153983;
      r[8]  =-0.9186520829307774; s[8]  =-0.9186520829307774; t[8]  =   0.7559562487923321; w[8]  =  0.01343628140709418;
      r[9]  =-0.9186520829307774; s[9]  = 0.7559562487923321; t[9]  =  -0.9186520829307774; w[9]  =  0.01343628140709418;
      r[10] = 0.7559562487923321; s[10] =-0.9186520829307774; t[10] =  -0.9186520829307774; w[10] =  0.01343628140709418;
      r[11] =-0.9186520829307774; s[11] =-0.9186520829307774; t[11] =  -0.9186520829307774; w[11] =  0.01343628140709418;
      r[12] = 0.2060113295832982; s[12] =-0.8726779962499649; t[12] =  -0.4606553370833684; w[12] =  0.06428571428571428;
      r[13] = 0.2060113295832982; s[13] =-0.8726779962499649; t[13] =  -0.8726779962499649; w[13] =  0.06428571428571428;
      r[14] =-0.8726779962499649; s[14] =-0.8726779962499649; t[14] =   0.2060113295832982; w[14] =  0.06428571428571428;
      r[15] =-0.4606553370833684; s[15] = 0.2060113295832982; t[15] =  -0.8726779962499649; w[15] =  0.06428571428571428;
      r[16] =-0.8726779962499649; s[16] =-0.4606553370833684; t[16] =   0.2060113295832982; w[16] =  0.06428571428571428;
      r[17] =-0.8726779962499649; s[17] = 0.2060113295832982; t[17] =  -0.8726779962499649; w[17] =  0.06428571428571428;
      r[18] =-0.4606553370833684; s[18] =-0.8726779962499649; t[18] =   0.2060113295832982; w[18] =  0.06428571428571428;
      r[19] =-0.8726779962499649; s[19] =-0.4606553370833684; t[19] =  -0.8726779962499649; w[19] =  0.06428571428571428;
      r[20] =-0.8726779962499649; s[20] =-0.8726779962499649; t[20] =  -0.4606553370833684; w[20] =  0.06428571428571428;
      r[21] =-0.8726779962499649; s[21] = 0.2060113295832982; t[21] =  -0.4606553370833684; w[21] =  0.06428571428571428;
      r[22] =-0.4606553370833684; s[22] =-0.8726779962499649; t[22] =  -0.8726779962499649; w[22] =  0.06428571428571428;
      r[23] = 0.2060113295832982; s[23] =-0.4606553370833684; t[23] =  -0.8726779962499649; w[23] =  0.06428571428571428;

      break;
    }
    case  7: {
      r[0]  =               -0.5; s[0]  =               -0.5; t[0]  =               -0.5; w[0]  =  0.1273137192855079;
      r[1]  =-0.3685977004435944; s[1]  =-0.3685977004435944; t[1]  =-0.8942068986692169; w[1]  = 0.05643944161328927;
      r[2]  =-0.3685977004435944; s[2]  =-0.8942068986692169; t[2]  =-0.3685977004435944; w[2]  = 0.05643944161328927;
      r[3]  =-0.8942068986692169; s[3]  =-0.3685977004435944; t[3]  =-0.3685977004435944; w[3]  = 0.05643944161328927;
      r[4]  =-0.3685977004435944; s[4]  =-0.3685977004435944; t[4]  =-0.3685977004435944; w[4]  = 0.05643944161328927;
      r[5]  =-0.1009796451967928; s[5]  =-0.8990203548032072; t[5]  =-0.8990203548032072; w[5]  = 0.04252923711047679;
      r[6]  =-0.8990203548032072; s[6]  =-0.1009796451967928; t[6]  =-0.8990203548032072; w[6]  = 0.04252923711047679;
      r[7]  =-0.1009796451967928; s[7]  =-0.1009796451967928; t[7]  =-0.8990203548032072; w[7]  = 0.04252923711047679;
      r[8]  =-0.1009796451967928; s[8]  =-0.8990203548032072; t[8]  =-0.1009796451967928; w[8]  = 0.04252923711047679;
      r[9]  =-0.8990203548032072; s[9]  =-0.1009796451967928; t[9]  =-0.1009796451967928; w[9]  = 0.04252923711047679;
      r[10] =-0.8990203548032072; s[10] =-0.8990203548032072; t[10] =-0.1009796451967928; w[10] = 0.04252923711047679;
      r[11] = 0.6216604821970968; s[11] =-0.9574690549170335; t[11] =-0.7067223723630298; w[11] = 0.01081436110653779;
      r[12] = 0.6216604821970968; s[12] =-0.9574690549170335; t[12] =-0.9574690549170335; w[12] = 0.01081436110653779;
      r[13] =-0.9574690549170335; s[13] =-0.9574690549170335; t[13] = 0.6216604821970968; w[13] = 0.01081436110653779;
      r[14] =-0.7067223723630298; s[14] = 0.6216604821970968; t[14] =-0.9574690549170335; w[14] = 0.01081436110653779;
      r[15] =-0.9574690549170335; s[15] =-0.7067223723630298; t[15] = 0.6216604821970968; w[15] = 0.01081436110653779;
      r[16] =-0.9574690549170335; s[16] = 0.6216604821970968; t[16] =-0.9574690549170335; w[16] = 0.01081436110653779;
      r[17] =-0.7067223723630298; s[17] =-0.9574690549170335; t[17] = 0.6216604821970968; w[17] = 0.01081436110653779;
      r[18] =-0.9574690549170335; s[18] =-0.7067223723630298; t[18] =-0.9574690549170335; w[18] = 0.01081436110653779;
      r[19] =-0.9574690549170335; s[19] =-0.9574690549170335; t[19] =-0.7067223723630298; w[19] = 0.01081436110653779;
      r[20] =-0.9574690549170335; s[20] = 0.6216604821970968; t[20] =-0.7067223723630298; w[20] = 0.01081436110653779;
      r[21] =-0.7067223723630298; s[21] =-0.9574690549170335; t[21] =-0.9574690549170335; w[21] = 0.01081436110653779;
      r[22] = 0.6216604821970968; s[22] =-0.7067223723630298; t[22] =-0.9574690549170335; w[22] = 0.01081436110653779;
      r[23] = 0.1503432751739998; s[23] =-0.6223323379479977; t[23] =-0.9056785992780041; w[23] = 0.04960950763777947;
      r[24] = 0.1503432751739997; s[24] =-0.6223323379479978; t[24] =-0.6223323379479978; w[24] = 0.04960950763777947;
      r[25] =-0.6223323379479977; s[25] =-0.6223323379479978; t[25] = 0.1503432751739997; w[25] = 0.04960950763777947;
      r[26] =-0.9056785992780041; s[26] = 0.1503432751739997; t[26] =-0.6223323379479978; w[26] = 0.04960950763777947;
      r[27] =-0.6223323379479978; s[27] =-0.9056785992780041; t[27] = 0.1503432751739997; w[27] = 0.04960950763777947;
      r[28] =-0.6223323379479977; s[28] = 0.1503432751739997; t[28] =-0.6223323379479978; w[28] = 0.04960950763777947;
      r[29] =-0.9056785992780041; s[29] =-0.6223323379479978; t[29] = 0.1503432751739997; w[29] = 0.04960950763777947;
      r[30] =-0.6223323379479978; s[30] =-0.9056785992780041; t[30] =-0.6223323379479978; w[30] = 0.04960950763777947;
      r[31] =-0.6223323379479978; s[31] =-0.6223323379479978; t[31] =-0.9056785992780041; w[31] = 0.04960950763777947;
      r[32] =-0.6223323379479978; s[32] = 0.1503432751739997; t[32] =-0.9056785992780041; w[32] = 0.04960950763777947;
      r[33] =-0.9056785992780041; s[33] =-0.6223323379479978; t[33] =-0.6223323379479978; w[33] = 0.04960950763777947;
      r[34] = 0.1503432751739998; s[34] =-0.9056785992780041; t[34] =-0.6223323379479977; w[34] = 0.04960950763777947;

      break;
    }
    case  8: {
      r[0]  =-0.6320721895815339; s[0]  =-0.6320721895815339; t[0]  =-0.1037834312553985; w[0]  =  0.07440058157405982;
      r[1]  =-0.6320721895815339; s[1]  =-0.1037834312553985; t[1]  =-0.6320721895815339; w[1]  =  0.07440058157405982;
      r[2]  =-0.1037834312553985; s[2]  =-0.6320721895815339; t[2]  =-0.6320721895815339; w[2]  =  0.07440058157405982;
      r[3]  =-0.6320721895815339; s[3]  =-0.6320721895815339; t[3]  =-0.6320721895815339; w[3]  =  0.07440058157405982;
      r[4]  = -0.923300199653955; s[4]  = -0.923300199653955; t[4]  = 0.7699005989618649; w[4]  = 0.007857765854252971;
      r[5]  = -0.923300199653955; s[5]  = 0.7699005989618649; t[5]  = -0.923300199653955; w[5]  = 0.007857765854252971;
      r[6]  =  0.769900598961865; s[6]  = -0.923300199653955; t[6]  = -0.923300199653955; w[6]  = 0.007857765854252971;
      r[7]  = -0.923300199653955; s[7]  = -0.923300199653955; t[7]  = -0.923300199653955; w[7]  = 0.007857765854252971;
      r[8]  =-0.8009956129757766; s[8]  =-0.8009956129757766; t[8]  = 0.4029868389273299; w[8]  =  0.03149884873824373;
      r[9]  =-0.8009956129757766; s[9]  = 0.4029868389273299; t[9]  =-0.8009956129757766; w[9]  =  0.03149884873824373;
      r[10] = 0.4029868389273298; s[10] =-0.8009956129757766; t[10] =-0.8009956129757766; w[10] =  0.03149884873824373;
      r[11] =-0.8009956129757766; s[11] =-0.8009956129757766; t[11] =-0.8009956129757766; w[11] =  0.03149884873824373;
      r[12] =-0.3707849034682734; s[12] =-0.3707849034682734; t[12] =-0.8876452895951799; w[12] =  0.05270161842714149;
      r[13] =-0.3707849034682734; s[13] =-0.8876452895951799; t[13] =-0.3707849034682734; w[13] =  0.05270161842714149;
      r[14] =-0.8876452895951799; s[14] =-0.3707849034682734; t[14] =-0.3707849034682734; w[14] =  0.05270161842714149;
      r[15] =-0.3707849034682734; s[15] =-0.3707849034682734; t[15] =-0.3707849034682734; w[15] =  0.05270161842714149;
      r[16] =-0.1260500583716688; s[16] =-0.8739499416283312; t[16] =-0.8739499416283312; w[16] =  0.04730584876233031;
      r[17] =-0.8739499416283312; s[17] =-0.1260500583716688; t[17] =-0.8739499416283312; w[17] =  0.04730584876233031;
      r[18] =-0.1260500583716688; s[18] =-0.1260500583716688; t[18] =-0.8739499416283312; w[18] =  0.04730584876233031;
      r[19] =-0.1260500583716688; s[19] =-0.8739499416283312; t[19] =-0.1260500583716688; w[19] =  0.04730584876233031;
      r[20] =-0.8739499416283312; s[20] =-0.1260500583716688; t[20] =-0.1260500583716688; w[20] =  0.04730584876233031;
      r[21] =-0.8739499416283312; s[21] =-0.8739499416283312; t[21] =-0.1260500583716688; w[21] =  0.04730584876233031;
      r[22] = 0.1583084315814429; s[22] =-0.5906770163268148; t[22] =-0.9769543989278133; w[22] =   0.0223681938462459;
      r[23] = 0.1583084315814429; s[23] =-0.5906770163268148; t[23] =-0.5906770163268148; w[23] =   0.0223681938462459;
      r[24] =-0.5906770163268148; s[24] =-0.5906770163268148; t[24] = 0.1583084315814429; w[24] =   0.0223681938462459;
      r[25] =-0.9769543989278133; s[25] = 0.1583084315814429; t[25] =-0.5906770163268148; w[25] =   0.0223681938462459;
      r[26] =-0.5906770163268148; s[26] =-0.9769543989278133; t[26] = 0.1583084315814429; w[26] =   0.0223681938462459;
      r[27] =-0.5906770163268148; s[27] = 0.1583084315814429; t[27] =-0.5906770163268148; w[27] =   0.0223681938462459;
      r[28] =-0.9769543989278133; s[28] =-0.5906770163268148; t[28] = 0.1583084315814429; w[28] =   0.0223681938462459;
      r[29] =-0.5906770163268148; s[29] =-0.9769543989278133; t[29] =-0.5906770163268148; w[29] =   0.0223681938462459;
      r[30] =-0.5906770163268148; s[30] =-0.5906770163268148; t[30] =-0.9769543989278133; w[30] =   0.0223681938462459;
      r[31] =-0.5906770163268148; s[31] = 0.1583084315814429; t[31] =-0.9769543989278133; w[31] =   0.0223681938462459;
      r[32] =-0.9769543989278133; s[32] =-0.5906770163268148; t[32] =-0.5906770163268148; w[32] =   0.0223681938462459;
      r[33] = 0.1583084315814429; s[33] =-0.9769543989278133; t[33] =-0.5906770163268148; w[33] =   0.0223681938462459;
      r[34] = 0.4419106530449903; s[34] =-0.9563482048314847; t[34] =-0.5292142433820209; w[34] = 0.009603721352467381;
      r[35] = 0.4419106530449903; s[35] =-0.9563482048314847; t[35] =-0.9563482048314847; w[35] = 0.009603721352467381;
      r[36] =-0.9563482048314847; s[36] =-0.9563482048314847; t[36] = 0.4419106530449902; w[36] = 0.009603721352467381;
      r[37] =-0.5292142433820209; s[37] = 0.4419106530449903; t[37] =-0.9563482048314847; w[37] = 0.009603721352467381;
      r[38] =-0.9563482048314847; s[38] =-0.5292142433820209; t[38] = 0.4419106530449903; w[38] = 0.009603721352467381;
      r[39] =-0.9563482048314847; s[39] = 0.4419106530449902; t[39] =-0.9563482048314847; w[39] = 0.009603721352467381;
      r[40] =-0.5292142433820209; s[40] =-0.9563482048314847; t[40] = 0.4419106530449903; w[40] = 0.009603721352467381;
      r[41] =-0.9563482048314849; s[41] =-0.5292142433820208; t[41] =-0.9563482048314847; w[41] = 0.009603721352467381;
      r[42] =-0.9563482048314849; s[42] =-0.9563482048314847; t[42] =-0.5292142433820208; w[42] = 0.009603721352467381;
      r[43] =-0.9563482048314847; s[43] = 0.4419106530449903; t[43] =-0.5292142433820209; w[43] = 0.009603721352467381;
      r[44] =-0.5292142433820208; s[44] =-0.9563482048314847; t[44] =-0.9563482048314847; w[44] = 0.009603721352467381;
      r[45] = 0.4419106530449903; s[45] =-0.5292142433820209; t[45] =-0.9563482048314847; w[45] = 0.009603721352467381;

      break;
    }
    case  9: {
      r[0]  =                -0.5; s[0]  =                -0.5; t[0]  =                -0.5; w[0]  =  0.07327761647615459;
      r[1]  = -0.8198449760801311; s[1]  = -0.8198449760801311; t[1]  =   0.459534928240393; w[1]  = 0.002848036499465571;
      r[2]  = -0.8198449760801311; s[2]  =   0.459534928240393; t[2]  = -0.8198449760801311; w[2]  = 0.002848036499465571;
      r[3]  =   0.459534928240393; s[3]  = -0.8198449760801311; t[3]  = -0.8198449760801311; w[3]  = 0.002848036499465571;
      r[4]  = -0.8198449760801311; s[4]  = -0.8198449760801311; t[4]  = -0.8198449760801311; w[4]  = 0.002848036499465571;
      r[5]  = -0.9160617305090143; s[5]  = -0.9160617305090143; t[5]  =  0.7481851915270431; w[5]  =  0.00957623445188957;
      r[6]  = -0.9160617305090143; s[6]  =  0.7481851915270431; t[6]  = -0.9160617305090143; w[6]  =  0.00957623445188957;
      r[7]  =  0.7481851915270431; s[7]  = -0.9160617305090143; t[7]  = -0.9160617305090143; w[7]  =  0.00957623445188957;
      r[8]  = -0.9160617305090143; s[8]  = -0.9160617305090143; t[8]  = -0.9160617305090143; w[8]  =  0.00957623445188957;
      r[9]  = -0.3560564433877504; s[9]  = -0.3560564433877504; t[9]  = -0.9318306698367487; w[9]  =  0.04059381195379974;
      r[10] = -0.3560564433877504; s[10] = -0.9318306698367487; t[10] = -0.3560564433877504; w[10] =  0.04059381195379974;
      r[11] = -0.9318306698367487; s[11] = -0.3560564433877504; t[11] = -0.3560564433877504; w[11] =  0.04059381195379974;
      r[12] = -0.3560564433877504; s[12] = -0.3560564433877504; t[12] = -0.3560564433877504; w[12] =  0.04059381195379974;
      r[13] =  -0.658427438836176; s[13] =  -0.658427438836176; t[13] =-0.02471768349147202; w[13] =  0.03397956979508031;
      r[14] =  -0.658427438836176; s[14] =-0.02471768349147202; t[14] =  -0.658427438836176; w[14] =  0.03397956979508031;
      r[15] =-0.02471768349147199; s[15] =  -0.658427438836176; t[15] =  -0.658427438836176; w[15] =  0.03397956979508031;
      r[16] =  -0.658427438836176; s[16] =  -0.658427438836176; t[16] =  -0.658427438836176; w[16] =  0.03397956979508031;
      r[17] = -0.2129618822280694; s[17] = -0.7870381177719306; t[17] = -0.7870381177719306; w[17] =  0.04875647105046832;
      r[18] = -0.7870381177719306; s[18] = -0.2129618822280694; t[18] = -0.7870381177719306; w[18] =  0.04875647105046832;
      r[19] = -0.2129618822280694; s[19] = -0.2129618822280694; t[19] = -0.7870381177719306; w[19] =  0.04875647105046832;
      r[20] = -0.2129618822280694; s[20] = -0.7870381177719306; t[20] = -0.2129618822280694; w[20] =  0.04875647105046832;
      r[21] = -0.7870381177719306; s[21] = -0.2129618822280694; t[21] = -0.2129618822280694; w[21] =  0.04875647105046832;
      r[22] = -0.7870381177719306; s[22] = -0.7870381177719306; t[22] = -0.2129618822280694; w[22] =  0.04875647105046832;
      r[23] =                -1.0; s[23] =-0.07906724107645879; t[23] = -0.8418655178470824; w[23] =  0.01022387218776507;
      r[24] =                -1.0; s[24] =-0.07906724107645879; t[24] =-0.07906724107645879; w[24] =  0.01022387218776507;
      r[25] =-0.07906724107645879; s[25] =-0.07906724107645879; t[25] =                -1.0; w[25] =  0.01022387218776507;
      r[26] = -0.8418655178470824; s[26] =                -1.0; t[26] =-0.07906724107645879; w[26] =  0.01022387218776507;
      r[27] =-0.07906724107645879; s[27] = -0.8418655178470824; t[27] =                -1.0; w[27] =  0.01022387218776507;
      r[28] =-0.07906724107645879; s[28] =                -1.0; t[28] =-0.07906724107645879; w[28] =  0.01022387218776507;
      r[29] = -0.8418655178470824; s[29] =-0.07906724107645879; t[29] =                -1.0; w[29] =  0.01022387218776507;
      r[30] =-0.07906724107645879; s[30] = -0.8418655178470824; t[30] =-0.07906724107645879; w[30] =  0.01022387218776507;
      r[31] =-0.07906724107645879; s[31] =-0.07906724107645879; t[31] = -0.8418655178470824; w[31] =  0.01022387218776507;
      r[32] =-0.07906724107645879; s[32] =                -1.0; t[32] = -0.8418655178470824; w[32] =  0.01022387218776507;
      r[33] = -0.8418655178470824; s[33] =-0.07906724107645879; t[33] =-0.07906724107645879; w[33] =  0.01022387218776507;
      r[34] =                -1.0; s[34] = -0.8418655178470824; t[34] =-0.07906724107645879; w[34] =  0.01022387218776507;
      r[35] =  0.1968616506697061; s[35] = -0.6334882372417514; t[35] = -0.9298851761862033; w[35] =  0.02791006448800061;
      r[36] =  0.1968616506697061; s[36] = -0.6334882372417514; t[36] = -0.6334882372417514; w[36] =  0.02791006448800061;
      r[37] = -0.6334882372417514; s[37] = -0.6334882372417514; t[37] =  0.1968616506697061; w[37] =  0.02791006448800061;
      r[38] = -0.9298851761862033; s[38] =  0.1968616506697061; t[38] = -0.6334882372417514; w[38] =  0.02791006448800061;
      r[39] = -0.6334882372417514; s[39] = -0.9298851761862033; t[39] =  0.1968616506697061; w[39] =  0.02791006448800061;
      r[40] = -0.6334882372417514; s[40] =  0.1968616506697061; t[40] = -0.6334882372417514; w[40] =  0.02791006448800061;
      r[41] = -0.9298851761862033; s[41] = -0.6334882372417514; t[41] =  0.1968616506697061; w[41] =  0.02791006448800061;
      r[42] = -0.6334882372417514; s[42] = -0.9298851761862033; t[42] = -0.6334882372417514; w[42] =  0.02791006448800061;
      r[43] = -0.6334882372417514; s[43] = -0.6334882372417514; t[43] = -0.9298851761862033; w[43] =  0.02791006448800061;
      r[44] = -0.6334882372417514; s[44] =  0.1968616506697061; t[44] = -0.9298851761862033; w[44] =  0.02791006448800061;
      r[45] = -0.9298851761862033; s[45] = -0.6334882372417514; t[45] = -0.6334882372417514; w[45] =  0.02791006448800061;
      r[46] =  0.1968616506697061; s[46] = -0.9298851761862033; t[46] = -0.6334882372417514; w[46] =  0.02791006448800061;
      r[47] =  0.4372341714004844; s[47] = -0.9329857011557074; t[47] = -0.5712627690890695; w[47] =  0.01349325330368666;
      r[48] =  0.4372341714004844; s[48] = -0.9329857011557074; t[48] = -0.9329857011557074; w[48] =  0.01349325330368666;
      r[49] = -0.9329857011557074; s[49] = -0.9329857011557074; t[49] =  0.4372341714004844; w[49] =  0.01349325330368666;
      r[50] = -0.5712627690890695; s[50] =  0.4372341714004844; t[50] = -0.9329857011557074; w[50] =  0.01349325330368666;
      r[51] = -0.9329857011557074; s[51] = -0.5712627690890695; t[51] =  0.4372341714004844; w[51] =  0.01349325330368666;
      r[52] = -0.9329857011557074; s[52] =  0.4372341714004844; t[52] = -0.9329857011557074; w[52] =  0.01349325330368666;
      r[53] = -0.5712627690890695; s[53] = -0.9329857011557074; t[53] =  0.4372341714004844; w[53] =  0.01349325330368666;
      r[54] = -0.9329857011557074; s[54] = -0.5712627690890695; t[54] = -0.9329857011557074; w[54] =  0.01349325330368666;
      r[55] = -0.9329857011557074; s[55] = -0.9329857011557074; t[55] = -0.5712627690890695; w[55] =  0.01349325330368666;
      r[56] = -0.9329857011557074; s[56] =  0.4372341714004844; t[56] = -0.5712627690890695; w[56] =  0.01349325330368666;
      r[57] = -0.5712627690890695; s[57] = -0.9329857011557074; t[57] = -0.9329857011557074; w[57] =  0.01349325330368666;
      r[58] =  0.4372341714004844; s[58] = -0.5712627690890695; t[58] = -0.9329857011557074; w[58] =  0.01349325330368666;

      break;
    }
    case 10: {
      r[0]  =                -0.5; s[0]  =                -0.5; t[0]  =                -0.5; w[0]  =  0.06111729903668942;
      r[1]  = -0.3902379693282467; s[1]  = -0.3902379693282467; t[1]  =   -0.82928609201526; w[1]  =  0.03648653859255376;
      r[2]  = -0.3902379693282467; s[2]  =   -0.82928609201526; t[2]  = -0.3902379693282467; w[2]  =  0.03648653859255376;
      r[3]  =   -0.82928609201526; s[3]  = -0.3902379693282467; t[3]  = -0.3902379693282467; w[3]  =  0.03648653859255376;
      r[4]  = -0.3902379693282467; s[4]  = -0.3902379693282467; t[4]  = -0.3902379693282467; w[4]  =  0.03648653859255376;
      r[5]  = -0.8646773501288808; s[5]  = -0.8646773501288808; t[5]  =  0.5940320503866422; w[5]  =  0.01772043697176758;
      r[6]  = -0.8646773501288808; s[6]  =  0.5940320503866422; t[6]  = -0.8646773501288808; w[6]  =  0.01772043697176758;
      r[7]  =  0.5940320503866422; s[7]  = -0.8646773501288808; t[7]  = -0.8646773501288808; w[7]  =  0.01772043697176758;
      r[8]  = -0.8646773501288808; s[8]  = -0.8646773501288808; t[8]  = -0.8646773501288808; w[8]  =  0.01772043697176758;
      r[9]  = -0.9473346828260848; s[9]  =-0.05266531717391518; t[9]  =-0.05266531717391518; w[9]  =  0.00811057586155056;
      r[10] =-0.05266531717391519; s[10] = -0.9473346828260848; t[10] =-0.05266531717391519; w[10] =  0.00811057586155056;
      r[11] = -0.9473346828260848; s[11] = -0.9473346828260848; t[11] = -0.0526653171739152; w[11] =  0.00811057586155056;
      r[12] = -0.9473346828260848; s[12] = -0.0526653171739152; t[12] = -0.9473346828260848; w[12] =  0.00811057586155056;
      r[13] =-0.05266531717391515; s[13] = -0.9473346828260848; t[13] = -0.9473346828260848; w[13] =  0.00811057586155056;
      r[14] =-0.05266531717391519; s[14] =-0.05266531717391519; t[14] = -0.9473346828260848; w[14] =  0.00811057586155056;
      r[15] = -0.2038789111878272; s[15] = -0.7961210888121728; t[15] = -0.7961210888121728; w[15] =  0.03066711559292157;
      r[16] = -0.7961210888121728; s[16] = -0.2038789111878272; t[16] = -0.7961210888121728; w[16] =  0.03066711559292157;
      r[17] = -0.2038789111878272; s[17] = -0.2038789111878272; t[17] = -0.7961210888121728; w[17] =  0.03066711559292157;
      r[18] = -0.2038789111878272; s[18] = -0.7961210888121728; t[18] = -0.2038789111878272; w[18] =  0.03066711559292157;
      r[19] = -0.7961210888121728; s[19] = -0.2038789111878272; t[19] = -0.2038789111878272; w[19] =  0.03066711559292157;
      r[20] = -0.7961210888121728; s[20] = -0.7961210888121728; t[20] = -0.2038789111878272; w[20] =  0.03066711559292157;
      r[21] = 0.06411454763598279; s[21] = -0.6048469883044114; t[21] =   -0.85442057102716; w[21] =  0.03277830640502917;
      r[22] = 0.06411454763598282; s[22] = -0.6048469883044114; t[22] = -0.6048469883044114; w[22] =  0.03277830640502917;
      r[23] = -0.6048469883044114; s[23] = -0.6048469883044114; t[23] = 0.06411454763598282; w[23] =  0.03277830640502917;
      r[24] =   -0.85442057102716; s[24] = 0.06411454763598279; t[24] = -0.6048469883044114; w[24] =  0.03277830640502917;
      r[25] = -0.6048469883044114; s[25] =   -0.85442057102716; t[25] = 0.06411454763598279; w[25] =  0.03277830640502917;
      r[26] = -0.6048469883044114; s[26] = 0.06411454763598282; t[26] = -0.6048469883044114; w[26] =  0.03277830640502917;
      r[27] =   -0.85442057102716; s[27] = -0.6048469883044114; t[27] = 0.06411454763598279; w[27] =  0.03277830640502917;
      r[28] = -0.6048469883044114; s[28] =   -0.85442057102716; t[28] = -0.6048469883044114; w[28] =  0.03277830640502917;
      r[29] = -0.6048469883044114; s[29] = -0.6048469883044114; t[29] =   -0.85442057102716; w[29] =  0.03277830640502917;
      r[30] = -0.6048469883044114; s[30] = 0.06411454763598279; t[30] =   -0.85442057102716; w[30] =  0.03277830640502917;
      r[31] =   -0.85442057102716; s[31] = -0.6048469883044114; t[31] = -0.6048469883044114; w[31] =  0.03277830640502917;
      r[32] = 0.06411454763598279; s[32] =   -0.85442057102716; t[32] = -0.6048469883044114; w[32] =  0.03277830640502917;
      r[33] =  0.8117588817115982; s[33] =                -1.0; t[33] = -0.8117588817115982; w[33] = 0.001276432944056112;
      r[34] =  0.8117588817115982; s[34] =                -1.0; t[34] =                -1.0; w[34] = 0.001276432944056112;
      r[35] =                -1.0; s[35] =                -1.0; t[35] =  0.8117588817115982; w[35] = 0.001276432944056112;
      r[36] = -0.8117588817115982; s[36] =  0.8117588817115982; t[36] =                -1.0; w[36] = 0.001276432944056112;
      r[37] =                -1.0; s[37] = -0.8117588817115982; t[37] =  0.8117588817115982; w[37] = 0.001276432944056112;
      r[38] =                -1.0; s[38] =  0.8117588817115982; t[38] =                -1.0; w[38] = 0.001276432944056112;
      r[39] = -0.8117588817115982; s[39] =                -1.0; t[39] =  0.8117588817115982; w[39] = 0.001276432944056112;
      r[40] =                -1.0; s[40] = -0.8117588817115982; t[40] =                -1.0; w[40] = 0.001276432944056112;
      r[41] =                -1.0; s[41] =                -1.0; t[41] = -0.8117588817115982; w[41] = 0.001276432944056112;
      r[42] =                -1.0; s[42] =  0.8117588817115982; t[42] = -0.8117588817115982; w[42] = 0.001276432944056112;
      r[43] = -0.8117588817115982; s[43] =                -1.0; t[43] =                -1.0; w[43] = 0.001276432944056112;
      r[44] =  0.8117588817115982; s[44] = -0.8117588817115982; t[44] =                -1.0; w[44] = 0.001276432944056112;
      r[45] =  0.3854615068170896; s[45] = -0.6928701162511337; t[45] = -0.9997212743148225; w[45] = 0.006560363793425748;
      r[46] =  0.3854615068170896; s[46] = -0.6928701162511336; t[46] = -0.6928701162511336; w[46] = 0.006560363793425748;
      r[47] = -0.6928701162511337; s[47] = -0.6928701162511335; t[47] =  0.3854615068170896; w[47] = 0.006560363793425748;
      r[48] = -0.9997212743148225; s[48] =  0.3854615068170896; t[48] = -0.6928701162511335; w[48] = 0.006560363793425748;
      r[49] = -0.6928701162511335; s[49] = -0.9997212743148225; t[49] =  0.3854615068170896; w[49] = 0.006560363793425748;
      r[50] = -0.6928701162511337; s[50] =  0.3854615068170896; t[50] = -0.6928701162511335; w[50] = 0.006560363793425748;
      r[51] = -0.9997212743148225; s[51] = -0.6928701162511335; t[51] =  0.3854615068170896; w[51] = 0.006560363793425748;
      r[52] = -0.6928701162511337; s[52] = -0.9997212743148225; t[52] = -0.6928701162511337; w[52] = 0.006560363793425748;
      r[53] = -0.6928701162511337; s[53] = -0.6928701162511337; t[53] = -0.9997212743148225; w[53] = 0.006560363793425748;
      r[54] = -0.6928701162511335; s[54] =  0.3854615068170896; t[54] = -0.9997212743148225; w[54] = 0.006560363793425748;
      r[55] = -0.9997212743148225; s[55] = -0.6928701162511336; t[55] = -0.6928701162511336; w[55] = 0.006560363793425748;
      r[56] =  0.3854615068170896; s[56] = -0.9997212743148225; t[56] = -0.6928701162511337; w[56] = 0.006560363793425748;
      r[57] = -0.6460142846738401; s[57] = -0.1823369737573523; t[57] = -0.9893117678114551; w[57] =  0.01289761909705229;
      r[58] = -0.6460142846738401; s[58] = -0.1823369737573524; t[58] = -0.1823369737573524; w[58] =  0.01289761909705229;
      r[59] = -0.1823369737573524; s[59] = -0.1823369737573524; t[59] = -0.6460142846738401; w[59] =  0.01289761909705229;
      r[60] = -0.9893117678114551; s[60] = -0.6460142846738401; t[60] = -0.1823369737573524; w[60] =  0.01289761909705229;
      r[61] = -0.1823369737573524; s[61] = -0.9893117678114551; t[61] = -0.6460142846738401; w[61] =  0.01289761909705229;
      r[62] = -0.1823369737573524; s[62] = -0.6460142846738401; t[62] = -0.1823369737573524; w[62] =  0.01289761909705229;
      r[63] = -0.9893117678114551; s[63] = -0.1823369737573524; t[63] = -0.6460142846738401; w[63] =  0.01289761909705229;
      r[64] = -0.1823369737573523; s[64] = -0.9893117678114551; t[64] = -0.1823369737573523; w[64] =  0.01289761909705229;
      r[65] = -0.1823369737573523; s[65] = -0.1823369737573523; t[65] = -0.9893117678114551; w[65] =  0.01289761909705229;
      r[66] = -0.1823369737573524; s[66] = -0.6460142846738401; t[66] = -0.9893117678114551; w[66] =  0.01289761909705229;
      r[67] = -0.9893117678114551; s[67] = -0.1823369737573524; t[67] = -0.1823369737573524; w[67] =  0.01289761909705229;
      r[68] = -0.6460142846738401; s[68] = -0.9893117678114551; t[68] = -0.1823369737573523; w[68] =  0.01289761909705229;
      r[69] =  0.3166885452555603; s[69] = -0.9215201459821879; t[69] = -0.4736482532911847; w[69] =  0.01504744303648042;
      r[70] =  0.3166885452555603; s[70] = -0.9215201459821878; t[70] = -0.9215201459821878; w[70] =  0.01504744303648042;
      r[71] = -0.9215201459821878; s[71] = -0.9215201459821878; t[71] =  0.3166885452555603; w[71] =  0.01504744303648042;
      r[72] = -0.4736482532911847; s[72] =  0.3166885452555603; t[72] = -0.9215201459821878; w[72] =  0.01504744303648042;
      r[73] = -0.9215201459821878; s[73] = -0.4736482532911847; t[73] =  0.3166885452555603; w[73] =  0.01504744303648042;
      r[74] = -0.9215201459821878; s[74] =  0.3166885452555603; t[74] = -0.9215201459821878; w[74] =  0.01504744303648042;
      r[75] = -0.4736482532911847; s[75] = -0.9215201459821878; t[75] =  0.3166885452555603; w[75] =  0.01504744303648042;
      r[76] = -0.9215201459821878; s[76] = -0.4736482532911847; t[76] = -0.9215201459821878; w[76] =  0.01504744303648042;
      r[77] = -0.9215201459821878; s[77] = -0.9215201459821878; t[77] = -0.4736482532911847; w[77] =  0.01504744303648042;
      r[78] = -0.9215201459821878; s[78] =  0.3166885452555603; t[78] = -0.4736482532911847; w[78] =  0.01504744303648042;
      r[79] = -0.4736482532911847; s[79] = -0.9215201459821878; t[79] = -0.9215201459821878; w[79] =  0.01504744303648042;
      r[80] =  0.3166885452555603; s[80] = -0.4736482532911847; t[80] = -0.9215201459821879; w[80] =  0.01504744303648042;

      break;
    }
    case 11: {
      r[0]   =                -0.5; s[0]   =                -0.5; t[0]   =                -0.5; w[0]   =  0.03784683898277254;
      r[1]   =  -0.397380485851065; s[1]   =  -0.397380485851065; t[1]   =  -0.807858542446805; w[1]   =  0.03166744105497066;
      r[2]   =  -0.397380485851065; s[2]   =  -0.807858542446805; t[2]   =  -0.397380485851065; w[2]   =  0.03166744105497066;
      r[3]   =  -0.807858542446805; s[3]   =  -0.397380485851065; t[3]   =  -0.397380485851065; w[3]   =  0.03166744105497066;
      r[4]   =  -0.397380485851065; s[4]   =  -0.397380485851065; t[4]   =  -0.397380485851065; w[4]   =  0.03166744105497066;
      r[5]   = -0.3448885483183461; s[5]   = -0.3448885483183461; t[5]   = -0.9653343550449618; w[5]   =  0.01522419400512382;
      r[6]   = -0.3448885483183461; s[6]   = -0.9653343550449618; t[6]   = -0.3448885483183461; w[6]   =  0.01522419400512382;
      r[7]   = -0.9653343550449618; s[7]   = -0.3448885483183461; t[7]   = -0.3448885483183461; w[7]   =  0.01522419400512382;
      r[8]   = -0.3448885483183461; s[8]   = -0.3448885483183461; t[8]   = -0.3448885483183461; w[8]   =  0.01522419400512382;
      r[9]   = -0.9428094459201933; s[9]   = -0.0571905540798067; t[9]   = -0.0571905540798067; w[9]   = 0.009963593188767755;
      r[10]  =-0.05719055407980668; s[10]  = -0.9428094459201934; t[10]  =-0.05719055407980668; w[10]  = 0.009963593188767755;
      r[11]  = -0.9428094459201933; s[11]  = -0.9428094459201933; t[11]  =-0.05719055407980667; w[11]  = 0.009963593188767755;
      r[12]  = -0.9428094459201933; s[12]  =-0.05719055407980667; t[12]  = -0.9428094459201933; w[12]  = 0.009963593188767755;
      r[13]  =-0.05719055407980672; s[13]  = -0.9428094459201933; t[13]  = -0.9428094459201933; w[13]  = 0.009963593188767755;
      r[14]  =-0.05719055407980668; s[14]  =-0.05719055407980668; t[14]  = -0.9428094459201934; w[14]  = 0.009963593188767755;
      r[15]  = -0.2544420245972532; s[15]  = -0.7455579754027468; t[15]  = -0.7455579754027468; w[15]  =  0.02286217079043096;
      r[16]  = -0.7455579754027468; s[16]  = -0.2544420245972532; t[16]  = -0.7455579754027468; w[16]  =  0.02286217079043096;
      r[17]  = -0.2544420245972532; s[17]  = -0.2544420245972532; t[17]  = -0.7455579754027468; w[17]  =  0.02286217079043096;
      r[18]  = -0.2544420245972532; s[18]  = -0.7455579754027468; t[18]  = -0.2544420245972532; w[18]  =  0.02286217079043096;
      r[19]  = -0.7455579754027468; s[19]  = -0.2544420245972532; t[19]  = -0.2544420245972532; w[19]  =  0.02286217079043096;
      r[20]  = -0.7455579754027468; s[20]  = -0.7455579754027468; t[20]  = -0.2544420245972532; w[20]  =  0.02286217079043096;
      r[21]  =  0.5227647114364038; s[21]  = -0.7755897474786564; t[21]  =  -0.971585216479091; w[21]  = 0.004639621843157262;
      r[22]  =  0.5227647114364038; s[22]  = -0.7755897474786564; t[22]  = -0.7755897474786564; w[22]  = 0.004639621843157262;
      r[23]  = -0.7755897474786564; s[23]  = -0.7755897474786564; t[23]  =  0.5227647114364038; w[23]  = 0.004639621843157262;
      r[24]  =  -0.971585216479091; s[24]  =  0.5227647114364038; t[24]  = -0.7755897474786564; w[24]  = 0.004639621843157262;
      r[25]  = -0.7755897474786564; s[25]  =  -0.971585216479091; t[25]  =  0.5227647114364038; w[25]  = 0.004639621843157262;
      r[26]  = -0.7755897474786564; s[26]  =  0.5227647114364038; t[26]  = -0.7755897474786564; w[26]  = 0.004639621843157262;
      r[27]  =  -0.971585216479091; s[27]  = -0.7755897474786564; t[27]  =  0.5227647114364038; w[27]  = 0.004639621843157262;
      r[28]  = -0.7755897474786564; s[28]  =  -0.971585216479091; t[28]  = -0.7755897474786564; w[28]  = 0.004639621843157262;
      r[29]  = -0.7755897474786564; s[29]  = -0.7755897474786564; t[29]  =  -0.971585216479091; w[29]  = 0.004639621843157262;
      r[30]  = -0.7755897474786564; s[30]  =  0.5227647114364038; t[30]  =  -0.971585216479091; w[30]  = 0.004639621843157262;
      r[31]  =  -0.971585216479091; s[31]  = -0.7755897474786564; t[31]  = -0.7755897474786564; w[31]  = 0.004639621843157262;
      r[32]  =  0.5227647114364038; s[32]  =  -0.971585216479091; t[32]  = -0.7755897474786564; w[32]  = 0.004639621843157262;
      r[33]  =  0.7831995665077649; s[33]  = -0.9669551040900525; t[33]  = -0.8492893583276598; w[33]  = 0.002112278397796427;
      r[34]  =  0.7831995665077649; s[34]  = -0.9669551040900525; t[34]  = -0.9669551040900525; w[34]  = 0.002112278397796427;
      r[35]  = -0.9669551040900525; s[35]  = -0.9669551040900525; t[35]  =   0.783199566507765; w[35]  = 0.002112278397796427;
      r[36]  = -0.8492893583276598; s[36]  =   0.783199566507765; t[36]  = -0.9669551040900525; w[36]  = 0.002112278397796427;
      r[37]  = -0.9669551040900525; s[37]  = -0.8492893583276598; t[37]  =   0.783199566507765; w[37]  = 0.002112278397796427;
      r[38]  = -0.9669551040900525; s[38]  =   0.783199566507765; t[38]  = -0.9669551040900525; w[38]  = 0.002112278397796427;
      r[39]  = -0.8492893583276598; s[39]  = -0.9669551040900525; t[39]  =   0.783199566507765; w[39]  = 0.002112278397796427;
      r[40]  = -0.9669551040900525; s[40]  = -0.8492893583276598; t[40]  = -0.9669551040900525; w[40]  = 0.002112278397796427;
      r[41]  = -0.9669551040900525; s[41]  = -0.9669551040900525; t[41]  = -0.8492893583276598; w[41]  = 0.002112278397796427;
      r[42]  = -0.9669551040900525; s[42]  =   0.783199566507765; t[42]  = -0.8492893583276598; w[42]  = 0.002112278397796427;
      r[43]  = -0.8492893583276598; s[43]  = -0.9669551040900525; t[43]  = -0.9669551040900525; w[43]  = 0.002112278397796427;
      r[44]  =  0.7831995665077649; s[44]  = -0.8492893583276598; t[44]  = -0.9669551040900525; w[44]  = 0.002112278397796427;
      r[45]  =  0.4152874915976885; s[45]  = -0.7537734583743153; t[45]  = -0.9077405748490578; w[45]  = 0.008720212321831882;
      r[46]  =  0.4152874915976886; s[46]  = -0.7537734583743154; t[46]  = -0.7537734583743154; w[46]  = 0.008720212321831882;
      r[47]  = -0.7537734583743154; s[47]  = -0.7537734583743154; t[47]  =  0.4152874915976886; w[47]  = 0.008720212321831882;
      r[48]  = -0.9077405748490577; s[48]  =  0.4152874915976885; t[48]  = -0.7537734583743154; w[48]  = 0.008720212321831882;
      r[49]  = -0.7537734583743154; s[49]  = -0.9077405748490577; t[49]  =  0.4152874915976885; w[49]  = 0.008720212321831882;
      r[50]  = -0.7537734583743154; s[50]  =  0.4152874915976886; t[50]  = -0.7537734583743154; w[50]  = 0.008720212321831882;
      r[51]  = -0.9077405748490577; s[51]  = -0.7537734583743154; t[51]  =  0.4152874915976885; w[51]  = 0.008720212321831882;
      r[52]  = -0.7537734583743154; s[52]  = -0.9077405748490577; t[52]  = -0.7537734583743154; w[52]  = 0.008720212321831882;
      r[53]  = -0.7537734583743154; s[53]  = -0.7537734583743154; t[53]  = -0.9077405748490577; w[53]  = 0.008720212321831882;
      r[54]  = -0.7537734583743154; s[54]  =  0.4152874915976885; t[54]  = -0.9077405748490577; w[54]  = 0.008720212321831882;
      r[55]  = -0.9077405748490578; s[55]  = -0.7537734583743154; t[55]  = -0.7537734583743154; w[55]  = 0.008720212321831882;
      r[56]  =  0.4152874915976885; s[56]  = -0.9077405748490578; t[56]  = -0.7537734583743153; w[56]  = 0.008720212321831882;
      r[57]  = 0.02508211404942551; s[57]  = -0.7430319617668262; t[57]  = -0.5390181905157732; w[57]  =  0.01923782271792721;
      r[58]  = 0.02508211404942554; s[58]  = -0.7430319617668262; t[58]  = -0.7430319617668262; w[58]  =  0.01923782271792721;
      r[59]  = -0.7430319617668262; s[59]  = -0.7430319617668262; t[59]  = 0.02508211404942554; w[59]  =  0.01923782271792721;
      r[60]  = -0.5390181905157732; s[60]  = 0.02508211404942554; t[60]  = -0.7430319617668262; w[60]  =  0.01923782271792721;
      r[61]  = -0.7430319617668262; s[61]  = -0.5390181905157732; t[61]  = 0.02508211404942551; w[61]  =  0.01923782271792721;
      r[62]  = -0.7430319617668262; s[62]  = 0.02508211404942554; t[62]  = -0.7430319617668262; w[62]  =  0.01923782271792721;
      r[63]  = -0.5390181905157732; s[63]  = -0.7430319617668262; t[63]  = 0.02508211404942554; w[63]  =  0.01923782271792721;
      r[64]  = -0.7430319617668262; s[64]  = -0.5390181905157732; t[64]  = -0.7430319617668262; w[64]  =  0.01923782271792721;
      r[65]  = -0.7430319617668262; s[65]  = -0.7430319617668262; t[65]  = -0.5390181905157732; w[65]  =  0.01923782271792721;
      r[66]  = -0.7430319617668262; s[66]  = 0.02508211404942551; t[66]  = -0.5390181905157732; w[66]  =  0.01923782271792721;
      r[67]  = -0.5390181905157732; s[67]  = -0.7430319617668262; t[67]  = -0.7430319617668262; w[67]  =  0.01923782271792721;
      r[68]  = 0.02508211404942551; s[68]  = -0.5390181905157732; t[68]  = -0.7430319617668262; w[68]  =  0.01923782271792721;
      r[69]  =  0.4045657634046023; s[69]  = -0.9481462410143495; t[69]  = -0.5082732813759032; w[69]  = 0.007336665507925923;
      r[70]  =  0.4045657634046023; s[70]  = -0.9481462410143495; t[70]  = -0.9481462410143495; w[70]  = 0.007336665507925923;
      r[71]  = -0.9481462410143495; s[71]  = -0.9481462410143495; t[71]  =  0.4045657634046023; w[71]  = 0.007336665507925923;
      r[72]  = -0.5082732813759032; s[72]  =  0.4045657634046023; t[72]  = -0.9481462410143495; w[72]  = 0.007336665507925923;
      r[73]  = -0.9481462410143495; s[73]  = -0.5082732813759032; t[73]  =  0.4045657634046023; w[73]  = 0.007336665507925923;
      r[74]  = -0.9481462410143495; s[74]  =  0.4045657634046023; t[74]  = -0.9481462410143495; w[74]  = 0.007336665507925923;
      r[75]  = -0.5082732813759032; s[75]  = -0.9481462410143495; t[75]  =  0.4045657634046023; w[75]  = 0.007336665507925923;
      r[76]  = -0.9481462410143495; s[76]  = -0.5082732813759032; t[76]  = -0.9481462410143495; w[76]  = 0.007336665507925923;
      r[77]  = -0.9481462410143495; s[77]  = -0.9481462410143495; t[77]  = -0.5082732813759032; w[77]  = 0.007336665507925923;
      r[78]  = -0.9481462410143495; s[78]  =  0.4045657634046023; t[78]  = -0.5082732813759032; w[78]  = 0.007336665507925923;
      r[79]  = -0.5082732813759032; s[79]  = -0.9481462410143495; t[79]  = -0.9481462410143495; w[79]  = 0.007336665507925923;
      r[80]  =  0.4045657634046023; s[80]  = -0.5082732813759032; t[80]  = -0.9481462410143495; w[80]  = 0.007336665507925923;
      r[81]  = 0.06079123224408745; s[81]  = -0.9501658362591827; t[81]  = -0.7183474030028485; w[81]  =  0.01693359003213852;
      r[82]  = -0.9501658362591827; s[82]  = -0.7183474030028485; t[82]  =  0.0607912322440875; w[82]  =  0.01693359003213852;
      r[83]  = -0.9501658362591827; s[83]  =  0.0607912322440875; t[83]  = -0.7183474030028485; w[83]  =  0.01693359003213852;
      r[84]  = -0.7183474030028484; s[84]  = 0.06079123224408749; t[84]  = -0.3922779929820564; w[84]  =  0.01693359003213852;
      r[85]  = -0.7183474030028485; s[85]  =  0.0607912322440875; t[85]  = -0.9501658362591828; w[85]  =  0.01693359003213852;
      r[86]  = -0.3922779929820563; s[86]  = 0.06079123224408749; t[86]  = -0.7183474030028484; w[86]  =  0.01693359003213852;
      r[87]  = -0.7183474030028485; s[87]  = -0.9501658362591827; t[87]  = -0.3922779929820564; w[87]  =  0.01693359003213852;
      r[88]  = -0.7183474030028484; s[88]  = -0.3922779929820564; t[88]  = 0.06079123224408749; w[88]  =  0.01693359003213852;
      r[89]  = -0.3922779929820565; s[89]  = -0.9501658362591827; t[89]  = -0.7183474030028485; w[89]  =  0.01693359003213852;
      r[90]  =  0.0607912322440875; s[90]  = -0.9501658362591827; t[90]  = -0.3922779929820564; w[90]  =  0.01693359003213852;
      r[91]  = -0.3922779929820564; s[91]  =  0.0607912322440875; t[91]  = -0.9501658362591827; w[91]  =  0.01693359003213852;
      r[92]  = -0.7183474030028485; s[92]  = -0.9501658362591828; t[92]  =  0.0607912322440875; w[92]  =  0.01693359003213852;
      r[93]  = -0.3922779929820564; s[93]  = -0.9501658362591827; t[93]  =  0.0607912322440875; w[93]  =  0.01693359003213852;
      r[94]  = 0.06079123224408749; s[94]  = -0.3922779929820564; t[94]  = -0.7183474030028484; w[94]  =  0.01693359003213852;
      r[95]  = -0.7183474030028485; s[95]  = -0.3922779929820564; t[95]  = -0.9501658362591827; w[95]  =  0.01693359003213852;
      r[96]  = -0.9501658362591827; s[96]  = -0.3922779929820565; t[96]  = -0.7183474030028485; w[96]  =  0.01693359003213852;
      r[97]  =  0.0607912322440875; s[97]  = -0.3922779929820564; t[97]  = -0.9501658362591827; w[97]  =  0.01693359003213852;
      r[98]  = -0.9501658362591827; s[98]  = -0.3922779929820564; t[98]  = 0.06079123224408747; w[98]  =  0.01693359003213852;
      r[99]  = -0.9501658362591827; s[99]  = 0.06079123224408747; t[99]  = -0.3922779929820564; w[99]  =  0.01693359003213852;
      r[100] = -0.3922779929820563; s[100] = -0.7183474030028484; t[100] = 0.06079123224408749; w[100] =  0.01693359003213852;
      r[101] = 0.06079123224408749; s[101] = -0.7183474030028484; t[101] = -0.3922779929820564; w[101] =  0.01693359003213852;
      r[102] = 0.06079123224408745; s[102] = -0.7183474030028485; t[102] = -0.9501658362591827; w[102] =  0.01693359003213852;
      r[103] = -0.9501658362591827; s[103] = -0.7183474030028485; t[103] = -0.3922779929820565; w[103] =  0.01693359003213852;
      r[104] = -0.3922779929820565; s[104] = -0.7183474030028485; t[104] = -0.9501658362591827; w[104] =  0.01693359003213852;

      break;
    }
    case 12: {
      r[0]   = -0.4191163350517226; s[0]   = -0.4191163350517226; t[0]   = -0.7426509948448321; w[0]   =   0.03722586152357994;
      r[1]   = -0.4191163350517226; s[1]   = -0.7426509948448321; t[1]   = -0.4191163350517226; w[1]   =   0.03722586152357994;
      r[2]   = -0.7426509948448321; s[2]   = -0.4191163350517226; t[2]   = -0.4191163350517226; w[2]   =   0.03722586152357994;
      r[3]   = -0.4191163350517226; s[3]   = -0.4191163350517226; t[3]   = -0.4191163350517226; w[3]   =   0.03722586152357994;
      r[4]   = -0.8117896607688906; s[4]   = -0.8117896607688906; t[4]   =  0.4353689823066721; w[4]   =   0.01079833060340124;
      r[5]   = -0.8117896607688906; s[5]   =  0.4353689823066721; t[5]   = -0.8117896607688906; w[5]   =   0.01079833060340124;
      r[6]   =   0.435368982306672; s[6]   = -0.8117896607688906; t[6]   = -0.8117896607688906; w[6]   =   0.01079833060340124;
      r[7]   = -0.8117896607688906; s[7]   = -0.8117896607688906; t[7]   = -0.8117896607688906; w[7]   =   0.01079833060340124;
      r[8]   = -0.6487514080793829; s[8]   = -0.6487514080793829; t[8]   =-0.05374577576185108; w[8]   =   0.02817918150655735;
      r[9]   = -0.6487514080793829; s[9]   =-0.05374577576185108; t[9]   = -0.6487514080793829; w[9]   =   0.02817918150655735;
      r[10]  =-0.05374577576185111; s[10]  = -0.6487514080793829; t[10]  = -0.6487514080793829; w[10]  =   0.02817918150655735;
      r[11]  = -0.6487514080793829; s[11]  = -0.6487514080793829; t[11]  = -0.6487514080793829; w[11]  =   0.02817918150655735;
      r[12]  = -0.2083097004570116; s[12]  = -0.7916902995429884; t[12]  = -0.7916902995429884; w[12]  =   0.02834338756489718;
      r[13]  = -0.7916902995429884; s[13]  = -0.2083097004570116; t[13]  = -0.7916902995429884; w[13]  =   0.02834338756489718;
      r[14]  = -0.2083097004570116; s[14]  = -0.2083097004570116; t[14]  = -0.7916902995429884; w[14]  =   0.02834338756489718;
      r[15]  = -0.2083097004570116; s[15]  = -0.7916902995429884; t[15]  = -0.2083097004570116; w[15]  =   0.02834338756489718;
      r[16]  = -0.7916902995429884; s[16]  = -0.2083097004570116; t[16]  = -0.2083097004570116; w[16]  =   0.02834338756489718;
      r[17]  = -0.7916902995429884; s[17]  = -0.7916902995429884; t[17]  = -0.2083097004570116; w[17]  =   0.02834338756489718;
      r[18]  = -0.9454374073009324; s[18]  =-0.05456259269906756; t[18]  =-0.05456259269906756; w[18]  =   0.00762096465940687;
      r[19]  =-0.05456259269906758; s[19]  = -0.9454374073009324; t[19]  =-0.05456259269906758; w[19]  =   0.00762096465940687;
      r[20]  = -0.9454374073009324; s[20]  = -0.9454374073009324; t[20]  = -0.0545625926990676; w[20]  =   0.00762096465940687;
      r[21]  = -0.9454374073009324; s[21]  = -0.0545625926990676; t[21]  = -0.9454374073009324; w[21]  =   0.00762096465940687;
      r[22]  =-0.05456259269906755; s[22]  = -0.9454374073009324; t[22]  = -0.9454374073009324; w[22]  =   0.00762096465940687;
      r[23]  =-0.05456259269906758; s[23]  =-0.05456259269906758; t[23]  = -0.9454374073009324; w[23]  =   0.00762096465940687;
      r[24]  = -0.5426252783987343; s[24]  = -0.2484192605435968; t[24]  = -0.9605362005140723; w[24]  =  0.009941214894086916;
      r[25]  = -0.5426252783987343; s[25]  = -0.2484192605435968; t[25]  = -0.2484192605435968; w[25]  =  0.009941214894086916;
      r[26]  = -0.2484192605435968; s[26]  = -0.2484192605435968; t[26]  = -0.5426252783987343; w[26]  =  0.009941214894086916;
      r[27]  = -0.9605362005140721; s[27]  = -0.5426252783987343; t[27]  = -0.2484192605435968; w[27]  =  0.009941214894086916;
      r[28]  = -0.2484192605435968; s[28]  = -0.9605362005140723; t[28]  = -0.5426252783987343; w[28]  =  0.009941214894086916;
      r[29]  = -0.2484192605435968; s[29]  = -0.5426252783987343; t[29]  = -0.2484192605435968; w[29]  =  0.009941214894086916;
      r[30]  = -0.9605362005140721; s[30]  = -0.2484192605435968; t[30]  = -0.5426252783987343; w[30]  =  0.009941214894086916;
      r[31]  = -0.2484192605435968; s[31]  = -0.9605362005140721; t[31]  = -0.2484192605435968; w[31]  =  0.009941214894086916;
      r[32]  = -0.2484192605435968; s[32]  = -0.2484192605435968; t[32]  = -0.9605362005140721; w[32]  =  0.009941214894086916;
      r[33]  = -0.2484192605435968; s[33]  = -0.5426252783987343; t[33]  = -0.9605362005140723; w[33]  =  0.009941214894086916;
      r[34]  = -0.9605362005140721; s[34]  = -0.2484192605435968; t[34]  = -0.2484192605435968; w[34]  =  0.009941214894086916;
      r[35]  = -0.5426252783987343; s[35]  = -0.9605362005140723; t[35]  = -0.2484192605435968; w[35]  =  0.009941214894086916;
      r[36]  = 0.02514623497982801; s[36]  = -0.5646581404659121; t[36]  = -0.8958299540480037; w[36]  =   0.01791292763274642;
      r[37]  =   0.025146234979828; s[37]  = -0.5646581404659121; t[37]  = -0.5646581404659121; w[37]  =   0.01791292763274642;
      r[38]  = -0.5646581404659121; s[38]  = -0.5646581404659121; t[38]  = 0.02514623497982797; w[38]  =   0.01791292763274642;
      r[39]  = -0.8958299540480037; s[39]  = 0.02514623497982796; t[39]  = -0.5646581404659121; w[39]  =   0.01791292763274642;
      r[40]  = -0.5646581404659121; s[40]  = -0.8958299540480037; t[40]  = 0.02514623497982799; w[40]  =   0.01791292763274642;
      r[41]  = -0.5646581404659121; s[41]  = 0.02514623497982797; t[41]  = -0.5646581404659121; w[41]  =   0.01791292763274642;
      r[42]  = -0.8958299540480037; s[42]  = -0.5646581404659121; t[42]  = 0.02514623497982796; w[42]  =   0.01791292763274642;
      r[43]  = -0.5646581404659121; s[43]  = -0.8958299540480037; t[43]  = -0.5646581404659121; w[43]  =   0.01791292763274642;
      r[44]  = -0.5646581404659121; s[44]  = -0.5646581404659121; t[44]  = -0.8958299540480037; w[44]  =   0.01791292763274642;
      r[45]  = -0.5646581404659121; s[45]  = 0.02514623497982799; t[45]  = -0.8958299540480037; w[45]  =   0.01791292763274642;
      r[46]  = -0.8958299540480037; s[46]  = -0.5646581404659121; t[46]  = -0.5646581404659121; w[46]  =   0.01791292763274642;
      r[47]  = 0.02514623497982801; s[47]  = -0.8958299540480037; t[47]  = -0.5646581404659121; w[47]  =   0.01791292763274642;
      r[48]  =  0.3827353334438319; s[48]  = -0.9651564497171329; t[48]  = -0.4524224340095662; w[48]  =  0.003515981327063368;
      r[49]  =  0.3827353334438319; s[49]  = -0.9651564497171329; t[49]  = -0.9651564497171329; w[49]  =  0.003515981327063368;
      r[50]  = -0.9651564497171329; s[50]  = -0.9651564497171329; t[50]  =   0.382735333443832; w[50]  =  0.003515981327063368;
      r[51]  = -0.4524224340095661; s[51]  =   0.382735333443832; t[51]  = -0.9651564497171329; w[51]  =  0.003515981327063368;
      r[52]  = -0.9651564497171329; s[52]  = -0.4524224340095662; t[52]  =  0.3827353334438319; w[52]  =  0.003515981327063368;
      r[53]  = -0.9651564497171329; s[53]  =   0.382735333443832; t[53]  = -0.9651564497171329; w[53]  =  0.003515981327063368;
      r[54]  = -0.4524224340095661; s[54]  = -0.9651564497171329; t[54]  =   0.382735333443832; w[54]  =  0.003515981327063368;
      r[55]  = -0.9651564497171329; s[55]  = -0.4524224340095662; t[55]  = -0.9651564497171329; w[55]  =  0.003515981327063368;
      r[56]  = -0.9651564497171329; s[56]  = -0.9651564497171329; t[56]  = -0.4524224340095662; w[56]  =  0.003515981327063368;
      r[57]  = -0.9651564497171329; s[57]  =  0.3827353334438319; t[57]  = -0.4524224340095662; w[57]  =  0.003515981327063368;
      r[58]  = -0.4524224340095661; s[58]  = -0.9651564497171329; t[58]  = -0.9651564497171329; w[58]  =  0.003515981327063368;
      r[59]  =  0.3827353334438319; s[59]  = -0.4524224340095662; t[59]  = -0.9651564497171329; w[59]  =  0.003515981327063368;
      r[60]  =  0.2761393513413102; s[60]  =  -0.860260573663048; t[60]  = -0.5556182040152142; w[60]  =   0.01316564240019808;
      r[61]  =  0.2761393513413101; s[61]  =  -0.860260573663048; t[61]  =  -0.860260573663048; w[61]  =   0.01316564240019808;
      r[62]  =  -0.860260573663048; s[62]  =  -0.860260573663048; t[62]  =  0.2761393513413102; w[62]  =   0.01316564240019808;
      r[63]  = -0.5556182040152142; s[63]  =  0.2761393513413102; t[63]  =  -0.860260573663048; w[63]  =   0.01316564240019808;
      r[64]  =  -0.860260573663048; s[64]  = -0.5556182040152142; t[64]  =  0.2761393513413102; w[64]  =   0.01316564240019808;
      r[65]  =  -0.860260573663048; s[65]  =  0.2761393513413102; t[65]  =  -0.860260573663048; w[65]  =   0.01316564240019808;
      r[66]  = -0.5556182040152142; s[66]  =  -0.860260573663048; t[66]  =  0.2761393513413102; w[66]  =   0.01316564240019808;
      r[67]  =  -0.860260573663048; s[67]  = -0.5556182040152142; t[67]  =  -0.860260573663048; w[67]  =   0.01316564240019808;
      r[68]  =  -0.860260573663048; s[68]  =  -0.860260573663048; t[68]  = -0.5556182040152142; w[68]  =   0.01316564240019808;
      r[69]  =  -0.860260573663048; s[69]  =  0.2761393513413102; t[69]  = -0.5556182040152142; w[69]  =   0.01316564240019808;
      r[70]  = -0.5556182040152142; s[70]  =  -0.860260573663048; t[70]  =  -0.860260573663048; w[70]  =   0.01316564240019808;
      r[71]  =  0.2761393513413102; s[71]  = -0.5556182040152142; t[71]  =  -0.860260573663048; w[71]  =   0.01316564240019808;
      r[72]  =  0.6774764765813509; s[72]  =  -0.948908768728992; t[72]  = -0.7796589391233667; w[72]  =   0.00421031530337885;
      r[73]  =  0.6774764765813508; s[73]  =  -0.948908768728992; t[73]  =  -0.948908768728992; w[73]  =   0.00421031530337885;
      r[74]  =  -0.948908768728992; s[74]  =  -0.948908768728992; t[74]  =  0.6774764765813509; w[74]  =   0.00421031530337885;
      r[75]  = -0.7796589391233667; s[75]  =  0.6774764765813509; t[75]  =  -0.948908768728992; w[75]  =   0.00421031530337885;
      r[76]  =  -0.948908768728992; s[76]  = -0.7796589391233667; t[76]  =  0.6774764765813509; w[76]  =   0.00421031530337885;
      r[77]  =  -0.948908768728992; s[77]  =  0.6774764765813509; t[77]  =  -0.948908768728992; w[77]  =   0.00421031530337885;
      r[78]  = -0.7796589391233667; s[78]  =  -0.948908768728992; t[78]  =  0.6774764765813509; w[78]  =   0.00421031530337885;
      r[79]  =  -0.948908768728992; s[79]  = -0.7796589391233668; t[79]  =  -0.948908768728992; w[79]  =   0.00421031530337885;
      r[80]  =  -0.948908768728992; s[80]  =  -0.948908768728992; t[80]  = -0.7796589391233668; w[80]  =   0.00421031530337885;
      r[81]  =  -0.948908768728992; s[81]  =  0.6774764765813509; t[81]  = -0.7796589391233667; w[81]  =   0.00421031530337885;
      r[82]  = -0.7796589391233668; s[82]  =  -0.948908768728992; t[82]  =  -0.948908768728992; w[82]  =   0.00421031530337885;
      r[83]  =  0.6774764765813509; s[83]  = -0.7796589391233667; t[83]  =  -0.948908768728992; w[83]  =   0.00421031530337885;
      r[84]  =  0.4516899137328004; s[84]  = -0.7307456223061921; t[84]  = -0.9901986691204159; w[84]  =  0.004750818519930626;
      r[85]  =  0.4516899137328003; s[85]  = -0.7307456223061922; t[85]  = -0.7307456223061922; w[85]  =  0.004750818519930626;
      r[86]  = -0.7307456223061921; s[86]  = -0.7307456223061923; t[86]  =  0.4516899137328004; w[86]  =  0.004750818519930626;
      r[87]  = -0.9901986691204159; s[87]  =  0.4516899137328004; t[87]  = -0.7307456223061923; w[87]  =  0.004750818519930626;
      r[88]  = -0.7307456223061923; s[88]  = -0.9901986691204159; t[88]  =  0.4516899137328004; w[88]  =  0.004750818519930626;
      r[89]  = -0.7307456223061921; s[89]  =  0.4516899137328004; t[89]  = -0.7307456223061923; w[89]  =  0.004750818519930626;
      r[90]  = -0.9901986691204159; s[90]  = -0.7307456223061923; t[90]  =  0.4516899137328004; w[90]  =  0.004750818519930626;
      r[91]  = -0.7307456223061922; s[91]  = -0.9901986691204159; t[91]  = -0.7307456223061922; w[91]  =  0.004750818519930626;
      r[92]  = -0.7307456223061922; s[92]  = -0.7307456223061922; t[92]  = -0.9901986691204159; w[92]  =  0.004750818519930626;
      r[93]  = -0.7307456223061923; s[93]  =  0.4516899137328004; t[93]  = -0.9901986691204159; w[93]  =  0.004750818519930626;
      r[94]  = -0.9901986691204159; s[94]  = -0.7307456223061922; t[94]  = -0.7307456223061922; w[94]  =  0.004750818519930626;
      r[95]  =  0.4516899137328004; s[95]  = -0.9901986691204159; t[95]  = -0.7307456223061921; w[95]  =  0.004750818519930626;
      r[96]  =  0.8845962806133972; s[96]  = -0.9719492361359849; t[96]  = -0.9406978083414277; w[96]  = 0.0003599363171224237;
      r[97]  =  0.8845962806133973; s[97]  = -0.9719492361359848; t[97]  = -0.9719492361359848; w[97]  = 0.0003599363171224237;
      r[98]  = -0.9719492361359848; s[98]  = -0.9719492361359848; t[98]  =  0.8845962806133972; w[98]  = 0.0003599363171224237;
      r[99]  = -0.9406978083414278; s[99]  =  0.8845962806133972; t[99]  = -0.9719492361359848; w[99]  = 0.0003599363171224237;
      r[100] = -0.9719492361359848; s[100] = -0.9406978083414278; t[100] =  0.8845962806133972; w[100] = 0.0003599363171224237;
      r[101] = -0.9719492361359848; s[101] =  0.8845962806133972; t[101] = -0.9719492361359848; w[101] = 0.0003599363171224237;
      r[102] = -0.9406978083414278; s[102] = -0.9719492361359848; t[102] =  0.8845962806133972; w[102] = 0.0003599363171224237;
      r[103] = -0.9719492361359848; s[103] = -0.9406978083414277; t[103] = -0.9719492361359848; w[103] = 0.0003599363171224237;
      r[104] = -0.9719492361359848; s[104] = -0.9719492361359848; t[104] = -0.9406978083414277; w[104] = 0.0003599363171224237;
      r[105] = -0.9719492361359848; s[105] =  0.8845962806133972; t[105] = -0.9406978083414278; w[105] = 0.0003599363171224237;
      r[106] = -0.9406978083414277; s[106] = -0.9719492361359848; t[106] = -0.9719492361359848; w[106] = 0.0003599363171224237;
      r[107] =  0.8845962806133972; s[107] = -0.9406978083414277; t[107] = -0.9719492361359849; w[107] = 0.0003599363171224237;
      r[108] =  0.1011535049643996; s[108] = -0.3535634288145484; t[108] = -0.9763709473096609; w[108] =  0.006935487029959765;
      r[109] = -0.3535634288145484; s[109] = -0.9763709473096609; t[109] =  0.1011535049643997; w[109] =  0.006935487029959765;
      r[110] = -0.3535634288145484; s[110] =  0.1011535049643997; t[110] = -0.9763709473096609; w[110] =  0.006935487029959765;
      r[111] = -0.9763709473096609; s[111] =  0.1011535049643996; t[111] = -0.7712191288401904; w[111] =  0.006935487029959765;
      r[112] = -0.9763709473096609; s[112] =  0.1011535049643996; t[112] = -0.3535634288145484; w[112] =  0.006935487029959765;
      r[113] = -0.7712191288401904; s[113] =  0.1011535049643997; t[113] = -0.9763709473096609; w[113] =  0.006935487029959765;
      r[114] = -0.9763709473096609; s[114] = -0.3535634288145484; t[114] = -0.7712191288401904; w[114] =  0.006935487029959765;
      r[115] = -0.9763709473096609; s[115] = -0.7712191288401904; t[115] =  0.1011535049643996; w[115] =  0.006935487029959765;
      r[116] = -0.7712191288401904; s[116] = -0.3535634288145484; t[116] = -0.9763709473096609; w[116] =  0.006935487029959765;
      r[117] =  0.1011535049643996; s[117] = -0.3535634288145484; t[117] = -0.7712191288401904; w[117] =  0.006935487029959765;
      r[118] = -0.7712191288401904; s[118] =  0.1011535049643996; t[118] = -0.3535634288145484; w[118] =  0.006935487029959765;
      r[119] = -0.9763709473096609; s[119] = -0.3535634288145484; t[119] =  0.1011535049643996; w[119] =  0.006935487029959765;
      r[120] = -0.7712191288401904; s[120] = -0.3535634288145484; t[120] =  0.1011535049643996; w[120] =  0.006935487029959765;
      r[121] =  0.1011535049643997; s[121] = -0.7712191288401904; t[121] = -0.9763709473096609; w[121] =  0.006935487029959765;
      r[122] = -0.9763709473096609; s[122] = -0.7712191288401904; t[122] = -0.3535634288145484; w[122] =  0.006935487029959765;
      r[123] = -0.3535634288145484; s[123] = -0.7712191288401904; t[123] = -0.9763709473096609; w[123] =  0.006935487029959765;
      r[124] =  0.1011535049643996; s[124] = -0.7712191288401904; t[124] = -0.3535634288145484; w[124] =  0.006935487029959765;
      r[125] = -0.3535634288145484; s[125] = -0.7712191288401904; t[125] =  0.1011535049643996; w[125] =  0.006935487029959765;
      r[126] = -0.3535634288145484; s[126] =  0.1011535049643996; t[126] = -0.7712191288401904; w[126] =  0.006935487029959765;
      r[127] = -0.7712191288401904; s[127] = -0.9763709473096609; t[127] =  0.1011535049643997; w[127] =  0.006935487029959765;
      r[128] =  0.1011535049643997; s[128] = -0.9763709473096609; t[128] = -0.7712191288401904; w[128] =  0.006935487029959765;
      r[129] =  0.1011535049643996; s[129] = -0.9763709473096609; t[129] = -0.3535634288145484; w[129] =  0.006935487029959765;
      r[130] = -0.3535634288145484; s[130] = -0.9763709473096609; t[130] = -0.7712191288401904; w[130] =  0.006935487029959765;
      r[131] = -0.7712191288401904; s[131] = -0.9763709473096609; t[131] = -0.3535634288145484; w[131] =  0.006935487029959765;

      break;
    }
  }
}

void FEMStandardElementClass::IntegrationPointsPyramid(void) {

  /*--- Set the number of integration points, depending on the order of
        polynomials that must be integrated exactly. ---*/
  switch( orderExact ) {
    case  0: nIntegration =   1; break;
    case  1: nIntegration =   1; break;
    case  2: nIntegration =   5; break;
    case  3: nIntegration =   6; break;
    case  4: nIntegration =  10; break;
    case  5: nIntegration =  15; break;
    case  6: nIntegration =  24; break;
    case  7: nIntegration =  31; break;
    case  8: nIntegration =  47; break;
    case  9: nIntegration =  62; break;
    default:
      cout << "FEMStandardElementClass::IntegrationPointsPyramid: Polynomial order not supported" << endl;
#ifndef HAVE_MPI
    exit(EXIT_FAILURE);
#else
    MPI_Abort(MPI_COMM_WORLD,1);
    MPI_Finalize();
#endif
  }

    /*--- Allocate the memory for the integration points and their weights. ---*/
  rIntegration.resize(nIntegration);
  sIntegration.resize(nIntegration);
  tIntegration.resize(nIntegration);
  wIntegration.resize(nIntegration);

  /*--- Set the pointers to the data arrays of the variables just allocated, such
        that the names are shorter. This is useful for the code below. ---*/
  su2double *r = rIntegration.data();
  su2double *s = sIntegration.data();
  su2double *t = tIntegration.data();
  su2double *w = wIntegration.data();

  /*--- Set the data for the integration points, depending on the order.
        These integration rules are obtained with the open source program
        Polyquad, developed by Freddie Witherden.           ---*/
  switch( orderExact ) {
    case  0:
    case  1: {
      r[0] = 0.0; s[0] = 0.0; t[0] = -0.5; w[0] = 2.666666666666667;
    
      break;
    }
    case  2: {
      r[0] =                0.0; s[0] =                0.0; t[0] = 0.5981411585404836; w[0] = 0.2950039595403516;
      r[1] = 0.4742123952580272; s[1] = 0.4742123952580272; t[1] =-0.6365944613162136; w[1] = 0.5929156767815785;
      r[2] = 0.4742123952580272; s[2] =-0.4742123952580272; t[2] =-0.6365944613162136; w[2] = 0.5929156767815785;
      r[3] =-0.4742123952580272; s[3] = 0.4742123952580272; t[3] =-0.6365944613162136; w[3] = 0.5929156767815785;
      r[4] =-0.4742123952580272; s[4] =-0.4742123952580272; t[4] =-0.6365944613162136; w[4] = 0.5929156767815785;

      break;
    }
    case  3: {
      r[0] =                0.0; s[0] =                0.0; t[0] = 0.1421672642720796; w[0] = 0.6739373843594919;
      r[1] =                0.0; s[1] =                0.0; t[1] =-0.9962540442494173; w[1] = 0.3054090848134727;
      r[2] = 0.5622126477146682; s[2] = 0.5622126477146682; t[2] =-0.6666666666666666; w[2] = 0.4218300493734253;
      r[3] = 0.5622126477146682; s[3] =-0.5622126477146682; t[3] =-0.6666666666666666; w[3] = 0.4218300493734253;
      r[4] =-0.5622126477146682; s[4] = 0.5622126477146682; t[4] =-0.6666666666666666; w[4] = 0.4218300493734253;
      r[5] =-0.5622126477146682; s[5] =-0.5622126477146682; t[5] =-0.6666666666666666; w[5] = 0.4218300493734253;

      break;
    }
    case  4: {
      r[0] =                0.0; s[0] =                0.0; t[0] =-0.7497260937825071; w[0] = 0.5516890735721394;
      r[1] =                0.0; s[1] =                0.0; t[1] = 0.3544655777722747; w[1] = 0.3033116884550452;
      r[2] = 0.6505815563982326; s[2] =                0.0; t[2] =-0.3552317008435727; w[2] = 0.2835322343715343;
      r[3] =                0.0; s[3] = 0.6505815563982326; t[3] =-0.3552317008435727; w[3] = 0.2835322343715343;
      r[4] =-0.6505815563982326; s[4] =                0.0; t[4] =-0.3552317008435727; w[4] = 0.2835322343715343;
      r[5] =                0.0; s[5] =-0.6505815563982326; t[5] =-0.3552317008435727; w[5] = 0.2835322343715343;
      r[6] =   0.65796699712169; s[6] =   0.65796699712169; t[6] =-0.9215034322023693; w[6] = 0.1693842417883357;
      r[7] =   0.65796699712169; s[7] =  -0.65796699712169; t[7] =-0.9215034322023693; w[7] = 0.1693842417883357;
      r[8] =  -0.65796699712169; s[8] =   0.65796699712169; t[8] =-0.9215034322023693; w[8] = 0.1693842417883357;
      r[9] =  -0.65796699712169; s[9] =  -0.65796699712169; t[9] =-0.9215034322023693; w[9] = 0.1693842417883357;

      break;
    }
    case  5: {
      r[0]  =                0.0; s[0]  =                0.0; t[0]  =-0.9948776933268603; w[0]  = 0.1628511313872155;
      r[1]  =                0.0; s[1]  =                0.0; t[1]  = 0.4597001791131636; w[1]  = 0.1824149770567491;
      r[2]  =                0.0; s[2]  =                0.0; t[2]  =-0.4137399867225705; w[2]  = 0.4561666547431934;
      r[3]  = 0.7161787206720601; s[3]  =                0.0; t[3]  =              -0.75; w[3]  = 0.1930739923649903;
      r[4]  =                0.0; s[4]  = 0.7161787206720601; t[4]  =              -0.75; w[4]  = 0.1930739923649903;
      r[5]  =-0.7161787206720601; s[5]  =                0.0; t[5]  =              -0.75; w[5]  = 0.1930739923649903;
      r[6]  =                0.0; s[6]  =-0.7161787206720601; t[6]  =              -0.75; w[6]  = 0.1930739923649903;
      r[7]  = 0.6978528932694547; s[7]  = 0.6978528932694547; t[7]  = -0.875569951422984; w[7]  = 0.1105728673872057;
      r[8]  = 0.6978528932694547; s[8]  =-0.6978528932694547; t[8]  = -0.875569951422984; w[8]  = 0.1105728673872057;
      r[9]  =-0.6978528932694547; s[9]  = 0.6978528932694547; t[9]  = -0.875569951422984; w[9]  = 0.1105728673872057;
      r[10] =-0.6978528932694547; s[10] =-0.6978528932694547; t[10] = -0.875569951422984; w[10] = 0.1105728673872057;
      r[11] =  0.429237163817977; s[11] =  0.429237163817977; t[11] =-0.1536304201868109; w[11] = 0.1626616161176814;
      r[12] =  0.429237163817977; s[12] = -0.429237163817977; t[12] =-0.1536304201868109; w[12] = 0.1626616161176814;
      r[13] = -0.429237163817977; s[13] =  0.429237163817977; t[13] =-0.1536304201868109; w[13] = 0.1626616161176814;
      r[14] = -0.429237163817977; s[14] = -0.429237163817977; t[14] =-0.1536304201868109; w[14] = 0.1626616161176814;

      break;
    }
    case  6: {
      r[0]  =                0.0; s[0]  =                0.0; t[0]  =-0.6788935976153283; w[0]  =  0.2831712115065466;
      r[1]  =                0.0; s[1]  =                0.0; t[1]  = 0.6180312275468769; w[1]  = 0.06690481255981294;
      r[2]  =                0.0; s[2]  =                0.0; t[2]  =-0.9276877344294127; w[2]  = 0.09145928219274543;
      r[3]  =                0.0; s[3]  =                0.0; t[3]  =-0.1422484612257312; w[3]  =  0.2609124805007913;
      r[4]  = 0.4339254093766991; s[4]  =                0.0; t[4]  = 0.1321491812466018; w[4]  = 0.08423354530179532;
      r[5]  =                0.0; s[5]  = 0.4339254093766991; t[5]  = 0.1321491812466018; w[5]  = 0.08423354530179532;
      r[6]  =-0.4339254093766991; s[6]  =                0.0; t[6]  = 0.1321491812466018; w[6]  = 0.08423354530179532;
      r[7]  =                0.0; s[7]  =-0.4339254093766991; t[7]  = 0.1321491812466018; w[7]  = 0.08423354530179532;
      r[8]  = 0.8345953511147083; s[8]  =                0.0; t[8]  = -0.805105317949076; w[8]  = 0.09853490116354532;
      r[9]  =                0.0; s[9]  = 0.8345953511147083; t[9]  = -0.805105317949076; w[9]  = 0.09853490116354532;
      r[10] =-0.8345953511147083; s[10] =                0.0; t[10] = -0.805105317949076; w[10] = 0.09853490116354532;
      r[11] =                0.0; s[11] =-0.8345953511147083; t[11] = -0.805105317949076; w[11] = 0.09853490116354532;
      r[12] = 0.5622642698597283; s[12] = 0.5622642698597283; t[12] =-0.9412721511677906; w[12] = 0.09909432092383597;
      r[13] = 0.5622642698597283; s[13] =-0.5622642698597283; t[13] =-0.9412721511677906; w[13] = 0.09909432092383597;
      r[14] =-0.5622642698597283; s[14] = 0.5622642698597283; t[14] =-0.9412721511677906; w[14] = 0.09909432092383597;
      r[15] =-0.5622642698597283; s[15] =-0.5622642698597283; t[15] =-0.9412721511677906; w[15] = 0.09909432092383597;
      r[16] = 0.4980978320663696; s[16] = 0.4980978320663696; t[16] =-0.4701693571867578; w[16] =  0.1970067584185277;
      r[17] = 0.4980978320663696; s[17] =-0.4980978320663696; t[17] =-0.4701693571867578; w[17] =  0.1970067584185277;
      r[18] =-0.4980978320663696; s[18] = 0.4980978320663696; t[18] =-0.4701693571867578; w[18] =  0.1970067584185277;
      r[19] =-0.4980978320663696; s[19] =-0.4980978320663696; t[19] =-0.4701693571867578; w[19] =  0.1970067584185277;
      r[20] = 0.9449464128891149; s[20] = 0.9449464128891149; t[20] =  -0.90429505212374; w[20] = 0.01218519416898826;
      r[21] = 0.9449464128891149; s[21] =-0.9449464128891149; t[21] =  -0.90429505212374; w[21] = 0.01218519416898826;
      r[22] =-0.9449464128891149; s[22] = 0.9449464128891149; t[22] =  -0.90429505212374; w[22] = 0.01218519416898826;
      r[23] =-0.9449464128891149; s[23] =-0.9449464128891149; t[23] =  -0.90429505212374; w[23] = 0.01218519416898826;

      break;
    }
    case  7: {
      r[0]  =                0.0; s[0]  =                0.0; t[0]  =-0.9999090434481293; w[0]  = 0.06666820462099099;
      r[1]  =                0.0; s[1]  =                0.0; t[1]  =  0.677365558291527; w[1]  = 0.04188580085475284;
      r[2]  =                0.0; s[2]  =                0.0; t[2]  = -0.212834236290852; w[2]  =   0.268080483940764;
      r[3]  = 0.8640987597877147; s[3]  =                0.0; t[3]  =-0.8666666666666667; w[3]  = 0.07117802478134111;
      r[4]  =                0.0; s[4]  = 0.8640987597877147; t[4]  =-0.8666666666666667; w[4]  = 0.07117802478134111;
      r[5]  =-0.8640987597877147; s[5]  =                0.0; t[5]  =-0.8666666666666667; w[5]  = 0.07117802478134111;
      r[6]  =                0.0; s[6]  =-0.8640987597877147; t[6]  =-0.8666666666666667; w[6]  = 0.07117802478134111;
      r[7]  = 0.6172133998483676; s[7]  =                0.0; t[7]  =-0.3333333333333334; w[7]  = 0.07656249999999998;
      r[8]  =                0.0; s[8]  = 0.6172133998483676; t[8]  =-0.3333333333333334; w[8]  = 0.07656249999999998;
      r[9]  =-0.6172133998483676; s[9]  =                0.0; t[9]  =-0.3333333333333334; w[9]  = 0.07656249999999998;
      r[10] =                0.0; s[10] =-0.6172133998483676; t[10] =-0.3333333333333334; w[10] = 0.07656249999999998;
      r[11] = 0.3541523732387406; s[11] = 0.3541523732387406; t[11] =-0.7414271760285904; w[11] =  0.1758938286698244;
      r[12] = 0.3541523732387406; s[12] =-0.3541523732387406; t[12] =-0.7414271760285904; w[12] =  0.1758938286698244;
      r[13] =-0.3541523732387406; s[13] = 0.3541523732387406; t[13] =-0.7414271760285904; w[13] =  0.1758938286698244;
      r[14] =-0.3541523732387406; s[14] =-0.3541523732387406; t[14] =-0.7414271760285904; w[14] =  0.1758938286698244;
      r[15] = 0.6143051052061707; s[15] = 0.6143051052061707; t[15] =-0.9999999846056514; w[15] = 0.03547692572372306;
      r[16] = 0.6143051052061707; s[16] =-0.6143051052061707; t[16] =-0.9999999846056514; w[16] = 0.03547692572372306;
      r[17] =-0.6143051052061707; s[17] = 0.6143051052061707; t[17] =-0.9999999846056514; w[17] = 0.03547692572372306;
      r[18] =-0.6143051052061707; s[18] =-0.6143051052061707; t[18] =-0.9999999846056514; w[18] = 0.03547692572372306;
      r[19] = 0.5248327023957946; s[19] = 0.5248327023957946; t[19] =-0.4189138623083963; w[19] = 0.09524944152202758;
      r[20] = 0.5248327023957946; s[20] =-0.5248327023957946; t[20] =-0.4189138623083963; w[20] = 0.09524944152202758;
      r[21] =-0.5248327023957946; s[21] = 0.5248327023957946; t[21] =-0.4189138623083963; w[21] = 0.09524944152202758;
      r[22] =-0.5248327023957946; s[22] =-0.5248327023957946; t[22] =-0.4189138623083963; w[22] = 0.09524944152202758;
      r[23] = 0.2541353514987326; s[23] = 0.2541353514987326; t[23] = 0.2110398539463917; w[23] = 0.07871762583939068;
      r[24] = 0.2541353514987326; s[24] =-0.2541353514987326; t[24] = 0.2110398539463917; w[24] = 0.07871762583939068;
      r[25] =-0.2541353514987326; s[25] = 0.2541353514987326; t[25] = 0.2110398539463917; w[25] = 0.07871762583939068;
      r[26] =-0.2541353514987326; s[26] =-0.2541353514987326; t[26] = 0.2110398539463917; w[26] = 0.07871762583939068;
      r[27] = 0.8027258573724061; s[27] = 0.8027258573724061; t[27] =-0.8397229410036666; w[27] = 0.03942969777623283;
      r[28] = 0.8027258573724061; s[28] =-0.8027258573724061; t[28] =-0.8397229410036666; w[28] = 0.03942969777623283;
      r[29] =-0.8027258573724061; s[29] = 0.8027258573724061; t[29] =-0.8397229410036666; w[29] = 0.03942969777623283;
      r[30] =-0.8027258573724061; s[30] =-0.8027258573724061; t[30] =-0.8397229410036666; w[30] = 0.03942969777623283;

      break;
    }
    case  8: {
      r[0]  =                0.0; s[0]  =                0.0; t[0]  =  -0.550300841839185; w[0]  =  0.2312967200769446;
      r[1]  =                0.0; s[1]  =                0.0; t[1]  =  0.3543960756117052; w[1]  =  0.0895296925217413;
      r[2]  =                0.0; s[2]  =                0.0; t[2]  =  0.7700578580753633; w[2]  = 0.01667753994764634;
      r[3]  = 0.7343313655310747; s[3]  =                0.0; t[3]  = -0.5812723517355866; w[3]  = 0.06765918408093817;
      r[4]  =                0.0; s[4]  = 0.7343313655310747; t[4]  = -0.5812723517355866; w[4]  = 0.06765918408093817;
      r[5]  =-0.7343313655310747; s[5]  =                0.0; t[5]  = -0.5812723517355866; w[5]  = 0.06765918408093817;
      r[6]  =                0.0; s[6]  =-0.7343313655310747; t[6]  = -0.5812723517355866; w[6]  = 0.06765918408093817;
      r[7]  = 0.4747235679824486; s[7]  =                0.0; t[7]  =-0.03583468456802791; w[7]  = 0.04230006651928538;
      r[8]  =                0.0; s[8]  = 0.4747235679824486; t[8]  =-0.03583468456802791; w[8]  = 0.04230006651928538;
      r[9]  =-0.4747235679824486; s[9]  =                0.0; t[9]  =-0.03583468456802791; w[9]  = 0.04230006651928538;
      r[10] =                0.0; s[10] =-0.4747235679824486; t[10] =-0.03583468456802791; w[10] = 0.04230006651928538;
      r[11] = 0.4915331508073653; s[11] =                0.0; t[11] = -0.9139201438073591; w[11] = 0.09881180834785219;
      r[12] =                0.0; s[12] = 0.4915331508073653; t[12] = -0.9139201438073591; w[12] = 0.09881180834785219;
      r[13] =-0.4915331508073653; s[13] =                0.0; t[13] = -0.9139201438073591; w[13] = 0.09881180834785219;
      r[14] =                0.0; s[14] =-0.4915331508073653; t[14] = -0.9139201438073591; w[14] = 0.09881180834785219;
      r[15] = 0.6014354459180917; s[15] = 0.6014354459180917; t[15] = -0.3292642989286104; w[15] = 0.02745707476038042;
      r[16] = 0.6014354459180917; s[16] =-0.6014354459180917; t[16] = -0.3292642989286104; w[16] = 0.02745707476038042;
      r[17] =-0.6014354459180917; s[17] = 0.6014354459180917; t[17] = -0.3292642989286104; w[17] = 0.02745707476038042;
      r[18] =-0.6014354459180917; s[18] =-0.6014354459180917; t[18] = -0.3292642989286104; w[18] = 0.02745707476038042;
      r[19] = 0.7336097872841104; s[19] = 0.7336097872841104; t[19] = -0.9999999921972655; w[19] = 0.01143429356535436;
      r[20] = 0.7336097872841104; s[20] =-0.7336097872841104; t[20] = -0.9999999921972655; w[20] = 0.01143429356535436;
      r[21] =-0.7336097872841104; s[21] = 0.7336097872841104; t[21] = -0.9999999921972655; w[21] = 0.01143429356535436;
      r[22] =-0.7336097872841104; s[22] =-0.7336097872841104; t[22] = -0.9999999921972655; w[22] = 0.01143429356535436;
      r[23] = 0.2811970261033676; s[23] = 0.2811970261033676; t[23] = -0.1471071439241793; w[23] =  0.1049900793442858;
      r[24] = 0.2811970261033676; s[24] =-0.2811970261033676; t[24] = -0.1471071439241793; w[24] =  0.1049900793442858;
      r[25] =-0.2811970261033676; s[25] = 0.2811970261033676; t[25] = -0.1471071439241793; w[25] =  0.1049900793442858;
      r[26] =-0.2811970261033676; s[26] =-0.2811970261033676; t[26] = -0.1471071439241793; w[26] =  0.1049900793442858;
      r[27] = 0.8386356598271838; s[27] = 0.8386356598271838; t[27] =   -0.82637681792063; w[27] = 0.01838174019191864;
      r[28] = 0.8386356598271838; s[28] =-0.8386356598271838; t[28] =   -0.82637681792063; w[28] = 0.01838174019191864;
      r[29] =-0.8386356598271838; s[29] = 0.8386356598271838; t[29] =   -0.82637681792063; w[29] = 0.01838174019191864;
      r[30] =-0.8386356598271838; s[30] =-0.8386356598271838; t[30] =   -0.82637681792063; w[30] = 0.01838174019191864;
      r[31] = 0.2519621813096428; s[31] = 0.2519621813096428; t[31] =  0.3917836190480479; w[31] = 0.02250417643355401;
      r[32] = 0.2519621813096428; s[32] =-0.2519621813096428; t[32] =  0.3917836190480479; w[32] = 0.02250417643355401;
      r[33] =-0.2519621813096428; s[33] = 0.2519621813096428; t[33] =  0.3917836190480479; w[33] = 0.02250417643355401;
      r[34] =-0.2519621813096428; s[34] =-0.2519621813096428; t[34] =  0.3917836190480479; w[34] = 0.02250417643355401;
      r[35] = 0.4943034129249919; s[35] = 0.4943034129249919; t[35] = -0.6479128212741551; w[35] =  0.1262187345108739;
      r[36] = 0.4943034129249919; s[36] =-0.4943034129249919; t[36] = -0.6479128212741551; w[36] =  0.1262187345108739;
      r[37] =-0.4943034129249919; s[37] = 0.4943034129249919; t[37] = -0.6479128212741551; w[37] =  0.1262187345108739;
      r[38] =-0.4943034129249919; s[38] =-0.4943034129249919; t[38] = -0.6479128212741551; w[38] =  0.1262187345108739;
      r[39] = 0.8747496344107996; s[39] = 0.3950260023113457; t[39] = -0.9184400985451371; w[39] = 0.03126676038782074;
      r[40] = 0.3950260023113457; s[40] = 0.8747496344107996; t[40] = -0.9184400985451371; w[40] = 0.03126676038782074;
      r[41] = 0.8747496344107996; s[41] =-0.3950260023113457; t[41] = -0.9184400985451371; w[41] = 0.03126676038782074;
      r[42] =-0.3950260023113457; s[42] = 0.8747496344107996; t[42] = -0.9184400985451371; w[42] = 0.03126676038782074;
      r[43] =-0.8747496344107996; s[43] = 0.3950260023113457; t[43] = -0.9184400985451371; w[43] = 0.03126676038782074;
      r[44] = 0.3950260023113457; s[44] =-0.8747496344107996; t[44] = -0.9184400985451371; w[44] = 0.03126676038782074;
      r[45] =-0.8747496344107996; s[45] =-0.3950260023113457; t[45] = -0.9184400985451371; w[45] = 0.03126676038782074;
      r[46] =-0.3950260023113457; s[46] =-0.8747496344107996; t[46] = -0.9184400985451371; w[46] = 0.03126676038782074;

      break;
    }
    case  9: {
      r[0]  =                0.0; s[0]  =                0.0; t[0]  =-0.03520741399018411; w[0]  =   0.1162269453717419;
      r[1]  =                0.0; s[1]  =                0.0; t[1]  = -0.5185054856130454; w[1]  =   0.1629908232113593;
      r[2]  = 0.5828238285288264; s[2]  =                0.0; t[2]  = -0.8333333333333333; w[2]  =  0.08418889450221589;
      r[3]  =                0.0; s[3]  = 0.5828238285288264; t[3]  = -0.8333333333333333; w[3]  =  0.08418889450221589;
      r[4]  =-0.5828238285288264; s[4]  =                0.0; t[4]  = -0.8333333333333333; w[4]  =  0.08418889450221589;
      r[5]  =                0.0; s[5]  =-0.5828238285288264; t[5]  = -0.8333333333333333; w[5]  =  0.08418889450221589;
      r[6]  = 0.6760969952248952; s[6]  =                0.0; t[6]  =  -0.460536437675081; w[6]  =  0.05945011023590639;
      r[7]  =                0.0; s[7]  = 0.6760969952248952; t[7]  =  -0.460536437675081; w[7]  =  0.05945011023590639;
      r[8]  =-0.6760969952248952; s[8]  =                0.0; t[8]  =  -0.460536437675081; w[8]  =  0.05945011023590639;
      r[9]  =                0.0; s[9]  =-0.6760969952248952; t[9]  =  -0.460536437675081; w[9]  =  0.05945011023590639;
      r[10] = 0.9258200932795072; s[10] =                0.0; t[10] = -0.9999999859734211; w[10] =  0.01112899365087397;
      r[11] =                0.0; s[11] = 0.9258200932795072; t[11] = -0.9999999859734211; w[11] =  0.01112899365087397;
      r[12] =-0.9258200932795072; s[12] =                0.0; t[12] = -0.9999999859734211; w[12] =  0.01112899365087397;
      r[13] =                0.0; s[13] =-0.9258200932795072; t[13] = -0.9999999859734211; w[13] =  0.01112899365087397;
      r[14] = 0.4319711731007538; s[14] =                0.0; t[14] = 0.06683561264898594; w[14] =   0.0307781514275427;
      r[15] =                0.0; s[15] = 0.4319711731007538; t[15] = 0.06683561264898594; w[15] =   0.0307781514275427;
      r[16] =-0.4319711731007538; s[16] =                0.0; t[16] = 0.06683561264898594; w[16] =   0.0307781514275427;
      r[17] =                0.0; s[17] =-0.4319711731007538; t[17] = 0.06683561264898594; w[17] =   0.0307781514275427;
      r[18] = 0.0799876304341426; s[18] = 0.0799876304341426; t[18] =  0.7450075164629527; w[18] = 0.005139815395587289;
      r[19] = 0.0799876304341426; s[19] =-0.0799876304341426; t[19] =  0.7450075164629527; w[19] = 0.005139815395587289;
      r[20] =-0.0799876304341426; s[20] = 0.0799876304341426; t[20] =  0.7450075164629527; w[20] = 0.005139815395587289;
      r[21] =-0.0799876304341426; s[21] =-0.0799876304341426; t[21] =  0.7450075164629527; w[21] = 0.005139815395587289;
      r[22] =  0.485798731208534; s[22] =  0.485798731208534; t[22] = 0.02829758509547473; w[22] = 0.007371458138614501;
      r[23] =  0.485798731208534; s[23] = -0.485798731208534; t[23] = 0.02829758509547473; w[23] = 0.007371458138614501;
      r[24] = -0.485798731208534; s[24] =  0.485798731208534; t[24] = 0.02829758509547473; w[24] = 0.007371458138614501;
      r[25] = -0.485798731208534; s[25] = -0.485798731208534; t[25] = 0.02829758509547473; w[25] = 0.007371458138614501;
      r[26] = 0.6809331572629768; s[26] = 0.6809331572629768; t[26] = -0.5523213926110991; w[26] =  0.03173539308745657;
      r[27] = 0.6809331572629768; s[27] =-0.6809331572629768; t[27] = -0.5523213926110991; w[27] =  0.03173539308745657;
      r[28] =-0.6809331572629768; s[28] = 0.6809331572629768; t[28] = -0.5523213926110991; w[28] =  0.03173539308745657;
      r[29] =-0.6809331572629768; s[29] =-0.6809331572629768; t[29] = -0.5523213926110991; w[29] =  0.03173539308745657;
      r[30] = 0.1640391767140107; s[30] = 0.1640391767140107; t[30] = -0.9293262207074541; w[30] =  0.02926621055687817;
      r[31] = 0.1640391767140107; s[31] =-0.1640391767140107; t[31] = -0.9293262207074541; w[31] =  0.02926621055687817;
      r[32] =-0.1640391767140107; s[32] = 0.1640391767140107; t[32] = -0.9293262207074541; w[32] =  0.02926621055687817;
      r[33] =-0.1640391767140107; s[33] =-0.1640391767140107; t[33] = -0.9293262207074541; w[33] =  0.02926621055687817;
      r[34] = 0.4057818773820147; s[34] = 0.4057818773820147; t[34] = -0.6132088798220949; w[34] =  0.09173809227488043;
      r[35] = 0.4057818773820147; s[35] =-0.4057818773820147; t[35] = -0.6132088798220949; w[35] =  0.09173809227488043;
      r[36] =-0.4057818773820147; s[36] = 0.4057818773820147; t[36] = -0.6132088798220949; w[36] =  0.09173809227488043;
      r[37] =-0.4057818773820147; s[37] =-0.4057818773820147; t[37] = -0.6132088798220949; w[37] =  0.09173809227488043;
      r[38] = 0.8626559970854881; s[38] = 0.8626559970854881; t[38] = -0.9421756515074574; w[38] =  0.01278107067534859;
      r[39] = 0.8626559970854881; s[39] =-0.8626559970854881; t[39] = -0.9421756515074574; w[39] =  0.01278107067534859;
      r[40] =-0.8626559970854881; s[40] = 0.8626559970854881; t[40] = -0.9421756515074574; w[40] =  0.01278107067534859;
      r[41] =-0.8626559970854881; s[41] =-0.8626559970854881; t[41] = -0.9421756515074574; w[41] =  0.01278107067534859;
      r[42] = 0.1718686868560564; s[42] = 0.1718686868560564; t[42] =  0.3817260245443895; w[42] =  0.03715178465040517;
      r[43] = 0.1718686868560564; s[43] =-0.1718686868560564; t[43] =  0.3817260245443895; w[43] =  0.03715178465040517;
      r[44] =-0.1718686868560564; s[44] = 0.1718686868560564; t[44] =  0.3817260245443895; w[44] =  0.03715178465040517;
      r[45] =-0.1718686868560564; s[45] =-0.1718686868560564; t[45] =  0.3817260245443895; w[45] =  0.03715178465040517;
      r[46] = 0.5481370836221994; s[46] = 0.5481370836221994; t[46] = -0.9827700059896685; w[46] =  0.02949590049620422;
      r[47] = 0.5481370836221994; s[47] =-0.5481370836221994; t[47] = -0.9827700059896685; w[47] =  0.02949590049620422;
      r[48] =-0.5481370836221994; s[48] = 0.5481370836221994; t[48] = -0.9827700059896685; w[48] =  0.02949590049620422;
      r[49] =-0.5481370836221994; s[49] =-0.5481370836221994; t[49] = -0.9827700059896685; w[49] =  0.02949590049620422;
      r[50] = 0.3583323806166026; s[50] = 0.3583323806166026; t[50] = -0.1987364322056734; w[50] =  0.09152296679602322;
      r[51] = 0.3583323806166026; s[51] =-0.3583323806166026; t[51] = -0.1987364322056734; w[51] =  0.09152296679602322;
      r[52] =-0.3583323806166026; s[52] = 0.3583323806166026; t[52] = -0.1987364322056734; w[52] =  0.09152296679602322;
      r[53] =-0.3583323806166026; s[53] =-0.3583323806166026; t[53] = -0.1987364322056734; w[53] =  0.09152296679602322;
      r[54] = 0.4485339792979572; s[54] = 0.8423279096089628; t[54] = -0.8333333333333333; w[54] =  0.03755669131647708;
      r[55] = 0.8423279096089628; s[55] = 0.4485339792979572; t[55] = -0.8333333333333333; w[55] =  0.03755669131647708;
      r[56] = 0.4485339792979572; s[56] =-0.8423279096089628; t[56] = -0.8333333333333333; w[56] =  0.03755669131647708;
      r[57] =-0.8423279096089628; s[57] = 0.4485339792979572; t[57] = -0.8333333333333333; w[57] =  0.03755669131647708;
      r[58] =-0.4485339792979572; s[58] = 0.8423279096089628; t[58] = -0.8333333333333333; w[58] =  0.03755669131647708;
      r[59] = 0.8423279096089628; s[59] =-0.4485339792979572; t[59] = -0.8333333333333333; w[59] =  0.03755669131647708;
      r[60] =-0.4485339792979572; s[60] =-0.8423279096089628; t[60] = -0.8333333333333333; w[60] =  0.03755669131647708;
      r[61] =-0.8423279096089628; s[61] =-0.4485339792979572; t[61] = -0.8333333333333333; w[61] =  0.03755669131647708;

      break;
    }
  }
}
