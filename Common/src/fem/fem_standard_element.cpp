/*!
 * \file fem_standard_element.cpp
 * \brief Functions for the FEM standard elements.
 * \author E. van der Weide
 * \version 8.0.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2023, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../include/fem/fem_standard_element.hpp"
#include "../../include/fem/fem_gauss_jacobi_quadrature.hpp"
#include "../../include/linear_algebra/blas_structure.hpp"

/*----------------------------------------------------------------------------------*/
/*          Public member functions of CFEMStandardElementBase.                     */
/*----------------------------------------------------------------------------------*/

unsigned short CFEMStandardElementBase::GetNDOFsStatic(unsigned short VTK_Type, unsigned short nPoly,
                                                       unsigned long typeErrorMessage) {
  unsigned short nDOFsEdge = nPoly + 1;
  unsigned short nDOFs = 0;  // To avoid a compiler warning.

  switch (VTK_Type) {
    case LINE:
      nDOFs = nDOFsEdge;
      break;

    case TRIANGLE:
      nDOFs = nDOFsEdge * (nDOFsEdge + 1) / 2;
      break;

    case QUADRILATERAL:
      nDOFs = nDOFsEdge * nDOFsEdge;
      break;

    case TETRAHEDRON:
      nDOFs = nDOFsEdge * (nDOFsEdge + 1) * (nDOFsEdge + 2) / 6;
      break;

    case HEXAHEDRON:
      nDOFs = nDOFsEdge * nDOFsEdge * nDOFsEdge;
      break;

    case PRISM:
      nDOFs = nDOFsEdge * nDOFsEdge * (nDOFsEdge + 1) / 2;
      break;

    case PYRAMID:
      nDOFs = nDOFsEdge * (nDOFsEdge + 1) * (2 * nDOFsEdge + 1) / 6;
      break;

    default:
      ostringstream message;
      message << "Unknown FEM element type, " << typeErrorMessage << ", encountered.";
      SU2_MPI::Error(message.str(), CURRENT_FUNCTION);
  }

  return nDOFs;
}

void CFEMStandardElementBase::InverseMatrix(unsigned short n, vector<su2double>& A) {
  /*--- Check the dimensions of A. ---*/
  unsigned long nEntities = n * n;
  if (A.size() != nEntities) SU2_MPI::Error("Wrong size of the A matrix in InverseMatrix", CURRENT_FUNCTION);

  /*--- Create a local matrix to carry out the actual inversion. ---*/
  vector<vector<su2double> > augmentedmatrix(n, vector<su2double>(2 * n));

  /*--- Copy the data from A into the first part of augmentedmatrix. Note
        that A is stored in column major order, such that also Lapack
        routines can be used to invert the matrix.       ---*/
  unsigned int ii = 0;
  for (unsigned short j = 0; j < n; ++j)
    for (unsigned short i = 0; i < n; ++i, ++ii) augmentedmatrix[i][j] = A[ii];

  /*--- Augmenting with identity matrix of similar dimensions ---*/
  for (unsigned short j = 0; j < n; ++j)
    for (unsigned short i = 0; i < n; ++i) augmentedmatrix[i][j + n] = i == j ? 1 : 0;

  /*--- Outer loop of the Gauss-Jordan elimination. ---*/
  for (unsigned short j = 0; j < n; ++j) {
    /*--- Find the pivot in the current column. ---*/
    unsigned short jj = j;
    su2double valMax = fabs(augmentedmatrix[j][j]);
    for (unsigned short i = j + 1; i < n; ++i) {
      su2double val = fabs(augmentedmatrix[i][j]);
      if (val > valMax) {
        jj = i;
        valMax = val;
      }
    }

    /* Swap the rows j and jj, if needed. */
    if (jj > j) {
      for (unsigned short k = j; k < 2 * n; ++k) {
        su2double valTmp = augmentedmatrix[j][k];
        augmentedmatrix[j][k] = augmentedmatrix[jj][k];
        augmentedmatrix[jj][k] = valTmp;
      }
    }

    /*--- Performing row operations to form required identity matrix out
          of the input matrix.              ---*/
    for (unsigned short i = 0; i < n; ++i) {
      if (i != j) {
        valMax = augmentedmatrix[i][j] / augmentedmatrix[j][j];
        for (unsigned short k = j; k < 2 * n; ++k) augmentedmatrix[i][k] -= valMax * augmentedmatrix[j][k];
      }
    }

    valMax = 1.0 / augmentedmatrix[j][j];
    for (unsigned short k = j; k < 2 * n; ++k) augmentedmatrix[j][k] *= valMax;
  }

  /*--- Store the inverse in A. Again column major order is used. ---*/
  ii = 0;
  for (unsigned short j = 0; j < n; ++j)
    for (unsigned short i = 0; i < n; ++i, ++ii) A[ii] = augmentedmatrix[i][j + n];
}

void CFEMStandardElementBase::Vandermonde1D(unsigned short nDOFs, const vector<su2double>& r, vector<su2double>& V) {
  /*--- Determine the number or rows of the Vandermonde matrix and check
        if the dimension of V is correct.     ---*/
  unsigned short nRows = r.size();
  unsigned long nEntities = nRows * nDOFs;
  if (V.size() != nEntities) SU2_MPI::Error("Wrong size of the V matrix", CURRENT_FUNCTION);

  /*--- Compute the Vandermonde matrix. ---*/
  unsigned int ii = 0;
  for (unsigned short i = 0; i < nDOFs; ++i) {
    for (unsigned short k = 0; k < nRows; ++k, ++ii) {
      V[ii] = NormJacobi(i, 0, 0, r[k]);
    }
  }
}

void CFEMStandardElementBase::GradVandermonde1D(unsigned short nDOFs, const vector<su2double>& r,
                                                vector<su2double>& VDr) {
  /*--- Determine the number or rows of the gradient of the Vandermonde matrix
        and check if the dimension of VDr is correct.     ---*/
  unsigned short nRows = r.size();
  unsigned long nEntities = nRows * nDOFs;
  if (VDr.size() != nEntities) SU2_MPI::Error("Wrong size of the VDr matrix", CURRENT_FUNCTION);

  /*--- Compute the gradient of the Vandermonde matrix. ---*/
  unsigned int ii = 0;
  for (unsigned short i = 0; i < nDOFs; ++i) {
    for (unsigned short k = 0; k < nRows; ++k, ++ii) {
      VDr[ii] = GradNormJacobi(i, 0, 0, r[k]);
    }
  }
}

/*----------------------------------------------------------------------------------*/
/*         Protected member functions of CFEMStandardElementBase.                   */
/*----------------------------------------------------------------------------------*/

CFEMStandardElementBase::CFEMStandardElementBase(unsigned short val_VTK_Type, unsigned short val_nPoly,
                                                 bool val_constJac, CConfig* config, unsigned short val_orderExact) {
  /*--- Copy the function arguments that must be stored in the member variables. ---*/
  VTK_Type = val_VTK_Type;
  constJacobian = val_constJac;

  /*--- Determine the polynomial degree that must be integrated exactly by the
        integration rule. If this degree is specified in the argument, just
        copy that value. Otherwise, determine its value. ---*/
  if (val_orderExact > 0) {
    orderExact = val_orderExact;
  } else {
    if (constJacobian)
      orderExact = (unsigned short)ceil(val_nPoly * config->GetQuadrature_Factor_Straight());
    else
      orderExact = (unsigned short)ceil(val_nPoly * config->GetQuadrature_Factor_Curved());
  }

  /*--- Determine the integration points. This depends on the element type. ---*/
  switch (VTK_Type) {
    case LINE:
      IntegrationPointsLine();
      break;
    case TRIANGLE:
      IntegrationPointsTriangle();
      break;
    case QUADRILATERAL:
      IntegrationPointsQuadrilateral();
      break;
    case TETRAHEDRON:
      IntegrationPointsTetrahedron();
      break;
    case PYRAMID:
      IntegrationPointsPyramid();
      break;
    case PRISM:
      IntegrationPointsPrism();
      break;
    case HEXAHEDRON:
      IntegrationPointsHexahedron();
      break;
  }
}

void CFEMStandardElementBase::CheckSumDerivativesLagrangianBasisFunctions(const unsigned short nPoints,
                                                                          const unsigned short nDOFs,
                                                                          const vector<su2double>& dLagBasisPoints) {
  /*--- Check for a zero sum of the derivatives in the given points. ---*/
  for (unsigned short j = 0; j < nPoints; ++j) {
    unsigned int jj = j * nDOFs;
    su2double val = 0.0;
    for (unsigned short i = 0; i < nDOFs; ++i) val += dLagBasisPoints[jj + i];

    if (fabs(val) > 1.e-6) SU2_MPI::Error("Difference is too large to be caused by roundoff", CURRENT_FUNCTION);
  }
}

void CFEMStandardElementBase::CheckSumLagrangianBasisFunctions(const unsigned short nPoints, const unsigned short nDOFs,
                                                               vector<su2double>& lagBasisPoints) {
  /*--- To reduce the error due to round off in the Lagrangian basis functions,
        make sure that the row sum is 1. Also check if the difference is not
        too large to be solely caused by roundoff.   ---*/
  for (unsigned short j = 0; j < nPoints; ++j) {
    unsigned int jj = j * nDOFs;
    su2double val = 0.0;
    for (unsigned short i = 0; i < nDOFs; ++i) val += lagBasisPoints[jj + i];

    if (fabs(val - 1.0) > 1.e-6) SU2_MPI::Error("Difference is too large to be caused by roundoff", CURRENT_FUNCTION);

    val = 1.0 / val;
    for (unsigned short i = 0; i < nDOFs; ++i) lagBasisPoints[jj + i] *= val;
  }
}

void CFEMStandardElementBase::Copy(const CFEMStandardElementBase& other) {
  VTK_Type = other.VTK_Type;
  nIntegration = other.nIntegration;
  constJacobian = other.constJacobian;
  orderExact = other.orderExact;

  rIntegration = other.rIntegration;
  sIntegration = other.sIntegration;
  tIntegration = other.tIntegration;
  wIntegration = other.wIntegration;
}

void CFEMStandardElementBase::DerivativesBasisFunctionsAdjacentElement(
    unsigned short VTK_TypeElem, unsigned short nPolyElem, const bool swapFaceInElement, unsigned short& nDOFsElem,
    vector<su2double>& drLagBasisIntegration, vector<su2double>& dsLagBasisIntegration,
    vector<su2double>& dtLagBasisIntegration) {
  /*--- Define a number of dummy variables, such that the general functions
        to compute the gradients of the basis functions can be used. ---*/
  vector<su2double> rDOFsDummy, sDOFsDummy, tDOFsDummy, matVandermondeInvDummy, lagBasisPointsDummy;

  /*--- Determine the type of the face. ---*/
  switch (VTK_Type) {
    case LINE: {
      /*--- The face element is a line. The adjacent element can be either a
            triangle or a quadrilateral. In order to use the member functions
            LagrangianBasisFunctionAndDerivativesTriangle and
            LagrangianBasisFunctionAndDerivativesQuadrilateral the parametric
            s-coordinates must be set to -1 for both the triangle and
            the quadrilateral. The convention used in this code is that the
            face is face 0 of standard element, which corresponds to s = -1. ---*/
      vector<su2double> sInt(nIntegration, -1.0);

      switch (VTK_TypeElem) {
        case TRIANGLE:
          LagrangianBasisFunctionAndDerivativesTriangle(nPolyElem, rIntegration, sInt, nDOFsElem, rDOFsDummy,
                                                        sDOFsDummy, matVandermondeInvDummy, lagBasisPointsDummy,
                                                        drLagBasisIntegration, dsLagBasisIntegration);
          break;

        case QUADRILATERAL:
          LagrangianBasisFunctionAndDerivativesQuadrilateral(nPolyElem, rIntegration, sInt, nDOFsElem, rDOFsDummy,
                                                             sDOFsDummy, matVandermondeInvDummy, lagBasisPointsDummy,
                                                             drLagBasisIntegration, dsLagBasisIntegration);
          break;
      }

      break;
    }

    case TRIANGLE: {
      /*---- The face element is a triangle. The adjacent element can be a tetrahedron,
             a pyramid or a prism. For the tetrahedron and the prism the convention is
             that the face is face 0 of these elements, which corresponds to a
             parametric coordinate t = -1. For the pyramid the situation is a bit more
             complicated. The face corresponds to face 1 of the pyramid. However, it
             is possible that the connectivity of triangle must be swapped compared to
             the connectivity of the corresponding face of the pyramid. If this is the
             case, also the parametric coordinates of the integration points must be
             swapped in order to get the correct behavior. ---*/

      switch (VTK_TypeElem) {
        case TETRAHEDRON: {
          vector<su2double> tInt(nIntegration, -1.0);

          LagrangianBasisFunctionAndDerivativesTetrahedron(nPolyElem, rIntegration, sIntegration, tInt, nDOFsElem,
                                                           rDOFsDummy, sDOFsDummy, tDOFsDummy, matVandermondeInvDummy,
                                                           lagBasisPointsDummy, drLagBasisIntegration,
                                                           dsLagBasisIntegration, dtLagBasisIntegration);
          break;
        }

        case PRISM: {
          vector<su2double> tInt(nIntegration, -1.0);

          LagrangianBasisFunctionAndDerivativesPrism(nPolyElem, rIntegration, sIntegration, tInt, nDOFsElem, rDOFsDummy,
                                                     sDOFsDummy, tDOFsDummy, matVandermondeInvDummy,
                                                     lagBasisPointsDummy, drLagBasisIntegration, dsLagBasisIntegration,
                                                     dtLagBasisIntegration);
          break;
        }

        case PYRAMID: {
          vector<su2double> rInt(nIntegration), sInt(nIntegration), tInt(nIntegration);
          if (swapFaceInElement) {
            for (unsigned short i = 0; i < nIntegration; ++i) {
              rInt[i] = sIntegration[i] + 0.5 * (rIntegration[i] + 1.0);
              sInt[i] = 0.5 * (rIntegration[i] - 1.0);
              tInt[i] = rIntegration[i];
            }
          } else {
            for (unsigned short i = 0; i < nIntegration; ++i) {
              rInt[i] = rIntegration[i] + 0.5 * (sIntegration[i] + 1.0);
              sInt[i] = 0.5 * (sIntegration[i] - 1.0);
              tInt[i] = sIntegration[i];
            }
          }

          LagrangianBasisFunctionAndDerivativesPyramid(
              nPolyElem, rInt, sInt, tInt, nDOFsElem, rDOFsDummy, sDOFsDummy, tDOFsDummy, matVandermondeInvDummy,
              lagBasisPointsDummy, drLagBasisIntegration, dsLagBasisIntegration, dtLagBasisIntegration);
          break;
        }
      }

      break;
    }

    case QUADRILATERAL: {
      /*--- The face element is a quadrilateral. The adjacent element can be a pyramid,
            a prism or a hexahedron. For the pyramid and the hexahedron the convention
            is that the face is face 0 of these elements, which corresponds to a
            parametric coordinate t = -1. For the prism, the convention is that the
            face is face 2 of the standard prism, which corresponds to a parametric
            coordinate s = -1. Furthermore, for the prism it could be necessary that a
            swap of the parametric coordinates of the face was needed in order to obtain
            the desired sequence. If that is the case, also the integration points of
            the quadrilateral face must be swapped. ---*/

      switch (VTK_TypeElem) {
        case HEXAHEDRON: {
          vector<su2double> tInt(nIntegration, -1.0);

          LagrangianBasisFunctionAndDerivativesHexahedron(nPolyElem, rIntegration, sIntegration, tInt, nDOFsElem,
                                                          rDOFsDummy, sDOFsDummy, tDOFsDummy, matVandermondeInvDummy,
                                                          lagBasisPointsDummy, drLagBasisIntegration,
                                                          dsLagBasisIntegration, dtLagBasisIntegration);
          break;
        }

        case PRISM: {
          vector<su2double> sInt(nIntegration, -1.0), rInt, tInt;

          if (swapFaceInElement) {
            rInt = sIntegration;
            tInt = rIntegration;
          } else {
            rInt = rIntegration;
            tInt = sIntegration;
          }

          LagrangianBasisFunctionAndDerivativesPrism(
              nPolyElem, rInt, sInt, tInt, nDOFsElem, rDOFsDummy, sDOFsDummy, tDOFsDummy, matVandermondeInvDummy,
              lagBasisPointsDummy, drLagBasisIntegration, dsLagBasisIntegration, dtLagBasisIntegration);
          break;
        }

        case PYRAMID: {
          vector<su2double> tInt(nIntegration, -1.0);

          LagrangianBasisFunctionAndDerivativesPyramid(nPolyElem, rIntegration, sIntegration, tInt, nDOFsElem,
                                                       rDOFsDummy, sDOFsDummy, tDOFsDummy, matVandermondeInvDummy,
                                                       lagBasisPointsDummy, drLagBasisIntegration,
                                                       dsLagBasisIntegration, dtLagBasisIntegration);
          break;
        }
      }

      break;
    }
  }
}

void CFEMStandardElementBase::LagrangianBasisFunctionAndDerivativesLine(
    const unsigned short nPoly, const vector<su2double>& rPoints, unsigned short& nDOFs, vector<su2double>& rDOFs,
    vector<su2double>& matVandermondeInv, vector<su2double>& lagBasisPoints, vector<su2double>& drLagBasisPoints) {
  /*--- Determine the number of points in which the functions must
        be determined. ---*/
  const unsigned short nPoints = rPoints.size();

  /*--- Determine the location of the DOFs of the line. ---*/
  nDOFs = nPoly + 1;
  rDOFs.resize(nDOFs);

  su2double dh = 2.0 / nPoly;
  for (unsigned i = 0; i < nDOFs; ++i) rDOFs[i] = -1.0 + i * dh;

  /*--- Compute the inverse of the Vandermonde matrix in the DOFs and
        compute the Vandermonde matrix in the given points. ---*/
  vector<su2double> V(nDOFs * nPoints);
  matVandermondeInv.resize(nDOFs * nDOFs);

  Vandermonde1D(nDOFs, rDOFs, matVandermondeInv);
  InverseMatrix(nDOFs, matVandermondeInv);

  Vandermonde1D(nDOFs, rPoints, V);

  /*--- Allocate the memory for lagBasisPoints and determine its values.
        The Lagrange basis functions in the points are equal to the
        interpolation coefficients from the DOFs to the points and are
        obtained from the matrix product V*Vinv. Note that the result is stored
        in row major order, because in this way the interpolation data for a
        point is contiguous in memory.        ---*/
  lagBasisPoints.resize(nDOFs * nPoints);
  MatMulRowMajor(nDOFs, nPoints, V, matVandermondeInv, lagBasisPoints);

  /*--- Compute the gradients of the 1D Vandermonde matrix in the
        points. The vector V can be used to store the data.     ---*/
  GradVandermonde1D(nDOFs, rPoints, V);

  /*--- Allocate the memory to store the derivatives in r-direction of the
        Lagrange basis functions in the points and determine them.
        The derivatives of the Lagrange basis functions in the points are
        obtained from the matrix product V*Vinv. Note that the result is stored
        in row major order, because in this way the gradient data for a point
        is contiguous in memory.  ---*/
  drLagBasisPoints.resize(nDOFs * nPoints);
  MatMulRowMajor(nDOFs, nPoints, V, matVandermondeInv, drLagBasisPoints);
}

void CFEMStandardElementBase::LagrangianBasisFunctionAndDerivativesTriangle(
    const unsigned short nPoly, const vector<su2double>& rPoints, const vector<su2double>& sPoints,
    unsigned short& nDOFs, vector<su2double>& rDOFs, vector<su2double>& sDOFs, vector<su2double>& matVandermondeInv,
    vector<su2double>& lagBasisPoints, vector<su2double>& drLagBasisPoints, vector<su2double>& dsLagBasisPoints) {
  /*--- Determine the number of points in which the functions must
        be determined. ---*/
  const unsigned short nPoints = rPoints.size();

  /*--- Determine the location of the DOFs of the standard triangle. ---*/
  nDOFs = (nPoly + 1) * (nPoly + 2) / 2;
  rDOFs.resize(nDOFs);
  sDOFs.resize(nDOFs);

  su2double dh = 2.0 / nPoly;

  unsigned int ii = 0;
  for (unsigned short j = 0; j <= nPoly; ++j) {
    su2double s = -1.0 + j * dh;
    unsigned short uppBoundI = nPoly - j;
    for (unsigned short i = 0; i <= uppBoundI; ++i, ++ii) {
      su2double r = -1.0 + i * dh;
      rDOFs[ii] = r;
      sDOFs[ii] = s;
    }
  }

  /*--- Compute the inverse of the Vandermonde matrix in the DOFs and
        compute the Vandermonde matrix in the points. ---*/
  vector<su2double> V(nDOFs * nPoints);
  matVandermondeInv.resize(nDOFs * nDOFs);

  Vandermonde2D_Triangle(nPoly, nDOFs, rDOFs, sDOFs, matVandermondeInv);
  InverseMatrix(nDOFs, matVandermondeInv);

  Vandermonde2D_Triangle(nPoly, nDOFs, rPoints, sPoints, V);

  /*--- Allocate the memory for lagBasisPoints and determine its values.
        The Lagrange basis functions in the points are equal to the interpolation
        coefficients from the DOFs to the points and are obtained from the matrix
        product V*Vinv. Note that the result is stored in row major order, because
        in this way the interpolation data for a point is contiguous in memory. ---*/
  lagBasisPoints.resize(nDOFs * nPoints);
  MatMulRowMajor(nDOFs, nPoints, V, matVandermondeInv, lagBasisPoints);

  /*--- Compute the gradients of the 2D Vandermonde matrix in the points. ---*/
  vector<su2double> VDr(nDOFs * nPoints), VDs(nDOFs * nPoints);
  GradVandermonde2D_Triangle(nPoly, nDOFs, rPoints, sPoints, VDr, VDs);

  /*--- Allocate the memory to store the derivatives in r- and s-direction of the
        Lagrange basis functions in the points and determine them. The derivatives
        of the Lagrange basis functions in the points are obtained from the matrix
        product VDr*Vinv and VDs*Vinv. Note that the result is stored in row major
        order, because in this way the gradient data for a point is contiguous
        in memory. ---*/
  drLagBasisPoints.resize(nDOFs * nPoints);
  dsLagBasisPoints.resize(nDOFs * nPoints);

  MatMulRowMajor(nDOFs, nPoints, VDr, matVandermondeInv, drLagBasisPoints);
  MatMulRowMajor(nDOFs, nPoints, VDs, matVandermondeInv, dsLagBasisPoints);
}

void CFEMStandardElementBase::LagrangianBasisFunctionAndDerivativesQuadrilateral(
    const unsigned short nPoly, const vector<su2double>& rPoints, const vector<su2double>& sPoints,
    unsigned short& nDOFs, vector<su2double>& rDOFs, vector<su2double>& sDOFs, vector<su2double>& matVandermondeInv,
    vector<su2double>& lagBasisPoints, vector<su2double>& drLagBasisPoints, vector<su2double>& dsLagBasisPoints) {
  /*--- Determine the number of points in which the functions must
        be determined. ---*/
  const unsigned short nPoints = rPoints.size();

  /*--- Determine the location of the DOFs of the standard quadrilateral. ---*/
  nDOFs = (nPoly + 1) * (nPoly + 1);
  rDOFs.resize(nDOFs);
  sDOFs.resize(nDOFs);

  su2double dh = 2.0 / nPoly;

  unsigned int ii = 0;
  for (unsigned short j = 0; j <= nPoly; ++j) {
    su2double s = -1.0 + j * dh;
    for (unsigned short i = 0; i <= nPoly; ++i, ++ii) {
      su2double r = -1.0 + i * dh;
      rDOFs[ii] = r;
      sDOFs[ii] = s;
    }
  }

  /*--- Compute the inverse of the Vandermonde matrix in the DOFs and
        compute the Vandermonde matrix in the points. ---*/
  vector<su2double> V(nDOFs * nPoints);
  matVandermondeInv.resize(nDOFs * nDOFs);

  Vandermonde2D_Quadrilateral(nPoly, nDOFs, rDOFs, sDOFs, matVandermondeInv);
  InverseMatrix(nDOFs, matVandermondeInv);

  Vandermonde2D_Quadrilateral(nPoly, nDOFs, rPoints, sPoints, V);

  /*--- Allocate the memory for lagBasisPoints and determine its values.
        The Lagrange basis functions in the points are equal to the interpolation
        coefficients from the DOFs to the points and are obtained from the matrix
        product V*Vinv. Note that the result is stored in row major order, because
        in this way the interpolation data for a point is contiguous in memory. ---*/
  lagBasisPoints.resize(nDOFs * nPoints);
  MatMulRowMajor(nDOFs, nPoints, V, matVandermondeInv, lagBasisPoints);

  /*--- Compute the gradients of the 2D Vandermonde matrix in the points. ---*/
  vector<su2double> VDr(nDOFs * nPoints), VDs(nDOFs * nPoints);
  GradVandermonde2D_Quadrilateral(nPoly, nDOFs, rPoints, sPoints, VDr, VDs);

  /*--- Allocate the memory to store the derivatives in r- and s-direction of the
        Lagrange basis functions in the points and determine them. The derivatives
        of the Lagrange basis functions in the points are obtained from the matrix
        product VDr*Vinv and VDr*Vinv. Note that the result is stored in row major
        order, because in this way the gradient data for a point is contiguous
        in memory. ---*/
  drLagBasisPoints.resize(nDOFs * nPoints);
  dsLagBasisPoints.resize(nDOFs * nPoints);

  MatMulRowMajor(nDOFs, nPoints, VDr, matVandermondeInv, drLagBasisPoints);
  MatMulRowMajor(nDOFs, nPoints, VDs, matVandermondeInv, dsLagBasisPoints);
}

void CFEMStandardElementBase::LagrangianBasisFunctionAndDerivativesTetrahedron(
    const unsigned short nPoly, const vector<su2double>& rPoints, const vector<su2double>& sPoints,
    const vector<su2double>& tPoints, unsigned short& nDOFs, vector<su2double>& rDOFs, vector<su2double>& sDOFs,
    vector<su2double>& tDOFs, vector<su2double>& matVandermondeInv, vector<su2double>& lagBasisPoints,
    vector<su2double>& drLagBasisPoints, vector<su2double>& dsLagBasisPoints, vector<su2double>& dtLagBasisPoints) {
  /*--- Determine the number of points in which the functions must
        be determined. ---*/
  const unsigned short nPoints = rPoints.size();

  /*--- Determine the location of the DOFs of the standard tetrahedron. ---*/
  nDOFs = (nPoly + 1) * (nPoly + 2) * (nPoly + 3) / 6;
  rDOFs.resize(nDOFs);
  sDOFs.resize(nDOFs);
  tDOFs.resize(nDOFs);

  su2double dh = 2.0 / nPoly;

  unsigned int ii = 0;
  for (unsigned short k = 0; k <= nPoly; ++k) {
    su2double t = -1.0 + k * dh;
    unsigned short uppBoundJ = nPoly - k;
    for (unsigned short j = 0; j <= uppBoundJ; ++j) {
      su2double s = -1.0 + j * dh;
      unsigned short uppBoundI = nPoly - k - j;
      for (unsigned short i = 0; i <= uppBoundI; ++i, ++ii) {
        su2double r = -1.0 + i * dh;
        rDOFs[ii] = r;
        sDOFs[ii] = s;
        tDOFs[ii] = t;
      }
    }
  }

  /*--- Compute the inverse of the Vandermonde matrix in the DOFs and
        compute the Vandermonde matrix in the points. ---*/
  vector<su2double> V(nDOFs * nPoints);
  matVandermondeInv.resize(nDOFs * nDOFs);

  Vandermonde3D_Tetrahedron(nPoly, nDOFs, rDOFs, sDOFs, tDOFs, matVandermondeInv);
  InverseMatrix(nDOFs, matVandermondeInv);

  Vandermonde3D_Tetrahedron(nPoly, nDOFs, rPoints, sPoints, tPoints, V);

  /*--- Allocate the memory for lagBasisPoints and determine its values.
        The Lagrange basis functions in the points are equal to the interpolation
        coefficients from the DOFs to the points and are obtained from the matrix
        product V*Vinv. Note that the result is stored in row major order, because
        in this way the interpolation data for a point is contiguous in memory. ---*/
  lagBasisPoints.resize(nDOFs * nPoints);
  MatMulRowMajor(nDOFs, nPoints, V, matVandermondeInv, lagBasisPoints);

  /*--- Compute the gradients of the 3D Vandermonde matrix in the points. ---*/
  vector<su2double> VDr(nDOFs * nPoints), VDs(nDOFs * nPoints), VDt(nDOFs * nPoints);
  GradVandermonde3D_Tetrahedron(nPoly, nDOFs, rPoints, sPoints, tPoints, VDr, VDs, VDt);

  /*--- Allocate the memory to store the derivatives in r-, s- and t-direction
        of the Lagrange basis functions in the points and determine them. The
        derivatives of the Lagrange basis functions in the points are obtained
        from the matrix product VDr*Vinv, VDr*Vinv and VDt*Vinv. Note that the
        result is stored in row major order, because in this way the gradient
        data for a point is contiguous in memory. ---*/
  drLagBasisPoints.resize(nDOFs * nPoints);
  dsLagBasisPoints.resize(nDOFs * nPoints);
  dtLagBasisPoints.resize(nDOFs * nPoints);

  MatMulRowMajor(nDOFs, nPoints, VDr, matVandermondeInv, drLagBasisPoints);
  MatMulRowMajor(nDOFs, nPoints, VDs, matVandermondeInv, dsLagBasisPoints);
  MatMulRowMajor(nDOFs, nPoints, VDt, matVandermondeInv, dtLagBasisPoints);
}

void CFEMStandardElementBase::LagrangianBasisFunctionAndDerivativesPyramid(
    const unsigned short nPoly, const vector<su2double>& rPoints, const vector<su2double>& sPoints,
    const vector<su2double>& tPoints, unsigned short& nDOFs, vector<su2double>& rDOFs, vector<su2double>& sDOFs,
    vector<su2double>& tDOFs, vector<su2double>& matVandermondeInv, vector<su2double>& lagBasisPoints,
    vector<su2double>& drLagBasisPoints, vector<su2double>& dsLagBasisPoints, vector<su2double>& dtLagBasisPoints) {
  /*--- Determine the number of points in which the functions must
        be determined. ---*/
  const unsigned short nPoints = rPoints.size();

  /*--- Allocate the memory for the DOFs of the standard pyramid. ---*/
  unsigned short nDOFsEdge = nPoly + 1;
  nDOFs = nDOFsEdge * (nDOFsEdge + 1) * (2 * nDOFsEdge + 1) / 6;
  rDOFs.resize(nDOFs);
  sDOFs.resize(nDOFs);
  tDOFs.resize(nDOFs);

  /*--- Determine the location of the DOFs of the standard pyramid.
        The outer loop is in the k-direction, which is from base to top. ---*/
  su2double dt = 2.0 / nPoly;
  unsigned short mPoly = nPoly;
  unsigned int ii = 0;

  for (unsigned short k = 0; k <= nPoly; ++k, --mPoly) {
    /*--- Determine the minimum and maximum value for r and s for this t-value. ---*/
    su2double t = -1.0 + k * dt;
    su2double rsMin = 0.5 * (t - 1.0);
    su2double rsMax = -rsMin;

    /*--- Determine the step size along the edges of the current quad.
          Take the exceptional situation mPoly == 0 into account to avoid a
          division by zero.     ---*/
    su2double dh = mPoly ? (rsMax - rsMin) / mPoly : su2double(0.0);

    /*--- Loop over the vertices of the current quadrilateral. ---*/
    for (unsigned short j = 0; j <= mPoly; ++j) {
      su2double s = rsMin + j * dh;
      for (unsigned short i = 0; i <= mPoly; ++i, ++ii) {
        su2double r = rsMin + i * dh;
        rDOFs[ii] = r;
        sDOFs[ii] = s;
        tDOFs[ii] = t;
      }
    }
  }

  /*--- Compute the inverse of the Vandermonde matrix in the DOFs and
        compute the Vandermonde matrix in the points. ---*/
  vector<su2double> V(nDOFs * nPoints);
  matVandermondeInv.resize(nDOFs * nDOFs);

  Vandermonde3D_Pyramid(nPoly, nDOFs, rDOFs, sDOFs, tDOFs, matVandermondeInv);
  InverseMatrix(nDOFs, matVandermondeInv);

  Vandermonde3D_Pyramid(nPoly, nDOFs, rPoints, sPoints, tPoints, V);

  /*--- Allocate the memory for lagBasisPoints and determine its values.
        The Lagrange basis functions in the points are equal to the interpolation
        coefficients from the DOFs to the points and are obtained from the matrix
        product V*Vinv. Note that the result is stored in row major order, because
        in this way the interpolation data for a point is contiguous in memory.  ---*/
  lagBasisPoints.resize(nDOFs * nPoints);
  MatMulRowMajor(nDOFs, nPoints, V, matVandermondeInv, lagBasisPoints);

  /*--- Compute the gradients of the 3D Vandermonde matrix in the points. ---*/
  vector<su2double> VDr(nDOFs * nPoints), VDs(nDOFs * nPoints), VDt(nDOFs * nPoints);
  GradVandermonde3D_Pyramid(nPoly, nDOFs, rPoints, sPoints, tPoints, VDr, VDs, VDt);

  /*--- Allocate the memory to store the derivatives in r-, s- and t-direction
        of the Lagrange basis functions in the points and determine them. The
        derivatives of the Lagrange basis functions in the points are obtained
        from the matrix product VDr*Vinv, VDr*Vinv and VDt*Vinv. Note that the
        result is stored in row major order, because in this way the gradient
        data for a point is contiguous in memory. ---*/
  drLagBasisPoints.resize(nDOFs * nPoints);
  dsLagBasisPoints.resize(nDOFs * nPoints);
  dtLagBasisPoints.resize(nDOFs * nPoints);

  MatMulRowMajor(nDOFs, nPoints, VDr, matVandermondeInv, drLagBasisPoints);
  MatMulRowMajor(nDOFs, nPoints, VDs, matVandermondeInv, dsLagBasisPoints);
  MatMulRowMajor(nDOFs, nPoints, VDt, matVandermondeInv, dtLagBasisPoints);
}

void CFEMStandardElementBase::LagrangianBasisFunctionAndDerivativesPrism(
    const unsigned short nPoly, const vector<su2double>& rPoints, const vector<su2double>& sPoints,
    const vector<su2double>& tPoints, unsigned short& nDOFs, vector<su2double>& rDOFs, vector<su2double>& sDOFs,
    vector<su2double>& tDOFs, vector<su2double>& matVandermondeInv, vector<su2double>& lagBasisPoints,
    vector<su2double>& drLagBasisPoints, vector<su2double>& dsLagBasisPoints, vector<su2double>& dtLagBasisPoints) {
  /*--- Determine the number of points in which the functions must
        be determined. ---*/
  const unsigned short nPoints = rPoints.size();

  /*--- Allocate the memory for the DOFs of the standard prism
        and determine its locations. ---*/
  unsigned short nDOFsEdge = nPoly + 1;
  nDOFs = nDOFsEdge * nDOFsEdge * (nDOFsEdge + 1) / 2;
  rDOFs.resize(nDOFs);
  sDOFs.resize(nDOFs);
  tDOFs.resize(nDOFs);

  su2double dh = 2.0 / nPoly;

  unsigned short ii = 0;
  for (unsigned short k = 0; k <= nPoly; ++k) {
    const su2double t = -1.0 + k * dh;

    for (unsigned short j = 0; j <= nPoly; ++j) {
      su2double s = -1.0 + j * dh;
      unsigned short uppBoundI = nPoly - j;
      for (unsigned short i = 0; i <= uppBoundI; ++i, ++ii) {
        su2double r = -1.0 + i * dh;
        rDOFs[ii] = r;
        sDOFs[ii] = s;
        tDOFs[ii] = t;
      }
    }
  }

  /*--- Compute the inverse of the Vandermonde matrix in the DOFs and
        compute the Vandermonde matrix in the points. ---*/
  vector<su2double> V(nDOFs * nPoints);
  matVandermondeInv.resize(nDOFs * nDOFs);

  Vandermonde3D_Prism(nPoly, nDOFs, rDOFs, sDOFs, tDOFs, matVandermondeInv);
  InverseMatrix(nDOFs, matVandermondeInv);

  Vandermonde3D_Prism(nPoly, nDOFs, rPoints, sPoints, tPoints, V);

  /*--- Allocate the memory for lagBasisPoints and determine its values.
        The Lagrange basis functions in the points are equal to the interpolation
        coefficients from the DOFs to the points and are obtained from the matrix
        product V*Vinv. Note that the result is stored in row major order, because
        in this way the interpolation data for a point is contiguous in memory.  ---*/
  lagBasisPoints.resize(nDOFs * nPoints);
  MatMulRowMajor(nDOFs, nPoints, V, matVandermondeInv, lagBasisPoints);

  /*--- Compute the gradients of the 3D Vandermonde matrix in the points. ---*/
  vector<su2double> VDr(nDOFs * nPoints), VDs(nDOFs * nPoints), VDt(nDOFs * nPoints);
  GradVandermonde3D_Prism(nPoly, nDOFs, rPoints, sPoints, tPoints, VDr, VDs, VDt);

  /*--- Allocate the memory to store the derivatives in r-, s- and t-direction
        of the Lagrange basis functions in the points and determine them. The
        derivatives of the Lagrange basis functions in the points are obtained
        from the matrix product VDr*Vinv, VDr*Vinv and VDt*Vinv. Note that the
        the result is stored in row major order, because in this way the gradient
        data for a point is contiguous in memory. ---*/
  drLagBasisPoints.resize(nDOFs * nPoints);
  dsLagBasisPoints.resize(nDOFs * nPoints);
  dtLagBasisPoints.resize(nDOFs * nPoints);

  MatMulRowMajor(nDOFs, nPoints, VDr, matVandermondeInv, drLagBasisPoints);
  MatMulRowMajor(nDOFs, nPoints, VDs, matVandermondeInv, dsLagBasisPoints);
  MatMulRowMajor(nDOFs, nPoints, VDt, matVandermondeInv, dtLagBasisPoints);
}

void CFEMStandardElementBase::LagrangianBasisFunctionAndDerivativesHexahedron(
    const unsigned short nPoly, const vector<su2double>& rPoints, const vector<su2double>& sPoints,
    const vector<su2double>& tPoints, unsigned short& nDOFs, vector<su2double>& rDOFs, vector<su2double>& sDOFs,
    vector<su2double>& tDOFs, vector<su2double>& matVandermondeInv, vector<su2double>& lagBasisPoints,
    vector<su2double>& drLagBasisPoints, vector<su2double>& dsLagBasisPoints, vector<su2double>& dtLagBasisPoints) {
  /*--- Determine the number of points in which the functions must
        be determined. ---*/
  const unsigned short nPoints = rPoints.size();

  /*--- Allocate the memory for the DOFs of the standard hexahedron
        and determine its locations. ---*/
  unsigned short nDOFsEdge = nPoly + 1;
  nDOFs = nDOFsEdge * nDOFsEdge * nDOFsEdge;
  rDOFs.resize(nDOFs);
  sDOFs.resize(nDOFs);
  tDOFs.resize(nDOFs);

  su2double dh = 2.0 / nPoly;

  unsigned short ii = 0;
  for (unsigned short k = 0; k <= nPoly; ++k) {
    su2double t = -1.0 + k * dh;
    for (unsigned short j = 0; j <= nPoly; ++j) {
      su2double s = -1.0 + j * dh;
      for (unsigned short i = 0; i <= nPoly; ++i, ++ii) {
        su2double r = -1.0 + i * dh;
        rDOFs[ii] = r;
        sDOFs[ii] = s;
        tDOFs[ii] = t;
      }
    }
  }

  /*--- Compute the inverse of the Vandermonde matrix in the DOFs and
        compute the Vandermonde matrix in the points. ---*/
  vector<su2double> V(nDOFs * nPoints);
  matVandermondeInv.resize(nDOFs * nDOFs);

  Vandermonde3D_Hexahedron(nPoly, nDOFs, rDOFs, sDOFs, tDOFs, matVandermondeInv);
  InverseMatrix(nDOFs, matVandermondeInv);

  Vandermonde3D_Hexahedron(nPoly, nDOFs, rPoints, sPoints, tPoints, V);

  /*--- Allocate the memory for lagBasisPoints and determine its values.
        The Lagrange basis functions in the points are equal to the interpolation
        coefficients from the DOFs to the points and are obtained from the matrix
        product V*Vinv. Note that the result is stored in row major order, because
        in this way the interpolation data for a point is contiguous in memory.  ---*/
  lagBasisPoints.resize(nDOFs * nPoints);
  MatMulRowMajor(nDOFs, nPoints, V, matVandermondeInv, lagBasisPoints);

  /*--- Compute the gradients of the 3D Vandermonde matrix in the points. ---*/
  vector<su2double> VDr(nDOFs * nPoints), VDs(nDOFs * nPoints), VDt(nDOFs * nPoints);
  GradVandermonde3D_Hexahedron(nPoly, nDOFs, rPoints, sPoints, tPoints, VDr, VDs, VDt);

  /*--- Allocate the memory to store the derivatives in r-, s- and t-direction
        of the Lagrange basis functions in the points and determine them. The
        derivatives of the Lagrange basis functions in the points are obtained
        from the matrix product VDr*Vinv, VDr*Vinv and VDt*Vinv. Note that the
        the result is stored in row major order, because in this way the gradient
        data for a point is contiguous in memory. ---*/
  drLagBasisPoints.resize(nDOFs * nPoints);
  dsLagBasisPoints.resize(nDOFs * nPoints);
  dtLagBasisPoints.resize(nDOFs * nPoints);

  MatMulRowMajor(nDOFs, nPoints, VDr, matVandermondeInv, drLagBasisPoints);
  MatMulRowMajor(nDOFs, nPoints, VDs, matVandermondeInv, dsLagBasisPoints);
  MatMulRowMajor(nDOFs, nPoints, VDt, matVandermondeInv, dtLagBasisPoints);
}

void CFEMStandardElementBase::MatMulRowMajor(const unsigned short nDOFs, const unsigned short nPoints,
                                             const vector<su2double>& A, const vector<su2double>& B,
                                             vector<su2double>& C) {
  /*--- Check if the dimensions of the matrices correspond to the
        assumptions made in this function.                    ---*/
  const unsigned int dimA = nDOFs * nPoints;
  const unsigned int dimB = nDOFs * nDOFs;

  if (A.size() != dimA || B.size() != dimB || C.size() != dimA)
    SU2_MPI::Error("Unexpected size of the matrices", CURRENT_FUNCTION);

  /*--- Carry out the actual matrix matrix multiplication and store the result
        in row major order (the matrices A and B are in column major order).  ---*/
  for (unsigned short j = 0; j < nDOFs; ++j) {
    for (unsigned short i = 0; i < nPoints; ++i) {
      const unsigned int ii = i * nDOFs + j;
      C[ii] = 0.0;

      for (unsigned short k = 0; k < nDOFs; ++k) {
        const unsigned int indA = k * nPoints + i;
        const unsigned int indB = j * nDOFs + k;

        C[ii] += A[indA] * B[indB];
      }
    }
  }
}

void CFEMStandardElementBase::SubConnForPlottingLine(const unsigned short nPoly, vector<unsigned short>& subConn) {
  /*--- Determine the local subconnectivity of the line element used for plotting
        purposes. This is rather trivial, because the line element is subdivided
        into nPoly linear line elements.                    ---*/
  unsigned short nnPoly = max(nPoly, (unsigned short)1);
  for (unsigned short i = 0; i < nnPoly; ++i) {
    subConn.push_back(i);
    subConn.push_back(i + 1);
  }
}

void CFEMStandardElementBase::SubConnForPlottingQuadrilateral(const unsigned short nPoly,
                                                              vector<unsigned short>& subConn) {
  /*--- Determine the local subconnectivity of the quadrilateral element used for
        plotting purposes. Note that the connectivity of the linear subelements
        obey the VTK connectivity rule of a quadrilateral, which is different
        from the connectivity for the high order quadrilateral. ---*/
  unsigned short nnPoly = max(nPoly, (unsigned short)1);
  for (unsigned short j = 0; j < nnPoly; ++j) {
    unsigned short jj = j * (nnPoly + 1);
    for (unsigned short i = 0; i < nnPoly; ++i) {
      const unsigned short n0 = jj + i;
      subConn.push_back(n0);
      const unsigned short n1 = n0 + 1;
      subConn.push_back(n1);
      const unsigned short n2 = n1 + nPoly + 1;
      subConn.push_back(n2);
      const unsigned short n3 = n2 - 1;
      subConn.push_back(n3);
    }
  }
}

void CFEMStandardElementBase::SubConnForPlottingTriangle(const unsigned short nPoly, vector<unsigned short>& subConn) {
  /*--- Determine the local subconnectivity of the triangular element used for
        plotting purposes. ---*/
  unsigned short jj = 0;

  /*--- Loop over subedges of the left boundary of the standard triangle. ---*/
  for (unsigned short j = 0; j < nPoly; ++j) {
    /*--- Check if the "down" elements must be written. ---*/
    if (j) {
      unsigned short kk = jj - (nPoly + 1 - j);  // Offset of the relevant DOF on the previous row.

      for (unsigned short i = 0; i < (nPoly - j); ++i) {
        unsigned short n0 = jj + i;
        unsigned short n1 = kk + i;
        unsigned short n2 = n0 + 1;

        subConn.push_back(n0);
        subConn.push_back(n1);
        subConn.push_back(n2);
      }
    }

    /*--- The "upp" elements must always be written.
          Determine the offset of the DOF on the next row. ---*/
    unsigned short kk = jj + (nPoly + 1 - j);

    for (unsigned short i = 0; i < (nPoly - j); ++i) {
      unsigned short n0 = jj + i;
      unsigned short n1 = n0 + 1;
      unsigned short n2 = kk + i;

      subConn.push_back(n0);
      subConn.push_back(n1);
      subConn.push_back(n2);
    }

    /*--- Set jj to kk for the next edge. ---*/
    jj = kk;
  }
}

void CFEMStandardElementBase::Vandermonde2D_Triangle(unsigned short nPoly, unsigned short nDOFs,
                                                     const vector<su2double>& r, const vector<su2double>& s,
                                                     vector<su2double>& V) {
  /*--- Determine the number or rows of the Vandermonde matrix and check
        if the dimension of V is correct.     ---*/
  unsigned short nRows = r.size();
  unsigned long nEntities = nRows * nDOFs;
  if (V.size() != nEntities) SU2_MPI::Error("Wrong size of the V matrix", CURRENT_FUNCTION);

  /*--- For a triangle the orthogonal basis for the reference element is obtained
        by a combination of a Jacobi polynomial and a Legendre polynomial. This
        is the result of the orthonormalization of the monomial basis. ---*/
  unsigned int ii = 0;
  for (unsigned short i = 0; i <= nPoly; ++i) {
    for (unsigned short j = 0; j <= (nPoly - i); ++j) {
      for (unsigned short k = 0; k < nRows; ++k, ++ii) {
        /*--- Determine the coefficients a and b. ---*/
        su2double a;
        if (fabs(s[k] - 1.0) < 1.e-8)
          a = -1.0;
        else
          a = 2.0 * (1.0 + r[k]) / (1.0 - s[k]) - 1.0;

        su2double b = s[k];

        /*--- Determine the value of the current basis function in this point. ---*/
        su2double tmp = 1.0;
        if (i) tmp = pow((1.0 - b), i);
        V[ii] = sqrt(2.0) * tmp * NormJacobi(i, 0, 0, a) * NormJacobi(j, 2 * i + 1, 0, b);
      }
    }
  }
}

void CFEMStandardElementBase::GradVandermonde2D_Triangle(unsigned short nPoly, unsigned short nDOFs,
                                                         const vector<su2double>& r, const vector<su2double>& s,
                                                         vector<su2double>& VDr, vector<su2double>& VDs) {
  /*--- Determine the number or rows of the gradient of the Vandermonde matrix
        and check if the dimensions of VDr and VDs are correct.     ---*/
  unsigned short nRows = r.size();
  unsigned long nEntities = nRows * nDOFs;
  if (VDr.size() != nEntities || VDs.size() != nEntities)
    SU2_MPI::Error("Wrong size of the VDr and/or VDs matrices", CURRENT_FUNCTION);

  /*--- For a triangle the orthogonal basis for the reference element is obtained
        by a combination of a Jacobi polynomial and a Legendre polynomial. This
        is the result of the orthonormalization of the monomial basis. ---*/
  unsigned int ii = 0;
  for (unsigned short i = 0; i <= nPoly; ++i) {
    for (unsigned short j = 0; j <= (nPoly - i); ++j) {
      for (unsigned k = 0; k < nRows; ++k, ++ii) {
        /*--- Determine the coefficients a and b. ---*/
        su2double a;
        if (fabs(s[k] - 1.0) < 1.e-8)
          a = -1.0;
        else
          a = 2.0 * (1.0 + r[k]) / (1.0 - s[k]) - 1.0;

        su2double b = s[k];

        /*--- Determine the value of the two 1D contributions to the 2D
              basis functions as well as the gradients of these basis
              functions w.r.t. to their arguments. ---*/
        su2double fa = NormJacobi(i, 0, 0, a);
        su2double gb = NormJacobi(j, 2 * i + 1, 0, b);
        su2double dfa = GradNormJacobi(i, 0, 0, a);
        su2double dgb = GradNormJacobi(j, 2 * i + 1, 0, b);

        /*--- Determine the gradients of the basis functions w.r.t. the
              coordinates r and s. The product rule must be used in order
              to change the derivative of a to the derivative of r and s. ---*/
        VDr[ii] = sqrt(2.0) * dfa * gb;
        VDs[ii] = VDr[ii];
        if (i > 0) {
          su2double tmp = 1.0;
          if (i - 1) tmp = pow((1.0 - b), (i - 1));

          VDr[ii] = 2.0 * tmp * VDr[ii];
          VDs[ii] = (a + 1.0) * tmp * VDs[ii] - i * tmp * sqrt(2.0) * fa * gb;
        }

        su2double tmp = 1.0;
        if (i) tmp = pow((1.0 - b), i);
        VDs[ii] += sqrt(2.0) * fa * dgb * tmp;
      }
    }
  }
}

void CFEMStandardElementBase::Vandermonde2D_Quadrilateral(unsigned short nPoly, unsigned short nDOFs,
                                                          const vector<su2double>& r, const vector<su2double>& s,
                                                          vector<su2double>& V) {
  /*--- Determine the number or rows of the Vandermonde matrix and check
        if the dimension of V is correct.     ---*/
  unsigned short nRows = r.size();
  unsigned long nEntities = nRows * nDOFs;
  if (V.size() != nEntities) SU2_MPI::Error("Wrong size of the V matrix", CURRENT_FUNCTION);

  /*--- For a quadrilateral the basis functions are the product of the 1D
        basis functions, which are the normalized Legendre polynomials.
        The Legendre polynomials are implemented via Jacobi polynomials. ---*/
  unsigned int ii = 0;
  for (unsigned short i = 0; i <= nPoly; ++i) {
    for (unsigned short j = 0; j <= nPoly; ++j) {
      for (unsigned short k = 0; k < nRows; ++k, ++ii) {
        V[ii] = NormJacobi(i, 0, 0, r[k]) * NormJacobi(j, 0, 0, s[k]);
      }
    }
  }
}

void CFEMStandardElementBase::GradVandermonde2D_Quadrilateral(unsigned short nPoly, unsigned short nDOFs,
                                                              const vector<su2double>& r, const vector<su2double>& s,
                                                              vector<su2double>& VDr, vector<su2double>& VDs) {
  /*--- Determine the number or rows of the gradient of the Vandermonde matrix
        and check if the dimensions of VDr and VDs are correct.     ---*/
  unsigned short nRows = r.size();
  unsigned long nEntities = nRows * nDOFs;
  if (VDr.size() != nEntities || VDs.size() != nEntities)
    SU2_MPI::Error("Wrong size of the VDr and/or VDs matrices", CURRENT_FUNCTION);

  /*--- For a quadrilateral the basis functions are the product of the 1D
        basis functions, which are the normalized Legendre polynomials.
        The Legendre polynomials are implemented via Jacobi polynomials.
        Hence the derivatives in r- and s-direction can be computed easily. ---*/
  unsigned int ii = 0;
  for (unsigned short i = 0; i <= nPoly; ++i) {
    for (unsigned short j = 0; j <= nPoly; ++j) {
      for (unsigned short k = 0; k < nRows; ++k, ++ii) {
        VDr[ii] = GradNormJacobi(i, 0, 0, r[k]) * NormJacobi(j, 0, 0, s[k]);
        VDs[ii] = GradNormJacobi(j, 0, 0, s[k]) * NormJacobi(i, 0, 0, r[k]);
      }
    }
  }
}

void CFEMStandardElementBase::Vandermonde3D_Tetrahedron(unsigned short nPoly, unsigned short nDOFs,
                                                        const vector<su2double>& r, const vector<su2double>& s,
                                                        const vector<su2double>& t, vector<su2double>& V) {
  /*--- Determine the number or rows of the Vandermonde matrix and check
        if the dimension of V is correct.     ---*/
  unsigned short nRows = r.size();
  unsigned long nEntities = nRows * nDOFs;
  if (V.size() != nEntities) SU2_MPI::Error("Wrong size of the V matrix", CURRENT_FUNCTION);

  /*--- For a tetrahedron the orthogonal basis for the reference element is obtained by a
        combination of Jacobi polynomials (of which the Legendre polynomials is a special
        case). This is the result of the orthonormalization of the monomial basis. ---*/
  unsigned int ii = 0;
  for (unsigned short i = 0; i <= nPoly; ++i) {
    for (unsigned short j = 0; j <= (nPoly - i); ++j) {
      for (unsigned short k = 0; k <= (nPoly - i - j); ++k) {
        for (unsigned short l = 0; l < nRows; ++l, ++ii) {
          /*--- Determine the coefficients a, b and c. ---*/
          su2double a, b;
          su2double tmp = s[l] + t[l];
          if (fabs(tmp) < 1.e-8)
            a = -1.0;
          else
            a = -1.0 - 2.0 * (1.0 + r[l]) / tmp;

          tmp = 1.0 - t[l];
          if (fabs(tmp) < 1.e-8)
            b = -1.0;
          else
            b = -1.0 + 2.0 * (1.0 + s[l]) / tmp;

          su2double c = t[l];

          /*--- Determine the value of the current basis function in this point. ---*/
          su2double tmpb = 1.0;
          if (i) tmpb = pow((1.0 - b), i);

          su2double tmpc = 1.0;
          if (i + j) tmpc = pow((1.0 - c), i + j);

          V[ii] = sqrt(8.0) * tmpb * tmpc * NormJacobi(i, 0, 0, a) * NormJacobi(j, 2 * i + 1, 0, b) *
                  NormJacobi(k, 2 * (i + j + 1), 0, c);
        }
      }
    }
  }
}

void CFEMStandardElementBase::GradVandermonde3D_Tetrahedron(unsigned short nPoly, unsigned short nDOFs,
                                                            const vector<su2double>& r, const vector<su2double>& s,
                                                            const vector<su2double>& t, vector<su2double>& VDr,
                                                            vector<su2double>& VDs, vector<su2double>& VDt) {
  /*--- Determine the number or rows of the gradient of the Vandermonde matrix
        and check if the dimensions of VDr, VDs and VDt are correct.     ---*/
  unsigned short nRows = r.size();
  unsigned long nEntities = nRows * nDOFs;
  if (VDr.size() != nEntities || VDs.size() != nEntities || VDt.size() != nEntities)
    SU2_MPI::Error("Wrong size of the VDr, VDs and VDt matrices", CURRENT_FUNCTION);

  /*--- For a tetrahedron the orthogonal basis for the reference element is obtained by a
        combination of Jacobi polynomials (of which the Legendre polynomials is a special
        case). This is the result of the orthonormalization of the monomial basis.
        Note that the sequence of the i, j and k loop must be identical to
        the evaluation of the Vandermonde matrix itself.                ---*/
  unsigned int ii = 0;
  for (unsigned short i = 0; i <= nPoly; ++i) {
    for (unsigned short j = 0; j <= (nPoly - i); ++j) {
      for (unsigned short k = 0; k <= (nPoly - i - j); ++k) {
        for (unsigned short l = 0; l < nRows; ++l, ++ii) {
          /*--- Determine the coefficients a, b and c. */
          su2double a, b;
          su2double tmp = s[l] + t[l];
          if (fabs(tmp) < 1.e-8)
            a = -1.0;
          else
            a = -1.0 - 2.0 * (1.0 + r[l]) / tmp;

          tmp = 1.0 - t[l];
          if (fabs(tmp) < 1.e-8)
            b = -1.0;
          else
            b = -1.0 + 2.0 * (1.0 + s[l]) / tmp;

          su2double c = t[l];

          /*--- Determine the value of the three 1D contributions to the 3D basis functions as
                well as the gradients of these basis functions w.r.t. to their arguments. ---*/
          su2double fa = NormJacobi(i, 0, 0, a);
          su2double gb = NormJacobi(j, 2 * i + 1, 0, b);
          su2double hc = NormJacobi(k, 2 * (i + j + 1), 0, c);
          su2double dfa = GradNormJacobi(i, 0, 0, a);
          su2double dgb = GradNormJacobi(j, 2 * i + 1, 0, b);
          su2double dhc = GradNormJacobi(k, 2 * (i + j + 1), 0, c);

          /*--- Compute the derivative of the basis function w.r.t. r. As r is only present in
                the parameter a the derivative of the basis function w.r.t. a is multiplied by
                dadr. Note that the implementation is such that all possible singularities are
                divided out of the expression.                                  ---*/
          VDr[ii] = sqrt(8.0) * dfa * gb * hc;
          if (i > 0) {
            VDr[ii] *= 4.0;
            if (i > 1) VDr[ii] *= pow((1.0 - b), (i - 1));
          }

          if (i + j > 1) VDr[ii] *= pow((1.0 - c), (i + j - 1));

          /*--- Compute the derivative of the basis function w.r.t. s. As s is present in both
                the parameters a and b, both variables must be taken into account when the
                derivative is computed. Note that the implementation is such that all possible
                singularities are divided out of the expression. The first part is the derivative
                of the basis function w.r.t. b multiplied by dbds. This value is stored, because
                it is needed later on to compute the derivative w.r.t. t.       ---*/
          VDs[ii] = dgb;
          if (i) VDs[ii] *= pow((1.0 - b), i);

          if (i > 0) {
            tmp = i * gb;
            if (i > 1) tmp *= pow((1.0 - b), (i - 1));
            VDs[ii] -= tmp;
          }

          if (i + j > 0) {
            VDs[ii] *= 2.0 * sqrt(8.0) * fa * hc;
            if (i + j > 1) VDs[ii] *= pow((1.0 - c), (i + j - 1));
          }

          su2double dPsidbXdbds = VDs[ii];

          /*--- Add the contribution from the derivative of the basis function
                w.r.t. a multiplied by dads.           ---*/
          VDs[ii] += 0.5 * (a + 1.0) * VDr[ii];

          /*--- Compute the derivative of the basis function w.r.t. t. As t is present in a, b and c,
                all parameters must be taken into account when the derivative is computed. Note that
                the implementation is such that all possible singularities are divided out of the
                expression. The first part is the derivative of the basis function w.r.t. c,
                which is equal to t.                                     ---*/
          VDt[ii] = dhc;
          if (i + j > 0) {
            VDt[ii] *= pow((1.0 - c), (i + j));

            tmp = (i + j) * hc;
            if (i + j > 1) tmp *= pow((1.0 - c), (i + j - 1));

            VDt[ii] -= tmp;
          }

          VDt[ii] *= sqrt(8.0) * fa * gb;
          if (i) VDt[ii] *= pow((1.0 - b), i);

          /*--- Add the contribution from the derivative of the basis function w.r.t. a multiplied
                by dadt and the derivative w.r.t. b multiplied by dbdt.           ---*/
          VDt[ii] += 0.5 * (a + 1.0) * VDr[ii] + 0.5 * (b + 1.0) * dPsidbXdbds;
        }
      }
    }
  }
}

void CFEMStandardElementBase::Vandermonde3D_Pyramid(unsigned short nPoly, unsigned short nDOFs,
                                                    const vector<su2double>& r, const vector<su2double>& s,
                                                    const vector<su2double>& t, vector<su2double>& V) {
  /*--- Determine the number or rows of the Vandermonde matrix and check
        if the dimension of V is correct.     ---*/
  unsigned short nRows = r.size();
  unsigned long nEntities = nRows * nDOFs;
  if (V.size() != nEntities) SU2_MPI::Error("Wrong size of the V matrix", CURRENT_FUNCTION);

  /*--- For a pyramid the orthogonal basis for the reference element is
        obtained by a combination of Jacobi polynomials (of which the Legendre
        polynomials is a special case). This is the result of the
        orthonormalization of the monomial basis. ---*/
  unsigned int ii = 0;
  for (unsigned short i = 0; i <= nPoly; ++i) {
    for (unsigned short j = 0; j <= nPoly; ++j) {
      unsigned short muij = max(i, j);
      const su2double scaleFact = pow(2, muij + 1);
      for (unsigned short k = 0; k <= (nPoly - muij); ++k) {
        for (unsigned short l = 0; l < nRows; ++l, ++ii) {
          /*--- Determine the coefficients a, b and c. ---*/
          su2double a, b;
          su2double tmp = 0.5 * (1.0 - t[l]);
          if (fabs(tmp) < 1.e-8)
            a = b = 0.0;
          else {
            a = r[l] / tmp;
            b = s[l] / tmp;
          }

          su2double c = t[l];

          /*--- Determine the value of the current basis function in this point.
                The multiplication with scaleFact is necessary, because this
                formulation is derived in literature for a t coordinate between
                0 and 1, while in this code t varies between -1 and 1. ---*/
          su2double tmpt = 1.0;
          if (muij) tmpt = pow(tmp, muij);

          V[ii] =
              scaleFact * tmpt * NormJacobi(i, 0, 0, a) * NormJacobi(j, 0, 0, b) * NormJacobi(k, 2 * (muij + 1), 0, c);
        }
      }
    }
  }
}

void CFEMStandardElementBase::GradVandermonde3D_Pyramid(unsigned short nPoly, unsigned short nDOFs,
                                                        const vector<su2double>& r, const vector<su2double>& s,
                                                        const vector<su2double>& t, vector<su2double>& VDr,
                                                        vector<su2double>& VDs, vector<su2double>& VDt) {
  /*--- Determine the number or rows of the gradient of the Vandermonde matrix
        and check if the dimensions of VDr, VDs and VDt are correct.     ---*/
  unsigned short nRows = r.size();
  unsigned long nEntities = nRows * nDOFs;
  if (VDr.size() != nEntities || VDs.size() != nEntities || VDt.size() != nEntities)
    SU2_MPI::Error("Wrong size of the VDr, VDs and VDt matrices", CURRENT_FUNCTION);

  /*--- For a pyramid the orthogonal basis for the reference element is
        obtained by a combination of Jacobi polynomials (of which the Legendre
        polynomials is a special case). This is the result of the
        orthonormalization of the monomial basis.
        Note that the sequence of the i, j and k loop must be identical to
        the evaluation of the Vandermonde matrix itself.  ---*/
  unsigned int ii = 0;
  for (unsigned short i = 0; i <= nPoly; ++i) {
    for (unsigned short j = 0; j <= nPoly; ++j) {
      unsigned short muij = max(i, j);
      const su2double scaleFact = pow(2, muij + 1);
      for (unsigned short k = 0; k <= (nPoly - muij); ++k) {
        for (unsigned short l = 0; l < nRows; ++l, ++ii) {
          /*--- Determine the coefficients a, b and c. ---*/
          su2double a, b;
          su2double tmp = 0.5 * (1.0 - t[l]);
          if (fabs(tmp) < 1.e-8)
            a = b = 0.0;
          else {
            a = r[l] / tmp;
            b = s[l] / tmp;
          }

          su2double c = t[l];

          /*--- Determine the value of the three 1D contributions to the 3D
                basis functions as well as the gradients of these basis
                functions w.r.t. to their arguments. ---*/
          su2double fa = NormJacobi(i, 0, 0, a);
          su2double gb = NormJacobi(j, 0, 0, b);
          su2double hc = NormJacobi(k, 2 * (muij + 1), 0, c);
          su2double dfa = GradNormJacobi(i, 0, 0, a);
          su2double dgb = GradNormJacobi(j, 0, 0, b);
          su2double dhc = GradNormJacobi(k, 2 * (muij + 1), 0, c);

          /*--- Compute the derivative of the basis function w.r.t. r and s.
                As r is only present in the parameter a the derivative of
                the basis function w.r.t. a is multiplied by dadr. A similar
                argument holds for s, which is only present in the parameter b.
                Note that the implementation is such that all possible
                singularities are divided out of the expression.  ---*/
          VDr[ii] = dfa * gb * hc;
          VDs[ii] = fa * dgb * hc;
          if (muij > 1) {
            su2double tmpt = pow(tmp, (muij - 1));
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
          VDt[ii] = dhc;
          if (muij > 0) VDt[ii] *= pow(tmp, muij);

          if (muij > 0) {
            su2double tmpt = 0.5 * muij * hc;
            if (muij > 1) tmpt *= pow(tmp, (muij - 1));
            VDt[ii] -= tmpt;
          }

          VDt[ii] *= fa * gb;

          /*--- Add the contribution from the derivative of the basis function
                w.r.t. a multiplied by dadt and the derivative w.r.t. b multiplied
                by dbdt.                      ---*/
          VDt[ii] += 0.5 * a * VDr[ii] + 0.5 * b * VDs[ii];

          /*--- Multiply the three derivatives with the scale factor to
                obtain the correct answers. See Vandermonde3D_Pyramid for
                the explanation. ---*/
          VDr[ii] *= scaleFact;
          VDs[ii] *= scaleFact;
          VDt[ii] *= scaleFact;
        }
      }
    }
  }
}

void CFEMStandardElementBase::Vandermonde3D_Prism(unsigned short nPoly, unsigned short nDOFs,
                                                  const vector<su2double>& r, const vector<su2double>& s,
                                                  const vector<su2double>& t, vector<su2double>& V) {
  /*--- Determine the number or rows of the Vandermonde matrix and check
        if the dimension of V is correct.     ---*/
  unsigned short nRows = r.size();
  unsigned long nEntities = nRows * nDOFs;
  if (V.size() != nEntities) SU2_MPI::Error("Wrong size of the V matrix", CURRENT_FUNCTION);

  /*--- For a prism the orthogonal basis for the reference element is a tensor
        product of the 1D basis functions in the structured direction of the prism
        and the basis functions of a triangle. For that triangle the orthogonal
        basis is obtained by a combination of a Jacobi polynomial and a Legendre
        polynomial. This is the result of the orthonormalization of the
        monomial basis.                   ---*/
  unsigned int ii = 0;
  for (unsigned short i = 0; i <= nPoly; ++i) {
    for (unsigned short j = 0; j <= (nPoly - i); ++j) {
      for (unsigned short k = 0; k <= nPoly; ++k) {
        for (unsigned short l = 0; l < nRows; ++l, ++ii) {
          /*--- Determine the coefficients a and b. ---*/
          su2double a;
          if (fabs(s[l] - 1.0) < 1.e-8)
            a = -1.0;
          else
            a = 2.0 * (1.0 + r[l]) / (1.0 - s[l]) - 1.0;

          su2double b = s[l];

          /*--- Determine the value of the current basis function in this point. ---*/
          su2double tmp = 1.0;
          if (i) tmp = pow((1.0 - b), i);

          V[ii] = sqrt(2.0) * tmp * NormJacobi(i, 0, 0, a) * NormJacobi(j, 2 * i + 1, 0, b) * NormJacobi(k, 0, 0, t[l]);
        }
      }
    }
  }
}

void CFEMStandardElementBase::GradVandermonde3D_Prism(unsigned short nPoly, unsigned short nDOFs,
                                                      const vector<su2double>& r, const vector<su2double>& s,
                                                      const vector<su2double>& t, vector<su2double>& VDr,
                                                      vector<su2double>& VDs, vector<su2double>& VDt) {
  /*--- Determine the number or rows of the gradient of the Vandermonde matrix
        and check if the dimensions of VDr, VDs and VDt are correct.     ---*/
  unsigned short nRows = r.size();
  unsigned long nEntities = nRows * nDOFs;
  if (VDr.size() != nEntities || VDs.size() != nEntities || VDt.size() != nEntities)
    SU2_MPI::Error("Wrong size of the VDr, VDs and VDt matrices", CURRENT_FUNCTION);

  /*--- For a prism the orthogonal basis for the reference element is a tensor
        product of the 1D basis functions in the structured direction of the prism
        and the basis functions of a triangle. Hence the derivative matrices also
        follows this tensor product rule.
        Note that the sequence of the i, j and k loop must be identical to
        the evaluation of the Vandermonde matrix itself.          ---*/
  unsigned int ii = 0;
  for (unsigned short i = 0; i <= nPoly; ++i) {
    for (unsigned short j = 0; j <= (nPoly - i); ++j) {
      for (unsigned short k = 0; k <= nPoly; ++k) {
        for (unsigned short l = 0; l < nRows; ++l, ++ii) {
          /*--- Determine the coefficients a and b. ---*/
          su2double a;
          if (fabs(s[l] - 1.0) < 1.e-8)
            a = -1.0;
          else
            a = 2.0 * (1.0 + r[l]) / (1.0 - s[l]) - 1.0;

          su2double b = s[l];

          /*--- Determine the value of the two 1D contributions to the 2D
                basis functions of the triangle as well as the gradients of
                these basis functions w.r.t. to its argument. ---*/
          su2double fa = NormJacobi(i, 0, 0, a);
          su2double gb = NormJacobi(j, 2 * i + 1, 0, b);
          su2double dfa = GradNormJacobi(i, 0, 0, a);
          su2double dgb = GradNormJacobi(j, 2 * i + 1, 0, b);

          /*--- Determine the gradients of the basis functions w.r.t. the
                coordinates r and s. The product rule must be used in order
                to change the derivative of a to the derivative of r and s. ---*/
          VDr[ii] = sqrt(2.0) * dfa * gb;
          VDs[ii] = VDr[ii];
          if (i > 0) {
            su2double tmp = 1.0;
            if (i > 1) tmp = pow((1.0 - b), (i - 1));

            VDr[ii] = 2.0 * tmp * VDr[ii];
            VDs[ii] = (a + 1.0) * tmp * VDs[ii] - i * tmp * sqrt(2.0) * fa * gb;
          }

          su2double tmp = 1.0;
          if (i > 0) tmp = pow((1.0 - b), i);

          VDs[ii] += sqrt(2.0) * fa * dgb * tmp;

          /*--- Multiply VDr and VDs with the contribution from the structured
                direction of the prism.                 ---*/
          VDr[ii] *= NormJacobi(k, 0, 0, t[l]);
          VDs[ii] *= NormJacobi(k, 0, 0, t[l]);

          /*--- Compute the derivative of the basis function in the t-direction,
                which is the structured direction.             ---*/
          VDt[ii] = sqrt(2.0) * tmp * fa * gb * GradNormJacobi(k, 0, 0, t[l]);
        }
      }
    }
  }
}

void CFEMStandardElementBase::Vandermonde3D_Hexahedron(unsigned short nPoly, unsigned short nDOFs,
                                                       const vector<su2double>& r, const vector<su2double>& s,
                                                       const vector<su2double>& t, vector<su2double>& V) {
  /*--- Determine the number or rows of the Vandermonde matrix and check
        if the dimension of V is correct.     ---*/
  unsigned short nRows = r.size();
  unsigned long nEntities = nRows * nDOFs;
  if (V.size() != nEntities) SU2_MPI::Error("Wrong size of the V matrix", CURRENT_FUNCTION);

  /*--- For a hexahedron the basis functions are the tensor product of the 1D
        basis functions, which are the normalized Legendre polynomials. Note
        that the Legendre polynomials are a special kind of Jacobi polynomials. ---*/
  unsigned int ii = 0;
  for (unsigned short i = 0; i <= nPoly; ++i) {
    for (unsigned short j = 0; j <= nPoly; ++j) {
      for (unsigned short k = 0; k <= nPoly; ++k) {
        for (unsigned short l = 0; l < nRows; ++l, ++ii) {
          V[ii] = NormJacobi(i, 0, 0, r[l]) * NormJacobi(j, 0, 0, s[l]) * NormJacobi(k, 0, 0, t[l]);
        }
      }
    }
  }
}

void CFEMStandardElementBase::GradVandermonde3D_Hexahedron(unsigned short nPoly, unsigned short nDOFs,
                                                           const vector<su2double>& r, const vector<su2double>& s,
                                                           const vector<su2double>& t, vector<su2double>& VDr,
                                                           vector<su2double>& VDs, vector<su2double>& VDt) {
  /*--- Determine the number or rows of the gradient of the Vandermonde matrix
        and check if the dimensions of VDr, VDs and VDt are correct.     ---*/
  unsigned short nRows = r.size();
  unsigned long nEntities = nRows * nDOFs;
  if (VDr.size() != nEntities || VDs.size() != nEntities || VDt.size() != nEntities)
    SU2_MPI::Error("Wrong size of the VDr, VDs and VDt matrices", CURRENT_FUNCTION);

  /*--- For a hexahedron the basis functions are the tensor product of the 1D
        basis functions, which are the normalized Legendre polynomials.
        The derivatives are therefore easy to compute.
        Note that the Legendre polynomials are a special kind of Jacobi polynomials.
        Also note that the sequence of the i, j and k loop must be identical to
        the evaluation of the Vandermonde matrix itself.      ---*/
  unsigned int ii = 0;
  for (unsigned short i = 0; i <= nPoly; ++i) {
    for (unsigned short j = 0; j <= nPoly; ++j) {
      for (unsigned short k = 0; k <= nPoly; ++k) {
        for (unsigned short l = 0; l < nRows; ++l, ++ii) {
          VDr[ii] = NormJacobi(j, 0, 0, s[l]) * NormJacobi(k, 0, 0, t[l]) * GradNormJacobi(i, 0, 0, r[l]);
          VDs[ii] = NormJacobi(i, 0, 0, r[l]) * NormJacobi(k, 0, 0, t[l]) * GradNormJacobi(j, 0, 0, s[l]);
          VDt[ii] = NormJacobi(i, 0, 0, r[l]) * NormJacobi(j, 0, 0, s[l]) * GradNormJacobi(k, 0, 0, t[l]);
        }
      }
    }
  }
}

su2double CFEMStandardElementBase::ViscousPenaltyParameter(const unsigned short VTK_TypeElem,
                                                           const unsigned short nPolyElem) const {
  /*--- Determine the element type and set the value of the penalty
        parameter accordingly. ---*/
  su2double penParam = 0;  // To avoid a compiler warning.

  switch (VTK_TypeElem) {
    case TRIANGLE:
      penParam = (nPolyElem + 1) * (nPolyElem + 2) * 0.5;
      break;
    case QUADRILATERAL:
      penParam = (nPolyElem + 1) * (nPolyElem + 1);
      break;
    case TETRAHEDRON:
      penParam = (nPolyElem + 1) * (nPolyElem + 3) / 3.0;
      break;
    case PYRAMID:
      if (VTK_Type == TRIANGLE)
        penParam = (nPolyElem + 1) * (2 * nPolyElem + 3) * 1.05 / 3.0;
      else
        penParam = (nPolyElem + 1) * (nPolyElem + 3) / 3.0;
      break;
    case PRISM:
      if (VTK_Type == TRIANGLE)
        penParam = (nPolyElem + 1) * (nPolyElem + 1);
      else
        penParam = (nPolyElem + 1) * (nPolyElem + 2) * 0.5;
      break;
    case HEXAHEDRON:
      penParam = (nPolyElem + 1) * (nPolyElem + 1);
      break;
  }

  /* Return the value of penParam. */
  return penParam;
}

/*----------------------------------------------------------------------------------*/
/*          Private member functions of CFEMStandardElementBase.                    */
/*----------------------------------------------------------------------------------*/

void CFEMStandardElementBase::GaussLegendrePoints1D(vector<su2double>& GLPoints, vector<su2double>& GLWeights) {
  /* The class used to determine the integration points operate on passivedoubles.
     Allocate the memory for the help vectors. */
  vector<passivedouble> GLPointsPas(GLPoints.size());
  vector<passivedouble> GLWeightsPas(GLWeights.size());

  /* Gauss Legendre quadrature is a special case of Gauss Jacobi integration.
     Determine the integration points for this case. */
  CGaussJacobiQuadrature GaussJacobi;
  GaussJacobi.GetQuadraturePoints(0.0, 0.0, -1.0, 1.0, GLPointsPas, GLWeightsPas);

  /* Copy the data back into GLPoints and GLWeights. */
  for (unsigned long i = 0; i < GLPoints.size(); ++i) {
    GLPoints[i] = GLPointsPas[i];
    GLWeights[i] = GLWeightsPas[i];
  }
}

su2double CFEMStandardElementBase::NormJacobi(unsigned short n, unsigned short alpha, unsigned short beta,
                                              su2double x) {
  /*--- Some abbreviations. ---*/
  su2double ap1 = alpha + 1;
  su2double bp1 = beta + 1;
  su2double apb = alpha + beta;
  su2double apbp1 = apb + 1;
  su2double apbp2 = apb + 2;
  su2double apbp3 = apb + 3;
  su2double b2ma2 = beta * beta - alpha * alpha;

  /*--- Determine the term, which involves the gamma function for P0. As the
        arguments are integers, this term can be computed easily, because
        Gamma(n+1) = n!. ---*/
  su2double Gamap1 = 1.0, Gambp1 = 1.0, Gamapbp2 = 1.0;
  for (unsigned short i = 2; i <= alpha; ++i) Gamap1 *= i;
  for (unsigned short i = 2; i <= beta; ++i) Gambp1 *= i;
  for (unsigned short i = 2; i <= (alpha + beta + 1); ++i) Gamapbp2 *= i;

  /*--- Initialize the normalized polynomials. ---*/
  su2double Pnm1 = sqrt(pow(0.5, apbp1) * Gamapbp2 / (Gamap1 * Gambp1));
  su2double Pn = 0.5 * Pnm1 * (apbp2 * x + alpha - beta) * sqrt(apbp3 / (ap1 * bp1));

  /*--- Take care of the special situation of n == 0. ---*/
  if (n == 0)
    Pn = Pnm1;
  else {
    /*--- The value of the normalized Jacobi polynomial must be obtained via recursion. ---*/
    for (unsigned short i = 2; i <= n; ++i) {
      /*--- Compute the coefficients a for i and i-1 and the coefficient bi. ---*/
      unsigned short j = i - 1;
      su2double tmp = 2 * j + apb;
      su2double aim1 = 2.0 * sqrt(j * (j + apb) * (j + alpha) * (j + beta) / ((tmp - 1.0) * (tmp + 1.0))) / tmp;

      su2double bi = b2ma2 / (tmp * (tmp + 2.0));

      tmp = 2 * i + apb;
      su2double ai = 2.0 * sqrt(i * (i + apb) * (i + alpha) * (i + beta) / ((tmp - 1.0) * (tmp + 1.0))) / tmp;

      /*--- Compute the new value of Pn and make sure to store Pnm1 correctly. ---*/
      tmp = Pnm1;
      Pnm1 = Pn;

      Pn = ((x - bi) * Pn - aim1 * tmp) / ai;
    }
  }

  /*--- Return Pn. ---*/
  return Pn;
}

su2double CFEMStandardElementBase::GradNormJacobi(unsigned short n, unsigned short alpha, unsigned short beta,
                                                  su2double x) {
  /*--- Make a distinction for n == 0 and n > 0. For n == 0 the derivative is
        zero, because the polynomial itself is constant. ---*/
  su2double grad;
  if (n == 0)
    grad = 0.0;
  else {
    su2double tmp = n * (n + alpha + beta + 1.0);
    grad = sqrt(tmp) * NormJacobi(n - 1, alpha + 1, beta + 1, x);
  }

  /*--- Return the gradient. ---*/
  return grad;
}

/*----------------------------------------------------------------------------------*/
/*           Public member functions of CFEMStandardElement.                        */
/*----------------------------------------------------------------------------------*/

CFEMStandardElement::CFEMStandardElement(unsigned short val_VTK_Type, unsigned short val_nPoly, bool val_constJac,
                                         CConfig* config, unsigned short val_orderExact,
                                         const vector<su2double>* rLocSolDOFs, const vector<su2double>* sLocSolDOFs,
                                         const vector<su2double>* tLocSolDOFs)

    : CFEMStandardElementBase(val_VTK_Type, val_nPoly, val_constJac, config, val_orderExact) {
  /*--- Copy the function arguments to the member variables. ---*/
  nPoly = val_nPoly;

  /*--- Determine the element type and compute the other member variables. ---*/
  switch (VTK_Type) {
    case LINE:
      DataStandardLine();
      break;
    case TRIANGLE:
      DataStandardTriangle();
      break;
    case QUADRILATERAL:
      DataStandardQuadrilateral();
      break;
    case TETRAHEDRON:
      DataStandardTetrahedron();
      break;
    case PYRAMID:
      DataStandardPyramid();
      break;
    case PRISM:
      DataStandardPrism();
      break;
    case HEXAHEDRON:
      DataStandardHexahedron();
      break;
  }

  /*--------------------------------------------------------------------------*/
  /*--- Create the data of the basis functions and its derivatives in the  ---*/
  /*--- integration points of the element.                                 ---*/
  /*--------------------------------------------------------------------------*/

  /*--- Check the sum of the Lagrangian basis functions. ---*/
  CheckSumLagrangianBasisFunctions(nIntegration, nDOFs, lagBasisIntegration);

  /*--- Check the sum of the derivatives of the Lagrangian basis functions. ---*/
  if (!drLagBasisIntegration.empty())
    CheckSumDerivativesLagrangianBasisFunctions(nIntegration, nDOFs, drLagBasisIntegration);
  if (!dsLagBasisIntegration.empty())
    CheckSumDerivativesLagrangianBasisFunctions(nIntegration, nDOFs, dsLagBasisIntegration);
  if (!dtLagBasisIntegration.empty())
    CheckSumDerivativesLagrangianBasisFunctions(nIntegration, nDOFs, dtLagBasisIntegration);

  /*--- Create the transpose of lagBasisIntegration. This is needed for the
        efficient residual computation for the ADER-DG predictor step. ---*/
  lagBasisIntegrationTrans.resize(lagBasisIntegration.size());

  unsigned int ii = 0;
  for (unsigned short j = 0; j < nDOFs; ++j) {
    for (unsigned short i = 0; i < nIntegration; ++i, ++ii) {
      const unsigned int ind = i * nDOFs + j;
      lagBasisIntegrationTrans[ii] = lagBasisIntegration[ind];
    }
  }

  /*--- For efficiency reasons it is beneficial to store lagBasisIntegration,
        drLagBasisIntegration, dsLagBasisIntegration and dtLagBasisIntegration in
        one array. ---*/
  const unsigned long sizeDerMat =
      drLagBasisIntegration.size() + dsLagBasisIntegration.size() + dtLagBasisIntegration.size();
  const unsigned long sizeMat = lagBasisIntegration.size() + sizeDerMat;

  matBasisIntegration.resize(sizeMat);
  matDerBasisIntTrans.resize(sizeDerMat);

  ii = 0;
  for (unsigned long i = 0; i < lagBasisIntegration.size(); ++i, ++ii) matBasisIntegration[ii] = lagBasisIntegration[i];

  for (unsigned long i = 0; i < drLagBasisIntegration.size(); ++i, ++ii)
    matBasisIntegration[ii] = drLagBasisIntegration[i];

  for (unsigned long i = 0; i < dsLagBasisIntegration.size(); ++i, ++ii)
    matBasisIntegration[ii] = dsLagBasisIntegration[i];

  for (unsigned long i = 0; i < dtLagBasisIntegration.size(); ++i, ++ii)
    matBasisIntegration[ii] = dtLagBasisIntegration[i];

  /* The transpose of the matrix containing the derivative information is needed
     for the volume integral. Create this info in such a way that the volume
     residual can be computed with one matrix multiplication when the fluxes
     are known. */
  unsigned short nDim;
  if (!dtLagBasisIntegration.empty())
    nDim = 3;
  else if (!dsLagBasisIntegration.empty())
    nDim = 2;
  else
    nDim = 1;

  ii = 0;
  for (unsigned short j = 0; j < nDOFs; ++j) {
    for (unsigned short i = 0; i < nIntegration; ++i) {
      for (unsigned short iDim = 0; iDim < nDim; ++iDim, ++ii) {
        const unsigned int ind = (iDim + 1) * nDOFs * nIntegration + i * nDOFs + j;
        matDerBasisIntTrans[ii] = matBasisIntegration[ind];
      }
    }
  }

  /*--------------------------------------------------------------------------*/
  /*--- Create the data of the derivatives of the basis functions in the   ---*/
  /*--- solution DOFs of the element.                                      ---*/
  /*--------------------------------------------------------------------------*/

  /* Store the location of the solution DOFs, depending whether or not these
     locations were specified when calling the constructor. */
  const vector<su2double> rSolDOFs = rLocSolDOFs ? *rLocSolDOFs : rDOFs;
  const vector<su2double> sSolDOFs = sLocSolDOFs ? *sLocSolDOFs : sDOFs;
  const vector<su2double> tSolDOFs = tLocSolDOFs ? *tLocSolDOFs : tDOFs;

  /* Create the basis functions and the matrix to store the derivatives in
     the basis functions in the solution DOFs. */
  CreateBasisFunctionsAndMatrixDerivatives(rSolDOFs, sSolDOFs, tSolDOFs, matVandermondeInv, lagBasisSolDOFs,
                                           matDerBasisSolDOFs);

  /*--------------------------------------------------------------------------*/
  /*--- Create the data of the derivatives of the basis functions in the   ---*/
  /*--- owned DOFs of the element. This data is needed for the computation ---*/
  /*--- of the derivatives of the metric terms in the integration points   ---*/
  /*--- when the grid DOFs do not coincide with the solution DOFs.         ---*/
  /*--------------------------------------------------------------------------*/

  vector<su2double> dummyLagBasis;
  vector<su2double> dummyMatVandermondeInv;
  CreateBasisFunctionsAndMatrixDerivatives(rDOFs, sDOFs, tDOFs, dummyMatVandermondeInv, dummyLagBasis,
                                           matDerBasisOwnDOFs);

  /*--------------------------------------------------------------------------*/
  /*--- Create the data of the second derivatives of the basis functions   ---*/
  /*--- in the integration points of the element. This data is needed for  ---*/
  /*--- the ADER preconditioning step when the Navier-Stokes equations are ---*/
  /*--- solved and can be obtained by multiplying matBasisIntegration and  ---*/
  /*--- matDerBasisSolDOFs.                                                ---*/
  /*--------------------------------------------------------------------------*/

  /* Definition of the object that is used to carry out the BLAS calls. */
  CBlasStructure blasFunctions;

  /* Easier storage of the offset between the derivatives for matBasisIntegration
     and matDerBasisSolDOFs. */
  const unsigned long offsetDerInt = nDOFs * nIntegration;
  const unsigned long offsetDerDOFs = nDOFs * nDOFs;

  /* Determine the size of the vector to store the second derivatives of the
     basis functions in the integration points. */
  const unsigned long sizeMat2ndDer = offsetDerInt * nDim * (nDim + 1) / 2;
  mat2ndDerBasisInt.resize(sizeMat2ndDer);

  /* Loop over the number of dimensions to carry out the matrix multiplications. */
  const su2double* matDerBasisInt = matBasisIntegration.data();
  su2double* mat2ndDerBasisIntPoint = mat2ndDerBasisInt.data();
  for (unsigned short jDim = 0; jDim < nDim; ++jDim) {
    /* Update the pointer matDerBasisInt to point to the correct derivative.
       As in matBasisIntegration also the interpolation data is stored, this
       update must happen before the matrix multiplication. */
    matDerBasisInt = matDerBasisInt + offsetDerInt;

    /* Second loop over the dimensions. Note that this loop only runs until
       iDim, because the Hessian is symmetric. */
    for (unsigned short iDim = 0; iDim <= jDim; iDim++) {
      /* Carry out the matrix multiplication. The last argument is NULL, such
         that this gemm call is ignored in the profiling. Replace by config if
         if should be included. */
      blasFunctions.gemm(nIntegration, nDOFs, nDOFs, matDerBasisInt, matDerBasisSolDOFs.data() + iDim * offsetDerDOFs,
                         mat2ndDerBasisIntPoint, nullptr);

      /* Update the pointer to the position where the next second
         derivative will be stored. */
      mat2ndDerBasisIntPoint = mat2ndDerBasisIntPoint + offsetDerInt;
    }
  }
}

void CFEMStandardElement::BasisFunctionsInPoint(const su2double* parCoor, vector<su2double>& lagBasis) {
  /* Determine the number of parametric dimensions, depending on the
     element type. */
  unsigned short nDimPar = 0;
  switch (VTK_Type) {
    case LINE:
      nDimPar = 1;
      break;
    case TRIANGLE:
    case QUADRILATERAL:
      nDimPar = 2;
      break;
    case TETRAHEDRON:
    case PYRAMID:
    case PRISM:
    case HEXAHEDRON:
      nDimPar = 3;
      break;
  }

  /* Allocate the memory for the help vectors for computing the Vandermonde
     matrix. */
  vector<vector<su2double> > rPoints(nDimPar, vector<su2double>(1));
  vector<su2double> V(nDOFs);

  /* Copy the parametric coordinates in rPoints, such that the functions to
     compute the Vandermonde matrices can be used. */
  for (unsigned long i = 0; i < rPoints.size(); ++i) rPoints[i][0] = parCoor[i];

  /* Determine the element type and call the appropriate function to compute
     the Vandermonde matrix. */
  switch (VTK_Type) {
    case LINE:
      Vandermonde1D(nDOFs, rPoints[0], V);
      break;

    case TRIANGLE:
      Vandermonde2D_Triangle(nPoly, nDOFs, rPoints[0], rPoints[1], V);
      break;

    case QUADRILATERAL:
      Vandermonde2D_Quadrilateral(nPoly, nDOFs, rPoints[0], rPoints[1], V);
      break;

    case TETRAHEDRON:
      Vandermonde3D_Tetrahedron(nPoly, nDOFs, rPoints[0], rPoints[1], rPoints[2], V);
      break;

    case PYRAMID:
      Vandermonde3D_Pyramid(nPoly, nDOFs, rPoints[0], rPoints[1], rPoints[2], V);
      break;

    case PRISM:
      Vandermonde3D_Prism(nPoly, nDOFs, rPoints[0], rPoints[1], rPoints[2], V);
      break;

    case HEXAHEDRON:
      Vandermonde3D_Hexahedron(nPoly, nDOFs, rPoints[0], rPoints[1], rPoints[2], V);
      break;
  }

  /* Carry out the matrix multiplication to obtain the values of the Lagrangian
     basis functions in the given parametric coordinate. */
  MatMulRowMajor(nDOFs, 1, V, matVandermondeInv, lagBasis);
}

void CFEMStandardElement::BasisFunctionsAndDerivativesInPoint(const su2double* parCoor, vector<su2double>& lagBasis,
                                                              vector<vector<su2double> >& dLagBasis) {
  /* Allocate the memory for the help vectors for computing the Vandermonde
     matrices and its derivatives. */
  vector<vector<su2double> > rPoints(dLagBasis.size(), vector<su2double>(1));
  vector<vector<su2double> > VDr(dLagBasis.size(), vector<su2double>(nDOFs));

  vector<su2double> V(nDOFs);

  /* Copy the parametric coordinates in rPoints, such that the functions to
     compute the Vandermonde matrices can be used. */
  for (unsigned long i = 0; i < rPoints.size(); ++i) rPoints[i][0] = parCoor[i];

  /* Determine the element type and call the appropriate function to compute
     the Vandermonde matrix and the derivative of the Vandermonde matrix. */
  switch (VTK_Type) {
    case LINE:
      Vandermonde1D(nDOFs, rPoints[0], V);
      GradVandermonde1D(nDOFs, rPoints[0], VDr[0]);
      break;

    case TRIANGLE:
      Vandermonde2D_Triangle(nPoly, nDOFs, rPoints[0], rPoints[1], V);
      GradVandermonde2D_Triangle(nPoly, nDOFs, rPoints[0], rPoints[1], VDr[0], VDr[1]);
      break;

    case QUADRILATERAL:
      Vandermonde2D_Quadrilateral(nPoly, nDOFs, rPoints[0], rPoints[1], V);
      GradVandermonde2D_Quadrilateral(nPoly, nDOFs, rPoints[0], rPoints[1], VDr[0], VDr[1]);
      break;

    case TETRAHEDRON:
      Vandermonde3D_Tetrahedron(nPoly, nDOFs, rPoints[0], rPoints[1], rPoints[2], V);
      GradVandermonde3D_Tetrahedron(nPoly, nDOFs, rPoints[0], rPoints[1], rPoints[2], VDr[0], VDr[1], VDr[2]);
      break;

    case PYRAMID:
      Vandermonde3D_Pyramid(nPoly, nDOFs, rPoints[0], rPoints[1], rPoints[2], V);
      GradVandermonde3D_Pyramid(nPoly, nDOFs, rPoints[0], rPoints[1], rPoints[2], VDr[0], VDr[1], VDr[2]);
      break;

    case PRISM:
      Vandermonde3D_Prism(nPoly, nDOFs, rPoints[0], rPoints[1], rPoints[2], V);
      GradVandermonde3D_Prism(nPoly, nDOFs, rPoints[0], rPoints[1], rPoints[2], VDr[0], VDr[1], VDr[2]);
      break;

    case HEXAHEDRON:
      Vandermonde3D_Hexahedron(nPoly, nDOFs, rPoints[0], rPoints[1], rPoints[2], V);
      GradVandermonde3D_Hexahedron(nPoly, nDOFs, rPoints[0], rPoints[1], rPoints[2], VDr[0], VDr[1], VDr[2]);
      break;
  }

  /* Carry out the matrix multiplication to obtain the values of the Lagrangian
     basis functions and its derivatives in the given parametric coordinate. */
  MatMulRowMajor(nDOFs, 1, V, matVandermondeInv, lagBasis);

  for (unsigned long i = 0; i < dLagBasis.size(); ++i)
    MatMulRowMajor(nDOFs, 1, VDr[i], matVandermondeInv, dLagBasis[i]);
}

bool CFEMStandardElement::SameStandardElement(unsigned short val_VTK_Type, unsigned short val_nPoly,
                                              bool val_constJac) {
  if (val_VTK_Type != VTK_Type) return false;
  if (val_nPoly != nPoly) return false;
  if (val_constJac != constJacobian) return false;

  return true;
}

/*----------------------------------------------------------------------------------*/
/*           Private member functions of CFEMStandardElement.                       */
/*----------------------------------------------------------------------------------*/

void CFEMStandardElement::Copy(const CFEMStandardElement& other) {
  CFEMStandardElementBase::Copy(other);

  nPoly = other.nPoly;
  nDOFs = other.nDOFs;

  VTK_Type1 = other.VTK_Type1;
  VTK_Type2 = other.VTK_Type2;

  rDOFs = other.rDOFs;
  sDOFs = other.sDOFs;
  tDOFs = other.tDOFs;

  lagBasisIntegration = other.lagBasisIntegration;
  lagBasisIntegrationTrans = other.lagBasisIntegrationTrans;
  lagBasisSolDOFs = other.lagBasisSolDOFs;

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

  matVandermondeInv = other.matVandermondeInv;
  matBasisIntegration = other.matBasisIntegration;
  matDerBasisIntTrans = other.matDerBasisIntTrans;
  matDerBasisSolDOFs = other.matDerBasisSolDOFs;
  matDerBasisOwnDOFs = other.matDerBasisOwnDOFs;
  mat2ndDerBasisInt = other.mat2ndDerBasisInt;
}

void CFEMStandardElement::CreateBasisFunctionsAndMatrixDerivatives(
    const vector<su2double>& rLoc, const vector<su2double>& sLoc, const vector<su2double>& tLoc,
    vector<su2double>& matVandermondeInv, vector<su2double>& lagBasis, vector<su2double>& matDerBasis) {
  /* Define the variables, such that the general functions to compute
     the gradients of the basis functions can be used. Note that some
     of these variables are dummy variables. */
  unsigned short nDOFsDummy;
  vector<su2double> rDOFsDummy, sDOFsDummy, tDOFsDummy;
  vector<su2double> drLagBasisLoc, dsLagBasisLoc, dtLagBasisLoc;

  /*--- Determine the element type and compute the basis functions and the
        gradients of the basis functions in the given DOFs. ---*/
  switch (VTK_Type) {
    case LINE:
      LagrangianBasisFunctionAndDerivativesLine(nPoly, rLoc, nDOFsDummy, rDOFsDummy, matVandermondeInv, lagBasis,
                                                drLagBasisLoc);
      break;

    case TRIANGLE:
      LagrangianBasisFunctionAndDerivativesTriangle(nPoly, rLoc, sLoc, nDOFsDummy, rDOFsDummy, sDOFsDummy,
                                                    matVandermondeInv, lagBasis, drLagBasisLoc, dsLagBasisLoc);
      break;

    case QUADRILATERAL:
      LagrangianBasisFunctionAndDerivativesQuadrilateral(nPoly, rLoc, sLoc, nDOFsDummy, rDOFsDummy, sDOFsDummy,
                                                         matVandermondeInv, lagBasis, drLagBasisLoc, dsLagBasisLoc);
      break;

    case TETRAHEDRON:
      LagrangianBasisFunctionAndDerivativesTetrahedron(nPoly, rLoc, sLoc, tLoc, nDOFsDummy, rDOFsDummy, sDOFsDummy,
                                                       tDOFsDummy, matVandermondeInv, lagBasis, drLagBasisLoc,
                                                       dsLagBasisLoc, dtLagBasisLoc);
      break;

    case PYRAMID:
      LagrangianBasisFunctionAndDerivativesPyramid(nPoly, rLoc, sLoc, tLoc, nDOFsDummy, rDOFsDummy, sDOFsDummy,
                                                   tDOFsDummy, matVandermondeInv, lagBasis, drLagBasisLoc,
                                                   dsLagBasisLoc, dtLagBasisLoc);
      break;

    case PRISM:
      LagrangianBasisFunctionAndDerivativesPrism(nPoly, rLoc, sLoc, tLoc, nDOFsDummy, rDOFsDummy, sDOFsDummy,
                                                 tDOFsDummy, matVandermondeInv, lagBasis, drLagBasisLoc, dsLagBasisLoc,
                                                 dtLagBasisLoc);
      break;

    case HEXAHEDRON:
      LagrangianBasisFunctionAndDerivativesHexahedron(nPoly, rLoc, sLoc, tLoc, nDOFsDummy, rDOFsDummy, sDOFsDummy,
                                                      tDOFsDummy, matVandermondeInv, lagBasis, drLagBasisLoc,
                                                      dsLagBasisLoc, dtLagBasisLoc);
      break;
  }

  /*--- Check the sum of the Lagrangian basis functions and the
        derivatives of the Lagrangian basis functions. ---*/
  const unsigned short nLoc = rLoc.size();

  CheckSumLagrangianBasisFunctions(nLoc, nDOFs, lagBasis);

  if (!drLagBasisLoc.empty()) CheckSumDerivativesLagrangianBasisFunctions(nLoc, nDOFs, drLagBasisLoc);
  if (!dsLagBasisLoc.empty()) CheckSumDerivativesLagrangianBasisFunctions(nLoc, nDOFs, dsLagBasisLoc);
  if (!dtLagBasisLoc.empty()) CheckSumDerivativesLagrangianBasisFunctions(nLoc, nDOFs, dtLagBasisLoc);

  /*--- For efficiency reasons drLagBasisLoc, dsLagBasisLoc and
        dtLagBasisLoc are stored in matDerBasis. ---*/
  const unsigned long sizeDerMat = drLagBasisLoc.size() + dsLagBasisLoc.size() + dtLagBasisLoc.size();

  matDerBasis.resize(sizeDerMat);

  unsigned int ii = 0;
  for (unsigned long i = 0; i < drLagBasisLoc.size(); ++i, ++ii) matDerBasis[ii] = drLagBasisLoc[i];

  for (unsigned long i = 0; i < dsLagBasisLoc.size(); ++i, ++ii) matDerBasis[ii] = dsLagBasisLoc[i];

  for (unsigned long i = 0; i < dtLagBasisLoc.size(); ++i, ++ii) matDerBasis[ii] = dtLagBasisLoc[i];
}

void CFEMStandardElement::DataStandardLine() {
  /*--- Determine the Lagrangian basis functions and its derivatives
        in the integration points. ---*/
  vector<su2double> matVandermondeInvDummy;
  LagrangianBasisFunctionAndDerivativesLine(nPoly, rIntegration, nDOFs, rDOFs, matVandermondeInvDummy,
                                            lagBasisIntegration, drLagBasisIntegration);

  /*--- Determine the local connectivity of the two "faces" of the line element.
        For a line element the faces are just points.   ---*/
  connFace0.reserve(1);
  connFace0.push_back(0);
  connFace1.reserve(1);
  connFace1.push_back(nPoly);

  /*--- Determine the local subconnectivity used for plotting purposes. ---*/
  SubConnForPlottingLine(nPoly, subConn1ForPlotting);

  /*--- Set the VTK_type(s) for this sub element. ---*/
  VTK_Type1 = LINE;
  VTK_Type2 = NONE;
}

void CFEMStandardElement::DataStandardTriangle() {
  /*--- Determine the Lagrangian basis functions and its derivatives
        in the integration points. ---*/
  vector<su2double> matVandermondeInvDummy;
  LagrangianBasisFunctionAndDerivativesTriangle(nPoly, rIntegration, sIntegration, nDOFs, rDOFs, sDOFs,
                                                matVandermondeInvDummy, lagBasisIntegration, drLagBasisIntegration,
                                                dsLagBasisIntegration);

  /*--- Determine the local connectivity of the three "faces" of the triangle.
        For a triangular element the faces are just lines. Make sure that the
        element is to the left of the face. ---*/
  connFace0.reserve(nPoly + 1);
  connFace1.reserve(nPoly + 1);
  connFace2.reserve(nPoly + 1);

  for (signed short i = 0; i <= nPoly; ++i) connFace0.push_back(i);
  for (signed short i = 0; i <= nPoly; ++i) connFace1.push_back((i + 1) * (nPoly + 1) - i * (i + 1) / 2 - 1);
  for (signed short i = nPoly; i >= 0; --i) connFace2.push_back(i * (nPoly + 1) - i * (i - 1) / 2);

  /*--- Determine the local subconnectivity of the triangular element used for
        plotting purposes. ---*/
  SubConnForPlottingTriangle(nPoly, subConn1ForPlotting);

  /*--- Set the VTK_type(s) for this sub element. ---*/
  VTK_Type1 = TRIANGLE;
  VTK_Type2 = NONE;
}

void CFEMStandardElement::DataStandardQuadrilateral() {
  /*--- Determine the Lagrangian basis functions and its derivatives
        in the integration points. ---*/
  vector<su2double> matVandermondeInvDummy;
  LagrangianBasisFunctionAndDerivativesQuadrilateral(nPoly, rIntegration, sIntegration, nDOFs, rDOFs, sDOFs,
                                                     matVandermondeInvDummy, lagBasisIntegration, drLagBasisIntegration,
                                                     dsLagBasisIntegration);

  /*--- Determine the local connectivity of the four "faces" of the quad element.
        For a quad element the faces are just lines. Make sure that the element
        is to the left of the face. ---*/
  connFace0.reserve(nPoly + 1);
  connFace1.reserve(nPoly + 1);
  connFace2.reserve(nPoly + 1);
  connFace3.reserve(nPoly + 1);

  unsigned short n0 = 0, n1 = nPoly, n2 = nDOFs - 1, n3 = nPoly * (nPoly + 1);

  for (signed short i = n0; i <= n1; ++i) connFace0.push_back(i);
  for (signed short i = n1; i <= n2; i += (nPoly + 1)) connFace1.push_back(i);
  for (signed short i = n2; i >= n3; --i) connFace2.push_back(i);
  for (signed short i = n3; i >= n0; i -= (nPoly + 1)) connFace3.push_back(i);

  /*--- Determine the local subconnectivity of the quadrilateral element used for
        plotting purposes. ---*/
  SubConnForPlottingQuadrilateral(nPoly, subConn1ForPlotting);

  /*--- Set the VTK_type(s) for this sub element. ---*/
  VTK_Type1 = QUADRILATERAL;
  VTK_Type2 = NONE;
}

void CFEMStandardElement::DataStandardTetrahedron() {
  /*--- Determine the Lagrangian basis functions and its derivatives
        in the integration points. ---*/
  vector<su2double> matVandermondeInvDummy;
  LagrangianBasisFunctionAndDerivativesTetrahedron(nPoly, rIntegration, sIntegration, tIntegration, nDOFs, rDOFs, sDOFs,
                                                   tDOFs, matVandermondeInvDummy, lagBasisIntegration,
                                                   drLagBasisIntegration, dsLagBasisIntegration, dtLagBasisIntegration);

  /*--- Determine the local connectivity of the four faces of the tetrahedron.
        For a tetrahedron the faces are triangles. ---*/
  unsigned short nDOFsTriangle = (nPoly + 1) * (nPoly + 2) / 2;
  connFace0.reserve(nDOFsTriangle);
  connFace1.reserve(nDOFsTriangle);
  connFace2.reserve(nDOFsTriangle);
  connFace3.reserve(nDOFsTriangle);

  unsigned int ii = 0;
  for (unsigned short k = 0; k <= nPoly; ++k) {
    unsigned short uppBoundJ = nPoly - k;
    for (unsigned short j = 0; j <= uppBoundJ; ++j) {
      unsigned short uppBoundI = nPoly - k - j;
      for (unsigned short i = 0; i <= uppBoundI; ++i, ++ii) {
        if (k == 0) connFace0.push_back(ii);
        if (j == 0) connFace1.push_back(ii);
        if (i == 0) connFace2.push_back(ii);
        if ((i + j + k) == nPoly) connFace3.push_back(ii);
      }
    }
  }

  /*--- Make sure that the element is to the left of the faces. ---*/
  unsigned short n0 = 0;
  unsigned short n1 = nPoly;
  unsigned short n2 = nDOFsTriangle - 1;
  unsigned short n3 = nDOFs - 1;

  ChangeDirectionTriangleConn(connFace0, n0, n1, n2);
  ChangeDirectionTriangleConn(connFace1, n0, n3, n1);
  ChangeDirectionTriangleConn(connFace2, n0, n2, n3);
  ChangeDirectionTriangleConn(connFace3, n1, n3, n2);

  /*--- Determine the local subconnectivity of the tetrahedron used for
        plotting purposes. The high order tetrahedron is split in several
        linear subtetrahedra.           ---*/
  SubConnTetrahedron();

  /*--- Set the VTK_type(s) for this sub element. ---*/
  VTK_Type1 = TETRAHEDRON;
  VTK_Type2 = NONE;
}

void CFEMStandardElement::DataStandardPyramid() {
  /*--- Determine the Lagrangian basis functions and its derivatives
        in the integration points. ---*/
  vector<su2double> matVandermondeInvDummy;
  LagrangianBasisFunctionAndDerivativesPyramid(nPoly, rIntegration, sIntegration, tIntegration, nDOFs, rDOFs, sDOFs,
                                               tDOFs, matVandermondeInvDummy, lagBasisIntegration,
                                               drLagBasisIntegration, dsLagBasisIntegration, dtLagBasisIntegration);

  /*--- Determine the local connectivity of the five faces of the pyramid.
        For a pyramid there are four triangular faces and one quadrilateral face. ---*/
  unsigned short nDOFsQuad = (nPoly + 1) * (nPoly + 1);
  unsigned short nDOFsTriangle = (nPoly + 1) * (nPoly + 2) / 2;

  connFace0.reserve(nDOFsQuad);
  connFace1.reserve(nDOFsTriangle);
  connFace2.reserve(nDOFsTriangle);
  connFace3.reserve(nDOFsTriangle);
  connFace4.reserve(nDOFsTriangle);

  unsigned short mPoly = nPoly;
  unsigned int ii = 0;
  for (unsigned short k = 0; k <= nPoly; ++k, --mPoly) {
    for (unsigned short j = 0; j <= mPoly; ++j) {
      for (unsigned short i = 0; i <= mPoly; ++i, ++ii) {
        if (k == 0) connFace0.push_back(ii);
        if (j == 0) connFace1.push_back(ii);
        if (j == mPoly) connFace2.push_back(ii);
        if (i == 0) connFace3.push_back(ii);
        if (i == mPoly) connFace4.push_back(ii);
      }
    }
  }

  /*--- Make sure that the element is to the left of the faces. ---*/
  unsigned short n0 = 0;
  unsigned short n1 = nPoly;
  unsigned short n2 = nDOFsQuad - 1;
  unsigned short n3 = n2 - nPoly;
  unsigned short n4 = nDOFs - 1;

  ChangeDirectionQuadConn(connFace0, n0, n1, n2, n3);
  ChangeDirectionTriangleConn(connFace1, n0, n4, n1);
  ChangeDirectionTriangleConn(connFace2, n3, n2, n4);
  ChangeDirectionTriangleConn(connFace3, n0, n3, n4);
  ChangeDirectionTriangleConn(connFace4, n1, n4, n2);

  /*--- Determine the local subconnectivity of the pyramid used for
        plotting purposes. The high order pyramid is split in several
        linear subpyramids and subtetrahedra, i.e. two element types. ---*/
  SubConnPyramid();

  /*--- Set the VTK_type(s) for this sub element. ---*/
  VTK_Type1 = PYRAMID;
  VTK_Type2 = TETRAHEDRON;
}

void CFEMStandardElement::DataStandardPrism() {
  /*--- Determine the Lagrangian basis functions and its derivatives
        in the integration points. ---*/
  vector<su2double> matVandermondeInvDummy;
  LagrangianBasisFunctionAndDerivativesPrism(nPoly, rIntegration, sIntegration, tIntegration, nDOFs, rDOFs, sDOFs,
                                             tDOFs, matVandermondeInvDummy, lagBasisIntegration, drLagBasisIntegration,
                                             dsLagBasisIntegration, dtLagBasisIntegration);

  /*--- Determine the local connectivity of the five faces of the prism.
        For a prism there are two triangular faces and three quadrilateral faces. ---*/
  unsigned short nDOFsQuad = (nPoly + 1) * (nPoly + 1);
  unsigned short nDOFsTriangle = (nPoly + 1) * (nPoly + 2) / 2;

  connFace0.reserve(nDOFsTriangle);
  connFace1.reserve(nDOFsTriangle);
  connFace2.reserve(nDOFsQuad);
  connFace3.reserve(nDOFsQuad);
  connFace4.reserve(nDOFsQuad);

  unsigned int ii = 0;
  for (unsigned short k = 0; k <= nPoly; ++k) {
    for (unsigned short j = 0; j <= nPoly; ++j) {
      unsigned short uppBoundI = nPoly - j;
      for (unsigned short i = 0; i <= uppBoundI; ++i, ++ii) {
        if (k == 0) connFace0.push_back(ii);
        if (k == nPoly) connFace1.push_back(ii);
        if (j == 0) connFace2.push_back(ii);
        if (i == 0) connFace3.push_back(ii);
        if ((i + j) == nPoly) connFace4.push_back(ii);
      }
    }
  }

  /*--- Make sure that the element is to the left of the faces. ---*/
  unsigned short n0 = 0;
  unsigned short n1 = nPoly;
  unsigned short n2 = nDOFsTriangle - 1;
  unsigned short n3 = n0 + nDOFsTriangle * nPoly;
  unsigned short n4 = n1 + nDOFsTriangle * nPoly;
  unsigned short n5 = n2 + nDOFsTriangle * nPoly;

  ChangeDirectionTriangleConn(connFace0, n0, n1, n2);
  ChangeDirectionTriangleConn(connFace1, n3, n5, n4);
  ChangeDirectionQuadConn(connFace2, n0, n3, n4, n1);
  ChangeDirectionQuadConn(connFace3, n0, n2, n5, n3);
  ChangeDirectionQuadConn(connFace4, n1, n4, n5, n2);

  /*--- Determine the local subconnectivity of the prism used for
        plotting purposes. The high order prism is split in several
        linear subprisms.                      ---*/
  SubConnPrism();

  /*--- Set the VTK_type(s) for this sub element. ---*/
  VTK_Type1 = PRISM;
  VTK_Type2 = NONE;
}

void CFEMStandardElement::DataStandardHexahedron() {
  /*--- Determine the Lagrangian basis functions and its derivatives
        in the integration points. ---*/
  vector<su2double> matVandermondeInvDummy;
  LagrangianBasisFunctionAndDerivativesHexahedron(nPoly, rIntegration, sIntegration, tIntegration, nDOFs, rDOFs, sDOFs,
                                                  tDOFs, matVandermondeInvDummy, lagBasisIntegration,
                                                  drLagBasisIntegration, dsLagBasisIntegration, dtLagBasisIntegration);

  /*--- Determine the local connectivity of the six faces of the hexahedron.
        For a hexahedron the faces are all quadrilateral faces. ---*/
  unsigned short nDOFsQuad = (nPoly + 1) * (nPoly + 1);

  connFace0.reserve(nDOFsQuad);
  connFace1.reserve(nDOFsQuad);
  connFace2.reserve(nDOFsQuad);
  connFace3.reserve(nDOFsQuad);
  connFace4.reserve(nDOFsQuad);
  connFace5.reserve(nDOFsQuad);

  unsigned int ii = 0;
  for (unsigned short k = 0; k <= nPoly; ++k) {
    for (unsigned short j = 0; j <= nPoly; ++j) {
      for (unsigned short i = 0; i <= nPoly; ++i, ++ii) {
        if (k == 0) connFace0.push_back(ii);
        if (k == nPoly) connFace1.push_back(ii);
        if (j == 0) connFace2.push_back(ii);
        if (j == nPoly) connFace3.push_back(ii);
        if (i == 0) connFace4.push_back(ii);
        if (i == nPoly) connFace5.push_back(ii);
      }
    }
  }

  /*--- Make sure that the element is to the left of the faces. ---*/
  unsigned short n0 = 0;
  unsigned short n1 = nPoly;
  unsigned short n2 = nDOFsQuad - 1;
  unsigned short n3 = n2 - nPoly;
  unsigned short n4 = n0 + nDOFsQuad * nPoly;
  unsigned short n5 = n1 + nDOFsQuad * nPoly;
  unsigned short n6 = n2 + nDOFsQuad * nPoly;
  unsigned short n7 = n3 + nDOFsQuad * nPoly;

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

  /*--- Set the VTK_type(s) for this sub element. ---*/
  VTK_Type1 = HEXAHEDRON;
  VTK_Type2 = NONE;
}

void CFEMStandardElement::SubConnTetrahedron() {
  /*--- Initialize the number of DOFs for the current edges to the number of
        DOFs of the edges present in the tetrahedron. Also initialize the
        current k offset to zero.    ---*/
  unsigned short nDOFsCurrentEdges = nPoly + 1;
  unsigned short offCurrentK = 0;

  /*--- Loop in the k-direction of the tetrahedron, which is along the edge
        from the first vertex to the last vertex of the tet.    ---*/
  for (unsigned short k = 0; k < nPoly; ++k) {
    /*--- Determine the offset for the next k. ---*/
    unsigned short offNextK = offCurrentK + nDOFsCurrentEdges * (nDOFsCurrentEdges + 1) / 2;

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

    for (unsigned short i = 0; i < (nDOFsCurrentEdges - 2); ++i) {
      /*--- Determine the lowest j-index on the current i-line that contributes
            to the subprism. Convert that index to the local vertex number in
            the tetrahedron.     ---*/
      unsigned short j = nDOFsCurrentEdges - 2 - i;
      n0 = j * nDOFsCurrentEdges + i - j * (j - 1) / 2 + offCurrentK;

      /*--- Increment the j index to obtain the n1 vertex of the subprism. ---*/
      ++j;
      n1 = j * nDOFsCurrentEdges + i - j * (j - 1) / 2 + offCurrentK;

      /*--- The n2 vertex of the prism is located on the next k-level.
            The i-index remains the same, but the j index must be adapted, because
            the number of DOFs on the edges on the next k-level is one less. ---*/
      j = nDOFsCurrentEdges - 2 - i;
      n2 = j * (nDOFsCurrentEdges - 1) + i - j * (j - 1) / 2 + offNextK;

      /*--- The n3 vertex is part of the upper triangle of the prism. Hence the
            i-index must be incremented by one, stored in ii. The j-index must
            be computed accordingly and converted to the 1D numbering. ---*/
      unsigned short ii = i + 1;
      j = nDOFsCurrentEdges - 2 - ii;
      n3 = j * nDOFsCurrentEdges + ii - j * (j - 1) / 2 + offCurrentK;

      /*--- Increment the j index to obtain the n4 vertex of the subprism. ---*/
      ++j;
      unsigned short n4 = j * nDOFsCurrentEdges + ii - j * (j - 1) / 2 + offCurrentK;

      /*--- The n5 vertex of the prism is located on the next k-level.
            The i-index remains the same, but the j index must be adapted, because
            the number of DOFs on the edges on the next k-level is one less. ---*/
      j = nDOFsCurrentEdges - 2 - ii;
      unsigned short n5 = j * (nDOFsCurrentEdges - 1) + ii - j * (j - 1) / 2 + offNextK;

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

    for (unsigned short i = 0; i < (nDOFsCurrentEdges - 2); ++i) {
      /*--- Define the variables n4 to n7, because these indices will be used
            again for the prism treated after the hexahedra. ---*/
      unsigned short n4, n5, n6, n7;

      /*--- Initialize n3, n2, n7 and n6 to the quad on the line j = 0. ---*/
      n3 = offCurrentK + i;
      n2 = n3 + 1;
      n7 = offNextK + i;
      n6 = n7 + 1;

      /*--- Loop in the j-direction for the hexahedra on this i-row. ---*/
      for (unsigned short j = 0; j < (nDOFsCurrentEdges - 3 - i); ++j) {
        /*--- Set the values of n0, n1, n4 and n5 from the previous hex. ---*/
        n0 = n3;
        n1 = n2;
        n4 = n7;
        n5 = n6;

        /*--- Nodes n3 and n7 are the j-neighbors of n0 and n4 respectively.
              Convert these (i,j,k) indices to the 1D index again. ---*/
        n3 = (j + 1) * nDOFsCurrentEdges + i - j * (j + 1) / 2 + offCurrentK;
        n7 = (j + 1) * (nDOFsCurrentEdges - 1) + i - j * (j + 1) / 2 + offNextK;

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
      unsigned short j = nDOFsCurrentEdges - 2 - i;

      n2 = j * nDOFsCurrentEdges + i - j * (j - 1) / 2 + offCurrentK;
      n5 = j * (nDOFsCurrentEdges - 1) + i - j * (j - 1) / 2 + offNextK;

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

void CFEMStandardElement::SubConnPyramid() {
  /*--- Initialize the number of DOFs for the current edges to the number of
        DOFs of the edges on the base of the pyramid. Also initialize the
        current k offset to zero.     ---*/
  unsigned short nDOFsCurrentEdges = nPoly + 1;
  unsigned short offCurrentK = 0;

  /*--- Loop in the k-direction of the pyramid. ---*/
  for (unsigned short k = 0; k < nPoly; ++k) {
    /*--------------------------------------------------------------------------
       Sub-pyramids in the same direction as the original pyramid.
      ------------------------------------------------------------------------*/

    /*--- Determine the index of the first vertex of the quadrilateral of the
          next k value.    ---*/
    unsigned short kk = offCurrentK + nDOFsCurrentEdges * nDOFsCurrentEdges;

    // Loop in j-direction of the current quad.
    for (unsigned short j = 0; j < (nDOFsCurrentEdges - 1); ++j) {
      /*--- Index of the first vertex along the j-row of the current quad. ---*/
      unsigned short jj = offCurrentK + j * nDOFsCurrentEdges;

      /*--- Loop in i-direction of the current quad. ---*/
      for (unsigned short i = 0; i < (nDOFsCurrentEdges - 1); ++i) {
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
    kk = offCurrentK + nDOFsCurrentEdges * nDOFsCurrentEdges;

    /*--- Loop in j-direction of the current quad. Note that the starting index
          of this loop is 1.                       ---*/
    for (unsigned short j = 1; j < (nDOFsCurrentEdges - 1); ++j) {
      /*--- Index of the first vertex along the j-row of the current quad. ---*/
      unsigned short jj = offCurrentK + j * nDOFsCurrentEdges;

      /*--- Loop in the i-direction of this quad. Again the starting index is 1. ---*/
      for (unsigned short i = 1; i < (nDOFsCurrentEdges - 1); ++i) {
        /*--- Determine the local indices of the corners of the quadrilateral
              of this subpyramid as well as the top of the subpyramid.  ---*/
        unsigned short n0 = kk + i - 1;
        unsigned short n1 = n0 + 1;
        unsigned short n2 = n1 + nDOFsCurrentEdges - 1;
        unsigned short n3 = n0 + nDOFsCurrentEdges - 1;
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
    kk = offCurrentK + nDOFsCurrentEdges * nDOFsCurrentEdges;

    /*--- Loop in the i-direction of the current quad. Note that the starting
          index must be 1.                         ---*/
    for (unsigned short i = 1; i < (nDOFsCurrentEdges - 1); ++i) {
      /*--- Loop in the j-direction of the current quad. This loop starts at 0. ---*/
      for (unsigned short j = 0; j < (nDOFsCurrentEdges - 1); ++j) {
        /*--- Determine the local indices of the 4 corner points of this tet.
              Its connectivity is stored in subConn2ForPlotting, because in
              subConn1ForPlotting the pyramids are stored. ---*/
        unsigned short n0 = kk + i - 1 + j * (nDOFsCurrentEdges - 1);
        unsigned short n1 = n0 + 1;
        unsigned short n2 = offCurrentK + i + j * nDOFsCurrentEdges;
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
    for (unsigned short j = 1; j < (nDOFsCurrentEdges - 1); ++j) {
      /*--- Index of the first vertex along the j-row of the current quad. ---*/
      unsigned short jj = offCurrentK + j * nDOFsCurrentEdges;

      /*--- Loop in the i-direction of the current quad. This loop starts at 0. ---*/
      for (unsigned short i = 0; i < (nDOFsCurrentEdges - 1); ++i) {
        /*--- Determine the local indices of the 4 corner points of this tet
              and store its connectivity in subConn2ForPlotting.   ---*/
        unsigned short n0 = kk + i;
        unsigned short n1 = jj + i;
        unsigned short n2 = n0 + nDOFsCurrentEdges - 1;
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
    offCurrentK += nDOFsCurrentEdges * nDOFsCurrentEdges;
    --nDOFsCurrentEdges;
  }
}

void CFEMStandardElement::SubConnPrism() {
  /*--- Determine the number of DOFs for a triangle. This is the offset in
        k-direction, the structured direction of a prisms.    ---*/
  unsigned short nDOFTria = (nPoly + 1) * (nPoly + 2) / 2;

  /*--- Loop in k-direction, which is the structured direction of the prism. ---*/
  for (unsigned short k = 0; k < nPoly; ++k) {
    /*--- Initialize the counter jj to the ID of the first vertex for this
          k-value. jj contains the ID of the first vertex on the "0-2 edge"
          of the current bottom triangle.                 ---*/
    unsigned short jj = k * nDOFTria;

    /*--- Loop over subedges of the left boundary of the standard triangle. ---*/
    for (unsigned short j = 0; j < nPoly; ++j) {
      /*--- Check if the "down" elements must be written. ---*/
      if (j) {
        /*--- Determine the offset of the relevant DOF on the previous row. ---*/
        unsigned short kk = jj - (nPoly + 1 - j);

        /*--- Loop over the edges of this row. ---*/
        for (unsigned short i = 0; i < (nPoly - j); ++i) {
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
      for (unsigned short i = 0; i < (nPoly - j); ++i) {
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

void CFEMStandardElement::SubConnHexahedron() {
  /*--- Determine the nodal offset in j- and k-direction. ---*/
  unsigned short jOff = nPoly + 1;
  unsigned short kOff = jOff * jOff;

  /*--- Loop over the subelements in k-direction. ---*/
  for (unsigned short k = 0; k < nPoly; ++k) {
    /*--- Abbreviate the offset in k-direction used in the connectivity. ---*/
    unsigned short kk = k * kOff;

    /*--- Loop over the subelements in j-direction. ---*/
    for (unsigned short j = 0; j < nPoly; ++j) {
      /*--- Abbreviate the offset in j-direction used in the connectivity. ---*/
      unsigned short jj = j * jOff;

      /*--- Loop over the subelements in i-direction. ---*/
      for (unsigned short i = 0; i < nPoly; ++i) {
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

unsigned short CFEMStandardElement::GetNDOFsPerSubElem(unsigned short val_VTK_Type) const {
  /*--- Distinguish between the possible element types for a linear
   sub-element and set the nDOFs accordingly. ---*/
  unsigned short nDOFsSubElem = 0;  // To avoid a compiler warning.
  switch (val_VTK_Type) {
    case NONE:
      nDOFsSubElem = 0;
      break;
    case LINE:
      nDOFsSubElem = 2;
      break;
    case TRIANGLE:
      nDOFsSubElem = 3;
      break;
    case QUADRILATERAL:
      nDOFsSubElem = 4;
      break;
    case TETRAHEDRON:
      nDOFsSubElem = 4;
      break;
    case PYRAMID:
      nDOFsSubElem = 5;
      break;
    case PRISM:
      nDOFsSubElem = 6;
      break;
    case HEXAHEDRON:
      nDOFsSubElem = 8;
      break;
    default:
      ostringstream message;
      message << "Impossible FEM sub element type, " << val_VTK_Type << ", encountered.";
      SU2_MPI::Error(message.str(), CURRENT_FUNCTION);
  }

  /* Return the number of DOFs for a subface. */
  return nDOFsSubElem;
}

void CFEMStandardElement::ChangeDirectionQuadConn(vector<unsigned short>& connQuad, unsigned short vert0,
                                                  unsigned short vert1, unsigned short vert2,
                                                  unsigned short vert3) const {
  /*--- Determine the indices of the 4 corner vertices of the quad. ---*/
  unsigned short ind0 = 0;
  unsigned short ind1 = nPoly;
  unsigned short ind2 = (nPoly + 1) * (nPoly + 1) - 1;
  unsigned short ind3 = ind2 - nPoly;

  /*--- There exists a linear mapping from the indices of the numbering used in the
        connectivity of this face to the indices of the target numbering. This
        mapping is of the form ii = a + b*i + c*j and jj = d + e*i + f*j, where
        ii,jj are the indices of the target numbering and i,j the indices of the
        numbering used for this face. The values of the coefficients a,b,c,d,e,f
        depend on how the corner points coincide with each other. This is
        determined below. The bool verticesDontMatch is there to check if vertices
        do not match. This should not happen, but it is checked for security. ---*/

  signed short a = 0, b = 0, c = 0, d = 0, e = 0, f = 0;  // Initialization to avoid a compiler warning.
  bool verticesDontMatch = false;

  if (vert0 == connQuad[ind0]) {
    /*--- Vert0 coincides with the first vertex of the face connectivity.
          Set the coefficients a and d accordingly.   ---*/
    a = d = 0;
    if (vert2 != connQuad[ind2]) verticesDontMatch = true;

    /*--- Check the situation for the neighboring vertices. ---*/
    if (vert1 == connQuad[ind1]) {
      /*--- The vertex numbering is the same for both faces. ---*/
      if (vert3 != connQuad[ind3]) verticesDontMatch = true;

      b = f = 1;
      c = e = 0;
    } else if (vert1 == connQuad[ind3]) {
      /*--- The i and j numbering are swapped. ---*/
      if (vert3 != connQuad[ind1]) verticesDontMatch = true;

      b = f = 0;
      c = e = 1;
    } else {
      verticesDontMatch = true;
    }
  } else if (vert1 == connQuad[ind0]) {
    /*--- Vert1 coincides with the first vertex of the face connectivity.
          Set the coefficients a and d accordingly.  ---*/
    a = nPoly;
    d = 0;
    if (vert3 != connQuad[ind2]) verticesDontMatch = true;

    /*--- Check the situation for the neighboring vertices. ---*/
    if (vert0 == connQuad[ind1]) {
      /*--- The i-direction is negated while the j-direction coincides. ---*/
      if (vert2 != connQuad[ind3]) verticesDontMatch = true;

      b = -1;
      f = 1;
      c = e = 0;
    } else if (vert0 == connQuad[ind3]) {
      /*--- The j-direction of the current face corresponds with the negative
            i-direction of the target, while the i-direction coincides with
            the j-direction of the target.     ---*/
      if (vert2 != connQuad[ind1]) verticesDontMatch = true;

      b = f = 0;
      c = -1;
      e = 1;
    } else {
      verticesDontMatch = true;
    }
  } else if (vert2 == connQuad[ind0]) {
    /*--- Vert2 coincides with the first vertex of the face connectivity.
          Set the coefficients a and d accordingly.  ---*/
    a = d = nPoly;
    if (vert0 != connQuad[ind2]) verticesDontMatch = true;

    /*--- Check the situation for the neighboring vertices. ---*/
    if (vert1 == connQuad[ind3]) {
      /*--- Both the i and j-direction are negated. ---*/
      if (vert3 != connQuad[ind1]) verticesDontMatch = true;

      b = f = -1;
      c = e = 0;
    } else if (vert1 == connQuad[ind1]) {
      /*--- The i and j-direction are negated and swapped. ---*/
      if (vert3 != connQuad[ind3]) verticesDontMatch = true;

      b = f = 0;
      c = e = -1;
    } else {
      verticesDontMatch = true;
    }
  } else if (vert3 == connQuad[ind0]) {
    /*--- Vert3 coincides with the first vertex of the face connectivity.
          Set the coefficients a and d accordingly.  ---*/
    a = 0;
    d = nPoly;
    if (vert1 != connQuad[ind2]) verticesDontMatch = true;

    /*--- Check the situation for the neighboring vertices. ---*/
    if (vert0 == connQuad[ind3]) {
      /*--- The i-direction coincides while the j-direction is negated. ---*/
      if (vert2 != connQuad[ind1]) verticesDontMatch = true;

      b = 1;
      f = -1;
      c = e = 0;
    } else if (vert0 == connQuad[ind1]) {
      /*--- The j-direction of the current face corresponds with the i-direction
            of the target, while the i-direction coincides with the negative
            j-direction of the target.    ---*/
      if (vert2 != connQuad[ind3]) verticesDontMatch = true;

      b = f = 0;
      c = 1;
      e = -1;
    } else {
      verticesDontMatch = true;
    }
  } else {
    verticesDontMatch = true;
  }

  /*--- If non-matching vertices have been found, terminate with an error message. ---*/
  if (verticesDontMatch) SU2_MPI::Error("Corner vertices do not match. This should not happen.", CURRENT_FUNCTION);

  /*--- Copy the connectivity, such that things works out correctly when carrying
        out the renumbering.      ---*/
  vector<unsigned short> connQuadOr = connQuad;

  /*--- Loop over the vertices of the original face to copy the connectivity data. ---*/
  unsigned short ind = 0;
  for (unsigned short j = 0; j <= nPoly; ++j) {
    for (unsigned short i = 0; i <= nPoly; ++i, ++ind) {
      /*--- Determine the ii and jj indices of the target, convert it to
            a 1D index and shore the modified index in connQuad. ---*/
      unsigned short ii = a + i * b + j * c;
      unsigned short jj = d + i * e + j * f;

      unsigned short iind = jj * (nPoly + 1) + ii;

      connQuad[iind] = connQuadOr[ind];
    }
  }
}

void CFEMStandardElement::ChangeDirectionTriangleConn(vector<unsigned short>& connTriangle, unsigned short vert0,
                                                      unsigned short vert1, unsigned short vert2) const {
  /*--- Determine the indices of the 3 corner vertices of the triangle. ---*/
  unsigned short ind0 = 0;
  unsigned short ind1 = nPoly;
  unsigned short ind2 = (nPoly + 1) * (nPoly + 2) / 2 - 1;

  /*--- There exists a linear mapping from the indices of the numbering used in the
        connectivity of this face to the indices of the target numbering. This
        mapping is of the form ii = a + b*i + c*j and jj = d + e*i + f*j, where
        ii,jj are the indices of the target numbering and i,j the indices of the
        numbering used for this face. The values of the coefficients a,b,c,d,e,f
        depend on how the corner points coincide with each other. This is
        determined below. The bool verticesDontMatch is there to check if vertices
        do not match. This should not happen, but it is checked for security. ---*/
  signed short a = 0, b = 0, c = 0, d = 0, e = 0, f = 0;
  bool verticesDontMatch = false;

  if (vert0 == connTriangle[ind0]) {
    /*--- Vert0 coincides with the first vertex of the face connectivity.
          Check the situation for the neighboring vertices.  ---*/
    if (vert1 == connTriangle[ind1]) {
      /*--- The vertex numbering is the same for both faces. ---*/
      if (vert2 != connTriangle[ind2]) verticesDontMatch = true;

      a = 0;
      b = 1;
      c = 0;
      d = 0;
      e = 0;
      f = 1;
    } else if (vert1 == connTriangle[ind2]) {
      /*--- The ii-index corresponds to the j-index and the
            jj-index corresponds to i-index.    ---*/
      if (vert2 != connTriangle[ind1]) verticesDontMatch = true;

      a = 0;
      b = 0;
      c = 1;
      d = 0;
      e = 1;
      f = 0;
    } else {
      verticesDontMatch = true;
    }
  } else if (vert0 == connTriangle[ind1]) {
    /*--- Vert0 coincides with the second vertex of the face connectivity.
          Check the situation for the neighboring vertices.    ---*/
    if (vert1 == connTriangle[ind0]) {
      /*--- The ii-index corresponds to a combination of the i and j index
            and the jj-index corresponds to the j-index.   ---*/
      if (vert2 != connTriangle[ind2]) verticesDontMatch = true;

      a = nPoly;
      b = -1;
      c = -1;
      d = 0;
      e = 0;
      f = 1;
    } else if (vert1 == connTriangle[ind2]) {
      /*--- The jj-index corresponds to a combination of the i and j index
            and the ii-index corresponds to the j-index.  ---*/
      if (vert2 != connTriangle[ind0]) verticesDontMatch = true;

      a = 0;
      b = 0;
      c = 1;
      d = nPoly;
      e = -1;
      f = -1;
    } else {
      verticesDontMatch = true;
    }
  } else if (vert0 == connTriangle[ind2]) {
    /*--- Vert0 coincides with the third vertex of the face connectivity.
          Check the situation for the neighboring vertices.    ---*/
    if (vert1 == connTriangle[ind0]) {
      /*--- The ii-index corresponds to a combination of the i and j index
            and the jj-index corresponds with the i-index.  ---*/
      if (vert2 != connTriangle[ind1]) verticesDontMatch = true;

      a = nPoly;
      b = -1;
      c = -1;
      d = 0;
      e = 1;
      f = 0;
    } else if (vert1 == connTriangle[ind1]) {
      /*--- The jj-index corresponds to a combination of the i and j index
            and the ii-index corresponds with the i-index.   ---*/
      if (vert2 != connTriangle[ind0]) verticesDontMatch = true;

      a = 0;
      b = 1;
      c = 0;
      d = nPoly;
      e = -1;
      f = -1;
    } else {
      verticesDontMatch = true;
    }
  } else {
    verticesDontMatch = true;
  }

  /*--- If non-matching vertices have been found, terminate with an error message. ---*/
  if (verticesDontMatch) SU2_MPI::Error("Corner vertices do not match. This should not happen.", CURRENT_FUNCTION);

  /*--- Copy the connectivity, such that things works out correctly when carrying
        out the renumbering.  ---*/
  vector<unsigned short> connTriangleOr = connTriangle;

  /*--- Loop over the vertices of the original face to copy the connectivity data. ---*/
  unsigned short ind = 0;
  for (unsigned short j = 0; j <= nPoly; ++j) {
    for (unsigned short i = 0; i <= (nPoly - j); ++i, ++ind) {
      /*--- Determine the ii and jj indices of the target, convert it to
            a 1D index and shore the modified index in connTriangle. ---*/
      unsigned short ii = a + i * b + j * c;
      unsigned short jj = d + i * e + j * f;

      unsigned short iind = jj * (nPoly + 1) + ii - jj * (jj - 1) / 2;

      connTriangle[iind] = connTriangleOr[ind];
    }
  }
}

/*----------------------------------------------------------------------------------*/
/*         Public member functions of CFEMStandardInternalFace.                     */
/*----------------------------------------------------------------------------------*/

CFEMStandardInternalFace::CFEMStandardInternalFace(unsigned short val_VTK_TypeFace, unsigned short val_VTK_TypeSide0,
                                                   unsigned short val_nPolySide0, unsigned short val_VTK_TypeSide1,
                                                   unsigned short val_nPolySide1, bool val_constJac,
                                                   bool val_swapFaceInElementSide0, bool val_swapFaceInElementSide1,
                                                   CConfig* config, unsigned short val_orderExact)

    : CFEMStandardElementBase(val_VTK_TypeFace, max(val_nPolySide0, val_nPolySide1), val_constJac, config,
                              val_orderExact) {
  /*--- Copy the function arguments to the member variables. ---*/
  nPolyElemSide0 = val_nPolySide0;
  VTK_TypeElemSide0 = val_VTK_TypeSide0;
  nPolyElemSide1 = val_nPolySide1;
  VTK_TypeElemSide1 = val_VTK_TypeSide1;
  swapFaceInElementSide0 = val_swapFaceInElementSide0;
  swapFaceInElementSide1 = val_swapFaceInElementSide1;

  /*--- Create dummy variable to store inverse of Vandermonde matrix temporarily. ---*/
  vector<su2double> matVandermondeInvDummy;

  /*--- Determine the constant in the penalty parameter in the viscous
        discretization. This constant depends on the element type on
        both sides of the face as well as their polynomial degree. ---*/
  const su2double penSide0 = ViscousPenaltyParameter(VTK_TypeElemSide0, nPolyElemSide0);
  const su2double penSide1 = ViscousPenaltyParameter(VTK_TypeElemSide1, nPolyElemSide1);

  penaltyConstantFace = max(penSide0, penSide1);

  /*--- Determine the Lagrangian basis functions and its gradients in the
        integration points for both sides of the face. ---*/
  switch (VTK_Type) {
    case LINE:
      LagrangianBasisFunctionAndDerivativesLine(nPolyElemSide0, rIntegration, nDOFsFaceSide0, rDOFsFaceSide0,
                                                matVandermondeInvDummy, lagBasisFaceIntegrationSide0,
                                                drLagBasisFaceIntegrationSide0);
      LagrangianBasisFunctionAndDerivativesLine(nPolyElemSide1, rIntegration, nDOFsFaceSide1, rDOFsFaceSide1,
                                                matVandermondeInvDummy, lagBasisFaceIntegrationSide1,
                                                drLagBasisFaceIntegrationSide1);
      break;

    case TRIANGLE:
      LagrangianBasisFunctionAndDerivativesTriangle(nPolyElemSide0, rIntegration, sIntegration, nDOFsFaceSide0,
                                                    rDOFsFaceSide0, sDOFsFaceSide0, matVandermondeInvDummy,
                                                    lagBasisFaceIntegrationSide0, drLagBasisFaceIntegrationSide0,
                                                    dsLagBasisFaceIntegrationSide0);
      LagrangianBasisFunctionAndDerivativesTriangle(nPolyElemSide1, rIntegration, sIntegration, nDOFsFaceSide1,
                                                    rDOFsFaceSide1, sDOFsFaceSide1, matVandermondeInvDummy,
                                                    lagBasisFaceIntegrationSide1, drLagBasisFaceIntegrationSide1,
                                                    dsLagBasisFaceIntegrationSide1);
      break;

    case QUADRILATERAL:
      LagrangianBasisFunctionAndDerivativesQuadrilateral(nPolyElemSide0, rIntegration, sIntegration, nDOFsFaceSide0,
                                                         rDOFsFaceSide0, sDOFsFaceSide0, matVandermondeInvDummy,
                                                         lagBasisFaceIntegrationSide0, drLagBasisFaceIntegrationSide0,
                                                         dsLagBasisFaceIntegrationSide0);
      LagrangianBasisFunctionAndDerivativesQuadrilateral(nPolyElemSide1, rIntegration, sIntegration, nDOFsFaceSide1,
                                                         rDOFsFaceSide0, sDOFsFaceSide1, matVandermondeInvDummy,
                                                         lagBasisFaceIntegrationSide1, drLagBasisFaceIntegrationSide1,
                                                         dsLagBasisFaceIntegrationSide1);
      break;
  }

  /*--- Check the sum of the Lagrangian basis functions. ---*/
  CheckSumLagrangianBasisFunctions(nIntegration, nDOFsFaceSide0, lagBasisFaceIntegrationSide0);
  CheckSumLagrangianBasisFunctions(nIntegration, nDOFsFaceSide1, lagBasisFaceIntegrationSide1);

  /*--- Create the transpose versions of lagBasisFaceIntegrationSide0 and
        lagBasisFaceIntegrationSide1. These are needed for an efficient computation
        of the residual of the faces. ---*/
  lagBasisFaceIntegrationTransposeSide0.resize(nIntegration * nDOFsFaceSide0);
  lagBasisFaceIntegrationTransposeSide1.resize(nIntegration * nDOFsFaceSide1);

  unsigned int ii = 0;
  for (unsigned short j = 0; j < nDOFsFaceSide0; ++j) {
    for (unsigned short i = 0; i < nIntegration; ++i, ++ii) {
      const unsigned int ind = i * nDOFsFaceSide0 + j;
      lagBasisFaceIntegrationTransposeSide0[ii] = lagBasisFaceIntegrationSide0[ind];
    }
  }

  ii = 0;
  for (unsigned short j = 0; j < nDOFsFaceSide1; ++j) {
    for (unsigned short i = 0; i < nIntegration; ++i, ++ii) {
      const unsigned int ind = i * nDOFsFaceSide1 + j;
      lagBasisFaceIntegrationTransposeSide1[ii] = lagBasisFaceIntegrationSide1[ind];
    }
  }

  /*--- Check the sum of the derivatives of the Lagrangian basis functions. ---*/
  if (!drLagBasisFaceIntegrationSide0.empty())
    CheckSumDerivativesLagrangianBasisFunctions(nIntegration, nDOFsFaceSide0, drLagBasisFaceIntegrationSide0);
  if (!dsLagBasisFaceIntegrationSide0.empty())
    CheckSumDerivativesLagrangianBasisFunctions(nIntegration, nDOFsFaceSide0, dsLagBasisFaceIntegrationSide0);
  if (!drLagBasisFaceIntegrationSide1.empty())
    CheckSumDerivativesLagrangianBasisFunctions(nIntegration, nDOFsFaceSide1, drLagBasisFaceIntegrationSide1);
  if (!dsLagBasisFaceIntegrationSide1.empty())
    CheckSumDerivativesLagrangianBasisFunctions(nIntegration, nDOFsFaceSide1, dsLagBasisFaceIntegrationSide1);

  /*--- Determine the derivatives of the Lagrangian basis functions of the
        adjacent elements in the integration points of the face for both sides
        of the face. ---*/
  DerivativesBasisFunctionsAdjacentElement(VTK_TypeElemSide0, nPolyElemSide0, swapFaceInElementSide0, nDOFsElemSide0,
                                           drLagBasisElemIntegrationSide0, dsLagBasisElemIntegrationSide0,
                                           dtLagBasisElemIntegrationSide0);

  DerivativesBasisFunctionsAdjacentElement(VTK_TypeElemSide1, nPolyElemSide1, swapFaceInElementSide1, nDOFsElemSide1,
                                           drLagBasisElemIntegrationSide1, dsLagBasisElemIntegrationSide1,
                                           dtLagBasisElemIntegrationSide1);

  /*--- Check the sum of the derivatives of the Lagrangian basis functions. ---*/
  if (!drLagBasisElemIntegrationSide0.empty())
    CheckSumDerivativesLagrangianBasisFunctions(nIntegration, nDOFsElemSide0, drLagBasisElemIntegrationSide0);
  if (!dsLagBasisElemIntegrationSide0.empty())
    CheckSumDerivativesLagrangianBasisFunctions(nIntegration, nDOFsElemSide0, dsLagBasisElemIntegrationSide0);
  if (!dtLagBasisElemIntegrationSide0.empty())
    CheckSumDerivativesLagrangianBasisFunctions(nIntegration, nDOFsElemSide0, dtLagBasisElemIntegrationSide0);

  if (!drLagBasisElemIntegrationSide1.empty())
    CheckSumDerivativesLagrangianBasisFunctions(nIntegration, nDOFsElemSide1, drLagBasisElemIntegrationSide1);
  if (!dsLagBasisElemIntegrationSide1.empty())
    CheckSumDerivativesLagrangianBasisFunctions(nIntegration, nDOFsElemSide1, dsLagBasisElemIntegrationSide1);
  if (!dtLagBasisElemIntegrationSide1.empty())
    CheckSumDerivativesLagrangianBasisFunctions(nIntegration, nDOFsElemSide1, dtLagBasisElemIntegrationSide1);

  /*--- If the BLAS routines are used it is beneficial to store all the
        derivatives of the element basis functions in one array. ---*/
  unsigned long sizeMatSide0 = drLagBasisElemIntegrationSide0.size() + dsLagBasisElemIntegrationSide0.size() +
                               dtLagBasisElemIntegrationSide0.size();
  unsigned long sizeMatSide1 = drLagBasisElemIntegrationSide1.size() + dsLagBasisElemIntegrationSide1.size() +
                               dtLagBasisElemIntegrationSide1.size();

  matDerBasisElemIntegrationSide0.resize(sizeMatSide0);
  matDerBasisElemIntegrationSide1.resize(sizeMatSide1);

  matDerBasisElemIntegrationTransposeSide0.resize(sizeMatSide0);
  matDerBasisElemIntegrationTransposeSide1.resize(sizeMatSide1);

  ii = 0;
  for (unsigned long i = 0; i < drLagBasisElemIntegrationSide0.size(); ++i, ++ii)
    matDerBasisElemIntegrationSide0[ii] = drLagBasisElemIntegrationSide0[i];

  for (unsigned long i = 0; i < dsLagBasisElemIntegrationSide0.size(); ++i, ++ii)
    matDerBasisElemIntegrationSide0[ii] = dsLagBasisElemIntegrationSide0[i];

  for (unsigned long i = 0; i < dtLagBasisElemIntegrationSide0.size(); ++i, ++ii)
    matDerBasisElemIntegrationSide0[ii] = dtLagBasisElemIntegrationSide0[i];

  ii = 0;
  for (unsigned long i = 0; i < drLagBasisElemIntegrationSide1.size(); ++i, ++ii)
    matDerBasisElemIntegrationSide1[ii] = drLagBasisElemIntegrationSide1[i];

  for (unsigned long i = 0; i < dsLagBasisElemIntegrationSide1.size(); ++i, ++ii)
    matDerBasisElemIntegrationSide1[ii] = dsLagBasisElemIntegrationSide1[i];

  for (unsigned long i = 0; i < dtLagBasisElemIntegrationSide1.size(); ++i, ++ii)
    matDerBasisElemIntegrationSide1[ii] = dtLagBasisElemIntegrationSide1[i];

  /*--- Also create the transpose of the matDerBasisElemIntegrationSide0 and
        matDerBasisElemIntegrationSide1, such that the residual of the
        symmetrizing terms can be computed efficiently. ---*/
  const unsigned short nDim = !dtLagBasisElemIntegrationSide0.empty() ? 3 : 2;

  ii = 0;
  for (unsigned short j = 0; j < nDOFsElemSide0; ++j) {
    for (unsigned short i = 0; i < nIntegration; ++i) {
      for (unsigned short iDim = 0; iDim < nDim; ++iDim, ++ii) {
        const unsigned int ind = iDim * nDOFsElemSide0 * nIntegration + i * nDOFsElemSide0 + j;
        matDerBasisElemIntegrationTransposeSide0[ii] = matDerBasisElemIntegrationSide0[ind];
      }
    }
  }

  ii = 0;
  for (unsigned short j = 0; j < nDOFsElemSide1; ++j) {
    for (unsigned short i = 0; i < nIntegration; ++i) {
      for (unsigned short iDim = 0; iDim < nDim; ++iDim, ++ii) {
        const unsigned int ind = iDim * nDOFsElemSide1 * nIntegration + i * nDOFsElemSide1 + j;
        matDerBasisElemIntegrationTransposeSide1[ii] = matDerBasisElemIntegrationSide1[ind];
      }
    }
  }
}

bool CFEMStandardInternalFace::SameStandardMatchingFace(unsigned short val_VTK_TypeFace, bool val_constJac,
                                                        unsigned short val_VTK_TypeSide0, unsigned short val_nPolySide0,
                                                        unsigned short val_VTK_TypeSide1, unsigned short val_nPolySide1,
                                                        bool val_swapFaceInElementSide0,
                                                        bool val_swapFaceInElementSide1) {
  if (val_VTK_TypeFace != VTK_Type) return false;
  if (val_constJac != constJacobian) return false;
  if (val_VTK_TypeSide0 != VTK_TypeElemSide0) return false;
  if (val_nPolySide0 != nPolyElemSide0) return false;
  if (val_VTK_TypeSide1 != VTK_TypeElemSide1) return false;
  if (val_nPolySide1 != nPolyElemSide1) return false;
  if (val_swapFaceInElementSide0 != swapFaceInElementSide0) return false;
  if (val_swapFaceInElementSide1 != swapFaceInElementSide1) return false;

  return true;
}

/*----------------------------------------------------------------------------------*/
/*        Private member functions of CFEMStandardInternalFace.                     */
/*----------------------------------------------------------------------------------*/

void CFEMStandardInternalFace::Copy(const CFEMStandardInternalFace& other) {
  CFEMStandardElementBase::Copy(other);

  nDOFsFaceSide0 = other.nDOFsFaceSide0;
  nDOFsFaceSide1 = other.nDOFsFaceSide1;

  nPolyElemSide0 = other.nPolyElemSide0;
  nPolyElemSide1 = other.nPolyElemSide1;
  nDOFsElemSide0 = other.nDOFsElemSide0;
  nDOFsElemSide1 = other.nDOFsElemSide1;

  VTK_TypeElemSide0 = other.VTK_TypeElemSide0;
  VTK_TypeElemSide1 = other.VTK_TypeElemSide1;

  swapFaceInElementSide0 = other.swapFaceInElementSide0;
  swapFaceInElementSide1 = other.swapFaceInElementSide1;

  penaltyConstantFace = other.penaltyConstantFace;

  rDOFsFaceSide0 = other.rDOFsFaceSide0;
  rDOFsFaceSide1 = other.rDOFsFaceSide1;
  sDOFsFaceSide0 = other.sDOFsFaceSide0;
  sDOFsFaceSide1 = other.sDOFsFaceSide1;

  lagBasisFaceIntegrationSide0 = other.lagBasisFaceIntegrationSide0;
  lagBasisFaceIntegrationSide1 = other.lagBasisFaceIntegrationSide1;

  lagBasisFaceIntegrationTransposeSide0 = other.lagBasisFaceIntegrationTransposeSide0;
  lagBasisFaceIntegrationTransposeSide1 = other.lagBasisFaceIntegrationTransposeSide1;

  drLagBasisFaceIntegrationSide0 = other.drLagBasisFaceIntegrationSide0;
  drLagBasisFaceIntegrationSide1 = other.drLagBasisFaceIntegrationSide1;
  dsLagBasisFaceIntegrationSide0 = other.dsLagBasisFaceIntegrationSide0;
  dsLagBasisFaceIntegrationSide1 = other.dsLagBasisFaceIntegrationSide1;

  drLagBasisElemIntegrationSide0 = other.drLagBasisElemIntegrationSide0;
  drLagBasisElemIntegrationSide1 = other.drLagBasisElemIntegrationSide1;
  dsLagBasisElemIntegrationSide0 = other.dsLagBasisElemIntegrationSide0;
  dsLagBasisElemIntegrationSide1 = other.dsLagBasisElemIntegrationSide1;
  dtLagBasisElemIntegrationSide0 = other.dtLagBasisElemIntegrationSide0;
  dtLagBasisElemIntegrationSide1 = other.dtLagBasisElemIntegrationSide1;

  matDerBasisElemIntegrationSide0 = other.matDerBasisElemIntegrationSide0;
  matDerBasisElemIntegrationTransposeSide0 = other.matDerBasisElemIntegrationTransposeSide0;

  matDerBasisElemIntegrationSide1 = other.matDerBasisElemIntegrationSide1;
  matDerBasisElemIntegrationTransposeSide1 = other.matDerBasisElemIntegrationTransposeSide1;
}

/*----------------------------------------------------------------------------------*/
/*         Public member functions of CFEMStandardBoundaryFace.                     */
/*----------------------------------------------------------------------------------*/

CFEMStandardBoundaryFace::CFEMStandardBoundaryFace(unsigned short val_VTK_TypeFace, unsigned short val_VTK_TypeElem,
                                                   unsigned short val_nPolyElem, bool val_constJac,
                                                   bool val_swapFaceInElement, CConfig* config,
                                                   unsigned short val_orderExact)

    : CFEMStandardElementBase(val_VTK_TypeFace, val_nPolyElem, val_constJac, config, val_orderExact) {
  /*--- Copy the function arguments to the member variables. ---*/
  nPolyElem = val_nPolyElem;
  VTK_TypeElem = val_VTK_TypeElem;
  swapFaceInElement = val_swapFaceInElement;

  /*--- Create dummy variable to store inverse of Vandermonde matrix temporarily. ---*/
  vector<su2double> matVandermondeInvDummy;

  /*--- Determine the constant in the penalty parameter in the viscous
        discretization. This constant depends on the element type as
        well as its polynomial degree. ---*/
  penaltyConstantFace = ViscousPenaltyParameter(VTK_TypeElem, nPolyElem);

  /*--- Determine the Lagrangian basis functions in the integration
        points of the face and the subconnectivity of the linear elements
        used for plotting. ---*/
  switch (VTK_Type) {
    case LINE:
      LagrangianBasisFunctionAndDerivativesLine(nPolyElem, rIntegration, nDOFsFace, rDOFsFace, matVandermondeInvDummy,
                                                lagBasisFaceIntegration, drLagBasisFaceIntegration);
      SubConnForPlottingLine(nPolyElem, subConnForPlotting);
      break;

    case TRIANGLE:
      LagrangianBasisFunctionAndDerivativesTriangle(nPolyElem, rIntegration, sIntegration, nDOFsFace, rDOFsFace,
                                                    sDOFsFace, matVandermondeInvDummy, lagBasisFaceIntegration,
                                                    drLagBasisFaceIntegration, dsLagBasisFaceIntegration);
      SubConnForPlottingTriangle(nPolyElem, subConnForPlotting);
      break;

    case QUADRILATERAL:
      LagrangianBasisFunctionAndDerivativesQuadrilateral(nPolyElem, rIntegration, sIntegration, nDOFsFace, rDOFsFace,
                                                         sDOFsFace, matVandermondeInvDummy, lagBasisFaceIntegration,
                                                         drLagBasisFaceIntegration, dsLagBasisFaceIntegration);
      SubConnForPlottingQuadrilateral(nPolyElem, subConnForPlotting);
      break;
  }

  /*--- Check the sum of the Lagrangian basis functions. ---*/
  CheckSumLagrangianBasisFunctions(nIntegration, nDOFsFace, lagBasisFaceIntegration);

  /*--- Create the transpose version of lagBasisFaceIntegration. This is needed
        for an efficient computation of the residual of the faces. ---*/
  lagBasisFaceIntegrationTranspose.resize(nIntegration * nDOFsFace);

  unsigned int ii = 0;
  for (unsigned short j = 0; j < nDOFsFace; ++j) {
    for (unsigned short i = 0; i < nIntegration; ++i, ++ii) {
      const unsigned int ind = i * nDOFsFace + j;
      lagBasisFaceIntegrationTranspose[ii] = lagBasisFaceIntegration[ind];
    }
  }

  /*--- Check the sum of the derivatives of the Lagrangian basis functions. ---*/
  if (!drLagBasisFaceIntegration.empty())
    CheckSumDerivativesLagrangianBasisFunctions(nIntegration, nDOFsFace, drLagBasisFaceIntegration);
  if (!dsLagBasisFaceIntegration.empty())
    CheckSumDerivativesLagrangianBasisFunctions(nIntegration, nDOFsFace, dsLagBasisFaceIntegration);

  /*--- Determine the derivatives of the Lagrangian basis functions of the
        adjacent element in the integration points of the face. ---*/
  DerivativesBasisFunctionsAdjacentElement(VTK_TypeElem, nPolyElem, swapFaceInElement, nDOFsElem,
                                           drLagBasisElemIntegration, dsLagBasisElemIntegration,
                                           dtLagBasisElemIntegration);

  /*--- Check the sum of the derivatives of the Lagrangian basis functions. ---*/
  if (!drLagBasisElemIntegration.empty())
    CheckSumDerivativesLagrangianBasisFunctions(nIntegration, nDOFsElem, drLagBasisElemIntegration);
  if (!dsLagBasisElemIntegration.empty())
    CheckSumDerivativesLagrangianBasisFunctions(nIntegration, nDOFsElem, dsLagBasisElemIntegration);
  if (!dtLagBasisElemIntegration.empty())
    CheckSumDerivativesLagrangianBasisFunctions(nIntegration, nDOFsElem, dtLagBasisElemIntegration);

  /*--- If the BLAS routines are used it is beneficial to store all the
        derivatives of the element basis functions in one array. ---*/
  unsigned long sizeMat =
      drLagBasisElemIntegration.size() + dsLagBasisElemIntegration.size() + dtLagBasisElemIntegration.size();

  matDerBasisElemIntegration.resize(sizeMat);
  matDerBasisElemIntegrationTranspose.resize(sizeMat);

  ii = 0;
  for (unsigned long i = 0; i < drLagBasisElemIntegration.size(); ++i, ++ii)
    matDerBasisElemIntegration[ii] = drLagBasisElemIntegration[i];

  for (unsigned long i = 0; i < dsLagBasisElemIntegration.size(); ++i, ++ii)
    matDerBasisElemIntegration[ii] = dsLagBasisElemIntegration[i];

  for (unsigned long i = 0; i < dtLagBasisElemIntegration.size(); ++i, ++ii)
    matDerBasisElemIntegration[ii] = dtLagBasisElemIntegration[i];

  /*--- Also create the transpose of the matDerBasisElemIntegration, such that
        the residual of the symmetrizing terms can be computed efficiently. ---*/
  const unsigned short nDim = !dtLagBasisElemIntegration.empty() ? 3 : 2;

  ii = 0;
  for (unsigned short j = 0; j < nDOFsElem; ++j) {
    for (unsigned short i = 0; i < nIntegration; ++i) {
      for (unsigned short iDim = 0; iDim < nDim; ++iDim, ++ii) {
        const unsigned int ind = iDim * nDOFsElem * nIntegration + i * nDOFsElem + j;
        matDerBasisElemIntegrationTranspose[ii] = matDerBasisElemIntegration[ind];
      }
    }
  }
}

unsigned short CFEMStandardBoundaryFace::GetNDOFsPerSubFace() const {
  /*--- Distinguish between the possible element types for a boundary surface
        and set the number nDOFs for a subface accordingly. ---*/
  unsigned short nDOFsSubface = 0;
  switch (VTK_Type) {
    case LINE:
      nDOFsSubface = 2;
      break;
    case TRIANGLE:
      nDOFsSubface = 3;
      break;
    case QUADRILATERAL:
      nDOFsSubface = 4;
      break;
    default:
      ostringstream message;
      message << "Impossible FEM surface element type, " << VTK_Type << ", encountered.";
      SU2_MPI::Error(message.str(), CURRENT_FUNCTION);
  }

  /* Return the number of DOFs for a subface. */
  return nDOFsSubface;
}

bool CFEMStandardBoundaryFace::SameStandardBoundaryFace(unsigned short val_VTK_TypeFace, bool val_constJac,
                                                        unsigned short val_VTK_TypeElem, unsigned short val_nPolyElem,
                                                        bool val_swapFaceInElem) {
  if (val_VTK_TypeFace != VTK_Type) return false;
  if (val_constJac != constJacobian) return false;
  if (val_VTK_TypeElem != VTK_TypeElem) return false;
  if (val_nPolyElem != nPolyElem) return false;
  if (val_swapFaceInElem != swapFaceInElement) return false;

  return true;
}

/*----------------------------------------------------------------------------------*/
/*        Private member functions of CFEMStandardBoundaryFace.                     */
/*----------------------------------------------------------------------------------*/

void CFEMStandardBoundaryFace::Copy(const CFEMStandardBoundaryFace& other) {
  CFEMStandardElementBase::Copy(other);

  nDOFsFace = other.nDOFsFace;
  nPolyElem = other.nPolyElem;
  nDOFsElem = other.nDOFsElem;
  VTK_TypeElem = other.VTK_TypeElem;
  swapFaceInElement = other.swapFaceInElement;

  penaltyConstantFace = other.penaltyConstantFace;

  rDOFsFace = other.rDOFsFace;
  sDOFsFace = other.sDOFsFace;

  lagBasisFaceIntegration = other.lagBasisFaceIntegration;
  lagBasisFaceIntegrationTranspose = other.lagBasisFaceIntegrationTranspose;
  drLagBasisFaceIntegration = other.drLagBasisFaceIntegration;
  dsLagBasisFaceIntegration = other.dsLagBasisFaceIntegration;

  drLagBasisElemIntegration = other.drLagBasisElemIntegration;
  dsLagBasisElemIntegration = other.dsLagBasisElemIntegration;
  dtLagBasisElemIntegration = other.dtLagBasisElemIntegration;

  matDerBasisElemIntegration = other.matDerBasisElemIntegration;
  matDerBasisElemIntegrationTranspose = other.matDerBasisElemIntegrationTranspose;

  subConnForPlotting = other.subConnForPlotting;
}
