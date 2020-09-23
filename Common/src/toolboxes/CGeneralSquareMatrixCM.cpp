/*!
 * \file CGeneralSquareMatrixCM.cpp
 * \brief Implementation of dense matrix helper class in Column Major order (see hpp).
 * \author Edwin van der Weide, Pedro Gomes.
 * \version 7.0.6 "Blackbird"
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

#include "../../include/toolboxes/CGeneralSquareMatrixCM.hpp"
#include "../../include/mpi_structure.hpp"

using namespace std;

#if defined(HAVE_MKL)
#include "mkl.h"
#ifndef HAVE_LAPACK
#define HAVE_LAPACK
#endif
#elif defined(HAVE_LAPACK)
/*--- Lapack / Blas routines used in CGeneralSquareMatrixCM. ---*/
extern "C" void dgetrf_(const int*, const int*, passivedouble*, const int*,
                        int*, int*);
extern "C" void dgetri_(const int*, passivedouble*, const int*, int*,
                        passivedouble*, const int*, int*);
extern "C" void dgemm_(char*, char*, const int*, const int*, const int*,
                       const passivedouble*, const passivedouble*,
                       const int *, const passivedouble*, const int*,
                       const passivedouble*, passivedouble*, const int*);
#define DGEMM dgemm_
#endif

void CGeneralSquareMatrixCM::Transpose() {

  for(int j=1; j<Size(); ++j)
    for(int i=0; i<j; ++i)
      swap(mat(i,j), mat(j,i));
}

void CGeneralSquareMatrixCM::Invert() {

#ifdef HAVE_LAPACK

  /*--- Computation of the inverse using the Lapack routines. ---*/
  int sz = Size();
  int info;
  vector<int> ipiv(sz);
  vector<passivedouble> work(sz);

  dgetrf_(&sz, &sz, mat.data(), &sz, ipiv.data(), &info);
  if(info != 0) SU2_MPI::Error(string("Matrix is singular"), CURRENT_FUNCTION);

  dgetri_(&sz, mat.data(), &sz, ipiv.data(), work.data(), &sz, &info);
  if(info != 0) SU2_MPI::Error(string("Matrix inversion failed"), CURRENT_FUNCTION);

#else

  /*--- Standard implementation without Lapack. Create a local augmented
        matrix to carry out the actual inversion. ---*/
  const int M = Size();
  ColMajorMatrix<passivedouble> augmentedmatrix(M, 2*M);

  /*--- Copy the data from A into the first part of augmentedmatrix and
        augmenting with identity matrix of similar dimensions. ---*/
  for(int j=0; j<M; ++j) {
    for(int i=0; i<M; ++i) {
      augmentedmatrix(i,j) = mat(i,j);
      augmentedmatrix(i,j+M) = (i == j) ? 1.0 : 0.0;
    }
  }

  /*--- Outer loop of the Gauss-Jordan elimination. ---*/
  for(int j=0; j<M; ++j) {

    /*--- Find the pivot in the current column. ---*/
    int jj = j;
    passivedouble  valMax = fabs(augmentedmatrix(j,j));
    for(int i=j+1; i<M; ++i) {
      passivedouble val = fabs(augmentedmatrix(i,j));
      if(val > valMax){
        jj = i;
        valMax = val;
      }
    }

    /*--- Swap the rows j and jj, if needed. ---*/
    if(jj > j)
      for(int k=j; k<2*M; ++k)
        swap(augmentedmatrix(j,k), augmentedmatrix(jj,k));

    /*--- Performing row operations to form required identity
          matrix out of the input matrix.  ---*/
    for(int i=0; i<M; ++i) {
      if(i != j) {
        valMax = augmentedmatrix(i,j)/augmentedmatrix(j,j);
        for(int k=j; k<2*M; ++k)
          augmentedmatrix(i,k) -= valMax*augmentedmatrix(j,k);
      }
    }

    valMax = 1.0/augmentedmatrix(j,j);
    for(int k=j; k<2*M; ++k)
      augmentedmatrix(j,k) *= valMax;
  }

  /*--- Store the inverse in mat. ---*/
  for(int j=0; j<M; ++j)
    for(int i=0; i<M; ++i)
      mat(i,j) = augmentedmatrix(i,j+M);

#endif
}

void CGeneralSquareMatrixCM::MatMatMult(const char                          side,
                                        const ColMajorMatrix<passivedouble> &mat_in,
                                        ColMajorMatrix<passivedouble>       &mat_out) const {

  /*--- Check the type of multiplication to be carried out. ---*/
  if (side == 'L' || side == 'l') {

    /*--- Left side: mat_out = this * mat_in. Set some sizes
          and allocate the memory for mat_out. ---*/
    const int M = Size(), N = mat_in.cols();
    assert(M == mat_in.rows());

    mat_out.resize(M,N);

#ifdef HAVE_LAPACK

    /*--- The Lapack/blas function dgemm is used to carry out
          the matrix matrix multiplication. ---*/
    passivedouble alpha = 1.0, beta = 0.0;
    char trans = 'N';

    DGEMM(&trans, &trans, &M, &N, &M, &alpha, mat.data(), &M,
          mat_in.data(), &M, &beta, mat_out.data(), &M);
#else
    /*--- Naive product. ---*/
    for (int i = 0; i < M; ++i) {
      for (int j = 0; j < N; ++j) {
        mat_out(i,j) = 0.0;
        for (int k = 0; k < M; ++k)
          mat_out(i,j) += Get(i,k) * mat_in(k,j);
      }
    }
#endif

  }
  else {

    /*--- Right_side: mat_out = mat_in * this. Set some sizes
          and allocate the memory for mat_out. ---*/
    const int M = mat_in.rows(), N = Size();
    assert(N == mat_in.cols());

    mat_out.resize(M,N);

#ifdef HAVE_LAPACK

    /*--- The Lapack/blas function dgemm is used to carry out
          the matrix matrix multiplication. ---*/
    passivedouble alpha = 1.0, beta = 0.0;
    char trans = 'N';

    DGEMM(&trans, &trans, &M, &N, &N, &alpha, mat_in.data(), &M,
          mat.data(), &N, &beta, mat_out.data(), &M);
#else
    /*--- Naive product. ---*/
    for (int i = 0; i < M; ++i) {
      for (int j = 0; j < N; ++j) {
        mat_out(i,j) = 0.0;
        for (int k = 0; k < N; ++k)
          mat_out(i,j) += mat_in(i,k) * Get(k,j);
      }
    }
#endif
  }
}
