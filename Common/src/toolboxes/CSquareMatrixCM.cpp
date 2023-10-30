/*!
 * \file CSquareMatrixCM.cpp
 * \brief Implementation of dense matrix helper class in Column Major order (see hpp).
 * \author Edwin van der Weide, Pedro Gomes.
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

#include "../../include/toolboxes/CSquareMatrixCM.hpp"
#include "../../include/parallelization/mpi_structure.hpp"
#include "../../include/linear_algebra/blas_structure.hpp"

using namespace std;

#if defined(HAVE_MKL)
#include "mkl.h"
#ifndef HAVE_LAPACK
#define HAVE_LAPACK
#endif
#elif defined(HAVE_LAPACK)
/*--- Lapack / Blas routines used in CSquareMatrixCM. ---*/
extern "C" void dgetrf_(const int*, const int*, passivedouble*, const int*, int*, int*);
extern "C" void dgetri_(const int*, passivedouble*, const int*, int*, passivedouble*, const int*, int*);
extern "C" void dgemm_(char*, char*, const int*, const int*, const int*, const passivedouble*, const passivedouble*,
                       const int*, const passivedouble*, const int*, const passivedouble*, passivedouble*, const int*);
#define DGEMM dgemm_
#endif

void CSquareMatrixCM::Transpose() {
  for (int j = 1; j < Size(); ++j)
    for (int i = 0; i < j; ++i) swap(mat(i, j), mat(j, i));
}

void CSquareMatrixCM::Invert() {
#ifdef HAVE_LAPACK

  /*--- Computation of the inverse using the Lapack routines. ---*/
  int sz = Size();
  int info;
  vector<int> ipiv(sz);
  vector<passivedouble> work(sz);

  dgetrf_(&sz, &sz, mat.data(), &sz, ipiv.data(), &info);
  if (info != 0) SU2_MPI::Error(string("Matrix is singular"), CURRENT_FUNCTION);

  dgetri_(&sz, mat.data(), &sz, ipiv.data(), work.data(), &sz, &info);
  if (info != 0) SU2_MPI::Error(string("Matrix inversion failed"), CURRENT_FUNCTION);

#else
  CBlasStructure::inverse(Size(), mat);
#endif
}

void CSquareMatrixCM::MatMatMult(const char side, const ColMajorMatrix<passivedouble>& mat_in,
                                 ColMajorMatrix<passivedouble>& mat_out) const {
  /*--- Check the type of multiplication to be carried out. ---*/
  if (side == 'L' || side == 'l') {
    /*--- Left side: mat_out = this * mat_in. Set some sizes
          and allocate the memory for mat_out. ---*/
    const int M = Size(), N = mat_in.cols();
    assert(M == static_cast<int>(mat_in.rows()));

    mat_out.resize(M, N);

#ifdef HAVE_LAPACK

    /*--- The Lapack/blas function dgemm is used to carry out
          the matrix matrix multiplication. ---*/
    passivedouble alpha = 1.0, beta = 0.0;
    char trans = 'N';

    DGEMM(&trans, &trans, &M, &N, &M, &alpha, mat.data(), &M, mat_in.data(), &M, &beta, mat_out.data(), &M);
#else
    /*--- Naive product. ---*/
    for (int i = 0; i < M; ++i) {
      for (int j = 0; j < N; ++j) {
        mat_out(i, j) = 0.0;
        for (int k = 0; k < M; ++k) mat_out(i, j) += mat(i, k) * mat_in(k, j);
      }
    }
#endif

  } else {
    /*--- Right_side: mat_out = mat_in * this. Set some sizes
          and allocate the memory for mat_out. ---*/
    const int M = mat_in.rows(), N = Size();
    assert(N == static_cast<int>(mat_in.cols()));

    mat_out.resize(M, N);

#ifdef HAVE_LAPACK

    /*--- The Lapack/blas function dgemm is used to carry out
          the matrix matrix multiplication. ---*/
    passivedouble alpha = 1.0, beta = 0.0;
    char trans = 'N';

    DGEMM(&trans, &trans, &M, &N, &N, &alpha, mat_in.data(), &M, mat.data(), &N, &beta, mat_out.data(), &M);
#else
    /*--- Naive product. ---*/
    for (int i = 0; i < M; ++i) {
      for (int j = 0; j < N; ++j) {
        mat_out(i, j) = 0.0;
        for (int k = 0; k < N; ++k) mat_out(i, j) += mat_in(i, k) * mat(k, j);
      }
    }
#endif
  }
}
