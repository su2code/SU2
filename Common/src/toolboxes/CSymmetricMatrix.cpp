/*!
 * \file CSymmetricMatrix.cpp
 * \brief Implementation of dense symmetric matrix helper class (see hpp).
 * \author Joel Ho, P. Gomes
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

#include "../../include/toolboxes/CSymmetricMatrix.hpp"
#include "../../include/parallelization/mpi_structure.hpp"
#include "../../include/linear_algebra/blas_structure.hpp"

using namespace std;

#if defined(HAVE_MKL)
#include "mkl.h"
#ifndef HAVE_LAPACK
#define HAVE_LAPACK
#endif
#elif defined(HAVE_LAPACK)
/*--- Lapack / Blas routines used in CSymmetricMatrix. ---*/
extern "C" void dsytrf_(const char*, const int*, passivedouble*, const int*, int*, passivedouble*, const int*, int*);
extern "C" void dsytri_(const char*, const int*, passivedouble*, const int*, const int*, passivedouble*, int*);
extern "C" void dpotrf_(const char*, const int*, passivedouble*, const int*, int*);
extern "C" void dpotri_(const char*, const int*, passivedouble*, const int*, int*);
extern "C" void dsymm_(const char*, const char*, const int*, const int*, const passivedouble*, const passivedouble*,
                       const int*, const passivedouble*, const int*, const passivedouble*, passivedouble*, const int*);
#define DSYMM dsymm_
#endif

void CSymmetricMatrix::Initialize(int N) { mat.resize(N, N); }

void CSymmetricMatrix::CholeskyDecompose() {
#ifndef HAVE_LAPACK
  int j;
  for (j = 0; j < Size(); ++j) {
    passivedouble sum = 0.0;
    for (int k = 0; k < j; ++k) sum -= pow(Get(j, k), 2);
    sum += Get(j, j);
    if (sum < 0.0) break;  // not SPD
    Set(j, j, sqrt(sum));

    for (int i = j + 1; i < Size(); ++i) {
      passivedouble sum = 0.0;
      for (int k = 0; k < j; ++k) sum -= Get(i, k) * Get(j, k);
      sum += Get(i, j);
      Set(i, j, sum / Get(j, j));
    }
  }
  if (j != Size()) SU2_MPI::Error("LLT factorization failed.", CURRENT_FUNCTION);
#endif
}

void CSymmetricMatrix::CalcInv(bool is_spd) {
#ifndef HAVE_LAPACK
  const int sz = Size();

  /*--- Compute inverse from decomposed matrices. ---*/
  if (is_spd) {
    CholeskyDecompose();

    /*--- Initialize inverse matrix. ---*/
    CSymmetricMatrix inv(sz);

    /*--- Compute L inverse. ---*/
    /*--- Solve smaller and smaller systems. ---*/
    for (int j = 0; j < sz; ++j) {
      /*--- Forward substitution. ---*/
      inv(j, j) = 1.0 / Get(j, j);

      for (int i = j + 1; i < sz; ++i) {
        passivedouble sum = 0.0;
        for (int k = j; k < i; ++k) sum -= Get(i, k) * inv(k, j);
        inv(i, j) = sum / Get(i, i);
      }
    }  // L inverse in inv

    /*--- Multiply inversed matrices overwrite mat. ---*/
    for (int j = 0; j < sz; ++j)
      for (int i = j; i < sz; ++i) {
        passivedouble sum = 0.0;
        for (int k = i; k < sz; ++k) sum += inv(k, i) * inv(k, j);
        Set(i, j, sum);
      }
  } else {
    auto inv = StealData();
    CBlasStructure::inverse(sz, inv);
    mat = move(inv);
  }
#endif
}

void CSymmetricMatrix::CalcInv_sytri() {
#ifdef HAVE_LAPACK
  const char uplo = 'L';
  const int sz = Size();
  int info;
  vector<int> ipiv(sz);

  /*--- Query the optimum work size. ---*/
  int query = -1;
  passivedouble tmp;
  dsytrf_(&uplo, &sz, mat.data(), &sz, ipiv.data(), &tmp, &query, &info);
  query = static_cast<int>(tmp);
  vector<passivedouble> work(query);

  /*--- Factorize and invert. ---*/
  dsytrf_(&uplo, &sz, mat.data(), &sz, ipiv.data(), work.data(), &query, &info);
  if (info != 0) SU2_MPI::Error("LDLT factorization failed.", CURRENT_FUNCTION);
  dsytri_(&uplo, &sz, mat.data(), &sz, ipiv.data(), work.data(), &info);
  if (info != 0) SU2_MPI::Error("Inversion with LDLT factorization failed.", CURRENT_FUNCTION);
#endif
}

void CSymmetricMatrix::CalcInv_potri() {
#ifdef HAVE_LAPACK
  const char uplo = 'L';
  const int sz = Size();
  int info;

  dpotrf_(&uplo, &sz, mat.data(), &sz, &info);
  if (info != 0) SU2_MPI::Error("LLT factorization failed.", CURRENT_FUNCTION);
  dpotri_(&uplo, &sz, mat.data(), &sz, &info);
  if (info != 0) SU2_MPI::Error("Inversion with LLT factorization failed.", CURRENT_FUNCTION);
#endif
}

void CSymmetricMatrix::Invert(const bool is_spd) {
#ifdef HAVE_LAPACK
  if (is_spd)
    CalcInv_potri();
  else
    CalcInv_sytri();
#else
  CalcInv(is_spd);
#endif
}

void CSymmetricMatrix::MatMatMult(const char side, const su2passivematrix& mat_in, su2passivematrix& mat_out) const {
  /*--- Left side: mat_out = this * mat_in. ---*/
  if (side == 'L' || side == 'l') {
    const int M = Size(), N = mat_in.cols();
    assert(M == static_cast<int>(mat_in.rows()));

    mat_out.resize(M, N);

#ifdef HAVE_LAPACK
    /*--- Right and lower because matrices are in row major order. ---*/
    const char side = 'R', uplo = 'L';
    const passivedouble alpha = 1.0, beta = 0.0;
    DSYMM(&side, &uplo, &N, &M, &alpha, mat.data(), &M, mat_in.data(), &N, &beta, mat_out.data(), &N);
#else  // Naive product
    for (int i = 0; i < M; ++i)
      for (int j = 0; j < N; ++j) {
        mat_out(i, j) = 0.0;
        for (int k = 0; k < M; ++k) mat_out(i, j) += Get(i, k) * mat_in(k, j);
      }
#endif
  }
  /*--- Right_side: mat_out = mat_in * this. ---*/
  else {
    const int M = mat_in.rows(), N = Size();
    assert(N == static_cast<int>(mat_in.cols()));

    mat_out.resize(M, N);

#ifdef HAVE_LAPACK
    /*--- Left and lower because matrices are in row major order. ---*/
    const char side = 'L', uplo = 'L';
    const passivedouble alpha = 1.0, beta = 0.0;
    DSYMM(&side, &uplo, &N, &M, &alpha, mat.data(), &N, mat_in.data(), &N, &beta, mat_out.data(), &N);
#else  // Naive product
    for (int i = 0; i < M; ++i)
      for (int j = 0; j < N; ++j) {
        mat_out(i, j) = 0.0;
        for (int k = 0; k < N; ++k) mat_out(i, j) += mat_in(i, k) * Get(k, j);
      }
#endif
  }
}

su2passivematrix CSymmetricMatrix::StealData() {
  /*--- Fill lower triangular part. ---*/
  for (int i = 1; i < Size(); ++i)
    for (int j = 0; j < i; ++j) mat(i, j) = mat(j, i);

  return move(mat);
}
