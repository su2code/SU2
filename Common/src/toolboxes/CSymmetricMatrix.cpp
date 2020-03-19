/*!
 * \file CSymmetricMatrix.cpp
 * \brief Implementation of dense symmetric matrix helper class (see hpp).
 * \author Joel Ho, P. Gomes
 * \version 7.0.2 "Blackbird"
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

#include "../../include/toolboxes/CSymmetricMatrix.hpp"
#include "../../include/mpi_structure.hpp"

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
#endif


void CSymmetricMatrix::Initialize(int N)
{
  sz = N;
  val_vec.resize(sz*sz);
  initialized = true;
}

void CSymmetricMatrix::CholeskyDecompose()
{
#ifndef HAVE_LAPACK
  assert(initialized && "Matrix not initialized.");

  for (int j = 0; j < sz; ++j) {
    passivedouble sum = 0.0;
    for (int k = 0; k < j; ++k) sum -= pow(Get(j,k), 2);
    sum += Get(j,j);
    assert(sum > 0.0 && "LLT failed, matrix is not SPD.");
    Set(j, j, sqrt(sum));

    for (int i = j+1; i < sz; ++i) {
      passivedouble sum = 0.0;
      for (int k = 0; k < j; ++k) sum -= Get(i,k) * Get(j,k);
      sum += Get(i,j);
      Set(i, j, sum / Get(j,j));
    }
  }
  decomposed = CHOLESKY;
#endif
}

void CSymmetricMatrix::LUDecompose()
{
#ifndef HAVE_LAPACK
  assert(initialized && "Matrix not initialized.");

  /*--- Copy matrix values to LU matrix, init permutation vec. ---*/
  decomp_vec.resize(sz*sz);
  perm_vec.resize(sz);
  for (int i = 0; i < sz; ++i) {
    perm_vec[i] = i;
    for (int j = i; j < sz; ++j) decomp(j,i) = decomp(i,j) = Get(i,j);
  }

  /*--- Decompose LU matrix. ---*/
  for (int j = 0; j < sz-1; ++j) {
    /*--- Search for maximum pivot and interchange rows. ---*/
    passivedouble pivot = decomp(j,j);
    int pivot_idx = j;
    for (int i = j+1; i < sz; ++i)
      if (abs(decomp(i,j)) > abs(pivot)) {
        pivot = decomp(i,j);
        pivot_idx = i;
      }

    if (pivot_idx != j) {
      swap(perm_vec[j], perm_vec[pivot_idx]);
      for (int k = 0; k < sz; ++k)
        swap(decomp(j,k), decomp(pivot_idx,k));
    }

    /*--- Perform elimination. ---*/
    for (int k = j+1; k < sz; ++k) decomp(k,j) /= pivot;

    for (int k = j+1; k < sz; ++k)
      for (int i = j+1; i < sz; ++i)
        decomp(i,k) -= decomp(j,k)*decomp(i,j);
  }

  decomposed = LU;
#endif
}

void CSymmetricMatrix::CalcInv()
{
#ifndef HAVE_LAPACK
  assert(initialized && "Matrix not initialized.");

  /*--- Decompose matrix if not done yet. ---*/
  if (decomposed == NONE) { LUDecompose(); }

  /*--- Compute inverse from decomposed matrices. ---*/
  switch (decomposed) {

    case CHOLESKY:
    {
      /*--- Initialize inverse matrix. ---*/
      vector<passivedouble> inv(sz*sz, 0.0);

      /*--- Compute L inverse. ---*/
      /*--- Solve smaller and smaller systems. ---*/
      for (int j = 0; j < sz; ++j) {
        /*--- Forward substitution. ---*/
        inv[IdxSym(j,j)] = 1.0 / Get(j,j);

        for (int i = j+1; i < sz; ++i) {
          passivedouble sum = 0.0;
          for (int k = j; k < i; ++k) sum -= Get(i,k) * inv[IdxSym(k,j)];
          inv[IdxSym(i,j)] = sum / Get(i,i);
        }
      } // L inverse in inv

      /*--- Multiply inversed matrices overwrite val_vec. ---*/
      for (int j = 0; j < sz; ++j)
        for (int i = j; i < sz; ++i) {
          passivedouble sum = 0.0;
          for (int k = i; k < sz; ++k) sum += inv[IdxSym(k,i)] * inv[IdxSym(k,j)];
          Set(i, j, sum);
        }

      break;
    }

    case LU:
    {
      /*--- Alias val_vec. ---*/
      vector<passivedouble>& inv = val_vec;

      /*--- Invert L and U matrices in place. ---*/
      for (int j = 0; j < sz; ++j) {
        inv[IdxFull(j,j)] = 1.0 / decomp(j,j);

        for (int i = j+1; i < sz; ++i) {
          inv[IdxFull(i,j)] = -decomp(i,j);
          inv[IdxFull(j,i)] = -decomp(j,i) * inv[IdxFull(j,j)];

          for (int k = j+1; k < i; ++k) {
            inv[IdxFull(i,j)] -= decomp(i,k) * inv[IdxFull(k,j)];
            inv[IdxFull(j,i)] -= decomp(k,i) * inv[IdxFull(j,k)];
          }
          if (j+1 <= i) inv[IdxFull(j,i)] /= decomp(i,i);
        }
      }
      // inverses in val_vec

      /*--- Multiply U_inv by L_inv, overwrite decomp_vec. ---*/
      for (int i = 0; i < sz; ++i)
        for (int j = 0; j < sz; ++j) {
          decomp(i,j) = 0.0;
          for (int k = max(i,j); k < sz; ++k)
            decomp(i,j) += inv[IdxFull(i,k)] * ((k==j)? 1.0 : inv[IdxFull(k,j)]);
        }

      /*--- Permute multiplied matrix to recover A_inv, overwrite val_vec. ---*/
      for (int j = 0; j < sz; ++j) {
        int k = perm_vec[j];
        for (int i = k; i < sz; ++i) Set(i, k, decomp(i,j));
      }

      /*--- Decomposition no longer needed. ---*/
      vector<passivedouble>().swap(decomp_vec);

      break;
    }
    default: assert(false && "Default (LU) decomposition failed.");
  }

  decomposed = NONE;
#endif
}

void CSymmetricMatrix::CalcInv_sytri()
{
#ifdef HAVE_LAPACK
  char uplo = 'L';
  int info;
  perm_vec.resize(sz); // ipiv array

  /*--- Query the optimum work size. ---*/
  int query = -1; passivedouble tmp;
  dsytrf_(&uplo, &sz, val_vec.data(), &sz, perm_vec.data(), &tmp, &query, &info);
  query = static_cast<int>(tmp);
  decomp_vec.resize(query); // work array

  /*--- Factorize and invert. ---*/
  dsytrf_(&uplo, &sz, val_vec.data(), &sz, perm_vec.data(), decomp_vec.data(), &query, &info);
  if (info!=0) SU2_MPI::Error("LDLT factorization failed.", CURRENT_FUNCTION);
  dsytri_(&uplo, &sz, val_vec.data(), &sz, perm_vec.data(), decomp_vec.data(), &info);
  if (info!=0) SU2_MPI::Error("Inversion with LDLT factorization failed.", CURRENT_FUNCTION);

  decomposed = NONE;
#endif
}

void CSymmetricMatrix::CalcInv_potri()
{
#ifdef HAVE_LAPACK
  char uplo = 'L';
  int info;

  dpotrf_(&uplo, &sz, val_vec.data(), &sz, &info);
  if (info!=0) SU2_MPI::Error("LLT factorization failed.", CURRENT_FUNCTION);
  dpotri_(&uplo, &sz, val_vec.data(), &sz, &info);
  if (info!=0) SU2_MPI::Error("Inversion with LLT factorization failed.", CURRENT_FUNCTION);

  decomposed = NONE;
#endif
}

void CSymmetricMatrix::Invert(const bool is_spd)
{
#ifdef HAVE_LAPACK
  if(is_spd) CalcInv_potri();
  else CalcInv_sytri();
#else
  if(!is_spd) LUDecompose();
  else CholeskyDecompose();
  CalcInv();
#endif
}

void CSymmetricMatrix::MatMatMult(const char side,
                                  const su2passivematrix& mat_in,
                                  su2passivematrix& mat_out) const
{
  /*--- Assumes row major storage of in/out matrices for LAPACK. ---*/
  static_assert(su2passivematrix::Storage == StorageType::RowMajor,"");

  /*--- Left side: mat_out = this * mat_in. ---*/
  if (side == 'L' || side == 'l') {
    int M = sz, N = mat_in.cols();
    assert(M == mat_in.rows());

    mat_out.resize(M,N);

#if !defined(HAVE_LAPACK) // Naive product
    for (int i = 0; i < M; ++i)
      for (int j = 0; j < N; ++j) {
        mat_out(i,j) = 0.0;
        for (int k = 0; k < M; ++k)
          mat_out(i,j) += Get(i,k) * mat_in(k,j);
      }
#elif defined(HAVE_MKL)
    passivedouble alpha = 1.0, beta = 0.0;
    cblas_dsymm(CblasRowMajor, CblasLeft, CblasUpper, M, N, alpha,
      val_vec.data(), M, mat_in.data(), N, beta, mat_out.data(), N);
#else // BLAS/LAPACK
    /*--- Right and lower because matrices are in row major order. ---*/
    const char side = 'R', uplo = 'L';
    const passivedouble alpha = 1.0, beta = 0.0;
    dsymm_(&side, &uplo, &N, &M, &alpha, val_vec.data(), &M,
           mat_in.data(), &N, &beta, mat_out.data(), &N);
#endif
  }
  /*--- Right_side: mat_out = mat_in * this. ---*/
  else {
    int M = mat_in.rows(), N = sz;
    assert(N == mat_in.cols());

    mat_out.resize(M,N);

#if !defined(HAVE_LAPACK) // Naive product
    for (int i = 0; i < M; ++i)
      for (int j = 0; j < N; ++j) {
        mat_out(i,j) = 0.0;
        for (int k = 0; k < N; ++k)
          mat_out(i,j) += mat_in(i,k) * Get(k,j);
      }
#elif defined(HAVE_MKL)
    passivedouble alpha = 1.0, beta = 0.0;
    cblas_dsymm(CblasRowMajor, CblasRight, CblasUpper, M, N, alpha,
      val_vec.data(), N, mat_in.data(), N, beta, mat_out.data(), N);
#else // BLAS/LAPACK
    /*--- Left and lower because matrices are in row major order. ---*/
    const char side = 'L', uplo = 'L';
    const passivedouble alpha = 1.0, beta = 0.0;
    dsymm_(&side, &uplo, &N, &M, &alpha, val_vec.data(), &N,
           mat_in.data(), &N, &beta, mat_out.data(), &N);
#endif
  }

}
