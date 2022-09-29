/*!
 * \file blas_structure.cpp
 * \brief Implementation of the functions that either simulate BLAS functionality
          or interface to an actual BLAS implementation.
 * \author E. van der Weide
 * \version 7.4.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2022, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../include/CConfig.hpp"
#include "../../include/linear_algebra/blas_structure.hpp"
#include "../../include/parallelization/omp_structure.hpp"
#include <cstring>

/* MKL or BLAS, if supported. */
#if defined(PRIMAL_SOLVER) && (defined (HAVE_MKL) || defined(HAVE_BLAS))

/* Function prototypes for the BLAS routines used. */
extern "C" void dgemm_(char*, char*, const int*, const int*, const int*,
                       const passivedouble*, const passivedouble*, const int*,
                       const passivedouble*, const int*,
                       const passivedouble*, passivedouble*, const int*);

extern "C" void dgemv_(char*, const int*, const int*, const passivedouble*,
                       const passivedouble*, const int*, const passivedouble*,
                       const int*, const passivedouble*, passivedouble*, const int*);

extern "C" double ddot_(int *n, double *x, int *incx, double *y, int *incy);

extern "C" void daxpy_(int *n, const double *a, const double *x,
                       const int *incx, double *y, const int *incy);
#endif

/* Constructor. Initialize the const member variables, if needed. */
CBlasStructure::CBlasStructure(void)
#if !defined(PRIMAL_SOLVER) || !(defined(HAVE_BLAS) || defined(HAVE_MKL))
  : mc (256), kc (128), nc (128)
#endif
{}

/* Dense matrix multiplication, gemm functionality. */
void CBlasStructure::gemm(const int M,         const int N,            const int K,
                          const int LDA,       const int LDB,          const int LDC,
                          const bool initZero, const passivedouble *A, const su2double *B,
                          su2double *C,        const CConfig *config) const {

  /* Initialize the variable for the timing, if profiling is active. */
#ifdef PROFILE
  double timeGemm;
  if( config ) config->GEMM_Tick(&timeGemm);
#endif

#if !defined(PRIMAL_SOLVER) || !(defined(HAVE_MKL) || defined(HAVE_BLAS))
  /* Native implementation of the matrix product. This implementation is based
     on https://github.com/flame/how-to-optimize-gemm. */
  gemm_imp(M, N, K, LDA, LDB, LDC, initZero, A, B, C);

#else // MKL and BLAS

  /* The standard blas routine dgemm is used for the multiplication. */
  passivedouble alpha = 1.0;
  passivedouble beta  = initZero ? 0.0 : 1.0;
  char trans = 'N';

  dgemm_(&trans, &trans, &M, &N, &K, &alpha, A, &LDA, B, &LDB, &beta, C, &LDC);

#endif

  /* Store the profiling information, if needed. */
#ifdef PROFILE
  if( config ) config->GEMM_Tock(timeGemm, M, N, K);
#endif
}

/* Dot product. */
su2double CBlasStructure::dot(int n, su2double *x, su2double *y) const {

#if defined(PRIMAL_SOLVER) && (defined (HAVE_BLAS) || defined(HAVE_MKL)) 

  /* The standard blas routine ddot is used. */
  int inc = 1;
  return ddot_(&n, x, &inc, y, &inc);
#else

  /* Native implementation of the dot product. */
  su2double result = 0.0;
  for(int i=0; i<n; ++i)
    result += x[i]*y[i];

  return result;

#endif
}

/* Function, to carry out the axpy operation, i.e y += a*x. */
void CBlasStructure::axpy(const int n,        const su2double a,
                          const su2double *x, su2double *y) const {

#if defined(PRIMAL_SOLVER) && (defined (HAVE_BLAS) || defined(HAVE_MKL)) 

  /* The standard blas routine daxpy is used. */
  int inc = 1;
  daxpy_(&n, &a, x, &inc, y, &inc);
    
#else

  /* Native implementation. */
  for(int i=0; i<n; ++i)
    y[i] += a*x[i];

#endif
}

/* Dense matrix vector multiplication, gemv functionality. */
void CBlasStructure::gemv(const int M,        const int N,   const passivedouble *A,
                          const su2double *x, su2double *y) const {

#if defined(PRIMAL_SOLVER) && (defined (HAVE_BLAS) || defined(HAVE_MKL))

  /* The standard blas routine dgemv is used for the multiplication. */
  passivedouble alpha = 1.0;
  passivedouble beta  = 0.0;
  int  inc   = 1;
  char trans = 'N';

  dgemv_(&trans, &M, &N, &alpha, A, &M, x, &inc, &beta, y, &inc);

#else

  /* Native implementation of the matix vector product.
     Initialize the elements of y to zero. */
  SU2_OMP_SIMD
  for(int i=0; i<M; ++i) y[i] = 0.0;

  /* Carry out the matrix vector product. */
  for(int l=0; l<N; ++l) {
    const passivedouble *AA = A + l*M;

    SU2_OMP_SIMD_IF_NOT_AD
    for(int k=0; k<M; ++k)
      y[k] += AA[k]*x[l];
  }

#endif
}

#if !defined(PRIMAL_SOLVER) || !(defined(HAVE_BLAS) || defined(HAVE_MKL))

/* Macros for accessing submatrices of a matmul using the leading dimension. */
#define A(i, j) a[(j)*lda + (i)]
#define B(i, j) b[(j)*ldb + (i)]
#define C(i, j) c[(j)*ldc + (i)]

/* Function, which performs the implementation of the gemm functionality.  */
void CBlasStructure::gemm_imp(const int m,         const int n,            const int k,
                              const int lda,       const int ldb,          const int ldc,
                              const bool initZero, const passivedouble *a, const su2double *b,
                              su2double *c) const {

  /* Initialize the elements of c to zero, if needed. */
  if( initZero ) {
    SU2_OMP_SIMD
    for(int i=0; i<(m*n); ++i) c[i] = 0.0;
  }

  /* The full matrix multiplication is split in several blocks.
     Loop over these blocks. */
  for(int p=0; p<k; p+=kc) {
    int pb = min(k-p, kc);
    for(int j=0; j<n; j+=nc) {
      int jb = min(n-j, nc);
      for(int i=0; i<m; i+=mc) {
        int ib = min(m-i, mc);

        /* Carry out the multiplication for this block. */
        gemm_inner(ib, jb, pb, &A(i, p), lda, &B(p, j), ldb, &C(i, j), ldc);
      }
    }
  }
}

/* Compute a portion of the c matrix one block at a time.
   Handle ragged edges with calls to a slow but general function. */
void CBlasStructure::gemm_inner(int m, int n, int k, const passivedouble *a, int lda,
                                const su2double *b, int ldb, su2double *c, int ldc) const {

  /* Carry out the multiplication for this block. At the
     moment simply a call to gemm_arbitrary. */
  gemm_arbitrary(m, n, k, a, lda, b, ldb, c, ldc);
}

/* Naive gemm implementation to handle arbitrary sized matrices. */
void CBlasStructure::gemm_arbitrary(int m, int n, int k, const passivedouble *a, int lda,
                                    const su2double *b, int ldb, su2double *c, int ldc) const {

  /* The order of these loops is tuned for column-major matrices. */
  for (int p = 0; p < k; p++) {
    for (int j = 0; j < n; j++) {
      SU2_OMP_SIMD_IF_NOT_AD
      for (int i = 0; i < m; i++) {
        C(i, j) += A(i, p) * B(p, j);
      }
    }
  }
}

#undef C
#undef B
#undef A

#endif
