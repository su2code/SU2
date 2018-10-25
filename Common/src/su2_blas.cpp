/*!
 * \file su2_blas.cpp
 * \brief Functions related to the BLAS functionality.
 * \author E. van der Weide
 * \version 6.0.1 "Cardinal"
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
 * Copyright (C) 2012-2016 SU2, the open-source CFD code.
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

#include "../include/su2_blas.hpp"
#include <cstring>

using namespace std;

#if !defined(HAVE_LIBXSMM) && !defined(HAVE_BLAS) && !defined(HAVE_MKL)

/*--- Create an unnamed namespace to keep the functions for the
      native implementation of the matrix product local. ---*/
namespace {

/* Macros for accessing submatrices of a matmul using the leading dimension. */
#define A(i, j) a[(j)*lda + (i)]
#define B(i, j) b[(j)*ldb + (i)]
#define C(i, j) c[(j)*ldc + (i)]

/* Naive gemm implementation to handle arbitrary sized matrices. */
void gemm_arbitrary(int m, int n, int k, const su2double *a, int lda,
                    const su2double *b, int ldb, su2double *c, int ldc) {

  /* The order of these loops is tuned for column-major matrices. */
  for (int p = 0; p < k; p++) {
    for (int j = 0; j < n; j++) {
      for (int i = 0; i < m; i++) {
        C(i, j) += A(i, p) * B(p, j);
      }
    }
  }
}

/* Blocking parameters for the outer kernel.  We multiply mc x kc blocks of A
   with kc x nc panels of B (this approach is referred to as `gebp` in the
   literature). */
const int mc = 256;
const int kc = 128;
const int nc = 128;

/* Compute a portion of C one block at a time.  Handle ragged edges with calls
   to a slow but general function. */
void gemm_inner(int m, int n, int k, const su2double *a, int lda,
                const su2double *b, int ldb, su2double *c, int ldc) {

  /* Carry out the multiplication for this block. At the
     moment simply a call to gemm_arbitrary. */
  gemm_arbitrary(m, n, k, a, lda, b, ldb, c, ldc);
}

/* Local implementation of the gemm functionality. */
void su2_gemm_imp(const int m,        const int n,        const int k,
                  const su2double *a, const su2double *b, su2double *c) {

  /* Initialize the elements of c to zero. */
  memset(c, 0, m*n*sizeof(su2double));

  /* Set the leading dimensions of the three matrices. */
  const int lda = m;
  const int ldb = k;
  const int ldc = m;

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

#undef C
#undef B
#undef A

} /* namespace */

#endif

/* Dense matrix multiplication, gemm functionality. */
void su2_gemm(const int M,        const int N,        const int K,
              const su2double *A, const su2double *B, su2double *C,
              CConfig *config) {

  /* Initialize the variable for the timing, if profiling is active. */
#ifdef PROFILE
  double timeGemm;
  if( config ) config->GEMM_Tick(&timeGemm);
#endif

#ifdef HAVE_LIBXSMM

  /* The gemm function of libxsmm is used to carry out the multiplication.
     Note that libxsmm_gemm expects the matrices in column major order. That's
     why the calling sequence is different from cblas_dgemm. */
  su2double alpha = 1.0;
  su2double beta  = 0.0;
  char trans = 'N';

  libxsmm_dgemm(&trans, &trans, &N, &M, &K, &alpha, B, &N, A, &K, &beta, C, &N);

#elif defined (HAVE_MKL) || defined(HAVE_BLAS)

  /* The standard blas routine dgemm is used for the multiplication.
     Call dgemm without transposing the matrices. In that case dgemm expects
     the matrices in column major order, see the comments for libxsmm. */
  su2double alpha = 1.0;
  su2double beta  = 0.0;
  char trans = 'N';

  dgemm_(&trans, &trans, &N, &M, &K, &alpha, B, &N, A, &K, &beta, C, &N);

#else

  /* Native implementation of the matrix product. This optimized implementation
     assumes that the matrices are in column major order. This can be
     accomplished by swapping N and M and A and B. This implementation is based
     on https://github.com/flame/how-to-optimize-gemm. */
  su2_gemm_imp(N, M, K, B, A, C);
  
#endif

  /* Store the profiling information, if needed. */
#ifdef PROFILE
  if( config ) config->GEMM_Tock(timeGemm, M, N, K);
#endif
}

/* Dense matrix vector multiplication, gemv functionality. */
void su2_gemv(const int M,        const int N,   const su2double *A,
              const su2double *x, su2double *y) {

#if defined (HAVE_BLAS) || defined(HAVE_MKL)

  /* The standard blas routine dgemv is used for the multiplication.
     Note that dgemv expects the matrices in column major order, while
     A is in row major order. This can be solved by using the transpose
     and switching M and N. */
     su2double alpha = 1.0;
     su2double beta  = 0.0;
     int       inc   = 1;
     char trans = 'T';

     dgemv_(&trans, &N, &M, &alpha, A, &N, x, &inc, &beta, y, &inc);

#else

  /* Native implementation of the matix vector product.
     Initialize the elements of y to zero. */
  memset(y, 0, M*sizeof(su2double));  

  /* Carry out the matrix vector product. */
  for(int k=0; k<M; ++k) {
    const su2double *AA = A + k*N;
    for(int l=0; l<N; ++l)
      y[k] += AA[l]*x[l];
  }

#endif
}
