/*!
 * \file dense_matrix_product.cpp
 * \brief Function used to carry out a dense matrix multiplication.
 * \author E. van der Weide
 * \version 4.3.0 "Cardinal"
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

#include "../include/dense_matrix_product.hpp"

void DenseMatrixProduct(const int M,        const int N,        const int K,
                        const su2double *A, const su2double *B, su2double *C) {

#ifdef HAVE_LIBXSMM

  /* The gemm function of libxsmm is used to carry out the multiplication.
     Note that libxsmm_gemm expects the matrices in column major order. That's
     why the calling sequence is different from cblas_dgemm. */
  libxsmm_gemm(NULL, NULL, N, M, K, NULL, B, NULL, A, NULL, NULL, C, NULL);

#elif defined (HAVE_CBLAS) || defined(HAVE_MKL)

  /* The standard blas routine dgemm is used for the multiplication. */
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K,
              1.0, A, K, B, N, 0.0, C, N);
#else

  /* Standard internal implementation of the matrix matrix multiplication. */
  for(int i=0; i<M; ++i) {
    const int jj = i*K;
    const int kk = i*N;
    for(int j=0; j<N; ++j) {
      const int ii = kk + j;
      C[ii] = 0.0;
      for(int k=0; k<K; ++k)
        C[ii] += A[jj+k]*B[k*N+j];
    }
  }

#endif
}
