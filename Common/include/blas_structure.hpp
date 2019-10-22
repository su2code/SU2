/*!
 * \file blas_structure.hpp
 * \brief Include files and headers of the functions for matrix and vector
          operations, which are typically found in the BLAS libraries.
          The functions are in the <i>blass_structure.cpp</i> file.
 * \author E. van der Weide
 * \version 7.0.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation 
 * (http://su2foundation.org)
 *
 * Copyright 2012-2019, SU2 Contributors (cf. AUTHORS.md)
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

#pragma once

#include "datatype_structure.hpp"
#include "config_structure.hpp"

/* LIBXSMM include files, if supported. */
#ifdef HAVE_LIBXSMM
#include "libxsmm.h"
#endif

/*!
 * \class CBlasStructure
 * \brief Class, which serves as an interface to the BLAS functionalities needed. 
 * \author: E. van der Weide
 * \version 7.0.0 "Blackbird"
 */
class CBlasStructure {
public:
  /*!
   * \brief Constructor of the class. Initialize the constant member variables.
   */
  CBlasStructure(void);

  /*!
   * \brief Constructor of the class. Nothing to be done.
   */
  ~CBlasStructure(void);

  /*!
   * \brief Function, which carries out a dense matrix product. It is a
            limited version of the BLAS gemm functionality..
   * \param[in]  M  - Number of rows of A and C.
   * \param[in]  N  - Number of columns of B and C.
   * \param[in]  K  - Number of columns of A and number of rows of B.
   * \param[in]  A  - Input matrix in the multiplication.
   * \param[in]  B  - Input matrix in the multiplication.
   * \param[out] C  - Result of the matrix product A*B.
   */
  void gemm(const int M,        const int N,        const int K,
            const su2double *A, const su2double *B, su2double *C,
            CConfig *config);

  /*!
   * \brief Function, which carries out a dense matrix vector product
            y = A x. It is a limited version of the BLAS gemv functionality.
   * \param[in]  M  - Number of rows of A and size of y.
   * \param[in]  N  - Number of columns of A and size of x.
   * \param[in]  A  - Input matrix in the multiplication, row major order.
   * \param[in]  x  - Input vector in the multiplication.
   * \param[out] y  - Result of the product A x.
   */
  void gemv(const int M,        const int N,   const su2double *A,
            const su2double *x, su2double *y);

  /*!
   * \brief Function, to carry out the axpy operation, i.e y += a*x.
   * \param[in]    n   - Number of elements in the vectors x and y.
   * \param[in]    a   - Specifies the scalar a.
   * \param[in]    x   - Array, which must be added, size of x must be
                         at least 1 + (n-1)*abs(incx).
   * param[in]    incx - Specifies the increment of x.
   * param[inout] y    - Vector to be updated, size of y must be
                         at least 1 + (n-1)*abs(incy).
   * param[in]    incy - Specifies the increment of y.
   */
  void axpy(const int n,    const su2double a,  const su2double *x,
            const int incx, su2double *y,       const int incy);

private:

#if !(defined(HAVE_LIBXSMM) || defined(HAVE_BLAS) || defined(HAVE_MKL)) || (defined(CODI_REVERSE_TYPE) || defined(CODI_FORWARD_TYPE))
    /* Blocking parameters for the outer kernel.  We multiply mc x kc blocks of
     the matrix A with kc x nc panels of the matrix B (this approach is referred
     to as `gebp` in the literature). */
  const int mc;
  const int kc;
  const int nc;

  /*!
   * \brief Function, which perform the implementation of the gemm functionality.
   * \param[in]  m  - Number of rows of a and c.
   * \param[in]  n  - Number of columns of b and c.
   * \param[in]  k  - Number of columns of a and number of rows of b.
   * \param[in]  a  - Input matrix in the multiplication.
   * \param[in]  b  - Input matrix in the multiplication.
   * \param[out] c  - Result of the matrix product a*b.
   */
  void gemm_imp(const int m,        const int n,        const int k,
                const su2double *a, const su2double *b, su2double *c);

  /*!
   * \brief Compute a portion of the c matrix one block at a time.
            Handle ragged edges with calls to a slow but general function.
   * \param[in]  m   - Number of rows of a and c.
   * \param[in]  n   - Number of columns of b and c.
   * \param[in]  k   - Number of columns of a and number of rows of b.
   * \param[in]  a   - Input matrix in the multiplication.
   * \param[in]  lda - Leading dimension of the matrix a.
   * \param[in]  b   - Input matrix in the multiplication.
   * \param[in]  ldb - Leading dimension of the matrix b.
   * \param[out] c   - Result of the matrix product a*b.
   * \param[in]  ldc - Leading dimension of the matrix c.
   */
  void gemm_inner(int m, int n, int k, const su2double *a, int lda,
                  const su2double *b, int ldb, su2double *c, int ldc);

  /*!
   * \brief Naive gemm implementation to handle arbitrary sized matrices.
   * \param[in]  m   - Number of rows of a and c.
   * \param[in]  n   - Number of columns of b and c.
   * \param[in]  k   - Number of columns of a and number of rows of b.
   * \param[in]  a   - Input matrix in the multiplication.
   * \param[in]  lda - Leading dimension of the matrix a.
   * \param[in]  b   - Input matrix in the multiplication.
   * \param[in]  ldb - Leading dimension of the matrix b.
   * \param[out] c   - Result of the matrix product a*b.
   * \param[in]  ldc - Leading dimension of the matrix c.
   */
  void gemm_arbitrary(int m, int n, int k, const su2double *a, int lda,
                      const su2double *b, int ldb, su2double *c, int ldc);
#endif
};
