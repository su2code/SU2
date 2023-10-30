/*!
 * \file blas_structure.hpp
 * \brief Include files and headers of the functions for matrix and vector
          operations, which are typically found in the BLAS libraries.
          The functions are in the <i>blass_structure.cpp</i> file.
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

#pragma once

/* LIBXSMM include files, if supported. */
#ifdef HAVE_LIBXSMM
#include "libxsmm.h"
#endif

class CConfig;

/*!
 * \class CBlasStructure
 * \ingroup BLAS
 * \brief Class, which serves as an interface to the BLAS functionalities needed.
 * \author: E. van der Weide
 * \version 8.0.0 "Harrier"
 */
class CBlasStructure {
 public:
  /*!
   * \brief Constructor of the class. Initialize the constant member variables.
   */
  CBlasStructure(void);

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
  void gemm(const int M, const int N, const int K, const su2double* A, const su2double* B, su2double* C,
            const CConfig* config);

  /*!
   * \brief Function, which carries out a dense matrix vector product
            y = A x. It is a limited version of the BLAS gemv functionality.
   * \param[in]  M  - Number of rows of A and size of y.
   * \param[in]  N  - Number of columns of A and size of x.
   * \param[in]  A  - Input matrix in the multiplication, row major order.
   * \param[in]  x  - Input vector in the multiplication.
   * \param[out] y  - Result of the product A x.
   */
  void gemv(const int M, const int N, const su2double* A, const su2double* x, su2double* y);

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
  void axpy(const int n, const su2double a, const su2double* x, const int incx, su2double* y, const int incy);

  /*!
   * \brief Invert a square matrix.
   * \param[in] M - Size.
   * \param[in,out] mat - Matrix, and inverse on exit.
   */
  template <class Mat>
  static void inverse(const int M, Mat& mat) {
    using Scalar = typename Mat::Scalar;

    /*--- Copy the data from A into the augmented matrix and initialize mat with the identity. ---*/
    Mat aug = mat;
    mat = Scalar(0);
    for (int j = 0; j < M; ++j) mat(j, j) = 1;

    /*--- Outer loop of the Gauss-Jordan elimination. ---*/
    for (int j = 0; j < M; ++j) {
      /*--- Find the pivot in the current column. ---*/
      int jj = j;
      Scalar valMax = fabs(aug(j, j));
      for (int i = j + 1; i < M; ++i) {
        Scalar val = fabs(aug(i, j));
        if (val > valMax) {
          jj = i;
          valMax = val;
        }
      }

      /*--- Swap the rows j and jj, if needed. ---*/
      if (jj > j) {
        for (int k = j; k < M; ++k) std::swap(aug(j, k), aug(jj, k));
        for (int k = 0; k < M; ++k) std::swap(mat(j, k), mat(jj, k));
      }

      /*--- Performing row operations to form required identity
            matrix out of the input matrix.  ---*/
      for (int i = 0; i < M; ++i) {
        if (i == j) continue;
        valMax = aug(i, j) / aug(j, j);
        for (int k = j; k < M; ++k) aug(i, k) -= valMax * aug(j, k);
        for (int k = 0; k < M; ++k) mat(i, k) -= valMax * mat(j, k);
      }

      valMax = 1.0 / aug(j, j);
      for (int k = j; k < M; ++k) aug(j, k) *= valMax;
      for (int k = 0; k < M; ++k) mat(j, k) *= valMax;
    }
  }

  /*!
   * \brief tred2
   * Author:
   *
   * Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
   * Klema, Moler.
   * C++ version by Aashwin Mishra and Jayant Mukhopadhaya.
   *
   * Reference:
   *
   * Martin, Reinsch, Wilkinson,
   * TRED2,
   * Numerische Mathematik,
   * Volume 11, pages 181-195, 1968.
   *
   * James Wilkinson, Christian Reinsch,
   * Handbook for Automatic Computation,
   * Volume II, Linear Algebra, Part 2,
   * Springer, 1971,
   * ISBN: 0387054146,
   * LC: QA251.W67.
   *
   * Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
   * Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
   * Matrix Eigensystem Routines, EISPACK Guide,
   * Lecture Notes in Computer Science, Volume 6,
   * Springer Verlag, 1976,
   * ISBN13: 978-3540075462,
   * LC: QA193.M37
   *
   * \param[in,out] V: matrix that needs to be decomposed
   * \param[in,out] d: holds eigenvalues
   * \param[in,out] e: work vector
   * \param[in] n: order of matrix V
   */
  template <class Mat, class Vec, class W>
  static void tred2(Mat& V, Vec& d, W& e, int n) {
    using Scalar = typename std::decay<decltype(e[0])>::type;

    int i, j, k;

    for (j = 0; j < n; j++) {
      d[j] = V[n - 1][j];
    }

    /* Householder reduction to tridiagonal form. */

    for (i = n - 1; i > 0; i--) {
      /* Scale to avoid under/overflow. */

      Scalar scale = 0.0;
      Scalar h = 0.0;
      for (k = 0; k < i; k++) {
        scale = scale + fabs(d[k]);
      }
      if (scale == 0.0) {
        e[i] = d[i - 1];
        for (j = 0; j < i; j++) {
          d[j] = V[i - 1][j];
          V[i][j] = 0.0;
          V[j][i] = 0.0;
        }
      } else {
        /* Generate Householder vector. */

        for (k = 0; k < i; k++) {
          d[k] /= scale;
          h += d[k] * d[k];
        }
        Scalar f = d[i - 1];
        Scalar g = sqrt(h);
        if (f > 0) {
          g = -g;
        }
        e[i] = scale * g;
        h = h - f * g;
        d[i - 1] = f - g;
        for (j = 0; j < i; j++) {
          e[j] = 0.0;
        }

        /* Apply similarity transformation to remaining columns. */

        for (j = 0; j < i; j++) {
          f = d[j];
          V[j][i] = f;
          g = e[j] + V[j][j] * f;
          for (k = j + 1; k <= i - 1; k++) {
            g += V[k][j] * d[k];
            e[k] += V[k][j] * f;
          }
          e[j] = g;
        }
        f = 0.0;
        for (j = 0; j < i; j++) {
          e[j] /= h;
          f += e[j] * d[j];
        }
        Scalar hh = f / (h + h);
        for (j = 0; j < i; j++) {
          e[j] -= hh * d[j];
        }
        for (j = 0; j < i; j++) {
          f = d[j];
          g = e[j];
          for (k = j; k <= i - 1; k++) {
            V[k][j] -= (f * e[k] + g * d[k]);
          }
          d[j] = V[i - 1][j];
          V[i][j] = 0.0;
        }
      }
      d[i] = h;
    }

    /* Accumulate transformations. */

    for (i = 0; i < n - 1; i++) {
      V[n - 1][i] = V[i][i];
      V[i][i] = 1.0;
      Scalar h = d[i + 1];
      if (h != 0.0) {
        for (k = 0; k <= i; k++) {
          d[k] = V[k][i + 1] / h;
        }
        for (j = 0; j <= i; j++) {
          Scalar g = 0.0;
          for (k = 0; k <= i; k++) {
            g += V[k][i + 1] * V[k][j];
          }
          for (k = 0; k <= i; k++) {
            V[k][j] -= g * d[k];
          }
        }
      }
      for (k = 0; k <= i; k++) {
        V[k][i + 1] = 0.0;
      }
    }
    for (j = 0; j < n; j++) {
      d[j] = V[n - 1][j];
      V[n - 1][j] = 0.0;
    }
    V[n - 1][n - 1] = 1.0;
    e[0] = 0.0;
  }

  /*!
   * \brief tql2
   * Author:
   *
   * Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
   * Klema, Moler.
   * C++ version by Aashwin Mishra and Jayant Mukhopadhaya.
   *
   * Reference:
   *
   * Bowdler, Martin, Reinsch, Wilkinson,
   * TQL2,
   * Numerische Mathematik,
   * Volume 11, pages 293-306, 1968.
   *
   * James Wilkinson, Christian Reinsch,
   * Handbook for Automatic Computation,
   * Volume II, Linear Algebra, Part 2,
   * Springer, 1971,
   * ISBN: 0387054146,
   * LC: QA251.W67.
   *
   * Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
   * Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
   * Matrix Eigensystem Routines, EISPACK Guide,
   * Lecture Notes in Computer Science, Volume 6,
   * Springer Verlag, 1976,
   * ISBN13: 978-3540075462,
   * LC: QA193.M37
   *
   * \param[in,out] V: matrix that will hold the eigenvectors
   * \param[in,out] d: array that will hold the ordered eigenvalues
   * \param[in,out] e: work vector
   * \param[in] n: order of matrix V
   */
  template <class Mat, class Vec, class W>
  static void tql2(Mat& V, Vec& d, W& e, int n) {
    using Scalar = typename std::decay<decltype(e[0])>::type;

    int i, j, k, l;
    for (i = 1; i < n; i++) {
      e[i - 1] = e[i];
    }
    e[n - 1] = 0.0;

    Scalar f = 0.0;
    Scalar tst1 = 0.0;
    Scalar eps = pow(2.0, -52.0);
    for (l = 0; l < n; l++) {
      /* Find small subdiagonal element */

      tst1 = max(tst1, (fabs(d[l]) + fabs(e[l])));
      int m = l;
      while (m < n) {
        if (fabs(e[m]) <= eps * tst1) {
          break;
        }
        m++;
      }

      /* If m == l, d[l] is an eigenvalue, */
      /* otherwise, iterate.               */

      if (m > l) {
        int iter = 0;
        do {
          iter = iter + 1; /* (Could check iteration count here.) */

          /* Compute implicit shift */

          Scalar g = d[l];
          Scalar p = (d[l + 1] - g) / (2.0 * e[l]);
          Scalar r = sqrt(p * p + 1.0);
          if (p < 0) {
            r = -r;
          }
          d[l] = e[l] / (p + r);
          d[l + 1] = e[l] * (p + r);
          Scalar dl1 = d[l + 1];
          Scalar h = g - d[l];
          for (i = l + 2; i < n; i++) {
            d[i] -= h;
          }
          f = f + h;

          /* Implicit QL transformation. */

          p = d[m];
          Scalar c = 1.0;
          Scalar c2 = c;
          Scalar c3 = c;
          Scalar el1 = e[l + 1];
          Scalar s = 0.0;
          Scalar s2 = 0.0;
          for (i = m - 1; i >= l; i--) {
            c3 = c2;
            c2 = c;
            s2 = s;
            g = c * e[i];
            h = c * p;
            r = sqrt(p * p + e[i] * e[i]);
            e[i + 1] = s * r;
            s = e[i] / r;
            c = p / r;
            p = c * d[i] - s * g;
            d[i + 1] = h + s * (c * g + s * d[i]);

            /* Accumulate transformation. */

            for (k = 0; k < n; k++) {
              h = V[k][i + 1];
              V[k][i + 1] = s * V[k][i] + c * h;
              V[k][i] = c * V[k][i] - s * h;
            }
          }
          p = -s * s2 * c3 * el1 * e[l] / dl1;
          e[l] = s * p;
          d[l] = c * p;

          /* Check for convergence. */

        } while (fabs(e[l]) > eps * tst1);
      }
      d[l] = d[l] + f;
      e[l] = 0.0;
    }

    /* Sort eigenvalues and corresponding vectors. */

    for (i = 0; i < n - 1; i++) {
      k = i;
      Scalar p = d[i];
      for (j = i + 1; j < n; j++) {
        if (d[j] < p) {
          k = j;
          p = d[j];
        }
      }
      if (k != i) {
        d[k] = d[i];
        d[i] = p;
        for (j = 0; j < n; j++) {
          p = V[j][i];
          V[j][i] = V[j][k];
          V[j][k] = p;
        }
      }
    }
  }

  /*!
   * \brief Decomposes the symmetric matrix A_ij, into eigenvectors and eigenvalues
   * \param[in] A_i: symmetric matrix to be decomposed
   * \param[in,out] Eig_Vec: stores the eigenvectors
   * \param[in,out] Eig_Val: stores the eigenvalues
   * \param[in] n: order of matrix A_ij
   * \param[in,out] e: work vector
   */
  template <class Mat, class Vec, class W>
  static void EigenDecomposition(const Mat& A_ij, Mat& Eig_Vec, Vec& Eig_Val, int n, W& e) {
    for (int iDim = 0; iDim < n; iDim++) {
      e[iDim] = 0.0;
      for (int jDim = 0; jDim < n; jDim++) {
        Eig_Vec[iDim][jDim] = A_ij[iDim][jDim];
      }
    }
    tred2(Eig_Vec, Eig_Val, e, n);
    tql2(Eig_Vec, Eig_Val, e, n);
  }

  /*!
   * \brief Recomposes the eigenvectors and eigenvalues into a matrix
   * \param[out] A_ij: recomposed matrix
   * \param[in] Eig_Vec: eigenvectors
   * \param[in] Eig_Val: eigenvalues
   * \param[in] n: order of matrix A_ij
   */
  template <class Mat, class Vec>
  static void EigenRecomposition(Mat& A_ij, const Mat& Eig_Vec, const Vec& Eig_Val, int n) {
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        A_ij[i][j] = 0.0;
        for (int k = 0; k < n; k++) A_ij[i][j] += Eig_Vec[i][k] * Eig_Val[k] * Eig_Vec[j][k];
      }
    }
  }

  /*!
   * \brief Algorithm to solve a linear system with a tridiagonal matrix.
   * \param[in] lower - lower diagonal
   * \param[in] main - main diagonal
   * \param[in,out] upper - upper diagonal (modified on exit)
   * \param[in,out] rhs - right hand side on entry, solution on exit
   * \note Same size for all vectors. Use row index for lower and upper vector (e.g. lower[0] does not matter).
   */
  template <class Vec, class Scalar = su2double>
  static void tdma(const Vec& lower, const Vec& main, Vec& upper, Vec& rhs) {
    const int N = main.size();

    upper[0] /= main[0];
    rhs[0] /= main[0];

    for (int i = 1; i < N; i++) {
      const Scalar denom = 1.0 / (main[i] - lower[i] * upper[i - 1]);
      upper[i] *= denom;
      rhs[i] = (rhs[i] - lower[i] * rhs[i - 1]) * denom;
    }

    for (int i = N - 2; i >= 0; i--) rhs[i] -= upper[i] * rhs[i + 1];
  }

 private:
#if !(defined(HAVE_LIBXSMM) || defined(HAVE_BLAS) || defined(HAVE_MKL)) || \
    (defined(CODI_REVERSE_TYPE) || defined(CODI_FORWARD_TYPE))
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
  void gemm_imp(const int m, const int n, const int k, const su2double* a, const su2double* b, su2double* c);

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
  void gemm_inner(int m, int n, int k, const su2double* a, int lda, const su2double* b, int ldb, su2double* c, int ldc);

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
  void gemm_arbitrary(int m, int n, int k, const su2double* a, int lda, const su2double* b, int ldb, su2double* c,
                      int ldc);
#endif
};
