/*!
 * \file CSysMatrix.inl
 * \brief Inline subroutines of the <i>CSysMatrix.hpp</i> file.
 * \note These are the "private" inlines, they are not needed outside
 *       of the .cpp file and so they are hidden to avoid triggering
 *       recompilation of other units when changes are made here.
 * \author F. Palacios, A. Bueno, T. Economon, P. Gomes
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

#include "CSysMatrix.hpp"

template <class ScalarType>
FORCEINLINE ScalarType* CSysMatrix<ScalarType>::GetBlock_ILUMatrix(unsigned long block_i, unsigned long block_j) {
  /*--- The position of the diagonal block is known which allows halving the search space. ---*/
  const auto end = (block_j < block_i) ? dia_ptr_ilu[block_i] : row_ptr_ilu[block_i + 1];
  for (auto index = (block_j < block_i) ? row_ptr_ilu[block_i] : dia_ptr_ilu[block_i]; index < end; ++index)
    if (col_ind_ilu[index] == block_j) return &ILU_matrix[index * nVar * nVar];
  return nullptr;
}

template <class ScalarType>
FORCEINLINE void CSysMatrix<ScalarType>::SetBlock_ILUMatrix(unsigned long block_i, unsigned long block_j,
                                                            ScalarType* val_block) {
  auto ilu_ij = GetBlock_ILUMatrix(block_i, block_j);
  if (!ilu_ij) return;
  MatrixCopy(val_block, ilu_ij);
}

namespace {

template <class T, bool alpha, bool beta, bool transp>
FORCEINLINE void gemv_impl(unsigned long n, unsigned long m, const T* a, const T* b, T* c) {
  /*---
   This is a templated version of GEMV with the constants as boolean
   template parameters so that they can be optimized away at compilation.
   This is still the traditional "row dot vector" method.
  ---*/
  if (!transp) {
    for (auto i = 0ul; i < n; i++) {
      if (!beta) c[i] = 0.0;
      for (auto j = 0ul; j < m; j++) c[i] += (alpha ? 1 : -1) * a[i * m + j] * b[j];
    }
  } else {
    if (!beta)
      for (auto j = 0ul; j < m; j++) c[j] = 0.0;
    for (auto i = 0ul; i < n; i++)
      for (auto j = 0ul; j < m; j++) c[j] += (alpha ? 1 : -1) * a[i * n + j] * b[i];
  }
}

template <class T>
FORCEINLINE void gemm_impl(unsigned long n, const T* a, const T* b, T* c) {
  /*--- Same deal as for GEMV but here only the type is templated. ---*/
  unsigned long i, j, k;
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      c[i * n + j] = 0.0;
      for (k = 0; k < n; k++) c[i * n + j] += a[i * n + k] * b[k * n + j];
    }
  }
}
}  // namespace

#define __MATVECPROD_SIGNATURE__(TYPE, NAME) \
  FORCEINLINE void CSysMatrix<TYPE>::NAME(const TYPE* matrix, const TYPE* vector, TYPE* product) const

#define MATVECPROD_SIGNATURE(NAME) \
  template <class ScalarType>      \
  __MATVECPROD_SIGNATURE__(ScalarType, NAME)

#if !defined(USE_MKL)
MATVECPROD_SIGNATURE(MatrixVectorProduct) {
  /*---
   Without MKL (default) picture copying the body of gemv_impl
   here and resolving the conditionals at compilation.
  ---*/
  gemv_impl<ScalarType, true, false, false>(nVar, nEqn, matrix, vector, product);
}

MATVECPROD_SIGNATURE(MatrixVectorProductAdd) {
  gemv_impl<ScalarType, true, true, false>(nVar, nEqn, matrix, vector, product);
}

MATVECPROD_SIGNATURE(MatrixVectorProductSub) {
  gemv_impl<ScalarType, false, true, false>(nVar, nEqn, matrix, vector, product);
}

template <class ScalarType>
FORCEINLINE void CSysMatrix<ScalarType>::MatrixMatrixProduct(const ScalarType* matrix_a, const ScalarType* matrix_b,
                                                             ScalarType* product) const {
  gemm_impl<ScalarType>(nVar, matrix_a, matrix_b, product);
}
#else
MATVECPROD_SIGNATURE(MatrixVectorProduct) {
  /*--- With MKL we use the just-in-time kernels instead of the naive implementation. ---*/
  MatrixVectorProductKernelBetaZero(MatrixVectorProductJitterBetaZero, const_cast<ScalarType*>(vector),
                                    const_cast<ScalarType*>(matrix), product);
}

MATVECPROD_SIGNATURE(MatrixVectorProductAdd) {
  MatrixVectorProductKernelBetaOne(MatrixVectorProductJitterBetaOne, const_cast<ScalarType*>(vector),
                                   const_cast<ScalarType*>(matrix), product);
}

MATVECPROD_SIGNATURE(MatrixVectorProductSub) {
  MatrixVectorProductKernelAlphaMinusOne(MatrixVectorProductJitterAlphaMinusOne, const_cast<ScalarType*>(vector),
                                         const_cast<ScalarType*>(matrix), product);
}

template <class ScalarType>
FORCEINLINE void CSysMatrix<ScalarType>::MatrixMatrixProduct(const ScalarType* matrix_a, const ScalarType* matrix_b,
                                                             ScalarType* product) const {
  MatrixMatrixProductKernel(MatrixMatrixProductJitter, const_cast<ScalarType*>(matrix_a),
                            const_cast<ScalarType*>(matrix_b), product);
}
#endif

#undef MATVECPROD_SIGNATURE
#undef __MATVECPROD_SIGNATURE__

template <class ScalarType>
FORCEINLINE void CSysMatrix<ScalarType>::Gauss_Elimination(unsigned long block_i, ScalarType* rhs) const {
  /*--- Copy block, as the algorithm modifies the matrix ---*/
  ScalarType block[MAXNVAR * MAXNVAR];
  MatrixCopy(&matrix[dia_ptr[block_i] * nVar * nVar], block);

  Gauss_Elimination(block, rhs);
}

template <class ScalarType>
FORCEINLINE void CSysMatrix<ScalarType>::InverseDiagonalBlock(unsigned long block_i, ScalarType* invBlock) const {
  /*--- Copy block, as the algorithm modifies the matrix ---*/
  ScalarType block[MAXNVAR * MAXNVAR];
  MatrixCopy(&matrix[dia_ptr[block_i] * nVar * nVar], block);

  MatrixInverse(block, invBlock);
}

template <class ScalarType>
FORCEINLINE void CSysMatrix<ScalarType>::InverseDiagonalBlock_ILUMatrix(unsigned long block_i,
                                                                        ScalarType* invBlock) const {
  /*--- Copy block, as the algorithm modifies the matrix ---*/
  ScalarType block[MAXNVAR * MAXNVAR];
  MatrixCopy(&ILU_matrix[dia_ptr_ilu[block_i] * nVar * nVar], block);

  MatrixInverse(block, invBlock);
}

template <class ScalarType>
FORCEINLINE void CSysMatrix<ScalarType>::RowProduct(const CSysVector<ScalarType>& vec, unsigned long row_i,
                                                    ScalarType* prod) const {
  for (auto iVar = 0ul; iVar < nVar; iVar++) prod[iVar] = 0.0;

  for (auto index = row_ptr[row_i]; index < row_ptr[row_i + 1]; index++) {
    auto col_j = col_ind[index];
    MatrixVectorProductAdd(&matrix[index * nVar * nEqn], &vec[col_j * nEqn], prod);
  }
}

template <class ScalarType>
FORCEINLINE void CSysMatrix<ScalarType>::UpperProduct(const CSysVector<ScalarType>& vec, unsigned long row_i,
                                                      unsigned long col_ub, ScalarType* prod) const {
  for (auto iVar = 0ul; iVar < nVar; iVar++) prod[iVar] = 0.0;

  for (auto index = dia_ptr[row_i] + 1; index < row_ptr[row_i + 1]; index++) {
    auto col_j = col_ind[index];
    /*--- Always include halos. ---*/
    if (col_j < col_ub || col_j >= nPointDomain)
      MatrixVectorProductAdd(&matrix[index * nVar * nEqn], &vec[col_j * nEqn], prod);
  }
}

template <class ScalarType>
FORCEINLINE void CSysMatrix<ScalarType>::LowerProduct(const CSysVector<ScalarType>& vec, unsigned long row_i,
                                                      unsigned long col_lb, ScalarType* prod) const {
  for (auto iVar = 0ul; iVar < nVar; iVar++) prod[iVar] = 0.0;

  for (auto index = row_ptr[row_i]; index < dia_ptr[row_i]; index++) {
    auto col_j = col_ind[index];
    if (col_j >= col_lb) MatrixVectorProductAdd(&matrix[index * nVar * nEqn], &vec[col_j * nEqn], prod);
  }
}

template <class ScalarType>
FORCEINLINE void CSysMatrix<ScalarType>::DiagonalProduct(const CSysVector<ScalarType>& vec, unsigned long row_i,
                                                         ScalarType* prod) const {
  MatrixVectorProduct(&matrix[dia_ptr[row_i] * nVar * nEqn], &vec[row_i * nEqn], prod);
}
