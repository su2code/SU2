/*!
 * \file CSysMatrix.inl
 * \brief Inline subroutines of the <i>CSysMatrix.hpp</i> file.
 * \note These are the "private" inlines, they are not needed outside
 *       of the .cpp file and so they are hidden to avoid triggering
 *       recompilation of other units when changes are made here.
 * \author F. Palacios, A. Bueno, T. Economon, P. Gomes
 * \version 8.5.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md)
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
/*--- Custom 8-bit float for per-row scales: layout [eeeee|sss], value = (1 + s/8) * 2^(-e).
 *    Exponent is always non-negative so all values are in (0, ~1.875].
 *    No sign, no zero — used to encode ratios in (0, 1] relative to the block peak. ---*/
inline uint8_t float8_encode(double v) {
  if (v <= 0.0) return static_cast<uint8_t>(31 << 3);  // clamp to minimum
  int fe;
  const double m = std::frexp(v, &fe);
  int e = 1 - fe;
  int s = static_cast<int>(std::round((2.0 * m - 1.0) * 8.0));
  if (s >= 8) {
    --e;
    s = 0;
  }  // carry: move to next higher power-of-two band (larger value → smaller e)
  return static_cast<uint8_t>((std::max(0, std::min(31, e)) << 3) | std::max(0, std::min(7, s)));
}

inline double float8_decode(uint8_t b) { return std::ldexp(1.0 + (b & 7) * 0.125, -(b >> 3)); }

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
  if (USE_QUANTIZATION && nVar > 1) {
    const auto k = dia_ptr[block_i];
    const ScalarType bscale = q_bscale[k];
    const uint8_t* qscale = &q_scale[k * nVar];
    const QuantType* q = &q_offdiag[k * nVar * nVar];
    for (auto r = 0ul; r < nVar; ++r) {
      const ScalarType row_scale = bscale * ScalarType(float8_decode(qscale[r]));
      for (auto c = 0ul; c < nVar; ++c) block[r * nVar + c] = row_scale * static_cast<ScalarType>(q[r * nVar + c]);
    }
  } else {
    MatrixCopy(&matrix[dia_ptr[block_i] * nVar * nVar], block);
  }
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
FORCEINLINE const ScalarType* CSysMatrix<ScalarType>::InvertDiagonalBlockILUMatrix(unsigned long block_i) {
  /*--- Copy block, as the algorithm modifies the matrix ---*/
  auto* Uii = &ILU_matrix[dia_ptr_ilu[block_i] * nVar * nVar];
  ScalarType block[MAXNVAR * MAXNVAR];
  MatrixCopy(Uii, block);
  MatrixInverse(block, Uii);
  return Uii;
}

template <class ScalarType>
FORCEINLINE void CSysMatrix<ScalarType>::QuantizedMatVecAdd(unsigned long k, const ScalarType* vec,
                                                            ScalarType* prod) const {
  const ScalarType bscale = q_bscale[k];
  const uint8_t* qscale = &q_scale[k * nVar];
  const QuantType* q = &q_offdiag[k * nVar * nVar];
  for (auto r = 0ul; r < nVar; ++r) {
    const ScalarType row_scale = bscale * ScalarType(float8_decode(qscale[r]));
    auto sum = ScalarType(0);
    for (auto c = 0ul; c < nVar; ++c) sum += static_cast<ScalarType>(q[r * nVar + c]) * vec[c];
    prod[r] += row_scale * sum;
  }
}

template <class ScalarType>
FORCEINLINE void CSysMatrix<ScalarType>::RowProduct(const CSysVector<ScalarType>& vec, unsigned long row_i,
                                                    ScalarType* prod) const {
  for (auto iVar = 0ul; iVar < nVar; iVar++) prod[iVar] = 0.0;

  for (auto index = row_ptr[row_i]; index < row_ptr[row_i + 1]; index++) {
    auto col_j = col_ind[index];
    if (USE_QUANTIZATION && nVar > 1) {
      QuantizedMatVecAdd(index, &vec[col_j * nEqn], prod);
    } else {
      MatrixVectorProductAdd(&matrix[index * nVar * nEqn], &vec[col_j * nEqn], prod);
    }
  }
}

template <class ScalarType>
FORCEINLINE void CSysMatrix<ScalarType>::UpperProduct(const CSysVector<ScalarType>& vec, unsigned long row_i,
                                                      unsigned long col_ub, ScalarType* prod) const {
  for (auto iVar = 0ul; iVar < nVar; iVar++) prod[iVar] = 0.0;

  for (auto index = dia_ptr[row_i] + 1; index < row_ptr[row_i + 1]; index++) {
    auto col_j = col_ind[index];
    /*--- Always include halos. ---*/
    if (col_j < col_ub || col_j >= nPointDomain) {
      if (USE_QUANTIZATION && nVar > 1) {
        QuantizedMatVecAdd(index, &vec[col_j * nEqn], prod);
      } else {
        MatrixVectorProductAdd(&matrix[index * nVar * nEqn], &vec[col_j * nEqn], prod);
      }
    }
  }
}

template <class ScalarType>
FORCEINLINE void CSysMatrix<ScalarType>::LowerProduct(const CSysVector<ScalarType>& vec, unsigned long row_i,
                                                      unsigned long col_lb, ScalarType* prod) const {
  for (auto iVar = 0ul; iVar < nVar; iVar++) prod[iVar] = 0.0;

  for (auto index = row_ptr[row_i]; index < dia_ptr[row_i]; index++) {
    auto col_j = col_ind[index];
    if (col_j >= col_lb) {
      if (USE_QUANTIZATION && nVar > 1) {
        QuantizedMatVecAdd(index, &vec[col_j * nEqn], prod);
      } else {
        MatrixVectorProductAdd(&matrix[index * nVar * nEqn], &vec[col_j * nEqn], prod);
      }
    }
  }
}

template <class ScalarType>
FORCEINLINE void CSysMatrix<ScalarType>::DiagonalProduct(const CSysVector<ScalarType>& vec, unsigned long row_i,
                                                         ScalarType* prod) const {
  if (USE_QUANTIZATION && nVar > 1) {
    for (auto iVar = 0ul; iVar < nVar; iVar++) prod[iVar] = 0.0;
    QuantizedMatVecAdd(dia_ptr[row_i], &vec[row_i * nEqn], prod);
  } else {
    MatrixVectorProduct(&matrix[dia_ptr[row_i] * nVar * nEqn], &vec[row_i * nEqn], prod);
  }
}
