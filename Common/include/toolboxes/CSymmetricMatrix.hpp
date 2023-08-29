/*!
 * \file CSymmetricMatrix.hpp
 * \brief Dense symmetric matrix, used for example in RBF interpolation.
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
#pragma once

#include <vector>
#include "../containers/C2DContainer.hpp"

/*!
 * \brief The matrix is symmetric but full storage is used as that gives much better
 * performance for some BLAS libraries (notably OpenBLAS). The code should be compiled
 * with LAPACK to use optimized matrix inversion and multiplication routines.
 * \ingroup BLAS
 */
class CSymmetricMatrix {
  static_assert(su2passivematrix::IsRowMajor, "Row major storage is assumed for LAPACK.");

 private:
  su2passivematrix mat;

  // Not optimized dense matrix factorization and inversion for portability.
  void CalcInv(bool is_spd);
  void CholeskyDecompose();
  // Matrix inversion using LAPACK routines (LDLT and LLT factorization).
  void CalcInv_sytri();
  void CalcInv_potri();

 public:
  CSymmetricMatrix() = default;
  CSymmetricMatrix(int N) { Initialize(N); }

  void Initialize(int N);

  inline int Size() const { return mat.rows(); }

  inline passivedouble Get(int i, int j) const { return mat(std::min(i, j), std::max(i, j)); }

  inline void Set(int i, int j, passivedouble val) { mat(std::min(i, j), std::max(i, j)) = val; }

  inline passivedouble& operator()(int i, int j) { return mat(std::min(i, j), std::max(i, j)); }

  inline const passivedouble& operator()(int i, int j) const { return mat(std::min(i, j), std::max(i, j)); }

  template <class ForwardIt>
  void MatVecMult(ForwardIt vec_in, ForwardIt vec_out) const {
    for (int i = 0; i < Size(); ++i) {
      *vec_out = 0.0;
      auto vec = vec_in;
      for (int k = 0; k < Size(); ++k) *vec_out += *(vec++) * Get(i, k);
      ++vec_out;
    }
  }

  void MatMatMult(const char side, const su2passivematrix& mat_in, su2passivematrix& mat_out) const;

  void Invert(bool is_spd = false);

  su2passivematrix StealData();
};
