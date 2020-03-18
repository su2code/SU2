/*!
 * \file CSymmetricMatrix.hpp
 * \brief Dense symmetric matrix, used for example in RBF interpolation.
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
#pragma once

#include <vector>
#include "C2DContainer.hpp"

using namespace std;

/*!
 * \brief The matrix is symmetric but full storage is used as that gives much better
 * performance for some BLAS libraries (notably OpenBLAS). The code should be compiled
 * with LAPACK to use optimized matrix inversion and multiplication routines.
 */
class CSymmetricMatrix {
private:
  enum DecompositionType { NONE, CHOLESKY, LU };

  vector<passivedouble> val_vec, decomp_vec;
  vector<int> perm_vec;
  int sz = 0;
  bool initialized = false;
  DecompositionType decomposed = NONE;

  inline void CheckBounds(int i, int j) const {
    assert(initialized && "Matrix not initialized.");
    assert(i>=0 && i<sz && j>=0 && j<sz && "Index to access matrix out of bounds.");
  }

  inline int IdxFull(int i, int j) const {CheckBounds(i,j); return i*sz + j;}

  inline int IdxSym(int i, int j) const {return IdxFull(min(i,j), max(i,j));}

  inline passivedouble& decomp(int i, int j) { return decomp_vec[IdxFull(i,j)]; }

  // Not optimized dense matrix factorization and inversion for portability.
  void CholeskyDecompose();
  void LUDecompose();
  void CalcInv();
  // Matrix inversion using LAPACK routines (LDLT and LLT factorization).
  void CalcInv_sytri();
  void CalcInv_potri();

public:
  CSymmetricMatrix() = default;
  CSymmetricMatrix(int N) {Initialize(N);}

  void Initialize(int N);

  inline int GetSize() const { return sz; }

  inline passivedouble Get(int i, int j) const { return val_vec[IdxSym(i,j)]; }

  inline void Set(int i, int j, passivedouble val) { val_vec[IdxSym(i,j)] = val; }

  inline passivedouble& operator() (int i, int j) { return val_vec[IdxSym(i,j)]; }

  inline const passivedouble& operator() (int i, int j) const { return val_vec[IdxSym(i,j)]; }

  void MatVecMult(passivedouble *v) const;

  void MatMatMult(const char side, su2passivematrix& mat_in, su2passivematrix& mat_out);

  void Invert(const bool is_spd);

};
