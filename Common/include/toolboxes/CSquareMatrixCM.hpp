/*!
 * \file CSquareMatrixCM.hpp
 * \brief Dense general square matrix, used for example in DG standard elements
 *        in Column Major order storage.
 * \author Edwin van der Weide, Pedro Gomes.
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
 * \brief Class to store a dense general square matrix that uses the Column
 *        Major order storage format. The code should be compiled with
 *        LAPACK to use optimized matrix inversion and multiplication routines.
 * \ingroup BLAS
 */
class CSquareMatrixCM {
  static_assert(ColMajorMatrix<passivedouble>::Storage == StorageType::ColumnMajor,
                "Column major storage is assumed for LAPACK.");

 private:
  ColMajorMatrix<passivedouble> mat; /*!< \brief Storage of the actual matrix. */

 public:
  /*!
   * \brief Default constructor. Nothing to be done.
   */
  CSquareMatrixCM() = default;

  /*!
   * \overload
   * \brief Overloaded constructor, which allocates the memory to store
   *        the matrix.
   * \param[in] N - Number of rows and colums of the matrix.
   */
  CSquareMatrixCM(int N) { Initialize(N); }

  /*!
   * \brief Operator, which makes available the given matrix element as a reference.
   * \param[in] i   - Row index of the matrix element.
   * \param[in] j   - Column index of the matrix element.
   * \return          Reference to element (i,j).
   */
  inline passivedouble& operator()(int i, int j) { return mat(i, j); }

  /*!
   * \brief Operator, which makes available the given matrix element as a const reference.
   * \param[in] i   - Row index of the matrix element.
   * \param[in] j   - Column index of the matrix element.
   * \return          Constant reference to element (i,j).
   */
  inline const passivedouble& operator()(int i, int j) const { return mat(i, j); }

  /*!
   * \brief Function, which makes available a reference to the actual matrix.
   * \return A reference to mat.
   */
  inline ColMajorMatrix<passivedouble>& GetMat() { return mat; }

  /*!
   * \brief Function, which makes available a const reference to the actual matrix.
   * \return A const reference to mat.
   */
  inline const ColMajorMatrix<passivedouble>& GetMat() const { return mat; }

  /*!
   * \brief Function, which allocates the memory for the matrix.
   * \param[in] N - Number of rows and colums of the matrix.
   */
  inline void Initialize(int N) { mat.resize(N, N); }

  /*!
   * \brief Function, which makes available the size of the matrix.
   * \return The number of rows, columns of the matrix.
   */
  inline int Size() const { return mat.rows(); }

  /*!
   * \brief Function, which carries out the matrix produc of the current matrix
   *        with mat_in and stores the result in mat_out.
   * \param[in]  side    - left: mat_out = this * mat_in, right: mat_out = mat_in * this
   * \param[in]  mat_in  - Matrix to be multiplied by the current matrix.
   * \param[out] mat_out - Matrix to store the result of the multiplication.
   */
  void MatMatMult(const char side, const ColMajorMatrix<passivedouble>& mat_in,
                  ColMajorMatrix<passivedouble>& mat_out) const;

  /*!
   * \brief Naive matrix-vector multiplication with general type.
   */
  template <class ForwardIt>
  void MatVecMult(ForwardIt vec_in, ForwardIt vec_out) const {
    for (int i = 0; i < Size(); ++i) {
      *vec_out = 0.0;
      auto vec = vec_in;
      for (int k = 0; k < Size(); ++k) *vec_out += *(vec++) * mat(i, k);
      ++vec_out;
    }
  }

  /*!
   * \brief Function, which inverts the matrix in-place.
   */
  void Invert();

  /*!
   * \brief Function, which transposes the matrix in-place.
   */
  void Transpose();
};
