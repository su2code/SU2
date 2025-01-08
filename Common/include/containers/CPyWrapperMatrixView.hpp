/*!
 * \file CPyWrapperMatrixView.hpp
 * \brief Simple matrix views to use with the python wrapper.
 * \author P. Gomes
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

#include <string>
#include <utility>
#include <vector>

#include "../geometry/dual_grid/CVertex.hpp"
#include "../parallelization/mpi_structure.hpp"
#include "C2DContainer.hpp"

/*!
 * \brief Python wrapper matrix interface
 * The accessors in this macro provide a passive double interface to matrix of
 * su2double. This can be extended to allow access to derivative information.
 * Classes that use this macro encapsulate the access to the raw data (su2double)
 * via the functions "Access(row, col) -> su2double&" (const and non-const versions).
 * We use a macro because it is diffult to use modern C++ idioms (e.g. CRTP) with
 * SWIG. In addition to "Access" classes must have member variables "rows_", "cols_",
 * and "read_only_".
 */
#define PY_WRAPPER_MATRIX_INTERFACE                                                                              \
  /*! \brief Returns the shape of the matrix. */                                                                 \
  std::pair<unsigned long, unsigned long> Shape() const { return std::make_pair(rows_, cols_); }                 \
                                                                                                                 \
  /*! \brief Returns whether the data is read-only [true] or if it can be modified [false]. */                   \
  bool IsReadOnly() const { return read_only_; }                                                                 \
                                                                                                                 \
  /*! \brief Gets the value for a (row, column) pair. */                                                         \
  passivedouble operator()(unsigned long row, unsigned long col) const { return Get(row, col); }                 \
                                                                                                                 \
  /*! \brief Gets the value for a (row, column) pair. */                                                         \
  passivedouble Get(unsigned long row, unsigned long col) const { return SU2_TYPE::GetValue(Access(row, col)); } \
                                                                                                                 \
  /*! \brief Gets the values for a row of the matrix. */                                                         \
  std::vector<passivedouble> Get(unsigned long row) const {                                                      \
    std::vector<passivedouble> vals(cols_);                                                                      \
    for (unsigned long j = 0; j < cols_; ++j) vals[j] = Get(row, j);                                             \
    return vals;                                                                                                 \
  }                                                                                                              \
  /*! \brief Sets the value for a (row, column) pair. This clears derivative information. */                     \
  void Set(unsigned long row, unsigned long col, passivedouble val) { Access(row, col) = val; }                  \
                                                                                                                 \
  /*! \brief Sets the values for a row of the matrix. */                                                         \
  void Set(unsigned long row, std::vector<passivedouble> vals) {                                                 \
    unsigned long j = 0;                                                                                         \
    for (const auto& val : vals) Set(row, j++, val);                                                             \
  }

/*!
 * \class CPyWrapperMatrixView
 * \ingroup PySU2
 * \brief This class wraps su2activematrix for the python wrapper matrix interface.
 * It is generaly used to wrap access to solver variables defined for the entire volume.
 */
class CPyWrapperMatrixView {
 protected:
  static_assert(su2activematrix::IsRowMajor, "");
  su2double* data_ = nullptr;
  unsigned long rows_ = 0, cols_ = 0;
  std::string name_;
  bool read_only_ = false;

  /*--- Define the functions required by the interface macro. ---*/
  inline const su2double& Access(unsigned long row, unsigned long col) const {
    if (row > rows_ || col > cols_) SU2_MPI::Error(name_ + " out of bounds", "CPyWrapperMatrixView");
    return data_[row * cols_ + col];
  }
  inline su2double& Access(unsigned long row, unsigned long col) {
    if (read_only_) SU2_MPI::Error(name_ + " is read-only", "CPyWrapperMatrixView");
    const auto& const_me = *this;
    return const_cast<su2double&>(const_me.Access(row, col));
  }

 public:
  CPyWrapperMatrixView() = default;

  /*!
   * \brief Construct the view of the matrix.
   * \note "name" should be set to the variable name being returned to give better information to users.
   * \note "read_only" can be set to true to prevent the data from being modified.
   */
  CPyWrapperMatrixView(su2activematrix& mat, const std::string& name, bool read_only)
      : data_(mat.data()), rows_(mat.rows()), cols_(mat.cols()), name_(name), read_only_(read_only) {}

  /*--- Use the macro to generate the interface. ---*/
  PY_WRAPPER_MATRIX_INTERFACE
};

/*!
 * \class CPyWrapperMarkerMatrixView
 * \ingroup PySU2
 * \brief This class wraps su2activematrix for the python wrapper matrix interface restricting it
 * to the vertices of a given marker.
 */
class CPyWrapperMarkerMatrixView {
 private:
  static_assert(su2activematrix::IsRowMajor, "");
  su2double* data_ = nullptr;
  const CVertex* const* vertices_ = nullptr;
  unsigned long rows_ = 0, cols_ = 0;
  std::string name_;
  bool read_only_ = false;

  /*--- Define the functions required by the interface macro. ---*/
  inline const su2double& Access(unsigned long row, unsigned long col) const {
    if (row > rows_ || col > cols_) SU2_MPI::Error(name_ + " out of bounds", "CPyWrapperMarkerMatrixView");
    return data_[vertices_[row]->GetNode() * cols_ + col];
  }
  inline su2double& Access(unsigned long row, unsigned long col) {
    if (read_only_) SU2_MPI::Error(name_ + " is read-only", "CPyWrapperMarkerMatrixView");
    const auto& const_me = *this;
    return const_cast<su2double&>(const_me.Access(row, col));
  }

 public:
  CPyWrapperMarkerMatrixView() = default;

  /*!
   * \brief Construct the view of the matrix.
   * \note "name" should be set to the variable name being returned to give better information to users.
   * \note "read_only" can be set to true to prevent the data from being modified.
   */
  CPyWrapperMarkerMatrixView(su2activematrix& mat, const CVertex* const* vertices, unsigned long n_vertices,
                             const std::string& name, bool read_only)
      : data_(mat.data()),
        vertices_(vertices),
        rows_(n_vertices),
        cols_(mat.cols()),
        name_(name),
        read_only_(read_only) {
    if (mat.rows() < n_vertices) {
      SU2_MPI::Error(name + " has fewer rows than the number of vertices in the marker.", "CPyWrapperMarkerMatrixView");
    }
  }

  /*--- Use the macro to generate the interface. ---*/
  PY_WRAPPER_MATRIX_INTERFACE
};
