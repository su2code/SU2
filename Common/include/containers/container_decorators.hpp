/*!
 * \file container_decorators.hpp
 * \brief Collection of small classes that decorate C2DContainer to
 * augment its functionality, e.g. give it extra dimensions.
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

#include "C2DContainer.hpp"

/// \addtogroup Containers
/// @{

/*!
 * \brief Class to represent a matrix (without owning the data, this just wraps a pointer).
 */
template <class T>
class CMatrixView {
 public:
  using Scalar = typename std::remove_const<T>::type;
  using Index = unsigned long;

 private:
  T* m_ptr;
  Index m_cols;

 public:
  CMatrixView(T* ptr = nullptr, Index cols = 0) : m_ptr(ptr), m_cols(cols) {}

  template <class U>
  friend class CMatrixView;
  template <class U>
  CMatrixView(const CMatrixView<U>& other) : m_ptr(other.m_ptr), m_cols(other.m_cols) {}

  explicit CMatrixView(su2matrix<Scalar>& mat) : m_ptr(mat.data()), m_cols(mat.cols()) {}

  template <class U = T, su2enable_if<std::is_const<U>::value> = 0>
  explicit CMatrixView(const su2matrix<Scalar>& mat) : m_ptr(mat.data()), m_cols(mat.cols()) {}

  const Scalar* operator[](Index i) const noexcept { return &m_ptr[i * m_cols]; }
  const Scalar& operator()(Index i, Index j) const noexcept { return m_ptr[i * m_cols + j]; }

  template <class U = T, su2enable_if<!std::is_const<U>::value> = 0>
  Scalar* operator[](Index i) noexcept {
    return &m_ptr[i * m_cols];
  }

  template <class U = T, su2enable_if<!std::is_const<U>::value> = 0>
  Scalar& operator()(Index i, Index j) noexcept {
    return m_ptr[i * m_cols + j];
  }

  friend CMatrixView operator+(CMatrixView mv, Index incr) { return CMatrixView(mv[incr], mv.m_cols); }
};

/*!
 * \class C3DContainerDecorator
 * \brief Decorate a matrix type (Storage) with 3 dimensions.
 */
template <class Storage>
class C3DContainerDecorator {
  static_assert(!Storage::IsVector, "Storage type must be a matrix.");
  static_assert(Storage::IsRowMajor, "Storage type must be row major.");

 public:
  using Scalar = typename Storage::Scalar;
  using Index = typename Storage::Index;
  static constexpr bool IsRowMajor = true;
  static constexpr bool IsColumnMajor = false;

  using Matrix = CMatrixView<Scalar>;
  using ConstMatrix = CMatrixView<const Scalar>;

  using CInnerIter = typename Storage::CInnerIter;
  template <class T, size_t N>
  using CInnerIterGather = typename Storage::template CInnerIterGather<simd::Array<T, N> >;

 private:
  Storage m_storage;
  Index m_innerSz;

 public:
  C3DContainerDecorator() = default;

  C3DContainerDecorator(Index length, Index rows, Index cols, Scalar value = 0) noexcept {
    resize(length, rows, cols, value);
  }

  void resize(Index length, Index rows, Index cols, Scalar value = 0) noexcept {
    m_innerSz = cols;
    m_storage.resize(length, rows * cols) = value;
  }

  /*!
   * \brief Container sizes.
   */
  Index size() const noexcept { return m_storage.size(); }
  Index length() const noexcept { return m_storage.rows(); }
  Index rows() const noexcept { return m_storage.cols() / m_innerSz; }
  Index cols() const noexcept { return m_innerSz; }

  /*!
   * \brief Element-wise access.
   */
  Scalar& operator()(Index i, Index j, Index k) noexcept { return m_storage(i, j * m_innerSz + k); }
  const Scalar& operator()(Index i, Index j, Index k) const noexcept { return m_storage(i, j * m_innerSz + k); }

  /*!
   * \brief Matrix access.
   */
  Matrix operator[](Index i) noexcept { return Matrix(m_storage[i], m_innerSz); }
  ConstMatrix operator[](Index i) const noexcept { return ConstMatrix(m_storage[i], m_innerSz); }

  /*!
   * \brief Matrix access with an offset.
   */
  Matrix operator()(Index i, Index j) noexcept { return Matrix(m_storage[i] + j * m_innerSz, m_innerSz); }
  ConstMatrix operator()(Index i, Index j) const noexcept {
    return ConstMatrix(m_storage[i] + j * m_innerSz, m_innerSz);
  }

  /*!
   * \brief Get a scalar iterator to the inner-most dimension of the container.
   */
  FORCEINLINE CInnerIter innerIter(Index i, Index j) const noexcept {
    return CInnerIter(&m_storage(i, j * m_innerSz), 1);
  }

  /*!
   * \brief Get a SIMD gather iterator to the inner-most dimension of the container.
   */
  template <class T, size_t N>
  FORCEINLINE CInnerIterGather<T, N> innerIter(simd::Array<T, N> i, Index j) const noexcept {
    return CInnerIterGather<T, N>(m_storage.data(), 1, i * m_storage.cols() + j * m_innerSz);
  }

  /*!
   * \brief Return copy of data in a static size container (see C2DContainer::get).
   * \param[in] i - Outer index.
   * \param[in] j - Starting middle index for the copy (amount determined by container size).
   */
  template <class StaticContainer, class Int>
  FORCEINLINE StaticContainer get(Int i, Index j = 0) const noexcept {
    return m_storage.template get<StaticContainer>(i, j * m_innerSz);
  }
};

/*!
 * \brief Some typedefs for 3D containers
 */
using C3DIntMatrix = C3DContainerDecorator<su2matrix<unsigned long> >;
using C3DDoubleMatrix = C3DContainerDecorator<su2activematrix>;
using CVectorOfMatrix = C3DDoubleMatrix;

/*!
 * \class C2DDummyLastView
 * \brief Helper class, adds dummy trailing dimension to a reference of a
 *        vector object making it a dummy matrix.
 * \note The constness of the object is derived from the template type, but
 *       we allways keep a reference, never a copy of the associated vector.
 */
template <class T>
struct C2DDummyLastView {
  static_assert(T::IsVector, "This class decorates vectors.");
  using Index = typename T::Index;
  using Scalar = typename T::Scalar;

  T& data;

  C2DDummyLastView() = delete;

  C2DDummyLastView(T& ref) : data(ref) {}

  template <class U = T, su2enable_if<!std::is_const<U>::value> = 0>
  Scalar& operator()(Index i, Index) noexcept {
    return data(i);
  }

  const Scalar& operator()(Index i, Index) const noexcept { return data(i); }
};

/*!
 * \class C3DDummyMiddleView
 * \brief Helper class, adds dummy middle dimension to a reference of a
 *        matrix object making it a dummy 3D array.
 * \note The constness of the object is derived from the template type, but
 *       we allways keep a reference, never a copy of the associated matrix.
 */
template <class T>
struct C3DDummyMiddleView {
  static_assert(!T::IsVector, "This class decorates matrices.");
  using Index = typename T::Index;
  using Scalar = typename T::Scalar;

  T& data;

  C3DDummyMiddleView() = delete;

  C3DDummyMiddleView(T& ref) : data(ref) {}

  template <class U = T, su2enable_if<!std::is_const<U>::value> = 0>
  Scalar& operator()(Index i, Index, Index k) noexcept {
    return data(i, k);
  }

  const Scalar& operator()(Index i, Index, Index k) const noexcept { return data(i, k); }
};

/*--- Helper functions to allocate containers of containers. ---*/

/*!
 * \brief Allocate a vector of varying-size vectors and initialize with some value.
 * \param[in] M - the first index is >=0 and < M
 * \param[in] N - the second index is >=0 and < N[first index]
 * \param[in,out] X - the vector of vectors
 * \param[in] val - the value for initialization, default is 0
 * \tparam IndexVector - type of N
 * \tparam VectorOfVector - type of X
 */
template <class IndexVector, class VectorOfVector, class Scalar = int>
inline void AllocVectorOfVectors(size_t M, const IndexVector& N, VectorOfVector& X, Scalar val = 0) {
  X.resize(M);
  for (size_t i = 0; i < M; ++i) {
    X[i].resize(N[i]);
    for (auto& x : X[i]) x = val;
  }
}

/*!
 * \overload Deduce outer size from index vector.
 */
template <class IndexVector, class VectorOfVector, class Scalar = int>
inline void AllocVectorOfVectors(const IndexVector& N, VectorOfVector& X, Scalar val = 0) {
  auto M = N.size();
  AllocVectorOfVectors(M, N, X, val);
}

/*!
 * \brief Allocate a vector of matrices with varying row count, and initialize with some value.
 * \param[in] M - the first index is >=0 and < M
 * \param[in] N - the second index is >=0 and < N[first index]
 * \param[in] P - the third index is >=0 and < P
 * \param[in,out] X - the vector of matrices
 * \param[in] val - the value for initialization, default is 0
 * \tparam IndexVector - type of N
 * \tparam VectorOfMatrix - type of X
 */
template <class IndexVector, class VectorOfMatrix, class Scalar = int>
inline void AllocVectorOfMatrices(size_t M, const IndexVector& N, size_t P, VectorOfMatrix& X, Scalar val = 0) {
  X.resize(M);
  for (size_t i = 0; i < M; ++i) {
    X[i].resize(N[i], P);
    for (auto& x : X[i]) x = val;
  }
}

/*!
 * \overload Deduce outer size from index vector.
 */
template <class IndexVector, class VectorOfMatrix, class Scalar = int>
inline void AllocVectorOfMatrices(const IndexVector& N, size_t P, VectorOfMatrix& X, Scalar val = 0) {
  auto M = N.size();
  AllocVectorOfMatrices(M, N, P, X, val);
}

/// @}
