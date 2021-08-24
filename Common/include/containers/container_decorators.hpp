/*!
 * \file container_decorators.hpp
 * \brief Collection of small classes that decorate C2DContainer to
 * augment its functionality, e.g. give it extra dimensions.
 * \author P. Gomes
 * \version 7.2.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2021, SU2 Contributors (cf. AUTHORS.md)
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

/*!
 * \class C3DContainerDecorator
 * \brief Decorate a matrix type (Storage) with 3 dimensions.
 */
template<class Storage>
class C3DContainerDecorator {
  static_assert(!Storage::IsVector, "Storage type must be a matrix.");
  static_assert(Storage::IsRowMajor, "Storage type must be row major.");
public:
  using Scalar = typename Storage::Scalar;
  using Index = typename Storage::Index;
  static constexpr bool IsRowMajor = true;
  static constexpr bool IsColumnMajor = false;

  using CInnerIter = typename Storage::CInnerIter;
  template<class T, size_t N>
  using CInnerIterGather = typename Storage::template CInnerIterGather<simd::Array<T,N> >;

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
    m_storage.resize(length, rows*cols) = value;
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
  Scalar& operator() (Index i, Index j, Index k) noexcept { return m_storage(i, j*m_innerSz + k); }
  const Scalar& operator() (Index i, Index j, Index k) const noexcept { return m_storage(i, j*m_innerSz + k); }

  /*!
   * \brief Get a scalar iterator to the inner-most dimension of the container.
   */
  FORCEINLINE CInnerIter innerIter(Index i, Index j) const noexcept {
    return CInnerIter(&m_storage(i, j*m_innerSz), 1);
  }

  /*!
   * \brief Get a SIMD gather iterator to the inner-most dimension of the container.
   */
  template<class T, size_t N>
  FORCEINLINE CInnerIterGather<T,N> innerIter(simd::Array<T,N> i, Index j) const noexcept {
    return CInnerIterGather<T,N>(m_storage.data(), 1, i*m_storage.cols() + j*m_innerSz);
  }

  /*!
   * \brief Return copy of data in a static size container (see C2DContainer::get).
   * \param[in] i - Outer index.
   * \param[in] j - Starting middle index for the copy (amount determined by container size).
   */
  template<class StaticContainer, class Int>
  FORCEINLINE StaticContainer get(Int i, Index j = 0) const noexcept {
    return m_storage.template get<StaticContainer>(i, j*m_innerSz);
  }
};

/*!
 * \brief Some typedefs for the
 */
using C3DIntMatrix = C3DContainerDecorator<su2matrix<unsigned long> >;
using C3DDoubleMatrix = C3DContainerDecorator<su2activematrix>;

/*!
 * \class CVectorOfMatrix
 * \brief This contrived container is used to store small matrices in a contiguous manner
 *        but still present the "su2double**" interface to the outside world.
 *        The "interface" part should be replaced by something more efficient, e.g. a "matrix view".
 */
class CVectorOfMatrix: public C3DDoubleMatrix {
private:
  su2matrix<Scalar*> interface;

public:
  CVectorOfMatrix() = default;

  CVectorOfMatrix(Index length, Index rows, Index cols, Scalar value = 0) noexcept {
    resize(length, rows, cols, value);
  }

  void resize(Index length, Index rows, Index cols, Scalar value = 0) noexcept {
    C3DDoubleMatrix::resize(length, rows, cols, value);
    interface.resize(length,rows);
    for(Index i=0; i<length; ++i)
      for(Index j=0; j<rows; ++j)
        interface(i,j) = &(*this)(i,j,0);
  }

  /*!
   * \brief Matrix-wise access.
   */
  Scalar** operator[] (Index i) noexcept { return interface[i]; }
  const Scalar* const* operator[] (Index i) const noexcept { return interface[i]; }
};

/*!
 * \class C2DDummyLastView
 * \brief Helper class, adds dummy trailing dimension to a reference of a
 *        vector object making it a dummy matrix.
 * \note The constness of the object is derived from the template type, but
 *       we allways keep a reference, never a copy of the associated vector.
 */
template<class T>
struct C2DDummyLastView
{
  static_assert(T::IsVector, "This class decorates vectors.");
  using Index = typename T::Index;
  using Scalar = typename T::Scalar;

  T& data;

  C2DDummyLastView() = delete;

  C2DDummyLastView(T& ref) : data(ref) {}

  template<class U = T, su2enable_if<!std::is_const<U>::value> = 0>
  Scalar& operator() (Index i, Index) noexcept
  {
    return data(i);
  }

  const Scalar& operator() (Index i, Index) const noexcept
  {
    return data(i);
  }
};

/*!
 * \class C3DDummyMiddleView
 * \brief Helper class, adds dummy middle dimension to a reference of a
 *        matrix object making it a dummy 3D array.
 * \note The constness of the object is derived from the template type, but
 *       we allways keep a reference, never a copy of the associated matrix.
 */
template<class T>
struct C3DDummyMiddleView
{
  static_assert(!T::IsVector, "This class decorates matrices.");
  using Index = typename T::Index;
  using Scalar = typename T::Scalar;

  T& data;

  C3DDummyMiddleView() = delete;

  C3DDummyMiddleView(T& ref) : data(ref) {}

  template<class U = T, su2enable_if<!std::is_const<U>::value> = 0>
  Scalar& operator() (Index i, Index, Index k) noexcept
  {
    return data(i,k);
  }

  const Scalar& operator() (Index i, Index, Index k) const noexcept
  {
    return data(i,k);
  }
};
