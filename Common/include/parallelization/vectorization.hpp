/*!
 * \file vectorization.hpp
 * \brief Implementation of a portable SIMD type.
 * \author P. Gomes
 * \version 7.0.5 "Blackbird"
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

#include "../basic_types/datatype_structure.hpp"
#include "../omp_structure.hpp"
#include <initializer_list>
#include <cmath>
#include <algorithm>

namespace simd {

/*--- Detect preferred SIMD size (bytes). ---*/
#if defined(__AVX512F__)
constexpr size_t SIMD_SIZE = 64;
#elif defined(__AVX__)
constexpr size_t SIMD_SIZE = 32;
#elif defined(__SSE__)
constexpr size_t SIMD_SIZE = 16;
#else
constexpr size_t SIMD_SIZE = 8;
#endif

/*!
 * \brief Convert the SIMD size (bytes) to a lenght (num elems).
 */
template<class T>
constexpr size_t simdLen() { return SIMD_SIZE / sizeof(T); }
template<>
constexpr size_t simdLen<su2double>() { return SIMD_SIZE / sizeof(passivedouble); }

/*!
 * \class Array
 * \brief A simple SIMD type relying mostly on implicit vectorization, i.e. done by
 * the compiler, as opposed to explicit (done via intrinsics or inline assembly).
 */
template<class Scalar_t, size_t N = simdLen<Scalar_t>()>
class Array {
#define FOREACH SU2_OMP_SIMD for(size_t k=0; k<N; ++k)
  static_assert(N > 0, "Invalid SIMD size");
public:
  using Scalar = Scalar_t;
  enum : size_t {Size = N};
  enum : size_t {Align = SIMD_SIZE};

private:
  alignas(Align) Scalar x_[N];

public:
  /*--- Constructors ---*/

  FORCEINLINE Array() = default;

  // broadcast
  FORCEINLINE Array(Scalar x) { bcast(x); }
  FORCEINLINE Array& operator= (Scalar x) { bcast(x); return *this; }

  // initialize with given values
  FORCEINLINE Array(std::initializer_list<Scalar> vals) {
    auto it = vals.begin(); FOREACH { x_[k] = *it; ++it; }
  }

  // linearly spaced
  FORCEINLINE Array(Scalar x0, Scalar dx) { FOREACH x_[k] = x0 + k*dx; }

  // load
  FORCEINLINE Array(const Scalar* ptr) { load(ptr); }

  // gather
  template<class T>
  FORCEINLINE Array(const Scalar* base_ptr, const T& offsets) { gather(base_ptr,offsets); }

  // copy / assign
  FORCEINLINE Array(const Array& other) { FOREACH x_[k] = other[k]; }

  FORCEINLINE Array& operator= (const Array& other) { FOREACH x_[k] = other[k]; return *this; }

  /*--- Accessors. ---*/

  FORCEINLINE Scalar& operator[] (size_t k) { return x_[k]; }

  FORCEINLINE const Scalar& operator[] (size_t k) const { return x_[k]; }

  FORCEINLINE void bcast(Scalar x) { FOREACH x_[k] = x; }

  FORCEINLINE void load(const Scalar* ptr) { FOREACH x_[k] = ptr[k]; }

  FORCEINLINE void store(Scalar* ptr) const { FOREACH ptr[k] = x_[k]; }

  template<class T>
  FORCEINLINE void gather(const Scalar* base_ptr, const T& offsets) { FOREACH x_[k] = base_ptr[offsets[k]]; }

  /*--- Compound math operators, "this" is not returned because it generates poor assembly. ---*/

#define MAKE_COMPOUND(OP)\
  FORCEINLINE void operator OP (Scalar x) { FOREACH x_[k] OP x; }\
  FORCEINLINE void operator OP (const Array& other) { FOREACH x_[k] OP other.x_[k]; }
  MAKE_COMPOUND(+=)
  MAKE_COMPOUND(-=)
  MAKE_COMPOUND(*=)
  MAKE_COMPOUND(/=)
#undef MAKE_COMPOUND

  /*--- Reductions. ---*/

  FORCEINLINE Scalar sum() const { Scalar s(0); FOREACH s += x_[k]; return s; }

  FORCEINLINE Scalar dot(const Array& other) const {
    Scalar s(0); FOREACH s += x_[k] * other[k]; return s;
  }
};

/*--- Math, logical, and relational operators, with arrays and scalars. ---*/

#define MAKE_OPERATOR(OP)\
template<class T, size_t N>\
FORCEINLINE Array<T,N> operator OP (const Array<T,N>& a, const Array<T,N>& b) {\
  Array<T,N> res; FOREACH res[k] = a[k] OP b[k]; return res;\
}\
template<class T, size_t N, class U>\
FORCEINLINE Array<T,N> operator OP (const Array<T,N>& a, U b) {\
  Array<T,N> res; FOREACH res[k] = a[k] OP b; return res;\
}\
template<class T, size_t N, class U>\
FORCEINLINE Array<T,N> operator OP (U b, const Array<T,N>& a) {\
  Array<T,N> res; FOREACH res[k] = b OP a[k]; return res;\
}

MAKE_OPERATOR(+)
MAKE_OPERATOR(-)
MAKE_OPERATOR(*)
MAKE_OPERATOR(/)
MAKE_OPERATOR(==)
MAKE_OPERATOR(!=)
MAKE_OPERATOR(>)
MAKE_OPERATOR(<=)
MAKE_OPERATOR(<)
MAKE_OPERATOR(>=)
MAKE_OPERATOR(&)
MAKE_OPERATOR(|)

#undef MAKE_OPERATOR

/*--- Functions of one argument, first macro param is the name of
 * the created function, IMPL is the scalar implementation. ---*/

#define MAKE_UNARY_FUN(NAME,IMPL)\
template<class T, size_t N>\
FORCEINLINE Array<T,N> NAME(const Array<T,N>& x) {\
  Array<T,N> res; FOREACH res[k] = IMPL(x[k]); return res;\
}

MAKE_UNARY_FUN(sqrt,::sqrt)
MAKE_UNARY_FUN(abs,std::abs)

#undef MAKE_UNARY_FUN

/*--- Functions of two arguments, with arrays and scalars. ---*/

#define MAKE_BINARY_FUN(NAME,IMPL)\
template<class T, size_t N>\
FORCEINLINE Array<T,N> NAME(const Array<T,N>& a, const Array<T,N>& b) {\
  Array<T,N> res; FOREACH res[k] = IMPL(a[k], b[k]); return res;\
}\
template<class T, size_t N, class U>\
FORCEINLINE Array<T,N> NAME(const Array<T,N>& a, U b) {\
  Array<T,N> res; FOREACH res[k] = IMPL(a[k], b); return res;\
}\
template<class T, size_t N, class U>\
FORCEINLINE Array<T,N> NAME(U b, const Array<T,N>& a) {\
  Array<T,N> res; FOREACH res[k] = IMPL(b, a[k]); return res;\
}

MAKE_BINARY_FUN(max,std::max)
MAKE_BINARY_FUN(min,std::min)
MAKE_BINARY_FUN(pow,::pow)

#undef MAKE_BINARY_FUN

/*--- AVX specializations. ---*/
#ifdef __AVX__
#include "x86intrin.h"

/*!
 * \brief Specialization for Array of 4 doubles, helps reducing register spills.
 */
template<>
class Array<double,4> {
public:
  using Scalar = double;
  enum : size_t {Size = 4};
  enum : size_t {N = Size};
  enum : size_t {Align = 32};

  union {
    __m256d ymm;
    Scalar x_[N];
  };

  /*--- Constructors ---*/

  FORCEINLINE Array() = default;

  // broadcast
  FORCEINLINE Array(Scalar x) { bcast(x); }
  FORCEINLINE Array& operator= (Scalar x) { bcast(x); return *this; }

  // initialize with given values
  FORCEINLINE Array(std::initializer_list<Scalar> vals) {
    auto it = vals.begin(); FOREACH { x_[k] = *it; ++it; }
  }

  // linearly spaced
  FORCEINLINE Array(Scalar x0, Scalar dx) { FOREACH x_[k] = x0 + k*dx; }

  // load
  FORCEINLINE Array(const Scalar* ptr) { load(ptr); }

  // gather
  template<class T>
  FORCEINLINE Array(const Scalar* base_ptr, const T& offsets) { gather(base_ptr,offsets); }

  // copy / assign
  FORCEINLINE Array(__m256d y) { ymm = y; }

  FORCEINLINE Array(const Array& other) { ymm = other.ymm; }

  FORCEINLINE Array& operator= (const Array& other) { ymm = other.ymm; return *this; }

  /*--- Accessors. ---*/

  FORCEINLINE Scalar& operator[] (size_t k) { return x_[k]; }

  FORCEINLINE const Scalar& operator[] (size_t k) const { return x_[k]; }

  FORCEINLINE void bcast(Scalar x) { ymm = _mm256_broadcast_sd(&x); }

  FORCEINLINE void load(const Scalar* ptr) { ymm = _mm256_loadu_pd(ptr); }

  FORCEINLINE void store(Scalar* ptr) const { _mm256_storeu_pd(ptr, ymm); }

  template<class T>
  FORCEINLINE void gather(const Scalar* base_ptr, const T& offsets) { FOREACH x_[k] = base_ptr[offsets[k]]; }

  /*--- Compound math operators, "this" is not returned because it generates poor assembly. ---*/

#define MAKE_COMPOUND(OP,IMPL)\
  FORCEINLINE void operator OP (Scalar x) { ymm = IMPL(ymm, _mm256_broadcast_sd(&x)); }\
  FORCEINLINE void operator OP (const Array& other) { ymm = IMPL(ymm, other.ymm); }
  MAKE_COMPOUND(+=, _mm256_add_pd)
  MAKE_COMPOUND(-=, _mm256_sub_pd)
  MAKE_COMPOUND(*=, _mm256_mul_pd)
  MAKE_COMPOUND(/=, _mm256_div_pd)
#undef MAKE_COMPOUND

  /*--- Reductions. ---*/

  FORCEINLINE Scalar sum() const { Scalar s(0); FOREACH s += x_[k]; return s; }

  FORCEINLINE Scalar dot(const Array& other) const {
    Scalar s(0); FOREACH s += x_[k] * other[k]; return s;
  }
};

#define MAKE_UNARY_FUN(NAME,IMPL)\
template<>\
FORCEINLINE Array<double,4> NAME(const Array<double,4>& x) {return IMPL(x.ymm);}

MAKE_UNARY_FUN(sqrt, _mm256_sqrt_pd)

#undef MAKE_UNARY_FUN

#define MAKE_BINARY_FUN(NAME,IMPL)\
template<>\
FORCEINLINE Array<double,4> NAME (const Array<double,4>& a, const Array<double,4>& b) {\
  return IMPL(a.ymm, b.ymm);\
}\
template<>\
FORCEINLINE Array<double,4> NAME (const Array<double,4>& a, double b) {\
  return IMPL(a.ymm, _mm256_broadcast_sd(&b));\
}\
template<>\
FORCEINLINE Array<double,4> NAME (double b, const Array<double,4>& a) {\
  return IMPL(_mm256_broadcast_sd(&b), a.ymm);\
}

MAKE_BINARY_FUN(operator+, _mm256_add_pd)
MAKE_BINARY_FUN(operator-, _mm256_sub_pd)
MAKE_BINARY_FUN(operator*, _mm256_mul_pd)
MAKE_BINARY_FUN(operator/, _mm256_div_pd)
MAKE_BINARY_FUN(max, _mm256_max_pd)
MAKE_BINARY_FUN(min, _mm256_min_pd)

#undef MAKE_BINARY_FUN

#endif // __AVX__

#undef FOREACH
} // namespace
