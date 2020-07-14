/*!
 * \file vectorization.hpp
 * \brief Implementation of a portable SIMD type.
 * \author P. Gomes
 * \version 7.0.6 "Blackbird"
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

#include "../linear_algebra/vector_expressions.hpp"
#include "../omp_structure.hpp"
#include <initializer_list>
#include <algorithm>
#include <cmath>

namespace simd {

using namespace VecExpr;

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
 * \brief A simple SIMD type relying on implicit vectorization, i.e. done by
 * the compiler, explicitly vectorized specializations are defined after.
 * \note This class gets its math operator overloads from CVecExpr, the
 * specializations do not use expression templates, IF YOU NEED A NEW FUNCTION,
 * define it both in vector_expressions.hpp and in special_vectorization.hpp.
 */
template<class Scalar_t, size_t N = simdLen<Scalar_t>()>
class Array : public CVecExpr<Array<Scalar_t,N>, Scalar_t> {
#define FOREACH SU2_OMP_SIMD for(size_t k=0; k<N; ++k)
  static_assert(N > 0, "Invalid SIMD size");
public:
  using Scalar = Scalar_t;
  enum : size_t {Size = N};
  enum : size_t {Align = SIMD_SIZE};

private:
  alignas(Align) Scalar x_[N];

public:
#define ARRAY_BOILERPLATE                                                     \
  /*!--- Access elements ---*/                                                \
  FORCEINLINE Scalar& operator[] (size_t k) { return x_[k]; }                 \
  FORCEINLINE const Scalar& operator[] (size_t k) const { return x_[k]; }     \
  /*!--- Constructors ---*/                                                   \
  FORCEINLINE Array() = default;                                              \
  FORCEINLINE Array(Scalar x) { bcast(x); }                                   \
  FORCEINLINE Array& operator= (Scalar x) { bcast(x); return *this; }         \
  FORCEINLINE Array(std::initializer_list<Scalar> vals) {                     \
    auto it = vals.begin(); FOREACH { x_[k] = *it; ++it; }                    \
  }                                                                           \
  FORCEINLINE Array(Scalar x0, Scalar dx) { FOREACH x_[k] = x0 + k*dx; }      \
  FORCEINLINE Array(const Scalar* ptr) { load(ptr); }                         \
  template<class T>                                                           \
  FORCEINLINE Array(const Scalar* beg, const T& off) { gather(beg,off); }     \
  /*!--- Reduction operations ---*/                                           \
  FORCEINLINE Scalar sum() const { Scalar s(0); FOREACH s+=x_[k]; return s; } \
  FORCEINLINE Scalar dot(const Array& other) const {                          \
    Scalar s(0); FOREACH s += x_[k] * other[k]; return s;                     \
  }

#if defined(CODI_REVERSE_TYPE) || defined(CODI_FORWARD_TYPE)
  template<class U = Scalar_t,
           typename std::enable_if<std::is_same<U,su2double>::value, bool>::type = 0>
  FORCEINLINE Array(passivedouble x) { bcast(x); }
  template<class U = Scalar_t,
           typename std::enable_if<std::is_same<U,su2double>::value, bool>::type = 0>
  FORCEINLINE Array& operator= (passivedouble x) { bcast(x); return *this; }
#endif

  ARRAY_BOILERPLATE

  /*! copy construct from expression. */
  template<class T, class U>
  FORCEINLINE Array(const CVecExpr<T,U>& expr) { FOREACH x_[k] = expr[k]; }

  /*! assign from expression. */
  template<class T, class U>
  FORCEINLINE Array& operator= (const CVecExpr<T,U>& expr) { FOREACH x_[k] = expr[k]; return *this; }

  /*--- Implementation of the construction primitives. ---*/

  FORCEINLINE void bcast(Scalar x) { FOREACH x_[k] = x; }
  FORCEINLINE void load(const Scalar* ptr) { FOREACH x_[k] = ptr[k]; }
  FORCEINLINE void loada(const Scalar* ptr) { load(ptr); }
  FORCEINLINE void store(Scalar* ptr) const { FOREACH ptr[k] = x_[k]; }
  FORCEINLINE void storea(Scalar* ptr) const { store(ptr); }
  FORCEINLINE void stream(Scalar* ptr) const { store(ptr); }
  template<class T>
  FORCEINLINE void gather(const Scalar* begin, const T& offsets) { FOREACH x_[k] = begin[offsets[k]]; }

  /*--- Compound math operators, "this" is not returned because it generates poor assembly. ---*/

#define MAKE_COMPOUND(OP)\
  FORCEINLINE void operator OP (Scalar x) { FOREACH x_[k] OP x; }\
  template<class T, class U>\
  FORCEINLINE void operator OP (const CVecExpr<T,U>& expr) { FOREACH x_[k] OP expr[k]; }
  MAKE_COMPOUND(+=)
  MAKE_COMPOUND(-=)
  MAKE_COMPOUND(*=)
  MAKE_COMPOUND(/=)
#undef MAKE_COMPOUND

#undef FOREACH
};

/*--- Explicit vectorization specializations. ---*/

#ifdef __AVX__
#include "x86intrin.h"

/*--- Size tags for overload resolution of some wrapper functions. ---*/
namespace SizeTag {
  struct FOUR {};
  struct EIGHT {};
  struct SIXTEEN {};
}

/*!
 * Create specialization for array of 4 doubles.
 */
#define ARRAY_T Array<double,4>
#define SCALAR_T double
#define REGISTER_T __m256d
#define SIZE_TAG SizeTag::FOUR()

/*--- abs forces the sign bit to 0 (& 0b0111...). ---*/
static const __m256d abs_mask_4d = _mm256_castsi256_pd(_mm256_set1_epi64x(0x7FFFFFFFFFFFFFFFL));
/*--- negation flips the sign bit (^ 0b1000...). ---*/
static const __m256d neg_mask_4d = _mm256_castsi256_pd(_mm256_set1_epi64x(0x8000000000000000L));

FORCEINLINE __m256d set1_p(SizeTag::FOUR, double p) { return _mm256_set1_pd(p); }
FORCEINLINE __m256d load_p(SizeTag::FOUR, const double* p) { return _mm256_load_pd(p); }
FORCEINLINE __m256d loadu_p(SizeTag::FOUR, const double* p) { return _mm256_loadu_pd(p); }
FORCEINLINE void store_p(double* p, __m256d x) { _mm256_store_pd(p,x); }
FORCEINLINE void storeu_p(double* p, __m256d x) { _mm256_storeu_pd(p,x); }
FORCEINLINE void stream_p(double* p, __m256d x) { _mm256_stream_pd(p,x); }

FORCEINLINE __m256d add_p(__m256d a, __m256d b) { return _mm256_add_pd(a,b); }
FORCEINLINE __m256d sub_p(__m256d a, __m256d b) { return _mm256_sub_pd(a,b); }
FORCEINLINE __m256d mul_p(__m256d a, __m256d b) { return _mm256_mul_pd(a,b); }
FORCEINLINE __m256d div_p(__m256d a, __m256d b) { return _mm256_div_pd(a,b); }
FORCEINLINE __m256d max_p(__m256d a, __m256d b) { return _mm256_max_pd(a,b); }
FORCEINLINE __m256d min_p(__m256d a, __m256d b) { return _mm256_min_pd(a,b); }

FORCEINLINE __m256d sqrt_p(__m256d x) { return _mm256_sqrt_pd(x); }
FORCEINLINE __m256d abs_p(__m256d x) { return _mm256_and_pd(x, abs_mask_4d); }
FORCEINLINE __m256d neg_p(__m256d x) { return _mm256_xor_pd(x, neg_mask_4d); }

#include "special_vectorization.hpp"

#endif // __AVX__

#ifdef __AVX512F__

/*!
 * Create specialization for array of 8 doubles.
 */
#define ARRAY_T Array<double,8>
#define SCALAR_T double
#define REGISTER_T __m512d
#define SIZE_TAG SizeTag::EIGHT()

static const __m512d abs_mask_8d = _mm512_castsi512_pd(_mm512_set1_epi64(0x7FFFFFFFFFFFFFFFL));
static const __m512d neg_mask_8d = _mm512_castsi512_pd(_mm512_set1_epi64(0x8000000000000000L));

FORCEINLINE __m512d set1_p(SizeTag::EIGHT, double p) { return _mm512_set1_pd(p); }
FORCEINLINE __m512d load_p(SizeTag::EIGHT, const double* p) { return _mm512_load_pd(p); }
FORCEINLINE __m512d loadu_p(SizeTag::EIGHT, const double* p) { return _mm512_loadu_pd(p); }
FORCEINLINE void store_p(double* p, __m512d x) { _mm512_store_pd(p,x); }
FORCEINLINE void storeu_p(double* p, __m512d x) { _mm512_storeu_pd(p,x); }
FORCEINLINE void stream_p(double* p, __m512d x) { _mm512_stream_pd(p,x); }

FORCEINLINE __m512d add_p(__m512d a, __m512d b) { return _mm512_add_pd(a,b); }
FORCEINLINE __m512d sub_p(__m512d a, __m512d b) { return _mm512_sub_pd(a,b); }
FORCEINLINE __m512d mul_p(__m512d a, __m512d b) { return _mm512_mul_pd(a,b); }
FORCEINLINE __m512d div_p(__m512d a, __m512d b) { return _mm512_div_pd(a,b); }
FORCEINLINE __m512d max_p(__m512d a, __m512d b) { return _mm512_max_pd(a,b); }
FORCEINLINE __m512d min_p(__m512d a, __m512d b) { return _mm512_min_pd(a,b); }

FORCEINLINE __m512d sqrt_p(__m512d x) { return _mm512_sqrt_pd(x); }
FORCEINLINE __m512d abs_p(__m512d x) { return _mm512_and_pd(x, abs_mask_8d); }
FORCEINLINE __m512d neg_p(__m512d x) { return _mm512_xor_pd(x, neg_mask_8d); }

#include "special_vectorization.hpp"

#endif // __AVX512F__

#undef ARRAY_BOILERPLATE

} // namespace
