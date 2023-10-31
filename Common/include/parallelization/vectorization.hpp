/*!
 * \file vectorization.hpp
 * \brief Implementation of a portable SIMD type.
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

#include "../linear_algebra/vector_expressions.hpp"
#include "../parallelization/omp_structure.hpp"
#include <initializer_list>
#include <algorithm>
#include <cmath>
#ifdef __SSE2__
#include "x86intrin.h"
#endif

namespace simd {
/// \addtogroup SIMD
/// @{

using namespace VecExpr;

/*--- Detect preferred SIMD size (bytes). This only covers x86 architectures. ---*/
#if defined(__AVX512F__)
constexpr size_t PREFERRED_SIZE = 64;
#elif defined(__AVX__)
constexpr size_t PREFERRED_SIZE = 32;
#elif defined(__SSE2__)
constexpr size_t PREFERRED_SIZE = 16;
#else
constexpr size_t PREFERRED_SIZE = 8;
#endif

/*!
 * \brief Convert the SIMD size (bytes) to a lenght (num elems).
 */
template <class T>
constexpr size_t preferredLen() {
  return PREFERRED_SIZE / sizeof(T);
}

template <>
constexpr size_t preferredLen<su2double>() {
#ifdef CODI_REVERSE_TYPE
  /*--- Use a SIMD size of 1 for reverse AD, larger sizes increase
   * the pre-accumulation time with no performance benefit. ---*/
  return 1;
#else
  /*--- For forward AD there is a performance benefit. This covers
   * forward AD and primal mode (su2double == passivedouble). ---*/
  return PREFERRED_SIZE / sizeof(passivedouble);
#endif
}

/*!
 * \class Array
 * \brief A simple SIMD type relying on implicit vectorization, i.e. done by
 * the compiler, explicitly vectorized specializations are defined after.
 * \note This class gets its math operator overloads from CVecExpr, the
 * specializations do not use expression templates, IF YOU NEED A NEW FUNCTION,
 * define it both in vector_expressions.hpp and in special_vectorization.hpp.
 */
template <class Scalar_t, size_t N = preferredLen<Scalar_t>()>
class Array : public CVecExpr<Array<Scalar_t, N>, Scalar_t> {
#define FOREACH for (size_t k = 0; k < N; ++k)
  static_assert(N > 0, "Invalid SIMD size");

 public:
  using Scalar = Scalar_t;
  enum : size_t { Size = N };
  enum : size_t { Align = Size * sizeof(Scalar) };
  static constexpr bool StoreAsRef = true;

 private:
  alignas(Size * sizeof(Scalar)) Scalar x_[N];

 public:
#define ARRAY_BOILERPLATE                                                  \
  /*!--- Access elements ---*/                                             \
  FORCEINLINE Scalar& operator[](size_t k) { return x_[k]; }               \
  FORCEINLINE const Scalar& operator[](size_t k) const { return x_[k]; }   \
  /*!--- Constructors ---*/                                                \
  FORCEINLINE Array() = default;                                           \
  FORCEINLINE Array(Scalar x) { bcast(x); }                                \
  FORCEINLINE Array(std::initializer_list<Scalar> vals) {                  \
    auto it = vals.begin();                                                \
    FOREACH {                                                              \
      x_[k] = *it;                                                         \
      ++it;                                                                \
    }                                                                      \
  }                                                                        \
  FORCEINLINE Array(Scalar x0, Scalar dx) { FOREACH x_[k] = x0 + k * dx; } \
  FORCEINLINE Array(const Scalar* ptr) { load(ptr); }                      \
  template <class T>                                                       \
  FORCEINLINE Array(const Scalar* beg, const T& off) {                     \
    gather(beg, off);                                                      \
  }                                                                        \
  /*!--- Reduction operations ---*/                                        \
  FORCEINLINE Scalar sum() const {                                         \
    Scalar s(0);                                                           \
    FOREACH { s += x_[k]; }                                                \
    return s;                                                              \
  }                                                                        \
  FORCEINLINE Scalar dot(const Array& other) const {                       \
    Scalar s(0);                                                           \
    FOREACH { s += x_[k] * other[k]; }                                     \
    return s;                                                              \
  }

#if defined(CODI_REVERSE_TYPE) || defined(CODI_FORWARD_TYPE)
  /*--- These are not very nice but without them it would not be
   * possible to assign literals to Arrays of active types. ---*/
  template <class U = Scalar, su2enable_if<std::is_same<U, su2double>::value> = 0>
  FORCEINLINE Array(passivedouble x) {
    bcast(x);
  }
  template <class U = Scalar, su2enable_if<std::is_same<U, su2double>::value> = 0>
  FORCEINLINE Array& operator=(passivedouble x) {
    bcast(x);
    return *this;
  }
#endif

  ARRAY_BOILERPLATE

  /*! \brief Copy construct from expression. */
  template <class U>
  FORCEINLINE Array(const CVecExpr<U, Scalar>& expr) {
    FOREACH x_[k] = expr.derived()[k];
  }

  /*--- Implementation of the construction primitives. ---*/

  FORCEINLINE void bcast(Scalar x) { FOREACH x_[k] = x; }
  FORCEINLINE void load(const Scalar* ptr) { FOREACH x_[k] = ptr[k]; }
  FORCEINLINE void loada(const Scalar* ptr) { load(ptr); }
  FORCEINLINE void store(Scalar* ptr) const { FOREACH ptr[k] = x_[k]; }
  FORCEINLINE void storea(Scalar* ptr) const { store(ptr); }
  FORCEINLINE void stream(Scalar* ptr) const { store(ptr); }
  template <class T>
  FORCEINLINE void gather(const Scalar* begin, const T& offsets) {
    FOREACH x_[k] = begin[offsets[k]];
  }

  /*--- Compound assignment operators. ---*/

#define MAKE_COMPOUND(OP)                                           \
  FORCEINLINE Array& operator OP(Scalar x) {                        \
    FOREACH { x_[k] OP x; }                                         \
    return *this;                                                   \
  }                                                                 \
  template <class U>                                                \
  FORCEINLINE Array& operator OP(const CVecExpr<U, Scalar>& expr) { \
    FOREACH { x_[k] OP expr.derived()[k]; }                         \
    return *this;                                                   \
  }
  MAKE_COMPOUND(=)
  MAKE_COMPOUND(+=)
  MAKE_COMPOUND(-=)
  MAKE_COMPOUND(*=)
  MAKE_COMPOUND(/=)
#undef MAKE_COMPOUND

#undef FOREACH
};

/*--- Explicit vectorization specializations, see e.g.
 * https://software.intel.com/sites/landingpage/IntrinsicsGuide/
 * for documentation on the "_mm*" functions. ---*/

/*--- Size tags for overload resolution of some wrapper functions. ---*/
namespace SizeTag {
struct TWO {};
struct FOUR {};
struct EIGHT {};
struct SIXTEEN {};
}  // namespace SizeTag

/*--- Constants for bitwise implementations. ---*/
/*--- abs forces the sign bit to 0 ("x" & 0b0111...). ---*/
constexpr auto abs_mask_d = 0x7FFFFFFFFFFFFFFFL;
/*--- negation flips the sign bit ("x" ^ 0b1000...). ---*/
constexpr auto sign_mask_d = 0x8000000000000000L;

#ifdef __SSE2__
/*!
 * Create specialization for array of 2 doubles (this should be always available).
 */
#define ARRAY_T Array<double, 2>
#define SCALAR_T double
#define REGISTER_T __m128d
#define SIZE_TAG SizeTag::TWO()

static const __m128d abs_mask_2d = _mm_castsi128_pd(_mm_set1_epi64x(abs_mask_d));
static const __m128d sign_mask_2d = _mm_castsi128_pd(_mm_set1_epi64x(sign_mask_d));
static const __m128d ones_2d = _mm_set1_pd(1);

FORCEINLINE __m128d set1_p(SizeTag::TWO, double p) { return _mm_set1_pd(p); }
FORCEINLINE __m128d load_p(SizeTag::TWO, const double* p) { return _mm_load_pd(p); }
FORCEINLINE __m128d loadu_p(SizeTag::TWO, const double* p) { return _mm_loadu_pd(p); }
FORCEINLINE void store_p(double* p, __m128d x) { _mm_store_pd(p, x); }
FORCEINLINE void storeu_p(double* p, __m128d x) { _mm_storeu_pd(p, x); }
FORCEINLINE void stream_p(double* p, __m128d x) { _mm_stream_pd(p, x); }

FORCEINLINE __m128d add_p(__m128d a, __m128d b) { return _mm_add_pd(a, b); }
FORCEINLINE __m128d sub_p(__m128d a, __m128d b) { return _mm_sub_pd(a, b); }
FORCEINLINE __m128d mul_p(__m128d a, __m128d b) { return _mm_mul_pd(a, b); }
FORCEINLINE __m128d div_p(__m128d a, __m128d b) { return _mm_div_pd(a, b); }
FORCEINLINE __m128d max_p(__m128d a, __m128d b) { return _mm_max_pd(a, b); }
FORCEINLINE __m128d min_p(__m128d a, __m128d b) { return _mm_min_pd(a, b); }

FORCEINLINE __m128d eq_p(__m128d a, __m128d b) { return _mm_and_pd(ones_2d, _mm_cmpeq_pd(a, b)); }
FORCEINLINE __m128d lt_p(__m128d a, __m128d b) { return _mm_and_pd(ones_2d, _mm_cmplt_pd(a, b)); }
FORCEINLINE __m128d le_p(__m128d a, __m128d b) { return _mm_and_pd(ones_2d, _mm_cmple_pd(a, b)); }
FORCEINLINE __m128d ne_p(__m128d a, __m128d b) { return _mm_and_pd(ones_2d, _mm_cmpneq_pd(a, b)); }
FORCEINLINE __m128d ge_p(__m128d a, __m128d b) { return _mm_and_pd(ones_2d, _mm_cmpge_pd(a, b)); }
FORCEINLINE __m128d gt_p(__m128d a, __m128d b) { return _mm_and_pd(ones_2d, _mm_cmpgt_pd(a, b)); }

FORCEINLINE __m128d sqrt_p(__m128d x) { return _mm_sqrt_pd(x); }
FORCEINLINE __m128d abs_p(__m128d x) { return _mm_and_pd(x, abs_mask_2d); }
FORCEINLINE __m128d neg_p(__m128d x) { return _mm_xor_pd(x, sign_mask_2d); }
FORCEINLINE __m128d sign_p(__m128d x) { return _mm_or_pd(ones_2d, _mm_and_pd(x, sign_mask_2d)); }

/*--- Generate specialization based on the defines
 * and functions above by including the header. ---*/

#include "special_vectorization.hpp"

#endif  // __SSE2__

#ifdef __AVX__
/*!
 * Create specialization for array of 4 doubles.
 */
#define ARRAY_T Array<double, 4>
#define SCALAR_T double
#define REGISTER_T __m256d
#define SIZE_TAG SizeTag::FOUR()

static const __m256d abs_mask_4d = _mm256_castsi256_pd(_mm256_set1_epi64x(abs_mask_d));
static const __m256d sign_mask_4d = _mm256_castsi256_pd(_mm256_set1_epi64x(sign_mask_d));
static const __m256d ones_4d = _mm256_set1_pd(1);

FORCEINLINE __m256d set1_p(SizeTag::FOUR, double p) { return _mm256_set1_pd(p); }
FORCEINLINE __m256d load_p(SizeTag::FOUR, const double* p) { return _mm256_load_pd(p); }
FORCEINLINE __m256d loadu_p(SizeTag::FOUR, const double* p) { return _mm256_loadu_pd(p); }
FORCEINLINE void store_p(double* p, __m256d x) { _mm256_store_pd(p, x); }
FORCEINLINE void storeu_p(double* p, __m256d x) { _mm256_storeu_pd(p, x); }
FORCEINLINE void stream_p(double* p, __m256d x) { _mm256_stream_pd(p, x); }

FORCEINLINE __m256d add_p(__m256d a, __m256d b) { return _mm256_add_pd(a, b); }
FORCEINLINE __m256d sub_p(__m256d a, __m256d b) { return _mm256_sub_pd(a, b); }
FORCEINLINE __m256d mul_p(__m256d a, __m256d b) { return _mm256_mul_pd(a, b); }
FORCEINLINE __m256d div_p(__m256d a, __m256d b) { return _mm256_div_pd(a, b); }
FORCEINLINE __m256d max_p(__m256d a, __m256d b) { return _mm256_max_pd(a, b); }
FORCEINLINE __m256d min_p(__m256d a, __m256d b) { return _mm256_min_pd(a, b); }

FORCEINLINE __m256d eq_p(__m256d a, __m256d b) { return _mm256_and_pd(ones_4d, _mm256_cmp_pd(a, b, 0)); }
FORCEINLINE __m256d lt_p(__m256d a, __m256d b) { return _mm256_and_pd(ones_4d, _mm256_cmp_pd(a, b, 1)); }
FORCEINLINE __m256d le_p(__m256d a, __m256d b) { return _mm256_and_pd(ones_4d, _mm256_cmp_pd(a, b, 2)); }
FORCEINLINE __m256d ne_p(__m256d a, __m256d b) { return _mm256_and_pd(ones_4d, _mm256_cmp_pd(a, b, 4)); }
FORCEINLINE __m256d ge_p(__m256d a, __m256d b) { return _mm256_and_pd(ones_4d, _mm256_cmp_pd(a, b, 13)); }
FORCEINLINE __m256d gt_p(__m256d a, __m256d b) { return _mm256_and_pd(ones_4d, _mm256_cmp_pd(a, b, 14)); }

FORCEINLINE __m256d sqrt_p(__m256d x) { return _mm256_sqrt_pd(x); }
FORCEINLINE __m256d abs_p(__m256d x) { return _mm256_and_pd(x, abs_mask_4d); }
FORCEINLINE __m256d neg_p(__m256d x) { return _mm256_xor_pd(x, sign_mask_4d); }
FORCEINLINE __m256d sign_p(__m256d x) { return _mm256_or_pd(ones_4d, _mm256_and_pd(x, sign_mask_4d)); }

#include "special_vectorization.hpp"

#endif  // __AVX__

#ifdef __AVX512F__
/*!
 * Create specialization for array of 8 doubles.
 */
#define ARRAY_T Array<double, 8>
#define SCALAR_T double
#define REGISTER_T __m512d
#define SIZE_TAG SizeTag::EIGHT()

static const __m512d abs_mask_8d = _mm512_castsi512_pd(_mm512_set1_epi64(abs_mask_d));
static const __m512d sign_mask_8d = _mm512_castsi512_pd(_mm512_set1_epi64(sign_mask_d));
static const __m512d ones_8d = _mm512_set1_pd(1);

FORCEINLINE __m512d set1_p(SizeTag::EIGHT, double p) { return _mm512_set1_pd(p); }
FORCEINLINE __m512d load_p(SizeTag::EIGHT, const double* p) { return _mm512_load_pd(p); }
FORCEINLINE __m512d loadu_p(SizeTag::EIGHT, const double* p) { return _mm512_loadu_pd(p); }
FORCEINLINE void store_p(double* p, __m512d x) { _mm512_store_pd(p, x); }
FORCEINLINE void storeu_p(double* p, __m512d x) { _mm512_storeu_pd(p, x); }
FORCEINLINE void stream_p(double* p, __m512d x) { _mm512_stream_pd(p, x); }

FORCEINLINE __m512d add_p(__m512d a, __m512d b) { return _mm512_add_pd(a, b); }
FORCEINLINE __m512d sub_p(__m512d a, __m512d b) { return _mm512_sub_pd(a, b); }
FORCEINLINE __m512d mul_p(__m512d a, __m512d b) { return _mm512_mul_pd(a, b); }
FORCEINLINE __m512d div_p(__m512d a, __m512d b) { return _mm512_div_pd(a, b); }
FORCEINLINE __m512d max_p(__m512d a, __m512d b) { return _mm512_max_pd(a, b); }
FORCEINLINE __m512d min_p(__m512d a, __m512d b) { return _mm512_min_pd(a, b); }

template <int opCode>
FORCEINLINE __m512d cmp_p(__m512d a, __m512d b) {
  return _mm512_mask_blend_pd(_mm512_cmp_pd_mask(a, b, opCode), _mm512_setzero_pd(), ones_8d);
}
FORCEINLINE __m512d eq_p(__m512d a, __m512d b) { return cmp_p<0>(a, b); }
FORCEINLINE __m512d lt_p(__m512d a, __m512d b) { return cmp_p<1>(a, b); }
FORCEINLINE __m512d le_p(__m512d a, __m512d b) { return cmp_p<2>(a, b); }
FORCEINLINE __m512d ne_p(__m512d a, __m512d b) { return cmp_p<4>(a, b); }
FORCEINLINE __m512d ge_p(__m512d a, __m512d b) { return cmp_p<13>(a, b); }
FORCEINLINE __m512d gt_p(__m512d a, __m512d b) { return cmp_p<14>(a, b); }

FORCEINLINE __m512d sqrt_p(__m512d x) { return _mm512_sqrt_pd(x); }
FORCEINLINE __m512d abs_p(__m512d x) { return _mm512_and_pd(x, abs_mask_8d); }
FORCEINLINE __m512d neg_p(__m512d x) { return _mm512_xor_pd(x, sign_mask_8d); }
FORCEINLINE __m512d sign_p(__m512d x) { return _mm512_or_pd(ones_8d, _mm512_and_pd(x, sign_mask_8d)); }

#include "special_vectorization.hpp"

#endif  // __AVX512F__

#undef ARRAY_BOILERPLATE

/// @}
}  // namespace simd
