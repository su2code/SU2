/*!
 * \file vectorization.hpp
 * \brief Implementation of a portable SIMD type.
 * \author P. Gomes
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

#include "../linear_algebra/vector_expressions.hpp"
#include "../parallelization/omp_structure.hpp"
#include <initializer_list>
#include <algorithm>
#include <cmath>
#ifdef __SSE2__
#include "x86intrin.h"
#endif
#if defined(__ARM_NEON) || defined(__ARM_NEON__)
#include <arm_neon.h>
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
#elif defined(__SSE2__) || defined(__ARM_NEON) || defined(__ARM_NEON__)
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
  enum : size_t { Align = Size * alignof(Scalar) };
  static constexpr bool StoreAsRef = true;

 private:
  alignas(Size * alignof(Scalar)) Scalar x_[N];

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
constexpr uint64_t abs_mask_d = 0x7FFFFFFFFFFFFFFFL;
constexpr uint32_t abs_mask_s = 0x7FFFFFFFU;
/*--- negation flips the sign bit ("x" ^ 0b1000...). ---*/
constexpr uint64_t sign_mask_d = 0x8000000000000000L;
constexpr uint32_t sign_mask_s = 0x80000000U;

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

/*!
 * Create specialization for array of 4 floats (this should be always available).
 */
#define ARRAY_T Array<float, 4>
#define SCALAR_T float
#define REGISTER_T __m128
#define SIZE_TAG SizeTag::FOUR()

static const __m128 abs_mask_4s = _mm_castsi128_ps(_mm_set1_epi32(abs_mask_s));
static const __m128 sign_mask_4s = _mm_castsi128_ps(_mm_set1_epi32(sign_mask_s));
static const __m128 ones_4s = _mm_set1_ps(1);

FORCEINLINE __m128 set1_p(SizeTag::FOUR, float p) { return _mm_set1_ps(p); }
FORCEINLINE __m128 load_p(SizeTag::FOUR, const float* p) { return _mm_load_ps(p); }
FORCEINLINE __m128 loadu_p(SizeTag::FOUR, const float* p) { return _mm_loadu_ps(p); }
FORCEINLINE void store_p(float* p, __m128 x) { _mm_store_ps(p, x); }
FORCEINLINE void storeu_p(float* p, __m128 x) { _mm_storeu_ps(p, x); }
FORCEINLINE void stream_p(float* p, __m128 x) { _mm_stream_ps(p, x); }

FORCEINLINE __m128 add_p(__m128 a, __m128 b) { return _mm_add_ps(a, b); }
FORCEINLINE __m128 sub_p(__m128 a, __m128 b) { return _mm_sub_ps(a, b); }
FORCEINLINE __m128 mul_p(__m128 a, __m128 b) { return _mm_mul_ps(a, b); }
FORCEINLINE __m128 div_p(__m128 a, __m128 b) { return _mm_div_ps(a, b); }
FORCEINLINE __m128 max_p(__m128 a, __m128 b) { return _mm_max_ps(a, b); }
FORCEINLINE __m128 min_p(__m128 a, __m128 b) { return _mm_min_ps(a, b); }

FORCEINLINE __m128 eq_p(__m128 a, __m128 b) { return _mm_and_ps(ones_4s, _mm_cmpeq_ps(a, b)); }
FORCEINLINE __m128 lt_p(__m128 a, __m128 b) { return _mm_and_ps(ones_4s, _mm_cmplt_ps(a, b)); }
FORCEINLINE __m128 le_p(__m128 a, __m128 b) { return _mm_and_ps(ones_4s, _mm_cmple_ps(a, b)); }
FORCEINLINE __m128 ne_p(__m128 a, __m128 b) { return _mm_and_ps(ones_4s, _mm_cmpneq_ps(a, b)); }
FORCEINLINE __m128 ge_p(__m128 a, __m128 b) { return _mm_and_ps(ones_4s, _mm_cmpge_ps(a, b)); }
FORCEINLINE __m128 gt_p(__m128 a, __m128 b) { return _mm_and_ps(ones_4s, _mm_cmpgt_ps(a, b)); }

FORCEINLINE __m128 sqrt_p(__m128 x) { return _mm_sqrt_ps(x); }
FORCEINLINE __m128 abs_p(__m128 x) { return _mm_and_ps(x, abs_mask_4s); }
FORCEINLINE __m128 neg_p(__m128 x) { return _mm_xor_ps(x, sign_mask_4s); }
FORCEINLINE __m128 sign_p(__m128 x) { return _mm_or_ps(ones_4s, _mm_and_ps(x, sign_mask_4s)); }

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

/*!
 * Create specialization for array of 8 floats.
 */
#define ARRAY_T Array<float, 8>
#define SCALAR_T float
#define REGISTER_T __m256
#define SIZE_TAG SizeTag::EIGHT()

static const __m256 abs_mask_8s = _mm256_castsi256_ps(_mm256_set1_epi32(abs_mask_s));
static const __m256 sign_mask_8s = _mm256_castsi256_ps(_mm256_set1_epi32(sign_mask_s));
static const __m256 ones_8s = _mm256_set1_ps(1);

FORCEINLINE __m256 set1_p(SizeTag::EIGHT, float p) { return _mm256_set1_ps(p); }
FORCEINLINE __m256 load_p(SizeTag::EIGHT, const float* p) { return _mm256_load_ps(p); }
FORCEINLINE __m256 loadu_p(SizeTag::EIGHT, const float* p) { return _mm256_loadu_ps(p); }
FORCEINLINE void store_p(float* p, __m256 x) { _mm256_store_ps(p, x); }
FORCEINLINE void storeu_p(float* p, __m256 x) { _mm256_storeu_ps(p, x); }
FORCEINLINE void stream_p(float* p, __m256 x) { _mm256_stream_ps(p, x); }

FORCEINLINE __m256 add_p(__m256 a, __m256 b) { return _mm256_add_ps(a, b); }
FORCEINLINE __m256 sub_p(__m256 a, __m256 b) { return _mm256_sub_ps(a, b); }
FORCEINLINE __m256 mul_p(__m256 a, __m256 b) { return _mm256_mul_ps(a, b); }
FORCEINLINE __m256 div_p(__m256 a, __m256 b) { return _mm256_div_ps(a, b); }
FORCEINLINE __m256 max_p(__m256 a, __m256 b) { return _mm256_max_ps(a, b); }
FORCEINLINE __m256 min_p(__m256 a, __m256 b) { return _mm256_min_ps(a, b); }

FORCEINLINE __m256 eq_p(__m256 a, __m256 b) { return _mm256_and_ps(ones_8s, _mm256_cmp_ps(a, b, 0)); }
FORCEINLINE __m256 lt_p(__m256 a, __m256 b) { return _mm256_and_ps(ones_8s, _mm256_cmp_ps(a, b, 1)); }
FORCEINLINE __m256 le_p(__m256 a, __m256 b) { return _mm256_and_ps(ones_8s, _mm256_cmp_ps(a, b, 2)); }
FORCEINLINE __m256 ne_p(__m256 a, __m256 b) { return _mm256_and_ps(ones_8s, _mm256_cmp_ps(a, b, 4)); }
FORCEINLINE __m256 ge_p(__m256 a, __m256 b) { return _mm256_and_ps(ones_8s, _mm256_cmp_ps(a, b, 13)); }
FORCEINLINE __m256 gt_p(__m256 a, __m256 b) { return _mm256_and_ps(ones_8s, _mm256_cmp_ps(a, b, 14)); }

FORCEINLINE __m256 sqrt_p(__m256 x) { return _mm256_sqrt_ps(x); }
FORCEINLINE __m256 abs_p(__m256 x) { return _mm256_and_ps(x, abs_mask_8s); }
FORCEINLINE __m256 neg_p(__m256 x) { return _mm256_xor_ps(x, sign_mask_8s); }
FORCEINLINE __m256 sign_p(__m256 x) { return _mm256_or_ps(ones_8s, _mm256_and_ps(x, sign_mask_8s)); }

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

/*!
 * Create specialization for array of 16 floats.
 */
#define ARRAY_T Array<float, 16>
#define SCALAR_T float
#define REGISTER_T __m512
#define SIZE_TAG SizeTag::SIXTEEN()

static const __m512 abs_mask_16s = _mm512_castsi512_ps(_mm512_set1_epi32(abs_mask_s));
static const __m512 sign_mask_16s = _mm512_castsi512_ps(_mm512_set1_epi32(sign_mask_s));
static const __m512 ones_16s = _mm512_set1_ps(1);

FORCEINLINE __m512 set1_p(SizeTag::SIXTEEN, float p) { return _mm512_set1_ps(p); }
FORCEINLINE __m512 load_p(SizeTag::SIXTEEN, const float* p) { return _mm512_load_ps(p); }
FORCEINLINE __m512 loadu_p(SizeTag::SIXTEEN, const float* p) { return _mm512_loadu_ps(p); }
FORCEINLINE void store_p(float* p, __m512 x) { _mm512_store_ps(p, x); }
FORCEINLINE void storeu_p(float* p, __m512 x) { _mm512_storeu_ps(p, x); }
FORCEINLINE void stream_p(float* p, __m512 x) { _mm512_stream_ps(p, x); }

FORCEINLINE __m512 add_p(__m512 a, __m512 b) { return _mm512_add_ps(a, b); }
FORCEINLINE __m512 sub_p(__m512 a, __m512 b) { return _mm512_sub_ps(a, b); }
FORCEINLINE __m512 mul_p(__m512 a, __m512 b) { return _mm512_mul_ps(a, b); }
FORCEINLINE __m512 div_p(__m512 a, __m512 b) { return _mm512_div_ps(a, b); }
FORCEINLINE __m512 max_p(__m512 a, __m512 b) { return _mm512_max_ps(a, b); }
FORCEINLINE __m512 min_p(__m512 a, __m512 b) { return _mm512_min_ps(a, b); }

template <int opCode>
FORCEINLINE __m512 cmp_p(__m512 a, __m512 b) {
  return _mm512_mask_blend_ps(_mm512_cmp_ps_mask(a, b, opCode), _mm512_setzero_ps(), ones_16s);
}
FORCEINLINE __m512 eq_p(__m512 a, __m512 b) { return cmp_p<0>(a, b); }
FORCEINLINE __m512 lt_p(__m512 a, __m512 b) { return cmp_p<1>(a, b); }
FORCEINLINE __m512 le_p(__m512 a, __m512 b) { return cmp_p<2>(a, b); }
FORCEINLINE __m512 ne_p(__m512 a, __m512 b) { return cmp_p<4>(a, b); }
FORCEINLINE __m512 ge_p(__m512 a, __m512 b) { return cmp_p<13>(a, b); }
FORCEINLINE __m512 gt_p(__m512 a, __m512 b) { return cmp_p<14>(a, b); }

FORCEINLINE __m512 sqrt_p(__m512 x) { return _mm512_sqrt_ps(x); }
FORCEINLINE __m512 abs_p(__m512 x) { return _mm512_and_ps(x, abs_mask_16s); }
FORCEINLINE __m512 neg_p(__m512 x) { return _mm512_xor_ps(x, sign_mask_16s); }
FORCEINLINE __m512 sign_p(__m512 x) { return _mm512_or_ps(ones_16s, _mm512_and_ps(x, sign_mask_16s)); }

#include "special_vectorization.hpp"

#endif  // __AVX512F__

#if defined(__ARM_NEON) || defined(__ARM_NEON__)
/*!
 * Create specialization for array of 2 doubles.
 */
#define ARRAY_T Array<double, 2>
#define SCALAR_T double
#define REGISTER_T float64x2_t
#define SIZE_TAG SizeTag::TWO()

static const uint64x2_t abs_mask_2d_u = vdupq_n_u64(abs_mask_d);
static const uint64x2_t sign_mask_2d_u = vdupq_n_u64(sign_mask_d);
static const uint64x2_t ones_2d_u = vreinterpretq_u64_f64(vdupq_n_f64(1.0));
static const uint64x2_t ones_2u = vdupq_n_u64(~0ULL);

FORCEINLINE float64x2_t set1_p(SizeTag::TWO, double p) { return vdupq_n_f64(p); }
FORCEINLINE float64x2_t load_p(SizeTag::TWO, const double* p) { return vld1q_f64(p); }
FORCEINLINE float64x2_t loadu_p(SizeTag::TWO, const double* p) { return vld1q_f64(p); }
FORCEINLINE void store_p(double* p, float64x2_t x) { vst1q_f64(p, x); }
FORCEINLINE void storeu_p(double* p, float64x2_t x) { vst1q_f64(p, x); }
/*--- No direct NEON equivalent to streaming stores. ---*/
FORCEINLINE void stream_p(double* p, float64x2_t x) { vst1q_f64(p, x); }

FORCEINLINE float64x2_t add_p(float64x2_t a, float64x2_t b) { return vaddq_f64(a, b); }
FORCEINLINE float64x2_t sub_p(float64x2_t a, float64x2_t b) { return vsubq_f64(a, b); }
FORCEINLINE float64x2_t mul_p(float64x2_t a, float64x2_t b) { return vmulq_f64(a, b); }
FORCEINLINE float64x2_t div_p(float64x2_t a, float64x2_t b) { return vdivq_f64(a, b); }
FORCEINLINE float64x2_t max_p(float64x2_t a, float64x2_t b) { return vmaxq_f64(a, b); }
FORCEINLINE float64x2_t min_p(float64x2_t a, float64x2_t b) { return vminq_f64(a, b); }

/*--- Comparisons return uint64x2_t masks. Convert to 0.0 / 1.0. ---*/
FORCEINLINE float64x2_t int2float(uint64x2_t a) { return vreinterpretq_f64_u64(a); }
FORCEINLINE float64x2_t cmp2float(uint64x2_t cmp) { return int2float(vandq_u64(ones_2d_u, cmp)); }
FORCEINLINE float64x2_t eq_p(float64x2_t a, float64x2_t b) { return cmp2float(vceqq_f64(a, b)); }
FORCEINLINE float64x2_t lt_p(float64x2_t a, float64x2_t b) { return cmp2float(vcltq_f64(a, b)); }
FORCEINLINE float64x2_t le_p(float64x2_t a, float64x2_t b) { return cmp2float(vcleq_f64(a, b)); }
FORCEINLINE float64x2_t ne_p(float64x2_t a, float64x2_t b) { return cmp2float(veorq_u64(ones_2u, vceqq_f64(a, b))); }
FORCEINLINE float64x2_t ge_p(float64x2_t a, float64x2_t b) { return cmp2float(vcgeq_f64(a, b)); }
FORCEINLINE float64x2_t gt_p(float64x2_t a, float64x2_t b) { return cmp2float(vcgtq_f64(a, b)); }

FORCEINLINE float64x2_t sqrt_p(float64x2_t x) { return vsqrtq_f64(x); }
FORCEINLINE float64x2_t abs_p(float64x2_t x) { return int2float(vandq_u64(vreinterpretq_u64_f64(x), abs_mask_2d_u)); }
FORCEINLINE float64x2_t neg_p(float64x2_t x) { return int2float(veorq_u64(vreinterpretq_u64_f64(x), sign_mask_2d_u)); }
FORCEINLINE float64x2_t sign_p(float64x2_t x) {
  return int2float(vorrq_u64(ones_2d_u, vandq_u64(vreinterpretq_u64_f64(x), sign_mask_2d_u)));
}

#include "special_vectorization.hpp"

/*!
 * Create specialization for array of 4 floats.
 */
#define ARRAY_T Array<float, 4>
#define SCALAR_T float
#define REGISTER_T float32x4_t
#define SIZE_TAG SizeTag::FOUR()

static const uint32x4_t abs_mask_4s_u = vdupq_n_u32(abs_mask_s);
static const uint32x4_t sign_mask_4s_u = vdupq_n_u32(sign_mask_s);
static const uint32x4_t ones_4s_u = vreinterpretq_u32_f32(vdupq_n_f32(1.0f));

FORCEINLINE float32x4_t set1_p(SizeTag::FOUR, float p) { return vdupq_n_f32(p); }
FORCEINLINE float32x4_t load_p(SizeTag::FOUR, const float* p) { return vld1q_f32(p); }
FORCEINLINE float32x4_t loadu_p(SizeTag::FOUR, const float* p) { return vld1q_f32(p); }
FORCEINLINE void store_p(float* p, float32x4_t x) { vst1q_f32(p, x); }
FORCEINLINE void storeu_p(float* p, float32x4_t x) { vst1q_f32(p, x); }
FORCEINLINE void stream_p(float* p, float32x4_t x) { vst1q_f32(p, x); }
FORCEINLINE float32x4_t add_p(float32x4_t a, float32x4_t b) { return vaddq_f32(a, b); }
FORCEINLINE float32x4_t sub_p(float32x4_t a, float32x4_t b) { return vsubq_f32(a, b); }
FORCEINLINE float32x4_t mul_p(float32x4_t a, float32x4_t b) { return vmulq_f32(a, b); }
FORCEINLINE float32x4_t div_p(float32x4_t a, float32x4_t b) { return vdivq_f32(a, b); }
FORCEINLINE float32x4_t max_p(float32x4_t a, float32x4_t b) { return vmaxq_f32(a, b); }
FORCEINLINE float32x4_t min_p(float32x4_t a, float32x4_t b) { return vminq_f32(a, b); }

FORCEINLINE float32x4_t int2float(uint32x4_t a) { return vreinterpretq_f32_u32(a); }
FORCEINLINE float32x4_t cmp2float(uint32x4_t cmp) { return int2float(vandq_u32(ones_4s_u, cmp)); }
FORCEINLINE float32x4_t eq_p(float32x4_t a, float32x4_t b) { return cmp2float(vceqq_f32(a, b)); }
FORCEINLINE float32x4_t lt_p(float32x4_t a, float32x4_t b) { return cmp2float(vcltq_f32(a, b)); }
FORCEINLINE float32x4_t le_p(float32x4_t a, float32x4_t b) { return cmp2float(vcleq_f32(a, b)); }
FORCEINLINE float32x4_t ne_p(float32x4_t a, float32x4_t b) { return cmp2float(vmvnq_u32(vceqq_f32(a, b))); }
FORCEINLINE float32x4_t ge_p(float32x4_t a, float32x4_t b) { return cmp2float(vcgeq_f32(a, b)); }
FORCEINLINE float32x4_t gt_p(float32x4_t a, float32x4_t b) { return cmp2float(vcgtq_f32(a, b)); }

FORCEINLINE float32x4_t sqrt_p(float32x4_t x) { return vsqrtq_f32(x); }
FORCEINLINE float32x4_t abs_p(float32x4_t x) { return int2float(vandq_u32(vreinterpretq_u32_f32(x), abs_mask_4s_u)); }
FORCEINLINE float32x4_t neg_p(float32x4_t x) { return int2float(veorq_u32(vreinterpretq_u32_f32(x), sign_mask_4s_u)); }
FORCEINLINE float32x4_t sign_p(float32x4_t x) {
  return int2float(vorrq_u32(ones_4s_u, vandq_u32(vreinterpretq_u32_f32(x), sign_mask_4s_u)));
}

#include "special_vectorization.hpp"

#endif  // __ARM_NEON__

#undef ARRAY_BOILERPLATE

/// @}
}  // namespace simd
