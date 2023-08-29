/*!
 * \file special_vectorization.hpp
 * \brief Code generator header to create specializations of simd::Array.
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

// no #pragma once, header needs to be included once per specialization.

/// \addtogroup SIMD
/// @{

/*!
 * \brief Symbols that need to be defined before including this header:
 * \param[in] ARRAY_T - The desired specialization of simd::Array.
 * \param[in] SCALAR_T - Scalar type.
 * \param[in] REGISTER_T - Intrinsic type.
 * \param[in] SIZE_TAG - Dummy object associated with the simd size.
 * \param[in] ARRAY_BOILERPLATE - Generates the general ctors, access, etc.
 * \note On top of that, the intrinsic functions must be wrapped to
 * strip their type / size characteristics, e.g. _mm256_add_pd -> add_p,
 * overload resolution will do the rest. The first four symbols are
 * undefined once we are done using them.
 */
template <>
class ARRAY_T {
#define FOREACH SU2_OMP_SIMD for (size_t k = 0; k < Size; ++k)
  template <class F, class S>
  FORCEINLINE static S second(F, S s) {
    return s;
  }

 public:
  using Scalar = SCALAR_T;
  using Register = REGISTER_T;
  enum : size_t { Align = alignof(Register) };
  enum : size_t { Size = sizeof(Register) / sizeof(Scalar) };

  /*--- The infamous union "hack", sue me. ---*/
  union {
    Register reg;
    Scalar x_[Size];
  };

  /*--- Same basic construction operations as the general case. ---*/

  ARRAY_BOILERPLATE

  /*--- Special construction using the "register type" directly. ---*/

  FORCEINLINE Array(Register y) { reg = y; }
  FORCEINLINE Array(const Array& other) noexcept { reg = other.reg; }

  /*--- Specialized construction primitives. ---*/

  FORCEINLINE void bcast(Scalar x) { reg = set1_p(SIZE_TAG, x); }
  FORCEINLINE void load(const Scalar* ptr) { reg = loadu_p(SIZE_TAG, ptr); }
  FORCEINLINE void loada(const Scalar* ptr) { reg = load_p(SIZE_TAG, ptr); }
  FORCEINLINE void store(Scalar* ptr) const { storeu_p(ptr, reg); }
  FORCEINLINE void storea(Scalar* ptr) const { store_p(ptr, reg); }
  FORCEINLINE void stream(Scalar* ptr) const { stream_p(ptr, reg); }
  template <class T>
  FORCEINLINE void gather(const Scalar* begin, const T& offsets) {
    FOREACH x_[k] = begin[offsets[k]];
  }

  /*--- Compound assignement operators. ---*/

#define MAKE_COMPOUND(OP, IMPL)                                 \
  FORCEINLINE Array& operator OP(Scalar x) {                    \
    reg = IMPL(reg, set1_p(SIZE_TAG, x));                       \
    return *this;                                               \
  }                                                             \
  FORCEINLINE Array& operator OP(const Array& other) noexcept { \
    reg = IMPL(reg, other.reg);                                 \
    return *this;                                               \
  }
  MAKE_COMPOUND(=, second)
  MAKE_COMPOUND(+=, add_p)
  MAKE_COMPOUND(-=, sub_p)
  MAKE_COMPOUND(*=, mul_p)
  MAKE_COMPOUND(/=, div_p)
#undef MAKE_COMPOUND

#undef FOREACH
};

/*!
 * SIMD overloads, NAME is the operator or function,
 * IMPL the intrinsic function that implements it.
 */
#define MAKE_UNARY_FUN(NAME, IMPL) \
  FORCEINLINE ARRAY_T NAME(const ARRAY_T& x) { return IMPL(x.reg); }

MAKE_UNARY_FUN(operator-, neg_p)
MAKE_UNARY_FUN(sqrt, sqrt_p)
MAKE_UNARY_FUN(abs, abs_p)
MAKE_UNARY_FUN(sign, sign_p)

#undef MAKE_UNARY_FUN

#define MAKE_BINARY_FUN(NAME, IMPL)                                                                   \
  FORCEINLINE ARRAY_T NAME(const ARRAY_T& a, const ARRAY_T& b) { return IMPL(a.reg, b.reg); }         \
  FORCEINLINE ARRAY_T NAME(const ARRAY_T& a, SCALAR_T b) { return IMPL(a.reg, set1_p(SIZE_TAG, b)); } \
  FORCEINLINE ARRAY_T NAME(SCALAR_T b, const ARRAY_T& a) { return IMPL(set1_p(SIZE_TAG, b), a.reg); }

MAKE_BINARY_FUN(operator+, add_p)
MAKE_BINARY_FUN(operator-, sub_p)
MAKE_BINARY_FUN(operator*, mul_p)
MAKE_BINARY_FUN(operator/, div_p)
MAKE_BINARY_FUN(operator<, lt_p)
MAKE_BINARY_FUN(operator>, gt_p)
MAKE_BINARY_FUN(operator==, eq_p)
MAKE_BINARY_FUN(operator!=, ne_p)
MAKE_BINARY_FUN(operator<=, le_p)
MAKE_BINARY_FUN(operator>=, ge_p)
MAKE_BINARY_FUN(fmax, max_p)
MAKE_BINARY_FUN(fmin, min_p)

#undef MAKE_BINARY_FUN

/*!
 * Compatibility mode overloads, element-wise implementation.
 */
#define FOREACH SU2_OMP_SIMD for (size_t k = 0; k < ARRAY_T::Size; ++k)

/*--- Functions of one (array) argument. ---*/

#define MAKE_UNARY_FUN(NAME, IMPL)             \
  FORCEINLINE ARRAY_T NAME(const ARRAY_T& x) { \
    ARRAY_T res;                               \
    FOREACH { res[k] = IMPL(x[k]); }           \
    return res;                                \
  }

#undef MAKE_UNARY_FUN

/*--- Functions of two arguments, with arrays and scalars. ---*/

#define MAKE_BINARY_FUN(NAME, IMPL)                              \
  FORCEINLINE ARRAY_T NAME(const ARRAY_T& a, const ARRAY_T& b) { \
    ARRAY_T res;                                                 \
    FOREACH { res[k] = IMPL(a[k], b[k]); }                       \
    return res;                                                  \
  }                                                              \
  FORCEINLINE ARRAY_T NAME(const ARRAY_T& a, SCALAR_T b) {       \
    ARRAY_T res;                                                 \
    FOREACH { res[k] = IMPL(a[k], b); }                          \
    return res;                                                  \
  }                                                              \
  FORCEINLINE ARRAY_T NAME(SCALAR_T b, const ARRAY_T& a) {       \
    ARRAY_T res;                                                 \
    FOREACH { res[k] = IMPL(b, a[k]); }                          \
    return res;                                                  \
  }

MAKE_BINARY_FUN(pow, ::pow)

#undef MAKE_BINARY_FUN

#undef FOREACH

#undef ARRAY_T
#undef SCALAR_T
#undef REGISTER_T
#undef SIZE_TAG

/// @}
