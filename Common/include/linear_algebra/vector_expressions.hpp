/*!
 * \file vector_expressions.hpp
 * \brief Expression templates for vector types with coefficient-wise operations.
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

#include "../basic_types/datatype_structure.hpp"
#include <type_traits>
#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <cmath>
#include <cstdint>

namespace VecExpr {
/// \addtogroup VecExpr
/// @{

/*!
 * \brief Base vector expression class.
 * \ingroup BLAS
 * \param[in] Derived - The class that inherits from this one to use the expressions.
 * \param[in] Scalar - Associated scalar type, prevents implicit conversions between exprs.
 * \note Derived classes must implement operator[], and at least operator= with
 * expressions (that is when they are evaluated). They must also contain a constexpr
 * boolean "StoreAsRef", indicating whether to store them by value (false) or by
 * reference (true), when composing expressions.
 * Expression classes need to be stored by value to allow nested expressions to
 * propagate correctly (i.e. if "Scalar" has its own expression templates).
 * Vector classes should be stored by reference to avoid copies, especially if they
 * allocate memory dynamically.
 */
template <class Derived, class Scalar>
class CVecExpr {
 public:
  /*!
   * \brief Cast the expression to Derived, usually to allow evaluation via operator[].
   */
  FORCEINLINE const Derived& derived() const { return static_cast<const Derived&>(*this); }

  // Allowed from C++14, allows nested expression propagation without
  // manually calling derived() on the expression being evaluated.
  // FORCEINLINE auto operator[] (size_t i) const { return derived()[i]; }
};

/*!
 * \brief Expression class to broadcast a scalar value. Allows implementing
 * "vector-scalar" operations re-using "vector-vector" expressions.
 */
template <class Scalar>
class Bcast : public CVecExpr<Bcast<Scalar>, Scalar> {
  Scalar x;

 public:
  static constexpr bool StoreAsRef = false;
  FORCEINLINE Bcast(const Scalar& x_) : x(x_) {}
  FORCEINLINE const Scalar& operator[](size_t) const { return x; }
};

/*!
 * \brief std::decay_t from C++14, used to allow implicit conversions
 * between scalar types, e.g. "CVecExpr<U,double>" + "int/double/etc.".
 */
template <class T>
using decay_t = typename std::decay<T>::type;

/*! \brief std::remove_reference_t from C++14, removes references from some type. */
template <class T>
using remove_reference_t = typename std::remove_reference<T>::type;

/*! \brief Mechanism to conditionally (based on "StoreAsRef") add lvalue reference to a type. */
template <class T, bool>
struct add_lref_if {
  using type = remove_reference_t<T>;
};
template <class T>
struct add_lref_if<T, true> {
  using type = remove_reference_t<T>&;
};
template <class T>
using store_t = typename add_lref_if<T, T::StoreAsRef>::type;

/*--- Namespace from which the math function implementations come. ---*/

#if defined(CODI_REVERSE_TYPE) || defined(CODI_FORWARD_TYPE)
namespace math = ::codi;
#else
namespace math = ::std;
#endif

/*--- Macro to simplify auto return type deduction in C++11, operator[] needs
 * it to allow inner expressions to propagate as the outer is evaluated.  ---*/

#define RETURNS(...) \
  ->decltype(__VA_ARGS__) { return __VA_ARGS__; }

/*--- Macro to create expression classes (EXPR) and overloads (FUN) for unary
 * functions, based on their coefficient-wise implementation (IMPL). ---*/

#define MAKE_UNARY_FUN(FUN, EXPR, IMPL)                             \
  /*!--- Expression class. ---*/                                    \
  template <class U, class Scalar>                                  \
  class EXPR : public CVecExpr<EXPR<U, Scalar>, Scalar> {           \
    store_t<const U> u;                                             \
                                                                    \
   public:                                                          \
    static constexpr bool StoreAsRef = false;                       \
    FORCEINLINE EXPR(const U& u_) : u(u_) {}                        \
    FORCEINLINE auto operator[](size_t i) const RETURNS(IMPL(u[i])) \
  };                                                                \
  /*!--- Function overload, returns an expression object. ---*/     \
  template <class U, class S>                                       \
  FORCEINLINE auto FUN(const CVecExpr<U, S>& u) RETURNS(EXPR<U, S>(u.derived()))

#define sign_impl(x) Scalar(1 - 2 * (x < 0))
MAKE_UNARY_FUN(operator-, minus_, -)
MAKE_UNARY_FUN(abs, abs_, math::abs)
MAKE_UNARY_FUN(sqrt, sqrt_, math::sqrt)
MAKE_UNARY_FUN(sign, sign_, sign_impl)
#undef sign_impl

#undef MAKE_UNARY_FUN

/*--- Macro to create expressions and overloads for binary functions. ---*/

// clang-format off
#define MAKE_BINARY_FUN(FUN, EXPR, IMPL)                                                                              \
  /*!--- Expression class. ---*/                                                                                      \
  template <class U, class V, class Scalar>                                                                           \
  class EXPR : public CVecExpr<EXPR<U, V, Scalar>, Scalar> {                                                          \
    store_t<const U> u;                                                                                               \
    store_t<const V> v;                                                                                               \
                                                                                                                      \
   public:                                                                                                            \
    static constexpr bool StoreAsRef = false;                                                                         \
    FORCEINLINE EXPR(const U& u_, const V& v_) : u(u_), v(v_) {}                                                      \
    FORCEINLINE auto operator[](size_t i) const RETURNS(IMPL(u[i], v[i]))                                             \
  };                                                                                                                  \
  /*!--- Vector with vector function overload. ---*/                                                                  \
  template <class U, class V, class S>                                                                                \
  FORCEINLINE auto FUN(const CVecExpr<U, S>& u, const CVecExpr<V, S>& v)                                              \
      RETURNS(EXPR<U, V, S>(u.derived(), v.derived()))                                                                \
  /*!--- Vector with scalar function overload. ---*/                                                                  \
  template <class U, class S>                                                                                         \
  FORCEINLINE auto FUN(const CVecExpr<U, S>& u, decay_t<S> v) RETURNS(EXPR<U, Bcast<S>, S>(u.derived(), Bcast<S>(v))) \
  /*!--- Scalar with vector function overload. ---*/                                                                  \
  template <class S, class V>                                                                                         \
  FORCEINLINE auto FUN(decay_t<S> u, const CVecExpr<V, S>& v) RETURNS(EXPR<Bcast<S>, V, S>(Bcast<S>(u), v.derived()))
// clang-format on

/*--- std::max/min have issues (because they return by reference).
 * fmin and fmax return by value and thus are fine, but they would force
 * conversions to double, to avoid that we provide integer overloads.
 * We use int32/64 instead of int/long to avoid issues with Windows,
 * where long is 32 bits (instead of 64 bits). ---*/

#define MAKE_FMINMAX_OVERLOADS(TYPE)                              \
  FORCEINLINE TYPE fmax(TYPE a, TYPE b) { return a < b ? b : a; } \
  FORCEINLINE TYPE fmin(TYPE a, TYPE b) { return a < b ? a : b; }
MAKE_FMINMAX_OVERLOADS(int32_t)
MAKE_FMINMAX_OVERLOADS(int64_t)
MAKE_FMINMAX_OVERLOADS(uint32_t)
MAKE_FMINMAX_OVERLOADS(uint64_t)
/*--- Make the float and double versions of fmin/max available in this
 * namespace to avoid ambiguous overloads. ---*/
using std::fmax;
using std::fmin;
#undef MAKE_FMINMAX_OVERLOADS

MAKE_BINARY_FUN(fmax, max_, fmax)
MAKE_BINARY_FUN(fmin, min_, fmin)
MAKE_BINARY_FUN(pow, pow_, math::pow)

/*--- sts::plus and co. were tried, the code was horrendous (due to the forced
 * conversion between different types) and creating functions for these ops
 * requires a lot of boilerplate (template args, auto return, etc.). ---*/

#define add_impl(a, b) a + b
#define sub_impl(a, b) a - b
#define mul_impl(a, b) a* b
#define div_impl(a, b) a / b
MAKE_BINARY_FUN(operator+, add_, add_impl)
MAKE_BINARY_FUN(operator-, sub_, sub_impl)
MAKE_BINARY_FUN(operator*, mul_, mul_impl)
MAKE_BINARY_FUN(operator/, div_, div_impl)
#undef add_impl
#undef sub_impl
#undef mul_impl
#undef div_impl

/*--- Relational operators need to be cast to the scalar type to allow vectorization.
 * TO_PASSIVE is used to convert active scalars to passive, which CoDi will then capture
 * by value in its expressions, and thus dangling references are avoided. No AD info
 * is lost since these operators are non-differentiable. ---*/

#define TO_PASSIVE(IMPL) SU2_TYPE::Passive<Scalar>::Value(IMPL)
#define le_impl(a, b) TO_PASSIVE(a <= b)
#define ge_impl(a, b) TO_PASSIVE(a >= b)
#define eq_impl(a, b) TO_PASSIVE(a == b)
#define ne_impl(a, b) TO_PASSIVE(a != b)
#define lt_impl(a, b) TO_PASSIVE(a < b)
#define gt_impl(a, b) TO_PASSIVE(a > b)
MAKE_BINARY_FUN(operator<=, le_, le_impl)
MAKE_BINARY_FUN(operator>=, ge_, ge_impl)
MAKE_BINARY_FUN(operator==, eq_, eq_impl)
MAKE_BINARY_FUN(operator!=, ne_, ne_impl)
MAKE_BINARY_FUN(operator<, lt_, lt_impl)
MAKE_BINARY_FUN(operator>, gt_, gt_impl)
#undef TO_PASSIVE
#undef le_impl
#undef ge_impl
#undef eq_impl
#undef ne_impl
#undef lt_impl
#undef gt_impl

#undef MAKE_BINARY_FUN

/// @}
}  // namespace VecExpr
