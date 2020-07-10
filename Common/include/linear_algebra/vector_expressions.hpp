/*!
 * \file vector_expressions.hpp
 * \brief Basic expression templates for vector types.
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

#include "../basic_types/datatype_structure.hpp"
#include <type_traits>
#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <cmath>

namespace VecExpr {

/*!
 * \brief Base vector expression class.
 * \param[in] Derived - The class that inherits (derives) from this one.
 * \param[in] Scalar_t - Return type of Derived::operator[] (which must be implemented).
 * \note Derived classes must be constructible from expressions (that is when the
 * expression tree is evaluated).
 */
template<class Derived, class Scalar_t>
class CVecExpr {
public:
  using Scalar = Scalar_t;
  FORCEINLINE Scalar operator[] (size_t i) const { return static_cast<const Derived&>(*this)[i]; }
};

template<class T>
using decay_t = typename std::decay<T>::type;

/*--- Macro to create expressions and overloads for unary functions. ---*/

#define MAKE_UNARY_FUN(FUN, NAME, IMPL)                                   \
template<class U>                                                         \
class NAME : public CVecExpr<NAME<U>, typename U::Scalar> {               \
  const U& u;                                                             \
public:                                                                   \
  using Scalar = typename U::Scalar;                                      \
  FORCEINLINE NAME(const U& u_) : u(u_) {}                                \
  FORCEINLINE Scalar operator[] (size_t i) const { return IMPL(u[i]); }   \
};                                                                        \
template<class U, class T>                                                \
FORCEINLINE NAME<U> FUN(const CVecExpr<U,T>& u) {                         \
  return NAME<U>(static_cast<const U&>(u));                               \
}

MAKE_UNARY_FUN(operator-, negate_v, -)
MAKE_UNARY_FUN(abs, abs_v, std::abs)
MAKE_UNARY_FUN(sqrt, sqrt_v, std::sqrt)

#undef MAKE_UNARY_FUN

/*--- Macro to create expressions and overloads for binary functions. ---*/

#define MAKE_BINARY_FUN(FUN, NAME1, NAME2, NAME3, IMPL)                       \
/*--- Vector with vector, expression (NAME1) and function overload. ---*/     \
template<class U, class V>                                                    \
class NAME1 : public CVecExpr<NAME1<U,V>, typename U::Scalar> {               \
  const U& u;                                                                 \
  const V& v;                                                                 \
public:                                                                       \
  using Scalar = typename U::Scalar;                                          \
  FORCEINLINE NAME1(const U& u_, const V& v_) : u(u_), v(v_) {}               \
  FORCEINLINE Scalar operator[] (size_t i) const { return IMPL(u[i], v[i]); } \
};                                                                            \
template<class U, class V, class T>                                           \
FORCEINLINE NAME1<U,V> FUN(const CVecExpr<U,T>& u, const CVecExpr<V,T>& v) {  \
  return NAME1<U,V>(static_cast<const U&>(u), static_cast<const V&>(v));      \
}                                                                             \
/*--- Vector with scalar, expression (NAME2) and function overload. ---*/     \
template<class U, class V>                                                    \
class NAME2 : public CVecExpr<NAME2<U,V>, typename U::Scalar> {               \
  const U& u;                                                                 \
  V v;                                                                        \
public:                                                                       \
  using Scalar = typename U::Scalar;                                          \
  FORCEINLINE NAME2(const U& u_, const V& v_) : u(u_), v(v_) {}               \
  FORCEINLINE Scalar operator[] (size_t i) const { return IMPL(u[i], v); }    \
};                                                                            \
template<class U, class T>                                                    \
FORCEINLINE NAME2<U,T> FUN(const CVecExpr<U,T>& u, decay_t<T> v) {            \
  return NAME2<U,T>(static_cast<const U&>(u), v);                             \
}                                                                             \
/*--- Scalar with vector, expression (NAME3) and function overload. ---*/     \
template<class U, class V>                                                    \
class NAME3 : public CVecExpr<NAME3<U,V>, typename U::Scalar> {               \
  const U& u;                                                                 \
  V v;                                                                        \
public:                                                                       \
  using Scalar = typename U::Scalar;                                          \
  FORCEINLINE NAME3(const U& u_, const V& v_) : u(u_), v(v_) {}               \
  FORCEINLINE Scalar operator[] (size_t i) const { return IMPL(v, u[i]); }    \
};                                                                            \
template<class U, class T>                                                    \
FORCEINLINE NAME3<U,T> FUN(decay_t<T> v, const CVecExpr<U,T>& u) {            \
  return NAME3<U,T>(static_cast<const U&>(u), v);                             \
}

/*--- Naming convention: v - vector (expression); s - scalar. ---*/

/*--- sts::plus and co. were tried, the code was horrendous
 * (due to the forced conversion between different types). ---*/

#define add_ss(a,b) a+b
#define sub_ss(a,b) a-b
#define mul_ss(a,b) a*b
#define div_ss(a,b) a/b
#define leq_ss(a,b) a<=b
#define geq_ss(a,b) a>=b
#define eq_ss(a,b) a==b
#define ne_ss(a,b) a!=b
#define lt_ss(a,b) a<b
#define gt_ss(a,b) a>b

MAKE_BINARY_FUN(operator+, add_vv, add_vs, add_sv, add_ss)
MAKE_BINARY_FUN(operator-, sub_vv, sub_vs, sub_sv, sub_ss)
MAKE_BINARY_FUN(operator*, mul_vv, mul_vs, mul_sv, mul_ss)
MAKE_BINARY_FUN(operator/, div_vv, div_vs, div_sv, div_ss)
MAKE_BINARY_FUN(operator<=, leq_vv, leq_vs, leq_sv, leq_ss)
MAKE_BINARY_FUN(operator>=, geq_vv, geq_vs, geq_sv, geq_ss)
MAKE_BINARY_FUN(operator==, eq_vv, eq_vs, eq_sv, eq_ss)
MAKE_BINARY_FUN(operator!=, ne_vv, ne_vs, ne_sv, ne_ss)
MAKE_BINARY_FUN(operator<, lt_vv, lt_vs, lt_sv, lt_ss)
MAKE_BINARY_FUN(operator>, gt_vv, gt_vs, gt_sv, gt_ss)
MAKE_BINARY_FUN(max, max_vv, max_vs, max_sv, std::max)
MAKE_BINARY_FUN(min, min_vv, min_vs, min_sv, std::min)
MAKE_BINARY_FUN(pow, pow_vv, pow_vs, pow_sv, std::pow)

#undef add_ss
#undef sub_ss
#undef mul_ss
#undef div_ss
#undef leq_ss
#undef geq_ss
#undef eq_ss
#undef ne_ss
#undef lt_ss
#undef gt_ss

#undef MAKE_BINARY_FUN

} // end namespace
