/*!
 * \file code_config.hpp
 * \brief Header file for collecting common macros, definitions and type configurations.
 * \author T. Albring, P. Gomes, J. Blühdorn
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

#include <type_traits>

#if defined(_MSC_VER)
#define FORCEINLINE __forceinline
#elif defined(__GNUC__) || defined(__clang__) || defined(__INTEL_COMPILER)
#define FORCEINLINE inline __attribute__((always_inline))
#else
#define FORCEINLINE inline
#endif

#if defined(__GNUC__) || defined(__clang__) || defined(__INTEL_COMPILER)
#define NEVERINLINE inline __attribute__((noinline))
#else
#define NEVERINLINE inline
#endif

#if defined(__INTEL_COMPILER)
/*--- Disable warnings related to inline attributes. ---*/
#pragma warning disable 2196
#pragma warning disable 3415
/*--- Disable warnings related to overloaded virtual. ---*/
#pragma warning disable 654
#pragma warning disable 1125
#if defined(CODI_FORWARD_TYPE) || defined(CODI_REVERSE_TYPE)
#pragma warning disable 1875
#endif
#endif

/*--- Convenience SFINAE typedef to conditionally
 * enable/disable function template overloads. ---*/
template <bool condition>
using su2enable_if = typename std::enable_if<condition, bool>::type;

/*--- Compile-time type selection. ---*/
template <bool B, class T, class F>
struct su2conditional {
  using type = T;
};
template <class T, class F>
struct su2conditional<false, T, F> {
  using type = F;
};

template <bool B, class T, class F>
using su2conditional_t = typename su2conditional<B, T, F>::type;

/*! \brief Static cast "In" to "Out", in debug builds a dynamic cast is used. */
template <class Out, class In>
FORCEINLINE Out su2staticcast_p(In ptr) {
  static_assert(std::is_pointer<In>::value, "This expects a pointer");
#ifndef NDEBUG
  return static_cast<Out>(ptr);
#else
  return dynamic_cast<Out>(ptr);
#endif
}

/*--- Detect compilation with OpenMP. ---*/
#if defined(_OPENMP)
#define HAVE_OMP
#endif

/*--- Depending on the datatype defined during the configuration,
 * include the correct definition, and create the main typedef. ---*/

#if defined(CODI_REVERSE_TYPE)  // reverse mode AD
#include "codi.hpp"
#include "codi/tools/data/externalFunctionUserData.hpp"

#if defined(HAVE_OMP)
using su2double = codi::RealReverseIndexOpenMPGen<double, double>;
#else
#if defined(CODI_JACOBIAN_LINEAR_TAPE)
using su2double = codi::RealReverse;
#elif defined(CODI_JACOBIAN_REUSE_TAPE)
using su2double = codi::RealReverseIndexGen<double, double, codi::ReuseIndexManager<int> >;
#elif defined(CODI_JACOBIAN_MULTIUSE_TAPE)
using su2double = codi::RealReverseIndex;
#elif defined(CODI_PRIMAL_LINEAR_TAPE)
using su2double = codi::RealReversePrimal;
#elif defined(CODI_PRIMAL_REUSE_TAPE)
using su2double = codi::RealReversePrimalIndexGen<double, double, codi::ReuseIndexManager<int> >;
#elif defined(CODI_PRIMAL_MULTIUSE_TAPE)
using su2double = codi::RealReversePrimalIndex;
#else
#error "Please define a CoDiPack tape."
#endif
#endif

#if defined(HAVE_OMP) || defined(CODI_JACOBIAN_REUSE_TAPE) || defined(CODI_JACOBIAN_MULTIUSE_TAPE) || \
    defined(CODI_PRIMAL_REUSE_TAPE) || defined(CODI_PRIMAL_MULTIUSE_TAPE)
#define CODI_INDEX_REUSE
#endif
#elif defined(CODI_FORWARD_TYPE)  // forward mode AD
#include "codi.hpp"
using su2double = codi::RealForward;
#else  // primal / direct / no AD
using su2double = double;
#endif

/*--- This type can be used for (rare) compatibility cases or for
 * computations that are intended to be (always) passive. ---*/
using passivedouble = double;

/*--- Define a type for potentially lower precision operations. ---*/
#ifdef USE_MIXED_PRECISION
using su2mixedfloat = float;
#else
using su2mixedfloat = passivedouble;
#endif

/*--- Detect if OpDiLib has to be used. ---*/
#if defined(HAVE_OMP) && defined(CODI_REVERSE_TYPE)
#ifndef __INTEL_COMPILER
#define HAVE_OPDI
#else
#warning Hybrid parallel reverse mode AD cannot be used with Intel compilers.
#endif

#if (_OPENMP >= 201811 && !defined(FORCE_OPDI_MACRO_BACKEND)) || defined(FORCE_OPDI_OMPT_BACKEND)
#define HAVE_OMPT
#endif
#endif
