/*!
 * \file code_config.hpp
 * \brief Header file for collecting common macros, definitions and type configurations.
 * \author T. Albring, P. Gomes, J. Blühdorn
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

#include <cstdint>
#include <type_traits>
#include <cmath>

#if defined(_MSC_VER)
#define PRAGMIZE(X) __pragma(X)
#else
#define PRAGMIZE(X) _Pragma(#X)
#endif

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

/*--- Marks a function callable from both host and device code, a no-op outside nvcc. ---*/
#ifdef __CUDACC__
#define SU2_CUDA_HOST_DEVICE __host__ __device__
#else
#define SU2_CUDA_HOST_DEVICE
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

/*--- Detect whether the CUDA kernels are part of this build. The .cu translation units
 * cannot be compiled with the CoDiPack defines (nvcc's device pass cannot parse the tape
 * machinery), and an object compiled with a different definition of su2double must not be
 * linked into an AD library. They are therefore only built into the primal libraries, and
 * all device dispatch has to be compiled out of the AD builds, which HAVE_CUDA alone does
 * not do because su2mixedfloat is a passive type there as well. ---*/
#if defined(HAVE_CUDA) && !defined(CODI_REVERSE_TYPE) && !defined(CODI_FORWARD_TYPE)
#define SU2_ENABLE_CUDA_KERNELS
#endif

/*--- No full single precision for AD builds. ---*/
#if (defined(CODI_REVERSE_TYPE) || defined(CODI_FORWARD_TYPE)) && defined(USE_SINGLE_PRECISION)
#undef USE_SINGLE_PRECISION
#endif

/*--- Default integer types. Currently used for rank-local sparse patterns. ---*/
using su2uint = uint32_t;
using su2int = int32_t;

/*--- This type can be used for (rare) compatibility cases or for
 * computations that are intended to be (always) passive. ---*/
#ifdef USE_SINGLE_PRECISION
using passivedouble = float;
#else
using passivedouble = double;
#endif

/*--- std::min/max do not compile if the arguments have inconsistent types, which
 * happens in single precision due to floating point literals (double by default).
 * These overloads delegate to fmin/fmax which do not have that problem. ---*/
#ifdef USE_SINGLE_PRECISION
namespace std {
FORCEINLINE float min(const float& a, const double& b) { return fmin(a, static_cast<float>(b)); }
FORCEINLINE float min(const double& b, const float& a) { return fmin(a, static_cast<float>(b)); }
FORCEINLINE float max(const float& a, const double& b) { return fmax(a, static_cast<float>(b)); }
FORCEINLINE float max(const double& b, const float& a) { return fmax(a, static_cast<float>(b)); }
}  // namespace std
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
#elif defined(CODI_TAG_TAPE)
using su2double = codi::RealReverseTag;
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
using su2double = passivedouble;
#endif

/*--- Define a type for potentially lower precision operations. ---*/
#ifndef CODI_FORWARD_TYPE
#ifdef USE_MIXED_PRECISION
using su2mixedfloat = float;
#else
using su2mixedfloat = passivedouble;
#endif
#else
/*--- There is no lower precision for forward AD so undefine the macro to simplify
 * the logic needed to deal with the multiple type configurations. ---*/
#ifdef USE_MIXED_PRECISION
#undef USE_MIXED_PRECISION
#endif
using su2mixedfloat = su2double;
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

#ifdef __GNUC__
#define SU2_IGNORE_WARNING(WARNING) \
  PRAGMIZE(GCC diagnostic push)     \
  PRAGMIZE(GCC diagnostic ignored WARNING)
#define SU2_RESTORE_WARNING PRAGMIZE(GCC diagnostic pop)
#else
#define SU2_IGNORE_WARNING(WARNING)
#define SU2_RESTORE_WARNING
#endif
