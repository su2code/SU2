/*!
 * \file omp_structure.hpp
 * \brief OpenMP interface header, provides compatibility functions
 *        if the code is built without OpenMP support.
 *        Parallel pragmas are defined here so that they can be
 *        completely "disabled" when compiling without OpenMP.
 * \note Do not include omp.h explicitly anywhere, use this header instead.
 * \note If you use an omp_*** function define a compatibility version here,
 *       if that is not practical use define "HAVE_OMP" to guard that function.
 * \note Always use the macro "SU2_OMP" to create OpenMP constructs, this is so
 *       we can disable pragmas. Other convenient pragmas are also defined here
 *       e.g. SU2_OMP_PARALLEL. Exotic pragmas of limited portability should be
 *       defined here with suitable fallback versions to limit the spread of
 *       compiler tricks in other areas of the code.
 * \author P. Gomes
 * \version 7.0.4 "Blackbird"
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

#include <type_traits>

#if defined(_MSC_VER)
#define PRAGMIZE(X) __pragma(X)
#else
#define PRAGMIZE(X) _Pragma(#X)
#endif

/*--- Detect compilation with OpenMP support, protect agaisnt
 *    using OpenMP with Reverse AD (not supported yet). ---*/
#if defined(_OPENMP) && !defined(CODI_REVERSE_TYPE)
#define HAVE_OMP
#include <omp.h>

/*--- The generic start of OpenMP constructs. ---*/
#define SU2_OMP(ARGS) PRAGMIZE(omp ARGS)

#else // Compile without OpenMP

/*--- Disable pragmas to quiet compilation warnings. ---*/
#define SU2_OMP(ARGS)

/*!
 * \brief Maximum number of threads available.
 */
inline constexpr int omp_get_max_threads(void) {return 1;}

/*!
 * \brief Number of threads in current team.
 */
inline constexpr int omp_get_num_threads(void) {return 1;}

/*!
 * \brief Set the maximum number of threads.
 */
inline void omp_set_num_threads(int) { }

/*!
 * \brief Index of current thread, akin to MPI rank.
 */
inline constexpr int omp_get_thread_num(void) {return 0;}

/*!
 * \brief Dummy lock type and associated functions.
 */
struct omp_lock_t {};
struct DummyVectorOfLocks {
  omp_lock_t l;
  inline omp_lock_t& operator[](int) {return l;}
};
inline void omp_init_lock(omp_lock_t*){}
inline void omp_set_lock(omp_lock_t*){}
inline void omp_unset_lock(omp_lock_t*){}
inline void omp_destroy_lock(omp_lock_t*){}

#endif

/*--- Convenience macros (do not use excessive nesting of macros). ---*/
#define SU2_OMP_SIMD SU2_OMP(simd)

#define SU2_OMP_MASTER SU2_OMP(master)
#define SU2_OMP_ATOMIC SU2_OMP(atomic)
#define SU2_OMP_BARRIER SU2_OMP(barrier)
#define SU2_OMP_CRITICAL SU2_OMP(critical)

#define SU2_OMP_PARALLEL SU2_OMP(parallel)
#define SU2_OMP_PARALLEL_(ARGS) SU2_OMP(parallel ARGS)
#define SU2_OMP_PARALLEL_ON(NTHREADS) SU2_OMP(parallel num_threads(NTHREADS))

#define SU2_OMP_FOR_DYN(CHUNK) SU2_OMP(for schedule(dynamic,CHUNK))
#define SU2_OMP_FOR_STAT(CHUNK) SU2_OMP(for schedule(static,CHUNK))


/*--- Convenience functions (e.g. to compute chunk sizes). ---*/

/*!
 * \brief Integer division rounding up.
 */
inline constexpr size_t roundUpDiv(size_t numerator, size_t denominator)
{
  return (numerator+denominator-1)/denominator;
}

/*!
 * \brief Round up to next multiple.
 */
inline constexpr size_t nextMultiple(size_t argument, size_t multiple)
{
  return roundUpDiv(argument, multiple) * multiple;
}

/*!
 * \brief Compute a chunk size based on totalWork and number of threads such that
 *        all threads get the same number of chunks (with limited size).
 * \param[in] totalWork - e.g. total number of loop iterations.
 * \param[in] numThreads - Number of threads that will share the work.
 * \param[in] maxChunkSize - Upper bound for chunk size.
 * \return The chunkSize.
 */
inline size_t computeStaticChunkSize(size_t totalWork,
                                     size_t numThreads,
                                     size_t maxChunkSize)
{
  if(!totalWork) return maxChunkSize;
  size_t workPerThread = roundUpDiv(totalWork, numThreads);
  size_t chunksPerThread = roundUpDiv(workPerThread, maxChunkSize);
  return roundUpDiv(workPerThread, chunksPerThread);
}

/*!
 * \brief Copy data from one array-like object to another in parallel.
 * \param[in] size - Number of elements.
 * \param[in] src - Source array.
 * \param[in] dst - Destination array.
 */
template<class T, class U>
void parallelCopy(size_t size, const T* src, U* dst)
{
  SU2_OMP_FOR_STAT(4196)
  for(size_t i=0; i<size; ++i) dst[i] = src[i];
}

/*!
 * \brief Set the entries of an array-like object to a constant value in parallel.
 * \param[in] size - Number of elements.
 * \param[in] val - Value to set.
 * \param[in] dst - Destination array.
 */
template<class T, class U>
void parallelSet(size_t size, T val, U* dst)
{
  SU2_OMP_FOR_STAT(4196)
  for(size_t i=0; i<size; ++i) dst[i] = val;
}

/*!
 * \brief Atomically update a (shared) lhs value with a (local) rhs value.
 * \note For types without atomic support (non-arithmetic) this is done via critical.
 * \param[in] rhs - Local variable being added to the shared one.
 * \param[in,out] lhs - Shared variable being updated.
 */
template<class T,
         typename std::enable_if<!std::is_arithmetic<T>::value,bool>::type = 0>
inline void atomicAdd(T rhs, T& lhs)
{
  SU2_OMP_CRITICAL
  lhs += rhs;
}
template<class T,
         typename std::enable_if<std::is_arithmetic<T>::value,bool>::type = 0>
inline void atomicAdd(T rhs, T& lhs)
{
  SU2_OMP_ATOMIC
  lhs += rhs;
}
