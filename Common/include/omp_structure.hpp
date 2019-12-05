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
 * \version 7.0.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2019, SU2 Contributors (cf. AUTHORS.md)
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

#if defined(_MSC_VER)
#define PRAGMIZE(X) __pragma(X)
#else
#define PRAGMIZE(X) _Pragma(#X)
#endif

/*--- Detect compilation with OpenMP support. ---*/
#ifdef _OPENMP
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
inline int omp_get_max_threads(void) {return 1;}

/*!
 * \brief Index of current thread, akin to MPI rank.
 */
inline int omp_get_thread_num(void) {return 0;}

#endif

/*--- Convenience macros (do not use excessive nesting of macros). ---*/
#define SU2_OMP_SIMD SU2_OMP(simd)

#define SU2_OMP_MASTER SU2_OMP(master)
#define SU2_OMP_BARRIER SU2_OMP(barrier)

#define SU2_OMP_PARALLEL SU2_OMP(parallel)
#define SU2_OMP_PARALLEL_(ARGS) SU2_OMP(parallel ARGS)
#define SU2_OMP_PARALLEL_ON(NTHREADS) SU2_OMP(parallel num_threads(NTHREADS))

#define SU2_OMP_FOR_DYN(CHUNK) SU2_OMP(for schedule(dynamic,CHUNK))
#define SU2_OMP_FOR_STAT(CHUNK) SU2_OMP(for schedule(static,CHUNK))

#define SU2_OMP_PAR_FOR_DYN(CHUNK) SU2_OMP(parallel for schedule(dynamic,CHUNK))
#define SU2_OMP_PAR_FOR_STAT(CHUNK) SU2_OMP(parallel for schedule(static,CHUNK))
