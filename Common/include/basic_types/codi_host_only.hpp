/*!
 * \file codi_host_only.hpp
 * \brief Keeps CoDiPack out of nvcc's device pass.
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

/*--- Must be the first include of every .cu, before anything can pull in codi.hpp.
 *
 * When __CUDA_ARCH__ is defined, i.e. only in the device pass, CoDiPack stamps
 * "__device__ __host__" on every one of its functions. nvcc then has to generate device
 * code for the entire tape, which does not work: the tape uses function-scope statics with
 * dynamic initializers (illegal in device code), takes the address of host statics, and
 * calls into std::map and std::bitset. The kernels never touch an active type, so none of
 * that is wanted in the first place.
 *
 * CODI_INLINE expands CODI_CUDAFunctionAttributes at each declaration, so emptying the
 * macro before CoDiPack is parsed leaves it host-only. Including the defining header first
 * makes the copy in config.h a no-op (it is "#pragma once"), so this definition is the one
 * that survives. The host pass is unaffected either way, it never defines __CUDA_ARCH__
 * and so already sees the same host-only declarations as the .cpp translation units. ---*/
#if defined(CODI_REVERSE_TYPE) || defined(CODI_FORWARD_TYPE)
#include "codi/tools/cuda/cudaFunctionAttributes.hpp"

#undef CODI_CUDAFunctionAttributes
#define CODI_CUDAFunctionAttributes
#endif
