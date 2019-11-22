/*!
 * \file allocation_toolbox.hpp
 * \brief Helper function and classes for memory allocation.
 *        Focus on portability across platforms.
 * \note  These are "kernel" functions, only to be used with good reason,
 *        always try to use higher level container classes.
 * \author P. Gomes, D. Kavolis
 * \version 6.2.0 "Falcon"
 *
 * The current SU2 release has been coordinated by the
 * SU2 International Developers Society <www.su2devsociety.org>
 * with selected contributions from the open-source community.
 *
 * The main research teams contributing to the current release are:
 *  - Prof. Juan J. Alonso's group at Stanford University.
 *  - Prof. Piero Colonna's group at Delft University of Technology.
 *  - Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *  - Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *  - Prof. Rafael Palacios' group at Imperial College London.
 *  - Prof. Vincent Terrapon's group at the University of Liege.
 *  - Prof. Edwin van der Weide's group at the University of Twente.
 *  - Lab. of New Concepts in Aeronautics at Tech. Institute of Aeronautics.
 *
 * Copyright 2012-2019, Francisco D. Palacios, Thomas D. Economon,
 *                      Tim Albring, and the SU2 contributors.
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

#if defined(_WIN32)
#include <malloc.h>
#else
#include <stdlib.h> 
#endif

#include <cassert>

namespace MemoryAllocation
{

inline constexpr bool is_power_of_two(size_t x)
{
  return x && !(x & (x-1));
}

/*!
 * \brief Aligned memory allocation compatible across platforms.
 * \param[in] alignment, in bytes, of the memory being allocated.
 * \param[in] size, also in bytes.
 * \return Pointer to memory, always use su2::aligned_free to deallocate.
 */
template<class T>
inline T* aligned_alloc(size_t alignment, size_t size) noexcept
{
  assert(is_power_of_two(alignment));

  if(alignment < alignof(void*)) alignment = alignof(void*);

  void* ptr = nullptr;

#if defined(__APPLE__)
  if(::posix_memalign(&ptr, alignment, size) != 0)
  {
    ptr = nullptr;
  }
#elif defined(_WIN32)
  ptr = _aligned_malloc(size, alignment);
#else
  ptr = ::aligned_alloc(alignment, size);
#endif
  return static_cast<T*>(ptr);
}

/*!
 * \brief Free memory allocated with su2::aligned_alloc.
 * \param[in] ptr, pointer to memory we want to release.
 */
template<class T>
inline void aligned_free(T* ptr) noexcept
{
#if defined(_WIN32)
  _aligned_free(ptr);
#else
  free(ptr);
#endif
}

} // namespace

