/*!
 * \file allocation_toolbox.hpp
 * \brief Helper function and classes for memory allocation.
 *        Focus on portability across platforms.
 * \note  These are "kernel" functions, only to be used with good reason,
 *        always try to use higher level container classes.
 * \author P. Gomes, D. Kavolis
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

#if defined(_WIN32)
#include <malloc.h>
#else
#include <stdlib.h>
#endif

#include <cstring>

#include <cassert>

namespace MemoryAllocation {

inline constexpr bool is_power_of_two(size_t x) { return x && !(x & (x - 1)); }

inline constexpr size_t round_up(size_t multiple, size_t x) { return ((x + multiple - 1) / multiple) * multiple; }

/*!
 * \brief Aligned memory allocation compatible across platforms.
 * \param[in] alignment, in bytes, of the memory being allocated.
 * \param[in] size, also in bytes.
 * \tparam ZeroInit, initialize memory to 0.
 * \return Pointer to memory, always use su2::aligned_free to deallocate.
 */
template <class T, bool ZeroInit = false>
inline T* aligned_alloc(size_t alignment, size_t size) noexcept {
  assert(is_power_of_two(alignment));

  if (alignment < alignof(void*)) alignment = alignof(void*);

  size = round_up(alignment, size);

  void* ptr = nullptr;

#if defined(__APPLE__)
  if (::posix_memalign(&ptr, alignment, size) != 0) {
    ptr = nullptr;
  }
#elif defined(_WIN32)
  ptr = _aligned_malloc(size, alignment);
#else
  ptr = ::aligned_alloc(alignment, size);
#endif
  if (ZeroInit) memset(ptr, 0, size);
  return static_cast<T*>(ptr);
}

/*!
 * \brief Free memory allocated with su2::aligned_alloc.
 * \param[in] ptr, pointer to memory we want to release.
 */
template <class T>
inline void aligned_free(T* ptr) noexcept {
#if defined(_WIN32)
  _aligned_free(ptr);
#else
  free(ptr);
#endif
}

}  // namespace MemoryAllocation
