/*!
 * \file classes_multiple_integers.hpp
 * \brief Header file for the classes that consists of multiple integer types.
 * \author E. van der Weide
 * \version 8.2.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2025, SU2 Contributors (cf. AUTHORS.md)
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

#include <iostream>

/*!
 * \struct CUnsignedLong2T
 * \brief Helper struct used to store two integral types as one entity.
 */
struct CUnsignedLong2T {
  unsigned long long0; /*!< \brief First integer to store in this class. */
  unsigned long long1; /*!< \brief Second integer to store in this class. */

  CUnsignedLong2T(unsigned long a = 0, unsigned long b = 0) : long0(a), long1(b) {}

  inline bool operator<(const CUnsignedLong2T& other) const {
    if (long0 != other.long0) return (long0 < other.long0);
    return (long1 < other.long1);
  }

  inline bool operator==(const CUnsignedLong2T& other) const {
    return (long0 == other.long0) && (long1 == other.long1);
  }
};

/*!
 * \struct CUnsignedShort2T
 * \brief Help struct used to store two integral types as one entity.
 */
struct CUnsignedShort2T {
  unsigned short short0; /*!< \brief First integer to store in this class. */
  unsigned short short1; /*!< \brief Second integer to store in this class. */

  CUnsignedShort2T(unsigned short a = 0, unsigned short b = 0) : short0(a), short1(b) {}

  inline bool operator<(const CUnsignedShort2T& other) const {
    if (short0 != other.short0) return (short0 < other.short0);
    return (short1 < other.short1);
  }

  inline bool operator==(const CUnsignedShort2T& other) const {
    return (short0 == other.short0) && (short1 == other.short1);
  }
};

/*!
 * \class CLong3T
 * \brief Help class used to store three longs as one entity.
 * \version 8.2.0 "Harrier"
 */
struct CLong3T {
  long long0 = 0; /*!< \brief First long to store in this class. */
  long long1 = 0; /*!< \brief Second long to store in this class. */
  long long2 = 0; /*!< \brief Third long to store in this class. */

  CLong3T() = default;

  CLong3T(const long a, const long b, const long c) {
    long0 = a;
    long1 = b;
    long2 = c;
  }

  bool operator<(const CLong3T& other) const {
    if (long0 != other.long0) return (long0 < other.long0);
    if (long1 != other.long1) return (long1 < other.long1);
    if (long2 != other.long2) return (long2 < other.long2);

    return false;
  }
};
