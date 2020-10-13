/*!
 * \file classes_multiple_integers.hpp
 * \brief Header file for the classes that consists of multiple integer types.
 * \author E. van der Weide
 * \version 7.0.7 "Blackbird"
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

#include <iostream>

/*!
 * \struct CUnsignedLong2T
 * \brief Helper struct used to store two integral types as one entity.
 */
struct CUnsignedLong2T {

  unsigned long long0;  /*!< \brief First integer to store in this class. */
  unsigned long long1;  /*!< \brief Second integer to store in this class. */

  /*!
   * \overload
   * \param[in] a - First element of the object.
   * \param[in] b - Second element of the object.
   */
  CUnsignedLong2T(unsigned long a = 0, unsigned long b = 0) : long0(a), long1(b) {}

  /*!
   * \brief Less than operator. Needed for the sorting and searching.
   * \param[in] - other   Object to be compared to
   * \return    - True if considered less and false otherwise.
   */
  inline bool operator<(const CUnsignedLong2T &other) const {
    if(long0 != other.long0)
      return (long0 < other.long0);
    return (long1 < other.long1);
  }

  /*!
   * \brief Equal operator. Needed for the searching.
   * \param[in] - other   Object to be compared to
   * \return    - True if considered equal and false otherwise.
   */
  inline bool operator==(const CUnsignedLong2T &other) const {
    return (long0 == other.long0) && (long1 == other.long1);
  }
};

/*!
 * \struct CUnsignedShort2T
 * \brief Help struct used to store two integral types as one entity.
 */
struct CUnsignedShort2T {

  unsigned short short0;  /*!< \brief First integer to store in this class. */
  unsigned short short1;  /*!< \brief Second integer to store in this class. */

  /*!
   * \overload
   * \param[in] a - First element of the object.
   * \param[in] b - Second element of the object.
   */
  CUnsignedShort2T(unsigned short a = 0, unsigned short b = 0) : short0(a), short1(b) {}

  /*!
   * \brief Less than operator. Needed for the sorting and searching.
   * \param[in] - other   Object to be compared to
   * \return    - True if considered less and false otherwise.
   */
  inline bool operator<(const CUnsignedShort2T &other) const {
    if(short0 != other.short0)
      return (short0 < other.short0);
    return (short1 < other.short1);
  }

  /*!
   * \brief Equal operator. Needed for the searching.
   * \param[in] - other   Object to be compared to
   * \return    - True if considered equal and false otherwise.
   */
  inline bool operator==(const CUnsignedShort2T &other) const {
    return (short0 == other.short0) && (short1 == other.short1);
  }
};

/*!
 * \struct CLong3T
 * \brief Help struct used to store three longs as one entity.
 * \version 7.0.7 "Blackbird"
 */
struct CLong3T {
  long long0 = 0;  /*!< \brief First long to store in this class. */
  long long1 = 0;  /*!< \brief Second long to store in this class. */
  long long2 = 0;  /*!< \brief Third long to store in this class. */

  /*!
   * \overload
   * \param[in] a - First element of the object.
   * \param[in] b - Second element of the object.
   * \param[in] c - Third element of the object
   */
  CLong3T(const long a, const long b, const long c) {long0 = a; long1 = b; long2 = c;}

  /*!
   * \brief Less than operator. Needed for the sorting and searching.
   * \param[in] - other   Object to be compared to
   * \return    - True if considered less and false otherwise.
   */
  inline bool operator<(const CLong3T &other) const {
    if(long0 != other.long0) return (long0 < other.long0);
    if(long1 != other.long1) return (long1 < other.long1);
    if(long2 != other.long2) return (long2 < other.long2);

    return false;
  }
};

/*!
 * \struct CUnsignedShort3T
 * \brief Help struct used to store three unsigned shorts as one entity.
 * \version 7.0.7 "Blackbird"
 */
struct CUnsignedShort3T {
  unsigned short short0;  /*!< \brief First short to store in this class. */
  unsigned short short1;  /*!< \brief Second short to store in this class. */
  unsigned short short2;  /*!< \brief Third short to store in this class. */

  /*!
   * \overload
   * \param[in] a - First element of the object.
   * \param[in] b - Second element of the object.
   * \param[in] c - Third element of the object
   */
  CUnsignedShort3T(const unsigned short a = 0,
                   const unsigned short b = 0,
                   const unsigned short c = 0) {
    short0 = a; short1 = b; short2 = c;
  }

  /*!
   * \brief Less than operator. Needed for the sorting and searching.
   * \param[in] - other   Object to be compared to
   * \return    - True if considered less and false otherwise.
   */
  inline bool operator<(const CUnsignedShort3T &other) const {
    if(short0 != other.short0) return (short0 < other.short0);
    if(short1 != other.short1) return (short1 < other.short1);
    if(short2 != other.short2) return (short2 < other.short2);

    return false;
  }

  /*!
   * \brief Equal operator. Needed for the searching.
   * \param[in] - other   Object to be compared to
   * \return    - True if considered equal and false otherwise.
   */
  inline bool operator==(const CUnsignedShort3T &other) const {
    return (short0 == other.short0) && (short1 == other.short1) &&
           (short2 == other.short2);
  }
};
