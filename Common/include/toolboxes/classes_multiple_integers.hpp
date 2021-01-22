/*!
 * \file classes_multiple_integers.hpp
 * \brief Header file for the classes that consists of multiple integer types.
 * \author E. van der Weide
 * \version 7.1.0 "Blackbird"
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

  unsigned long long0 = 0;  /*!< \brief First integer to store in this class. */
  unsigned long long1 = 0;  /*!< \brief Second integer to store in this class. */

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

  unsigned short short0 = 0;  /*!< \brief First integer to store in this class. */
  unsigned short short1 = 0;  /*!< \brief Second integer to store in this class. */

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
 * \version 7.1.0 "Blackbird"
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
  CLong3T(const long a = 0, const long b = 0, const long c = 0)
    : long0(a), long1(b), long2(c) {}

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
 * \version 7.1.0 "Blackbird"
 */
struct CUnsignedShort3T {
  unsigned short short0 = 0;  /*!< \brief First short to store in this class. */
  unsigned short short1 = 0;  /*!< \brief Second short to store in this class. */
  unsigned short short2 = 0;  /*!< \brief Third short to store in this class. */

  /*!
   * \overload
   * \param[in] a - First element of the object.
   * \param[in] b - Second element of the object.
   * \param[in] c - Third element of the object
   */
  CUnsignedShort3T(const unsigned short a = 0,
                   const unsigned short b = 0,
                   const unsigned short c = 0)
    : short0(a), short1(b), short2(c) {}

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

/*!
 * \struct CUnsignedShort4T
 * \brief Help struct used to store four unsigned shorts as one entity.
 * \version 7.1.0 "Blackbird"
 */
struct CUnsignedShort4T {
  unsigned short short0 = 0;  /*!< \brief First short to store in this class. */
  unsigned short short1 = 0;  /*!< \brief Second short to store in this class. */
  unsigned short short2 = 0;  /*!< \brief Third short to store in this class. */
  unsigned short short3 = 0;  /*!< \brief Fourth short to store in this class. */

  /*!
   * \overload
   * \param[in] a - First element of the object.
   * \param[in] b - Second element of the object.
   * \param[in] c - Third element of the object.
   * \param[in] d - Fourth element of the object.
   */
  CUnsignedShort4T(const unsigned short a = 0,
                   const unsigned short b = 0,
                   const unsigned short c = 0,
                   const unsigned short d = 0)
    : short0(a), short1(b), short2(c), short3(d) {}

  /*!
   * \brief Less than operator. Needed for the sorting and searching.
   * \param[in] - other   Object to be compared to
   * \return    - True if considered less and false otherwise.
   */
  inline bool operator<(const CUnsignedShort4T &other) const {
    if(short0 != other.short0) return (short0 < other.short0);
    if(short1 != other.short1) return (short1 < other.short1);
    if(short2 != other.short2) return (short2 < other.short2);
    if(short3 != other.short3) return (short3 < other.short3);

    return false;
  }

  /*!
   * \brief Equal operator. Needed for the searching.
   * \param[in] - other   Object to be compared to
   * \return    - True if considered equal and false otherwise.
   */
  inline bool operator==(const CUnsignedShort4T &other) const {
    return (short0 == other.short0) && (short1 == other.short1) &&
           (short2 == other.short2) && (short3 == other.short3);
  }
};

/*!
 * \struct CUnsignedShort8T
 * \brief Help struct used to store eight unsigned shorts as one entity.
 * \version 7.1.0 "Blackbird"
 */
struct CUnsignedShort8T {
  unsigned short short0 = 0;  /*!< \brief First short to store in this class. */
  unsigned short short1 = 0;  /*!< \brief Second short to store in this class. */
  unsigned short short2 = 0;  /*!< \brief Third short to store in this class. */
  unsigned short short3 = 0;  /*!< \brief Fourth short to store in this class. */
  unsigned short short4 = 0;  /*!< \brief Fifth short to store in this class. */
  unsigned short short5 = 0;  /*!< \brief Sixth short to store in this class. */
  unsigned short short6 = 0;  /*!< \brief Seventh short to store in this class. */
  unsigned short short7 = 0;  /*!< \brief Eighth short to store in this class. */

  /*!
   * \overload
   * \param[in] a - First element of the object.
   * \param[in] b - Second element of the object.
   * \param[in] c - Third element of the object.
   * \param[in] d - Fourth element of the object.
   * \param[in] e - Fifth element of the object.
   * \param[in] f - Sixth element of the object.
   * \param[in] g - Seventh element of the object.
   * \param[in] h - Eighth element of the object.
   */
  CUnsignedShort8T(const unsigned short a = 0,
                   const unsigned short b = 0,
                   const unsigned short c = 0,
                   const unsigned short d = 0,
                   const unsigned short e = 0,
                   const unsigned short f = 0,
                   const unsigned short g = 0,
                   const unsigned short h = 0)
    : short0(a), short1(b), short2(c), short3(d),
      short4(e), short5(f), short6(g), short7(h) {}

  /*!
   * \brief Less than operator. Needed for the sorting and searching.
   * \param[in] - other   Object to be compared to
   * \return    - True if considered less and false otherwise.
   */
  inline bool operator<(const CUnsignedShort8T &other) const {
    if(short0 != other.short0) return (short0 < other.short0);
    if(short1 != other.short1) return (short1 < other.short1);
    if(short2 != other.short2) return (short2 < other.short2);
    if(short3 != other.short3) return (short3 < other.short3);
    if(short4 != other.short4) return (short4 < other.short4);
    if(short5 != other.short5) return (short5 < other.short5);
    if(short6 != other.short6) return (short6 < other.short6);
    if(short7 != other.short7) return (short7 < other.short7);

    return false;
  }

  /*!
   * \brief Equal operator. Needed for the searching.
   * \param[in] - other   Object to be compared to
   * \return    - True if considered equal and false otherwise.
   */
  inline bool operator==(const CUnsignedShort8T &other) const {
    return (short0 == other.short0) && (short1 == other.short1) &&
           (short2 == other.short2) && (short3 == other.short3) &&
           (short4 == other.short4) && (short5 == other.short5) &&
           (short6 == other.short6) && (short7 == other.short7);
  }
};

/*!
 * \struct CUnsignedShort11T
 * \brief Help struct used to store eleven unsigned shorts as one entity.
 * \version 7.1.0 "Blackbird"
 */
struct CUnsignedShort11T {
  unsigned short short0  = 0;  /*!< \brief First short to store in this class. */
  unsigned short short1  = 0;  /*!< \brief Second short to store in this class. */
  unsigned short short2  = 0;  /*!< \brief Third short to store in this class. */
  unsigned short short3  = 0;  /*!< \brief Fourth short to store in this class. */
  unsigned short short4  = 0;  /*!< \brief Fifth short to store in this class. */
  unsigned short short5  = 0;  /*!< \brief Sixth short to store in this class. */
  unsigned short short6  = 0;  /*!< \brief Seventh short to store in this class. */
  unsigned short short7  = 0;  /*!< \brief Eighth short to store in this class. */
  unsigned short short8  = 0;  /*!< \brief Nineth short to store in this class. */
  unsigned short short9  = 0;  /*!< \brief Tenth short to store in this class. */
  unsigned short short10 = 0;  /*!< \brief Eleventh short to store in this class. */

  /*!
   * \brief Less than operator. Needed for the sorting and searching.
   * \param[in] - other   Object to be compared to
   * \return    - True if considered less and false otherwise.
   */
  inline bool operator<(const CUnsignedShort11T &other) const {
    if(short0  != other.short0)  return (short0  < other.short0);
    if(short1  != other.short1)  return (short1  < other.short1);
    if(short2  != other.short2)  return (short2  < other.short2);
    if(short3  != other.short3)  return (short3  < other.short3);
    if(short4  != other.short4)  return (short4  < other.short4);
    if(short5  != other.short5)  return (short5  < other.short5);
    if(short6  != other.short6)  return (short6  < other.short6);
    if(short7  != other.short7)  return (short7  < other.short7);
    if(short8  != other.short8)  return (short8  < other.short8);
    if(short9  != other.short9)  return (short9  < other.short9);
    if(short10 != other.short10) return (short10 < other.short10);

    return false;
  }

  /*!
   * \brief Equal operator. Needed for the searching.
   * \param[in] - other   Object to be compared to
   * \return    - True if considered equal and false otherwise.
   */
  inline bool operator==(const CUnsignedShort11T &other) const {
    return (short0  == other.short0) && (short1 == other.short1) &&
           (short2  == other.short2) && (short3 == other.short3) &&
           (short4  == other.short4) && (short5 == other.short5) &&
           (short6  == other.short6) && (short7 == other.short7) &&
           (short8  == other.short8) && (short9 == other.short9) &&
           (short10 == other.short10);
  }
};
