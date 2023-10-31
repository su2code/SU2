/*!
 * \file CBBoxTargetClass.hpp
 * \brief Class for storing the information of a possible bounding box candidate
           during a minimum distance search.
 * \author E. van der Weide
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
#include "../basic_types/datatype_structure.hpp"

/*!
 * \class CBBoxTargetClass
 * \ingroup ADT
 * \brief  Class for storing the information of a possible bounding box candidate
           during a minimum distance search.
 * \author E. van der Weide
 * \version 8.0.0 "Harrier"
 */
struct CBBoxTargetClass {
  unsigned long boundingBoxID;  /*!< \brief Corresponding bounding box ID. */
  su2double possibleMinDist2;   /*!< \brief Possible minimimum distance squared to the
                                            given coordinate. */
  su2double guaranteedMinDist2; /*!< \brief Guaranteed minimum distance squared to the
                                            given coordinate. */

  /*!
   * \brief Constructor of the class. Nothing to be done.
   */
  CBBoxTargetClass() = default;

  /*!
   * \brief Constructor of the class, which initializes the member variables.
   * \param[in] val_BBoxID    - ID of the bounding box to be stored.
   * \param[in] val_posDist2  - Possible minimum distance squared to the target
                                for this bounding box.
   * \param[in] val_guarDist2 - Guaranteed minimum distance squared to the target
                                for this bounding box.
   */
  inline CBBoxTargetClass(const unsigned long val_BBoxID, const su2double val_posDist2, const su2double val_guarDist2)
      : boundingBoxID(val_BBoxID), possibleMinDist2(val_posDist2), guaranteedMinDist2(val_guarDist2) {}

  /*!
   * \brief Less than operator. Needed for the sorting of the candidates.
   * \param[in] other  Object to which the current object must be compared.
   */
  inline bool operator<(const CBBoxTargetClass& other) const {
    /* Make sure that the bounding boxes with the smallest possible distances
       are stored first. */
    if (possibleMinDist2 < other.possibleMinDist2) return true;
    if (possibleMinDist2 > other.possibleMinDist2) return false;
    if (guaranteedMinDist2 < other.guaranteedMinDist2) return true;
    if (guaranteedMinDist2 > other.guaranteedMinDist2) return false;
    if (boundingBoxID < other.boundingBoxID) return true;
    return false;
  }
};
