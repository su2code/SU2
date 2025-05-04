/*!
 * \file CMatchingFace.hpp
 * \brief Header file for the class CMatchingFace.
 *        The implementations are in the <i>CMatchingFace.cpp</i> file.
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

#include "../../code_config.hpp"
#include <iostream>
#include <algorithm>
#include <cmath>

using namespace std;

/*!
 * \class CMatchingFace
 * \brief Help class used to determine whether or not (periodic) faces match.
 */
class CMatchingFace {
 public:
  unsigned short nCornerPoints; /*!< \brief Number of corner points of the face. */
  unsigned short nDim;          /*!< \brief Number of spatial dimensions. */
  unsigned short nPoly;         /*!< \brief Polynomial degree of the face. */
  unsigned short nDOFsElem;     /*!< \brief Number of DOFs of the relevant adjacent element. */
  unsigned short elemType;      /*!< \brief Type of the adjacent element. */
  unsigned long elemID;         /*!< \brief The relevant adjacent element ID. */
  su2double cornerCoor[4][3];   /*!< \brief Coordinates of the corner points of the face. */
  su2double tolForMatching;     /*!< \brief Tolerance for this face for matching points. */

  /* Constructor. Initialize the member variables to zero. */
  CMatchingFace();

  /* Copy constructor and assignment operator. */
  inline CMatchingFace(const CMatchingFace& other) { Copy(other); }

  inline CMatchingFace& operator=(const CMatchingFace& other) {
    Copy(other);
    return (*this);
  }

  /* Less than operator. Needed for the sorting and searching. */
  bool operator<(const CMatchingFace& other) const;

  /*--- Member function, which sorts the coordinates of the face. ---*/
  void SortFaceCoordinates();

 private:
  /*--- Copy function, which copies the data of the given object into the current object. ---*/
  void Copy(const CMatchingFace& other);
};
