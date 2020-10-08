/*!
 * \file CMatchingFace.hpp
 * \brief Header file for the class CMatchingFace.
 *        The implementations are in the <i>CMatchingFace.cpp</i> file.
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

#include "../../mpi_structure.hpp"
#include "../../option_structure.hpp"

#include <vector>
#include <algorithm>
#include <climits>

using namespace std;

/*!
 * \class CMatchingFace
 * \brief Help class used to determine whether or not (periodic) faces match.
 */
class CMatchingFace {
public:
  unsigned short nCornerPoints;          /*!< \brief Number of corner points of the face. */
  unsigned short nDim;                   /*!< \brief Number of spatial dimensions. */
  unsigned short nPoly;                  /*!< \brief Polynomial degree of the face. */
  unsigned short nDOFsElem;              /*!< \brief Number of DOFs of the relevant adjacent element. */
  unsigned short elemType;               /*!< \brief Type of the adjacent element. */
  unsigned long  elemID;                 /*!< \brief The relevant adjacent element ID. */
  su2double cornerCoor[4][3];            /*!< \brief Coordinates of the corner points of the face. */
  su2double tolForMatching;              /*!< \brief Tolerance for this face for matching points. */

  /*!
   * \brief Default constructor.
   */
  CMatchingFace();

  /*!
   * \brief Copy contructor.
   * \param[in] - other   Object to be copied.
   */
  inline CMatchingFace(const CMatchingFace &other) { Copy(other); }

  /*!
   * \brief Assignment operator.
   * \param[in] - other   Object to be copied.
   * \return    - Reference to the current object.
   */
  inline CMatchingFace& operator=(const CMatchingFace &other) { Copy(other); return (*this); }

  /*!
   * \brief Less than operator. Needed for the sorting and searching.
   * \param[in] - other   Object to be compared to
   * \return    - True if considered less and false otherwise.
   */
  bool operator<(const CMatchingFace &other) const;

  /*!
   * \brief Function, which sorts the coordinates of the face.
   */
  void SortFaceCoordinates(void);

private:
  /*!
   * \brief Function, which copies the data of the given object into the current object
   * \param[in] - other   Object to be copied.
   */
  void Copy(const CMatchingFace &other);
};
