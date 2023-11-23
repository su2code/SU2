/*!
 * \file CADTComparePointClass.hpp
 * \brief subroutines for comparing two points in an alternating digital tree (ADT).
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
 * \class CADTComparePointClass
 * \ingroup ADT
 * \brief  Functor, used for the sorting of the points when building an ADT.
 * \author E. van der Weide
 */
class CADTComparePointClass {
 private:
  const su2double* pointCoor;          /*!< \brief Pointer to the coordinates of the points. */
  const unsigned short splitDirection; /*!< \brief Split direction used in the sorting. */
  const unsigned short nDim;           /*!< \brief Number of spatial dimensions stored in the coordinates. */

 public:
  /*!
   * \brief Constructor of the class. The member variables are initialized.
   * \param[in] coor      Pointer to the coordinates of the points.
   * \param[in] splitDir  Direction that must be used to sort the coordinates.
   * \param[in] nDimADT   Number of spatial dimensions of the ADT and coordinates.
   */
  CADTComparePointClass(const su2double* coor, const unsigned short splitDir, const unsigned short nDimADT)
      : pointCoor(coor), splitDirection(splitDir), nDim(nDimADT) {}

  /*!
   * \brief Operator used for the sorting of the points.
   * \param[in] p0  Index of the first point to be compared.
   * \param[in] p1  Index of the second point to be compared.
   */
  inline bool operator()(const unsigned long p0, const unsigned long p1) const {
    return pointCoor[nDim * p0 + splitDirection] < pointCoor[nDim * p1 + splitDirection];
  }

  /*!
   * \brief Default constructor of the class, disabled.
   */
  CADTComparePointClass() = delete;
};
