/*!
 * \file CMeshElement.hpp
 * \brief Declaration and inlines of the class
 *        to define the variables of the mesh movement.
 * \author Ruben Sanchez
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

#include <cstdlib>

#include "../../../Common/include/CConfig.hpp"

class CMeshElement {
protected:

  su2double Ref_Volume;       /*!< \brief Store the reference coordinates of the mesh. */
  su2double Curr_Volume;      /*!< \brief Store the current coordinates of the mesh. */

  su2double WallDistance;     /*!< \brief Store the distance of the center of the element to the wall. */

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] val_coor - Values of the coordinates (initialization value).
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CMeshElement(void);

  /*!
   * \brief Destructor of the class.
   */
  ~CMeshElement(void) = default;

  /*!
   * \brief Get the value of the element volume with undeformed coordinates.
   * \return Value of the element volume with reference coordinates.
   */
  inline su2double GetRef_Volume(void) const { return Ref_Volume; }

  /*!
   * \brief Set the value of the element volume with undeformed coordinates.
   * \param[in] val_volume - Value of the reference volume.
   */
  inline void SetRef_Volume(su2double val_volume) { Ref_Volume = val_volume; }

  /*!
   * \brief Get the value of the element volume with deformed coordinates.
   * \return Value of the element volume with deformed coordinates.
   */
  inline su2double GetCurr_Volume(void) const { return Curr_Volume; }

  /*!
   * \brief Set the value of the element distance to the nearest wall with deformed coordinates.
   * \param[in] val_volume - Value of the element distance to the nearest wall.
   */
  inline void SetCurr_Volume(su2double val_volume) { Curr_Volume = val_volume; }

  /*!
   * \brief Get the value of the element distance to the nearest wall with undeformed coordinates.
   * \return Value of the element distance to the nearest wall with reference coordinates.
   */
  inline su2double GetWallDistance(void) const { return WallDistance; }

  /*!
   * \brief Set the value of the element distance to the nearest wall with undeformed coordinates.
   * \param[in] val_volume - Value of the element distance to the nearest wall.
   */
  inline void SetWallDistance(su2double val_volume) { WallDistance = val_volume; }

};
