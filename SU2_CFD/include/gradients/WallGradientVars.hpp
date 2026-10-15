/*!
 * \file WallGradientVars.hpp
 * \brief Description of the variables with known derivatives at no-slip walls.
 * \version 8.5.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md)
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

/*!
 * \brief Variables of a gradient array that have a known tangential derivative at no-slip walls.
 * \ingroup FvmAlgos
 */
struct WallGradientVars {
  int idxVel = -1;          /*!< \brief First velocity component (zero at the wall), -1 if none. */
  bool divFreeVel = false;  /*!< \brief Also impose div(u) = 0 at the wall. */
  int idxTemp = -1;         /*!< \brief Temperature (uniform at isothermal walls), -1 if none. */
  int zeroBegin = 0;        /*!< \brief First variable that is zero at smooth walls. */
  int zeroEnd = 0;          /*!< \brief One past the last variable that is zero at smooth walls. */

  bool Active() const { return idxVel >= 0 || idxTemp >= 0 || zeroEnd > zeroBegin; }
};
