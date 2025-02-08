/*!
 * \file CSortBoundaryFaces.hpp
 * \brief Header file for the class CSortBoundaryFaces.
 *        The implementations are in the <i>CSortBoundaryFaces.cpp</i> file.
 * \author E. van der Weide
 * \version 8.1.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2024, SU2 Contributors (cf. AUTHORS.md)
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
 * \class CSortBoundaryFaces
 * \brief Functor, used for a different sorting of the faces than the < operator
 *        of CSurfaceElementFEM.
 * \author E. van der Weide
 * \version 8.1.0 "Harrier"
 */
struct CSurfaceElementFEM;  // Forward declaration to avoid problems.
struct CSortBoundaryFaces {
  /*!
   * \brief Operator used for the comparison.
   * \param[in] f0 - First boundary face in the comparison.
   * \param[in] f1 - Second boundary face in the comparison.
   */
  bool operator()(const CSurfaceElementFEM& f0, const CSurfaceElementFEM& f1);
};
