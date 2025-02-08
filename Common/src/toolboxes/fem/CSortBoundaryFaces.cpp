/*!
 * \file CSortFaces.cpp
 * \brief Functor, used for a different sorting of the surface elements
 *        than the < operator of CSurfaceElementFEM.
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

#include "../../../include/toolboxes/fem/CSortBoundaryFaces.hpp"
#include "../../../include/fem/fem_geometry_structure.hpp"

bool CSortBoundaryFaces::operator()(const CSurfaceElementFEM& f0, const CSurfaceElementFEM& f1) {
  /* First sorting criterion is the index of the standard element. The
     boundary faces should be sorted per standard element. Note that the
     time level is not taken into account here, because it is assumed that
     the surface elements to be sorted belong to one time level. */
  if (f0.indStandardElement != f1.indStandardElement) return f0.indStandardElement < f1.indStandardElement;

  /* The standard elements are the same. The second criterion is the
     corresponding volume IDs of the surface elements. */
  return f0.volElemID < f1.volElemID;
}
