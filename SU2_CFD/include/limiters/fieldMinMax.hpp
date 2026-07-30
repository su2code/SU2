/*!
 * \file fieldMinMax.hpp
 * \brief Small helpers to track the min/max of a field over the neighbors of a point.
 * \note These allow the extrema to be gathered by whichever neighbor loop is convenient,
 *       the one of the limiter or the one of the Green-Gauss gradient.
 * \author P. Gomes
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

#include "../../../Common/include/basic_types/datatype_structure.hpp"

/*!
 * \brief Tracks the min/max of a field over the neighbors of one point, in the rows of the
 *        min/max matrices of a solver.
 * \ingroup FvmAlgos
 */
template <class FieldType>
struct CRowMinMax {
  FieldType& fieldMin;
  FieldType& fieldMax;
  size_t iPoint;
  FORCEINLINE void update(size_t iVar, const su2double& val) {
    fieldMin(iPoint, iVar) = fmin(fieldMin(iPoint, iVar), val);
    fieldMax(iPoint, iVar) = fmax(fieldMax(iPoint, iVar), val);
  }
};

/*!
 * \brief As above but in plain per-variable arrays, for when the extrema are not stored.
 * \ingroup FvmAlgos
 */
struct CPointMinMax {
  su2double* minVals;
  su2double* maxVals;
  FORCEINLINE void update(size_t iVar, const su2double& val) {
    minVals[iVar] = fmin(minVals[iVar], val);
    maxVals[iVar] = fmax(maxVals[iVar], val);
  }
};

/*!
 * \brief Stand-in for when the extrema of the field have already been determined elsewhere.
 * \ingroup FvmAlgos
 */
struct CNoMinMax {
  FORCEINLINE void update(size_t, const su2double&) {}
};
