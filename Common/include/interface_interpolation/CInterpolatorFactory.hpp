/*!
 * \file CInterpolatorFactory.hpp
 * \brief Factory to generate interpolator objects.
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

class CConfig;
class CGeometry;
class CInterpolator;

namespace CInterpolatorFactory {
/*!
 * \brief Factory method for CInterpolator objects.
 * \ingroup Interfaces
 * \param[in] geometry_container - Geometrical definition of the problem.
 * \param[in] config - Definition of the particular problem.
 * \param[in] transpInterpolator - Transpose interpolator.
 * \param[in] iZone - Index of the donor zone.
 * \param[in] jZone - Index of the target zone.
 * \param[in] verbose - If true, print information to screen.
 * \return Pointer to interpolator on the heap, caller is responsible for deletion.
 */
CInterpolator* CreateInterpolator(CGeometry**** geometry_container, const CConfig* const* config,
                                  const CInterpolator* transpInterpolator, unsigned iZone, unsigned jZone,
                                  bool verbose = true);
}  // namespace CInterpolatorFactory
