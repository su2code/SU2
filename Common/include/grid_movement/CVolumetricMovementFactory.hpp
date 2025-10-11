/*!
 * \file CVolumetricMovementFactory.hpp
 * \brief Factory to generate volumetric mover objects.
 * \version 8.3.0 "Harrier"
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

class CConfig;
class CGeometry;
class CVolumetricMovement;

namespace CVolumetricMovementFactory {
/*!
 * \brief Factory method for CVolumetricMovement objects.
 * \param[in] geometry_container - Geometrical definition of the problem.
 * \param[in] config - Definition of the particular problem.
 * \return Pointer to the allocated volumetric mover.
 */

CVolumetricMovement* CreateCVolumetricMovement(CGeometry* geometry, CConfig* config);
}  // namespace CVolumetricMovementFactory
