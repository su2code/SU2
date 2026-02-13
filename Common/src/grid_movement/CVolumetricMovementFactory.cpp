/*!
 * \file CVolumetricMovementFactory.cpp
 * \brief Factory to generate volumetric mover objects.
 * \version 8.4.0 "Harrier"
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

#include "../../include/CConfig.hpp"

#include "../../include/grid_movement/CVolumetricMovementFactory.hpp"
#include "../../include/grid_movement/CRadialBasisFunctionInterpolation.hpp"
#include "../../include/grid_movement/CLinearElasticity.hpp"

namespace CVolumetricMovementFactory {

CVolumetricMovement* CreateCVolumetricMovement(CGeometry* geometry, CConfig* config) {
  const auto type = config->GetDeform_Kind();
  CVolumetricMovement* VolumetricMovement = nullptr;
  switch (type) {
    case DEFORM_KIND::RBF:
      VolumetricMovement = new CRadialBasisFunctionInterpolation(geometry, config);
      break;
    case DEFORM_KIND::ELASTIC:
      VolumetricMovement = new CLinearElasticity(geometry, config);
  }
  return VolumetricMovement;
}

}  // namespace CVolumetricMovementFactory
