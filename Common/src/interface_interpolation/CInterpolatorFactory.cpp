/*!
 * \file CInterpolatorFactory.cpp
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

#include "../../include/CConfig.hpp"

#include "../../include/interface_interpolation/CInterpolatorFactory.hpp"
#include "../../include/interface_interpolation/CIsoparametric.hpp"
#include "../../include/interface_interpolation/CMirror.hpp"
#include "../../include/interface_interpolation/CNearestNeighbor.hpp"
#include "../../include/interface_interpolation/CRadialBasisFunction.hpp"
#include "../../include/interface_interpolation/CSlidingMesh.hpp"

namespace CInterpolatorFactory {
CInterpolator* CreateInterpolator(CGeometry**** geometry_container, const CConfig* const* config,
                                  const CInterpolator* transpInterpolator, unsigned iZone, unsigned jZone,
                                  bool verbose) {
  CInterpolator* interpolator = nullptr;

  /*--- Only print information on master node. ---*/
  verbose &= (SU2_MPI::GetRank() == MASTER_NODE);

  /*--- Type of interpolation defined by donor. ---*/
  const auto type = config[iZone]->GetKindInterpolation();

  if (verbose) cout << " Setting coupling ";

  /*--- Conservative interpolation is not applicable to the sliding
   *    mesh approach so that case is handled first. Then we either
   *    return a CMirror if the target requires conservative inter-
   *    polation, or the type of interpolator defined by "type". ---*/

  if (type == INTERFACE_INTERPOLATOR::WEIGHTED_AVERAGE) {
    if (verbose) cout << "using a sliding mesh approach." << endl;
    interpolator = new CSlidingMesh(geometry_container, config, iZone, jZone);
  } else if (config[jZone]->GetConservativeInterpolation()) {
    if (verbose) cout << "using the mirror approach, \"transposing\" coefficients from opposite mesh." << endl;
    interpolator = new CMirror(geometry_container, config, transpInterpolator, iZone, jZone);
  } else {
    switch (type) {
      case INTERFACE_INTERPOLATOR::ISOPARAMETRIC:
        if (verbose) cout << "using the isoparametric approach." << endl;
        interpolator = new CIsoparametric(geometry_container, config, iZone, jZone);
        break;

      case INTERFACE_INTERPOLATOR::NEAREST_NEIGHBOR:
        if (verbose) cout << "using a nearest neighbor approach." << endl;
        interpolator = new CNearestNeighbor(geometry_container, config, iZone, jZone);
        break;

      case INTERFACE_INTERPOLATOR::RADIAL_BASIS_FUNCTION:
        if (verbose) cout << "using a radial basis function approach." << endl;
        interpolator = new CRadialBasisFunction(geometry_container, config, iZone, jZone);
        break;

      default:
        SU2_MPI::Error("Unknown type of interpolation.", CURRENT_FUNCTION);
    }
  }

  if (verbose) interpolator->PrintStatistics();

  return interpolator;
}
}  // namespace CInterpolatorFactory
