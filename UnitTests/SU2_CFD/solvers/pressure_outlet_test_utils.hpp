/*!
 * \file pressure_outlet_test_utils.hpp
 * \brief Independent wave increments for pressure-outlet matrix tests.
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

#include "../../../SU2_CFD/include/fluid/CFluidModel.hpp"
#include "../../../Common/include/toolboxes/geometry_toolbox.hpp"

// Independent primitive-variable wave decomposition for checking the implicit outlet projection.
inline void OutgoingOutletIncrement(CFluidModel& gas, unsigned short dim, const su2double* state,
                                    const su2double* normal, su2double gridNormal, su2double turbulentK,
                                    unsigned short column, su2double* increment) {
  const auto rho = state[0];
  su2double velocity[3] = {}, dv[3] = {};
  const su2double drho = column == 0 ? 1.0 : 0.0;
  for (auto d = 0u; d < dim; ++d) {
    velocity[d] = state[d + 1] / rho;
    dv[d] = ((column == d + 1 ? 1.0 : 0.0) - velocity[d] * drho) / rho;
  }
  const auto energy = state[dim + 1] / rho;
  gas.SetTDState_rhoe(rho, energy - .5 * GeometryToolbox::SquaredNorm(dim, velocity) - turbulentK);
  const auto sound = gas.GetSoundSpeed(), dpde = gas.GetdPde_rho(), dpdrho = gas.GetdPdrho_e();
  const su2double de =
      ((column == dim + 1 ? 1.0 : 0.0) - energy * drho) / rho - GeometryToolbox::DotProduct(dim, velocity, dv);
  const su2double dp = dpdrho * drho + dpde * de;
  const su2double q = GeometryToolbox::DotProduct(dim, velocity, normal) - gridNormal;
  const su2double dq = GeometryToolbox::DotProduct(dim, dv, normal);
  su2double projectedRho = 0, projectedP = 0, projectedVelocity[3] = {};
  if (q >= 0) {
    projectedRho = drho - dp / (sound * sound);
    for (auto d = 0u; d < dim; ++d) projectedVelocity[d] = dv[d] - dq * normal[d];
  }
  for (const auto sign : {-1, 1}) {
    if (q + sign * sound < 0) continue;
    const su2double amplitude = .5 * (dp + sign * rho * sound * dq);
    projectedRho += amplitude / (sound * sound);
    projectedP += amplitude;
    for (auto d = 0u; d < dim; ++d) projectedVelocity[d] += sign * amplitude / (rho * sound) * normal[d];
  }
  increment[0] = projectedRho;
  for (auto d = 0u; d < dim; ++d) increment[d + 1] = velocity[d] * projectedRho + rho * projectedVelocity[d];
  increment[dim + 1] = energy * projectedRho + rho * (projectedP - dpdrho * projectedRho) / dpde +
                       rho * GeometryToolbox::DotProduct(dim, velocity, projectedVelocity);
}

// Boundary state that CEulerSolver's static-pressure outlets should use, assembled from the fluid-model maps.
// reservoir = {T_total, P_total, direction_x, direction_y, direction_z}, or nullptr without MARKER_OUTLET_BACKFLOW.
inline CFluidModel::PressureOutletStatus ExpectedOutletState(CFluidModel& gas, unsigned short dim,
                                                             const su2double* state, const su2double* normal,
                                                             su2double gridNormal, su2double pressure,
                                                             su2double* boundary, su2double turbulentK, su2double gamma,
                                                             const su2double* reservoir, bool steady) {
  auto status = gas.ComputePressureOutletState(dim, state, normal, gridNormal, pressure, boundary, turbulentK, gamma);
  if (reservoir != nullptr && status != CFluidModel::PressureOutletStatus::InvalidState) {
    // Inward flow relative to the grid comes from the reservoir.
    if (status == CFluidModel::PressureOutletStatus::Backflow ||
        GeometryToolbox::DotProduct(dim, boundary + 1, normal) / boundary[0] < gridNormal)
      // A zero total pressure means the imposed static pressure.
      status =
          gas.ComputeReservoirInflowState(dim, state, normal, gridNormal, reservoir[1] > 0 ? reservoir[1] : pressure,
                                          reservoir[0], reservoir + 2, boundary, turbulentK, gamma);
  } else if (status == CFluidModel::PressureOutletStatus::Backflow && steady) {
    status = gas.ComputeSonicBackflowState(dim, state, normal, gridNormal, pressure, boundary, turbulentK, gamma);
  }
  return status;
}
