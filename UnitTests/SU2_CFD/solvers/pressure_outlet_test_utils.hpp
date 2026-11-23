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

#include <iostream>
#include <memory>
#include <sstream>
#include "../../../SU2_CFD/include/solvers/CEulerSolver.hpp"
#include "../../../SU2_CFD/include/fluid/CIdealGas.hpp"
#include "../../../SU2_CFD/include/fluid/CPengRobinson.hpp"
#include "../../../SU2_CFD/include/fluid/CVanDerWaalsGas.hpp"
#include "../../../Common/include/toolboxes/geometry_toolbox.hpp"

// Gases of the outlet tests: 0 ideal (gamma 1.4), 1 Peng-Robinson and 2 van der Waals (siloxane MDM).
inline std::unique_ptr<CFluidModel> MakeOutletTestGas(unsigned short kind) {
  if (kind == 0) return std::unique_ptr<CFluidModel>(new CIdealGas(1.4, 35.152));
  if (kind == 1) return std::unique_ptr<CFluidModel>(new CPengRobinson(1.0125, 35.152, 1415200, 564.1, .529));
  return std::unique_ptr<CFluidModel>(new CVanDerWaalsGas(1.0125, 35.152, 1415200, 564.1));
}

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

// Boundary state that the static-pressure outlets of a steady run use (CEulerSolver::BuildPressureOutletState):
// the outlet state, the sonic backflow state for flow entering through every characteristic, or else the interior.
// Returns true if a fallback was used.
inline bool ExpectedSteadyOutletState(CFluidModel& gas, unsigned short dim, const su2double* state,
                                      const su2double* normal, su2double gridNormal, su2double pressure,
                                      su2double turbulentK, bool idealGas, su2double gamma, su2double* boundary) {
  using Status = CEulerSolver::PressureOutletStatus;
  auto status = CEulerSolver::ComputePressureOutletState(&gas, dim, state, normal, gridNormal, pressure, turbulentK,
                                                         idealGas, gamma, boundary);
  if (status == Status::Valid) return false;
  if (status == Status::Backflow)
    status = CEulerSolver::ComputeSonicBackflowState(&gas, dim, state, normal, gridNormal, pressure, turbulentK,
                                                     idealGas, gamma, boundary);
  if (status != Status::Valid)
    for (auto k = 0u; k < dim + 2u; ++k) boundary[k] = state[k];
  return true;
}

// Fallback count of the last boundary call: the monitor sums and resets it (its warnings are discarded here).
inline unsigned long CountOutletFallbacks(CSolver* solver, const CConfig* config) {
  std::ostringstream quiet;
  auto* buffer = std::cout.rdbuf(quiet.rdbuf());
  auto* euler = static_cast<CEulerSolver*>(solver);
  euler->MonitorPressureOutletFallbacks(config);
  std::cout.rdbuf(buffer);
  return euler->GetPressureOutletFallbacks();
}
