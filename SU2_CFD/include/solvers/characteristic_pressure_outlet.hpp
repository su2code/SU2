/*!
 * \file characteristic_pressure_outlet.hpp
 * \brief Nonlinear characteristic pressure-outlet reconstruction.
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

#include "../fluid/CFluidModel.hpp"
#include "../../../Common/include/toolboxes/geometry_toolbox.hpp"
#include <array>
#include <cmath>

namespace PressureOutlet {
enum class Status { Valid, InvalidState, Backflow };

/*--- Integrate the outgoing acoustic invariant on the interior isentrope. A
 * pressure below the sonic limit cannot impose an incoming characteristic. ---*/
inline Status State(CFluidModel& fluid, unsigned short nDim, const su2double* interior, const su2double* normal,
                    su2double gridNormal, su2double targetPressure, su2double* boundary, su2double turbulentK = 0) {
  for (auto v = 0u; v < nDim + 2; ++v)
    if (!std::isfinite(SU2_TYPE::GetValue(interior[v]))) return Status::InvalidState;
  const auto rho = interior[0];
  if (!(rho > 0) || !(targetPressure > 0) || !std::isfinite(SU2_TYPE::GetValue(targetPressure)))
    return Status::InvalidState;
  su2double velocity[3] = {};
  for (auto d = 0u; d < nDim; ++d) velocity[d] = interior[d + 1] / rho;
  const auto energy = interior[nDim + 1] / rho - .5 * GeometryToolbox::SquaredNorm(nDim, velocity) - turbulentK;
  fluid.SetTDState_rhoe(rho, energy);
  const auto pressure = fluid.GetPressure(), sound2 = fluid.GetSoundSpeed2();
  auto admissible = [&]() {
    const auto p = SU2_TYPE::GetValue(fluid.GetPressure()), t = SU2_TYPE::GetValue(fluid.GetTemperature());
    const auto a2 = SU2_TYPE::GetValue(fluid.GetSoundSpeed2());
    return p > 0 && t > 0 && a2 > 0 && std::isfinite(p + t + a2);
  };
  if (!admissible()) return Status::InvalidState;
  const auto sound = sqrt(sound2);
  const auto normalVelocity = GeometryToolbox::DotProduct(nDim, velocity, normal) - gridNormal;
  if (normalVelocity < 0) return Status::Backflow;
  if (normalVelocity >= sound) {
    for (auto v = 0u; v < nDim + 2; ++v) boundary[v] = interior[v];
    return Status::Valid;
  }

  using Pair = std::array<su2double, 2>;
  Pair state = {rho, normalVelocity};
  auto evaluate = [&](su2double logPressure, const Pair& y, Pair& derivative) {
    if (!(y[0] > 0) || !std::isfinite(SU2_TYPE::GetValue(y[0] + y[1]))) return false;
    const auto p = exp(logPressure);
    fluid.SetTDState_Prho(p, y[0]);
    if (!admissible() || fabs(fluid.GetPressure() / p - 1) > 1e-6) return false;
    const auto a2 = fluid.GetSoundSpeed2();
    derivative = {p / a2, -p / (y[0] * sqrt(a2))};
    return true;
  };
  auto step = [&](su2double x, const Pair& y, su2double h, Pair& out) {
    Pair k1, k2, k3, k4, trial;
    if (!evaluate(x, y, k1)) return false;
    for (auto k = 0u; k < 2; ++k) trial[k] = y[k] + .5 * h * k1[k];
    if (!evaluate(x + .5 * h, trial, k2)) return false;
    for (auto k = 0u; k < 2; ++k) trial[k] = y[k] + .5 * h * k2[k];
    if (!evaluate(x + .5 * h, trial, k3)) return false;
    for (auto k = 0u; k < 2; ++k) trial[k] = y[k] + h * k3[k];
    if (!evaluate(x + h, trial, k4)) return false;
    for (auto k = 0u; k < 2; ++k) out[k] = y[k] + h * (k1[k] + 2 * k2[k] + 2 * k3[k] + k4[k]) / 6;
    return evaluate(x + h, out, trial);
  };
  su2double x = log(pressure), end = log(targetPressure);
  const auto direction = end >= x ? 1.0 : -1.0;
  su2double h = direction * .05;
  bool finished = false;
  for (auto iteration = 0u; iteration < 10000; ++iteration) {
    if (fabs(end - x) < 1e-13) {
      finished = true;
      break;
    }
    if (fabs(h) > fabs(end - x)) h = end - x;
    Pair coarse, half, fine;
    bool valid = step(x, state, h, coarse) && step(x, state, .5 * h, half) && step(x + .5 * h, half, .5 * h, fine);
    su2double error = 1;
    if (valid)
      error =
          fmax(fabs(fine[0] - coarse[0]) / (rho + fabs(fine[0])), fabs(fine[1] - coarse[1]) / (sound + fabs(fine[1])));
    if (!valid || error > 1e-10) {
      h *= .5;
      if (fabs(h) < 1e-12) return Status::InvalidState;
      continue;
    }
    const auto nextSound = sqrt(fluid.GetSoundSpeed2());
    if (fine[1] >= nextSound) {
      su2double lower = 0, upper = 1;
      for (auto k = 0u; k < 40; ++k) {
        const auto fraction = .5 * (lower + upper);
        if (!step(x, state, fraction * h * .5, half) || !step(x + fraction * h * .5, half, fraction * h * .5, fine))
          return Status::InvalidState;
        if (fine[1] >= sqrt(fluid.GetSoundSpeed2()))
          upper = fraction;
        else
          lower = fraction;
      }
      x += .5 * (lower + upper) * h;
      state = fine;
      finished = true;
      break;
    }
    x += h;
    state = fine;
    if (error < 1e-12) h *= 2;
  }
  if (!finished) return Status::InvalidState;
  if (state[1] < 0) return Status::Backflow;
  fluid.SetTDState_Prho(exp(x), state[0]);
  if (!admissible()) return Status::InvalidState;
  boundary[0] = state[0];
  for (auto d = 0u; d < nDim; ++d) {
    velocity[d] += (state[1] - normalVelocity) * normal[d];
    boundary[d + 1] = state[0] * velocity[d];
  }
  boundary[nDim + 1] =
      state[0] * (fluid.GetStaticEnergy() + .5 * GeometryToolbox::SquaredNorm(nDim, velocity) + turbulentK);
  return Status::Valid;
}
}  // namespace PressureOutlet
