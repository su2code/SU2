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
                    su2double gridNormal, su2double targetPressure, su2double* boundary, su2double turbulentK = 0, su2double idealGamma = 0) {
  for (auto v = 0u; v < nDim + 2; ++v)
    if (!std::isfinite(SU2_TYPE::GetValue(interior[v]))) return Status::InvalidState;
  const auto rho = interior[0];
  if (!(rho > 0) || !(targetPressure > 0) || !std::isfinite(SU2_TYPE::GetValue(targetPressure)))
    return Status::InvalidState;
  su2double velocity[3] = {};
  for (auto d = 0u; d < nDim; ++d) velocity[d] = interior[d + 1] / rho;
  // SST total energy contains frozen turbulent kinetic energy; the EOS needs thermal energy only.
  const su2double energy = interior[nDim + 1] / rho - .5 * GeometryToolbox::SquaredNorm(nDim, velocity) - turbulentK;
  fluid.SetTDState_rhoe(rho, energy);
  const auto pressure = fluid.GetPressure(), sound2 = fluid.GetSoundSpeed2();
  auto admissible = [&]() {
    const auto p = SU2_TYPE::GetValue(fluid.GetPressure()), t = SU2_TYPE::GetValue(fluid.GetTemperature());
    const auto a2 = SU2_TYPE::GetValue(fluid.GetSoundSpeed2());
    return p > 0 && t > 0 && a2 > 0 && std::isfinite(p + t + a2);
  };
  if (!admissible()) return Status::InvalidState;
  const su2double sound = sqrt(sound2);
  const su2double normalVelocity = GeometryToolbox::DotProduct(nDim, velocity, normal) - gridNormal;
  // Subsonic reverse flow retains one outgoing acoustic wave. Extrapolate entropy and tangential
  // velocity from the interior, as in the existing pressure-outlet closure; fully incoming flow needs inlet data.
  if (normalVelocity <= -sound) return Status::Backflow;
  if (normalVelocity >= sound) {
    for (auto v = 0u; v < nDim + 2; ++v) boundary[v] = interior[v];
    return Status::Valid;
  }

  // The closed-form ideal-gas invariant avoids numerical integration in ordinary outlets.
  if (idealGamma > 1) {
    const su2double gm1 = idealGamma - 1;
    const su2double sonicSound = (gm1 * normalVelocity + 2 * sound) / (idealGamma + 1);
    const su2double targetSound = sound * pow(targetPressure / pressure, gm1 / (2 * idealGamma));
    const su2double boundarySound = fmax(targetSound, sonicSound);
    const su2double density = rho * pow(boundarySound / sound, 2 / gm1);
    const su2double q = normalVelocity + 2 * (sound - boundarySound) / gm1;
    if (q <= -boundarySound) return Status::Backflow;
    const su2double p = pressure * pow(density / rho, idealGamma);
    fluid.SetTDState_Prho(p, density);
    if (!admissible()) return Status::InvalidState;
    boundary[0] = density;
    for (auto d = 0u; d < nDim; ++d) {
      velocity[d] += (q - normalVelocity) * normal[d];
      boundary[d + 1] = density * velocity[d];
    }
    boundary[nDim + 1] =
        density * (fluid.GetStaticEnergy() + .5 * GeometryToolbox::SquaredNorm(nDim, velocity) + turbulentK);
    return Status::Valid;
  }

  // Along an isentrope: d(rho)/dp=1/a^2, de/dp=p/(rho^2*a^2), dq/dp=-1/(rho*a).
  // Advancing thermal energy avoids relying on an iterative pressure/density inversion at every stage.
  // Log pressure makes the step size independent of units.
  using StateVector = std::array<su2double, 3>;  // Density, thermal energy, relative normal velocity.
  StateVector state = {rho, energy, normalVelocity};
  auto evaluate = [&](su2double logPressure, const StateVector& y, StateVector& derivative) {
    if (!(y[0] > 0) || !std::isfinite(SU2_TYPE::GetValue(y[0] + y[1] + y[2]))) return false;
    const su2double p = exp(logPressure);
    fluid.SetTDState_rhoe(y[0], y[1]);
    // RK stage predictors need not lie exactly on the pressure curve. Reject EOS clipping instead.
    if (!admissible() || fabs(fluid.GetDensity() / y[0] - 1) > 1e-6 ||
        fabs(fluid.GetStaticEnergy() - y[1]) > 1e-6 * (sound2 + fabs(y[1]))) return false;
    const auto a2 = fluid.GetSoundSpeed2();
    derivative = {p / a2, p * p / (y[0] * y[0] * a2), -p / (y[0] * sqrt(a2))};
    return true;
  };
  auto step = [&](su2double x, const StateVector& y, su2double h, StateVector& out) {
    StateVector k1, k2, k3, k4, trial;
    if (!evaluate(x, y, k1)) return false;
    for (auto k = 0u; k < 3; ++k) trial[k] = y[k] + .5 * h * k1[k];
    if (!evaluate(x + .5 * h, trial, k2)) return false;
    for (auto k = 0u; k < 3; ++k) trial[k] = y[k] + .5 * h * k2[k];
    if (!evaluate(x + .5 * h, trial, k3)) return false;
    for (auto k = 0u; k < 3; ++k) trial[k] = y[k] + h * k3[k];
    if (!evaluate(x + h, trial, k4)) return false;
    for (auto k = 0u; k < 3; ++k) out[k] = y[k] + h * (k1[k] + 2 * k2[k] + 2 * k3[k] + k4[k]) / 6;
    return evaluate(x + h, out, trial);
  };
  su2double x = log(pressure), end = log(targetPressure);
  const auto direction = end >= x ? 1.0 : -1.0;
  su2double h = direction * .05;
  bool finished = false;
  for (auto iteration = 0u; iteration < 10000; ++iteration) {
    if (fabs(end - x) < 1e-13) {
      // A negligible pressure interval can still carry an AD sensitivity.
      StateVector derivative;
      if (!evaluate(x, state, derivative)) return Status::InvalidState;
      for (auto k = 0u; k < 3; ++k) state[k] += (end - x) * derivative[k];
      x = end;
      finished = true;
      break;
    }
    if (fabs(h) > fabs(end - x)) h = end - x;
    // Step doubling controls the nonlinear integration error without assuming a constant gamma.
    StateVector coarse, half, fine;
    bool valid = step(x, state, h, coarse) && step(x, state, .5 * h, half) && step(x + .5 * h, half, .5 * h, fine);
    su2double error = 1;
    if (valid)
      error =
          fmax(fabs(fine[0] - coarse[0]) / (rho + fabs(fine[0])), fmax(fabs(fine[1] - coarse[1]) / (fabs(energy) + sound2 + fabs(fine[1])),
               fabs(fine[2] - coarse[2]) / (sound + fabs(fine[2]))));
    if (!valid || error > 1e-10) {
      h *= .5;
      if (fabs(h) < 1e-12) return Status::InvalidState;
      continue;
    }
    const su2double nextSound = sqrt(fluid.GetSoundSpeed2());
    // Stop expansion where q-a changes sign: prescribed pressure can no longer enter the domain.
    if (fine[2] >= nextSound) {
      su2double lower = 0, upper = 1;
      for (auto k = 0u; k < 40; ++k) {
        const su2double fraction = .5 * (lower + upper);
        if (!step(x, state, fraction * h * .5, half) || !step(x + fraction * h * .5, half, fraction * h * .5, fine))
          return Status::InvalidState;
        if (fine[2] >= sqrt(fluid.GetSoundSpeed2()))
          upper = fraction;
        else
          lower = fraction;
      }
      x += .5 * (lower + upper) * h;
      state = fine;
      // Bisection locates the sonic state but does not differentiate its moving root.
      StateVector derivative;
      if (!evaluate(x, state, derivative)) return Status::InvalidState;
      const su2double sonicResidual = state[2] - sqrt(fluid.GetSoundSpeed2());
      constexpr double delta = 1e-5;
      StateVector shifted = state, unused;
      for (auto k = 0u; k < 3; ++k) shifted[k] += delta * derivative[k];
      if (!evaluate(x + delta, shifted, unused)) return Status::InvalidState;
      const double soundPlus = SU2_TYPE::GetValue(fluid.GetSoundSpeed());
      shifted = state;
      for (auto k = 0u; k < 3; ++k) shifted[k] -= delta * derivative[k];
      if (!evaluate(x - delta, shifted, unused)) return Status::InvalidState;
      const double soundMinus = SU2_TYPE::GetValue(fluid.GetSoundSpeed());
      const double slope = SU2_TYPE::GetValue(derivative[2]) - (soundPlus - soundMinus) / (2 * delta);
      if (!std::isfinite(slope) || fabs(slope) < 1e-12) return Status::InvalidState;
      const su2double correction = -sonicResidual / slope;
      for (auto k = 0u; k < 3; ++k) state[k] += correction * derivative[k];
      x += correction;
      finished = true;
      break;
    }
    x += h;
    state = fine;
    if (error < 1e-12) h *= 2;
  }
  if (!finished) return Status::InvalidState;
  fluid.SetTDState_rhoe(state[0], state[1]);
  if (!admissible()) return Status::InvalidState;
  if (fabs(fluid.GetPressure() / exp(x) - 1) > 1e-6) return Status::InvalidState;
  if (state[2] <= -fluid.GetSoundSpeed()) return Status::Backflow;
  boundary[0] = state[0];
  for (auto d = 0u; d < nDim; ++d) {
    velocity[d] += (state[2] - normalVelocity) * normal[d];
    boundary[d + 1] = state[0] * velocity[d];
  }
  boundary[nDim + 1] =
      state[0] * (fluid.GetStaticEnergy() + .5 * GeometryToolbox::SquaredNorm(nDim, velocity) + turbulentK);
  return Status::Valid;
}
/*--- Differentiate the state map, rather than freezing the prescribed-pressure ghost state. ---*/
inline Status Jacobian(CFluidModel& fluid, unsigned short nDim, const su2double* interior, const su2double* normal,
                       su2double gridNormal, su2double pressure, const su2double* boundary, su2double** jacobian,
                       su2double turbulentK = 0, su2double idealGamma = 0) {
  const auto nVar = nDim + 2;
  const su2double rho = interior[0];
  su2double velocity[3] = {};
  for (auto d = 0u; d < nDim; ++d) velocity[d] = interior[d + 1] / rho;
  fluid.SetTDState_rhoe(rho, interior[nDim + 1] / rho - .5 * GeometryToolbox::SquaredNorm(nDim, velocity) - turbulentK);
  const su2double sound = fluid.GetSoundSpeed();
  if (GeometryToolbox::DotProduct(nDim, velocity, normal) - gridNormal >= sound) {
    for (auto row = 0u; row < nVar; ++row)
      for (auto column = 0u; column < nVar; ++column) jacobian[row][column] = row == column ? 1 : 0;
    return Status::Valid;
  }
  // Differentiate the complete EOS reconstruction, including the moving sonic limit. The one-sided
  // fallback handles perturbations that leave the admissible EOS/characteristic range.
  for (auto column = 0u; column < nVar; ++column) {
    su2double plus[5], minus[5], bPlus[5], bMinus[5];
    for (auto k = 0u; k < nVar; ++k) plus[k] = minus[k] = interior[k];
    double scale = SU2_TYPE::GetValue(rho);
    if (column > 0) scale *= SU2_TYPE::GetValue(sound);
    if (column == nVar - 1) scale *= SU2_TYPE::GetValue(sound);
    const double delta = 1e-5 * fmax(fabs(SU2_TYPE::GetValue(interior[column])), scale);
    plus[column] += delta;
    minus[column] -= delta;
    const bool validPlus = State(fluid, nDim, plus, normal, gridNormal, pressure, bPlus, turbulentK, idealGamma) == Status::Valid;
    const bool validMinus = State(fluid, nDim, minus, normal, gridNormal, pressure, bMinus, turbulentK, idealGamma) == Status::Valid;
    if (!validPlus && !validMinus) return Status::InvalidState;
    for (auto row = 0u; row < nVar; ++row) {
      if (validPlus && validMinus) jacobian[row][column] = (bPlus[row] - bMinus[row]) / (2 * delta);
      else if (validPlus) jacobian[row][column] = (bPlus[row] - boundary[row]) / delta;
      else jacobian[row][column] = (boundary[row] - bMinus[row]) / delta;
    }
  }
  return Status::Valid;
}
}  // namespace PressureOutlet
