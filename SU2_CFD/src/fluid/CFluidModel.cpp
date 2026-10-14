/*!
 * \file CFluidModel.cpp
 * \brief Source of the fluid model base class containing thermo-physical subroutines.
 * \author S.Vitale, M.Pini, G.Gori, A.Guardone, P.Colonna, T. Economon
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

#include "../../include/fluid/CFluidModel.hpp"

#include <utility>
#include <array>

#include "../../../Common/include/toolboxes/geometry_toolbox.hpp"

#include "../../include/fluid/CConstantConductivity.hpp"
#include "../../include/fluid/CConstantConductivityRANS.hpp"
#include "../../include/fluid/CConstantDiffusivity.hpp"
#include "../../include/fluid/CConstantPrandtl.hpp"
#include "../../include/fluid/CConstantPrandtlRANS.hpp"
#include "../../include/fluid/CConstantSchmidt.hpp"
#include "../../include/fluid/CConstantViscosity.hpp"
#include "../../include/fluid/CFluidScalar.hpp"
#include "../../include/fluid/CPolynomialConductivity.hpp"
#include "../../include/fluid/CPolynomialConductivityRANS.hpp"
#include "../../include/fluid/CPolynomialViscosity.hpp"
#include "../../include/fluid/CSutherland.hpp"
#include "../../include/fluid/CFluidFlamelet.hpp"
#include "../../include/fluid/CCoolPropViscosity.hpp"
#include "../../include/fluid/CConstantLewisDiffusivity.hpp"
#include "../../include/fluid/CCoolPropConductivity.hpp"

unique_ptr<CViscosityModel> CFluidModel::MakeLaminarViscosityModel(const CConfig* config, unsigned short iSpecies) {
  switch (config->GetKind_ViscosityModel()) {
    case VISCOSITYMODEL::CONSTANT:
      return unique_ptr<CConstantViscosity>(new CConstantViscosity(config->GetMu_ConstantND(iSpecies)));
    case VISCOSITYMODEL::COOLPROP:
      return unique_ptr<CCoolPropViscosity>(new CCoolPropViscosity(config->GetFluid_Name()));
    case VISCOSITYMODEL::SUTHERLAND:
      return unique_ptr<CSutherland>(new CSutherland(config->GetMu_RefND(iSpecies),
                                                     config->GetMu_Temperature_RefND(iSpecies),
                                                     config->GetMu_SND(iSpecies)));
    case VISCOSITYMODEL::POLYNOMIAL:
      return unique_ptr<CPolynomialViscosity<N_POLY_COEFFS>>(
          new CPolynomialViscosity<N_POLY_COEFFS>(config->GetMu_PolyCoeffND()));
    case VISCOSITYMODEL::FLAMELET:
      /*--- Viscosity is obtained from the LUT ---*/
      return nullptr;
    default:
      SU2_MPI::Error("Viscosity model not available.", CURRENT_FUNCTION);
      return nullptr;
  }
}

void CFluidModel::SetLaminarViscosityModel(const CConfig* config) {
  LaminarViscosity = MakeLaminarViscosityModel(config, 0);
}

unique_ptr<CConductivityModel> CFluidModel::MakeThermalConductivityModel(const CConfig* config,
                                                                         unsigned short iSpecies) {
  switch (config->GetKind_ConductivityModel()) {
    case CONDUCTIVITYMODEL::CONSTANT:
      if (config->GetKind_ConductivityModel_Turb() == CONDUCTIVITYMODEL_TURB::CONSTANT_PRANDTL) {
        return unique_ptr<CConstantConductivityRANS>(
            new CConstantConductivityRANS(config->GetThermal_Conductivity_ConstantND(iSpecies),
                                          config->GetPrandtl_Turb(iSpecies)));
      } else {
        return unique_ptr<CConstantConductivity>(
            new CConstantConductivity(config->GetThermal_Conductivity_ConstantND(iSpecies)));
      }
      break;
    case CONDUCTIVITYMODEL::COOLPROP:
      if (config->GetKind_ConductivityModel_Turb() == CONDUCTIVITYMODEL_TURB::CONSTANT_PRANDTL) {
        return unique_ptr<CConstantConductivityRANS>(
            new CConstantConductivityRANS(config->GetThermal_Conductivity_ConstantND(iSpecies),
                                          config->GetPrandtl_Turb(iSpecies)));
      } else {
        return unique_ptr<CCoolPropConductivity>(
            new CCoolPropConductivity(config->GetFluid_Name()));
      }
      break;
    case CONDUCTIVITYMODEL::CONSTANT_PRANDTL:
      if (config->GetKind_ConductivityModel_Turb() == CONDUCTIVITYMODEL_TURB::CONSTANT_PRANDTL) {
        return unique_ptr<CConstantPrandtlRANS>(
            new CConstantPrandtlRANS(config->GetPrandtl_Lam(iSpecies), config->GetPrandtl_Turb(iSpecies)));
      } else {
        return unique_ptr<CConstantPrandtl>(new CConstantPrandtl(config->GetPrandtl_Lam(iSpecies)));
      }
      break;
    case CONDUCTIVITYMODEL::POLYNOMIAL:
      if (config->GetKind_ConductivityModel_Turb() == CONDUCTIVITYMODEL_TURB::CONSTANT_PRANDTL) {
        return unique_ptr<CPolynomialConductivityRANS<N_POLY_COEFFS>>(
            new CPolynomialConductivityRANS<N_POLY_COEFFS>(config->GetKt_PolyCoeffND(), config->GetPrandtl_Turb()));
      } else {
        return unique_ptr<CPolynomialConductivity<N_POLY_COEFFS>>(
            new CPolynomialConductivity<N_POLY_COEFFS>(config->GetKt_PolyCoeffND()));
      }
      break;
    case CONDUCTIVITYMODEL::FLAMELET:
      /*--- Conductivity is obtained from the LUT ---*/
      return nullptr;
    default:
      SU2_MPI::Error("Conductivity model not available.", CURRENT_FUNCTION);
      return nullptr;
  }
}

void CFluidModel::SetThermalConductivityModel(const CConfig* config) {
  ThermalConductivity = MakeThermalConductivityModel(config, 0);
}

unique_ptr<CDiffusivityModel> CFluidModel::MakeMassDiffusivityModel(const CConfig* config, unsigned short iSpecies) {
  switch (config->GetKind_Diffusivity_Model()) {
    case DIFFUSIVITYMODEL::CONSTANT_DIFFUSIVITY:
      return unique_ptr<CConstantDiffusivity>(new CConstantDiffusivity(config->GetDiffusivity_ConstantND()));
      break;
    case DIFFUSIVITYMODEL::CONSTANT_SCHMIDT:
      return unique_ptr<CConstantSchmidt>(new CConstantSchmidt(config->GetSchmidt_Number_Laminar()));
      break;
    case DIFFUSIVITYMODEL::UNITY_LEWIS:
      return unique_ptr<CConstantLewisDiffusivity>(new CConstantLewisDiffusivity(1.0));
      break;
    case DIFFUSIVITYMODEL::CONSTANT_LEWIS:
      return unique_ptr<CConstantLewisDiffusivity>(
          new CConstantLewisDiffusivity(config->GetConstant_Lewis_Number(iSpecies)));
      break;
    case DIFFUSIVITYMODEL::FLAMELET:
      /*--- Diffusivity is obtained from the LUT ---*/
      return nullptr;
      break;
    default:
      SU2_MPI::Error("Diffusivity model not available.", CURRENT_FUNCTION);
      return nullptr;
  }
}

void CFluidModel::SetMassDiffusivityModel(const CConfig* config) {
  MassDiffusivity = MakeMassDiffusivityModel(config, 0);
}

/*!
 * \brief Build a pressure-outlet state while preserving interior entropy and tangential velocity.
 *
 * The acoustic wave speeds are q-a and q+a, where q is outward normal velocity relative to the
 * grid and a is sound speed. For 0 < q < a, only q-a enters the domain, so pressure supplies the
 * one missing condition. The outgoing acoustic relation is dq = -dp/(rho*a) at constant entropy.
 * Expansion stops at q=a: lower external pressure can no longer send information upstream.
 *
 * Characteristic boundary theory: K. W. Thompson, "Time dependent boundary conditions for
 * hyperbolic systems", J. Comput. Phys. 68 (1987), 1-24, https://doi.org/10.1016/0021-9991(87)90041-6;
 * T. J. Poinsot and S. K. Lele, "Boundary conditions for direct simulations of compressible viscous
 * flows", J. Comput. Phys. 101 (1992), 104-129, https://doi.org/10.1016/0021-9991(92)90046-2.
 * These references cover characteristic wave counting. The finite-change integration along the
 * equation of state (EOS) isentrope and the sonic-root derivative below are implementation choices.
 *
 * Evaluations overwrite this fluid model's thermodynamic state.
 * \param[in] nDim - Spatial dimension (2 or 3).
 * \param[in] interior - Conservative state [rho, rho*u_1, ..., rho*u_nDim, rho*E].
 * \param[in] normal - Outward unit normal; tangential velocity is unchanged.
 * \param[in] gridNormal - Grid velocity dotted with normal, in the same units as interior velocity.
 * \param[in] targetPressure - Positive static thermodynamic pressure in the interior state's units.
 * \param[out] boundary - Conservative outlet state; use only when the returned status is Valid.
 * \param[in] turbulentK - Specific SST turbulent kinetic energy, held fixed during reconstruction.
 * \param[in] idealGamma - Constant ideal-gas gamma; pass zero to integrate the general EOS.
 * \return Valid, InvalidState (invalid input, EOS state or failed integration), or Backflow.
 */
CFluidModel::PressureOutletStatus CFluidModel::ComputePressureOutletState(
    unsigned short nDim, const su2double* interior, const su2double* normal, su2double gridNormal,
    su2double targetPressure, su2double* boundary, su2double turbulentK, su2double idealGamma) {
  // Reject non-finite inputs and nonpositive density or pressure before using the EOS.
  for (auto v = 0u; v < nDim + 2; ++v)
    if (!std::isfinite(SU2_TYPE::GetValue(interior[v]))) return PressureOutletStatus::InvalidState;
  const auto rho = interior[0];
  if (!(rho > 0) || !(targetPressure > 0) || !std::isfinite(SU2_TYPE::GetValue(targetPressure)))
    return PressureOutletStatus::InvalidState;
  // Recover velocity and thermal energy from the interior conservative state.
  su2double velocity[3] = {};
  for (auto d = 0u; d < nDim; ++d) velocity[d] = interior[d + 1] / rho;
  // SST total energy contains frozen turbulent kinetic energy; the EOS needs thermal energy only.
  const su2double energy = interior[nDim + 1] / rho - .5 * GeometryToolbox::SquaredNorm(nDim, velocity) - turbulentK;
  SetTDState_rhoe(rho, energy);
  const auto pressure = GetPressure(), sound2 = GetSoundSpeed2();
  // Reuse this check whenever an EOS evaluation produces a new thermodynamic state.
  auto admissible = [&]() {
    const auto p = SU2_TYPE::GetValue(GetPressure()), t = SU2_TYPE::GetValue(GetTemperature());
    const auto a2 = SU2_TYPE::GetValue(GetSoundSpeed2());
    return p > 0 && t > 0 && a2 > 0 && std::isfinite(p + t + a2);
  };
  if (!admissible()) return PressureOutletStatus::InvalidState;
  const su2double sound = sqrt(sound2);
  // Wave direction depends on normal speed relative to the boundary, not total flow speed.
  const su2double normalVelocity = GeometryToolbox::DotProduct(nDim, velocity, normal) - gridNormal;
  // Subsonic reverse flow uses interior entropy and tangential velocity as inflow assumptions.
  // Fully incoming flow needs inlet data that a pressure-only boundary cannot provide.
  if (normalVelocity <= -sound) return PressureOutletStatus::Backflow;
  if (normalVelocity >= sound) {
    // No wave enters the domain, so the prescribed pressure cannot change this state.
    for (auto v = 0u; v < nDim + 2; ++v) boundary[v] = interior[v];
    return PressureOutletStatus::Valid;
  }

  // Ideal gases have a closed-form solution, so no numerical integration is needed.
  if (idealGamma > 1) {
    const su2double gm1 = idealGamma - 1;
    // Compute the sound speed at the sonic limit and at the requested pressure.
    const su2double sonicSound = (gm1 * normalVelocity + 2 * sound) / (idealGamma + 1);
    const su2double targetSound = sound * pow(targetPressure / pressure, gm1 / (2 * idealGamma));
    // Limit expansion to sonic speed, then recover density and normal velocity.
    const su2double boundarySound = fmax(targetSound, sonicSound);
    const su2double density = rho * pow(boundarySound / sound, 2 / gm1);
    const su2double q = normalVelocity + 2 * (sound - boundarySound) / gm1;
    if (q <= -boundarySound) return PressureOutletStatus::Backflow;
    // Evaluate the EOS at the actual outlet pressure, which may be limited by the sonic state.
    const su2double p = pressure * pow(density / rho, idealGamma);
    SetTDState_Prho(p, density);
    if (!admissible()) return PressureOutletStatus::InvalidState;
    // Change only normal velocity and convert back to conservative variables.
    boundary[0] = density;
    for (auto d = 0u; d < nDim; ++d) {
      velocity[d] += (q - normalVelocity) * normal[d];
      boundary[d + 1] = density * velocity[d];
    }
    // Rebuild total energy, including the unchanged SST turbulent kinetic energy.
    boundary[nDim + 1] =
        density * (GetStaticEnergy() + .5 * GeometryToolbox::SquaredNorm(nDim, velocity) + turbulentK);
    return PressureOutletStatus::Valid;
  }

  // For a general EOS, follow the constant-entropy path by integrating density, energy and velocity.
  // Use log pressure so each step represents a relative pressure change.
  using StateVector = std::array<su2double, 3>;  // Density, thermal energy, relative normal velocity.
  StateVector state = {rho, energy, normalVelocity};
  // Return how the three state variables change with log pressure at this integration stage.
  auto evaluate = [&](su2double logPressure, const StateVector& y, StateVector& derivative) {
    if (!(y[0] > 0) || !std::isfinite(SU2_TYPE::GetValue(y[0] + y[1] + y[2]))) return false;
    const su2double p = exp(logPressure);
    SetTDState_rhoe(y[0], y[1]);
    // Intermediate stages can have pressure error, but the EOS must not silently change their density or energy.
    if (!admissible() || fabs(GetDensity() / y[0] - 1) > 1e-6 ||
        fabs(GetStaticEnergy() - y[1]) > 1e-6 * (sound2 + fabs(y[1]))) return false;
    const auto a2 = GetSoundSpeed2();
    // Constant-entropy derivatives: d(rho)/dp=1/a^2, de/dp=p/(rho^2*a^2), dq/dp=-1/(rho*a).
    derivative = {p / a2, p * p / (y[0] * y[0] * a2), -p / (y[0] * sqrt(a2))};
    return true;
  };
  // Classical fourth-order Runge-Kutta; each stage must remain inside the EOS validity range.
  auto step = [&](su2double x, const StateVector& y, su2double h, StateVector& out) {
    StateVector k1, k2, k3, k4, trial;
    // Sample the slope at the start and at two estimates of the midpoint.
    if (!evaluate(x, y, k1)) return false;
    for (auto k = 0u; k < 3; ++k) trial[k] = y[k] + .5 * h * k1[k];
    if (!evaluate(x + .5 * h, trial, k2)) return false;
    for (auto k = 0u; k < 3; ++k) trial[k] = y[k] + .5 * h * k2[k];
    if (!evaluate(x + .5 * h, trial, k3)) return false;
    // Sample the endpoint, then combine all four slopes to advance the state.
    for (auto k = 0u; k < 3; ++k) trial[k] = y[k] + h * k3[k];
    if (!evaluate(x + h, trial, k4)) return false;
    for (auto k = 0u; k < 3; ++k) out[k] = y[k] + h * (k1[k] + 2 * k2[k] + 2 * k3[k] + k4[k]) / 6;
    // Validate the combined state and leave the EOS at that endpoint.
    return evaluate(x + h, out, trial);
  };
  // Start at the interior pressure and choose the direction toward the requested pressure.
  su2double x = log(pressure), end = log(targetPressure);
  const auto direction = end >= x ? 1.0 : -1.0;
  su2double h = direction * .05;
  bool finished = false;
  for (auto iteration = 0u; iteration < 10000; ++iteration) {
    if (fabs(end - x) < 1e-13) {
      // Finish the tiny remaining pressure change while retaining its sensitivity to the inputs.
      StateVector derivative;
      if (!evaluate(x, state, derivative)) return PressureOutletStatus::InvalidState;
      for (auto k = 0u; k < 3; ++k) state[k] += (end - x) * derivative[k];
      x = end;
      finished = true;
      break;
    }
    // Shorten the last step so it cannot pass the requested pressure.
    if (fabs(h) > fabs(end - x)) h = end - x;
    // Compare one full step with two half steps and retain the finer solution when accepted.
    StateVector coarse, half, fine;
    bool valid = step(x, state, h, coarse) && step(x, state, .5 * h, half) && step(x + .5 * h, half, .5 * h, fine);
    su2double error = 1;
    // Scale each error by a typical density, energy or speed so their sizes can be compared.
    if (valid)
      error =
          fmax(fabs(fine[0] - coarse[0]) / (rho + fabs(fine[0])), fmax(fabs(fine[1] - coarse[1]) / (fabs(energy) + sound2 + fabs(fine[1])),
               fabs(fine[2] - coarse[2]) / (sound + fabs(fine[2]))));
    // Retry with a smaller step if the EOS rejects a stage or the integration error is too large.
    if (!valid || error > 1e-10) {
      h *= .5;
      if (fabs(h) < 1e-12) return PressureOutletStatus::InvalidState;
      continue;
    }
    const su2double nextSound = sqrt(GetSoundSpeed2());
    // If this step reaches sonic speed, stop expansion at the crossing instead of at the target pressure.
    if (fine[2] >= nextSound) {
      su2double lower = 0, upper = 1;
      // Bisect the fraction of this step, reintegrating from the last accepted state each time.
      for (auto k = 0u; k < 40; ++k) {
        const su2double fraction = .5 * (lower + upper);
        if (!step(x, state, fraction * h * .5, half) || !step(x + fraction * h * .5, half, fraction * h * .5, fine))
          return PressureOutletStatus::InvalidState;
        if (fine[2] >= sqrt(GetSoundSpeed2()))
          upper = fraction;
        else
          lower = fraction;
      }
      // Keep the state and pressure at the sonic crossing.
      x += .5 * (lower + upper) * h;
      state = fine;
      // Bisection finds the sonic point; a Newton correction recovers how it moves with the inputs.
      StateVector derivative;
      if (!evaluate(x, state, derivative)) return PressureOutletStatus::InvalidState;
      const su2double sonicResidual = state[2] - sqrt(GetSoundSpeed2());
      // Sample sound speed on either side to estimate how q-a changes along the pressure path.
      constexpr double delta = 1e-5;
      StateVector shifted = state, unused;
      for (auto k = 0u; k < 3; ++k) shifted[k] += delta * derivative[k];
      if (!evaluate(x + delta, shifted, unused)) return PressureOutletStatus::InvalidState;
      const double soundPlus = SU2_TYPE::GetValue(GetSoundSpeed());
      shifted = state;
      for (auto k = 0u; k < 3; ++k) shifted[k] -= delta * derivative[k];
      if (!evaluate(x - delta, shifted, unused)) return PressureOutletStatus::InvalidState;
      const double soundMinus = SU2_TYPE::GetValue(GetSoundSpeed());
      // Only the residual needs automatic differentiation here; the slope's derivative has negligible weight.
      const double slope = SU2_TYPE::GetValue(derivative[2]) - (soundPlus - soundMinus) / (2 * delta);
      if (!std::isfinite(slope) || fabs(slope) < 1e-12) return PressureOutletStatus::InvalidState;
      // Apply the correction to both the state and its stopping pressure.
      const su2double correction = -sonicResidual / slope;
      for (auto k = 0u; k < 3; ++k) state[k] += correction * derivative[k];
      x += correction;
      finished = true;
      break;
    }
    // Accept the two-half-step result; grow the next step if the error was comfortably small.
    x += h;
    state = fine;
    if (error < 1e-12) h *= 2;
  }
  // Reject an unfinished integration and check the final EOS state before building the output.
  if (!finished) return PressureOutletStatus::InvalidState;
  SetTDState_rhoe(state[0], state[1]);
  if (!admissible()) return PressureOutletStatus::InvalidState;
  // Check the accumulated pressure error at the target or sonic stop, not at RK predictors.
  if (fabs(GetPressure() / exp(x) - 1) > 1e-6) return PressureOutletStatus::InvalidState;
  if (state[2] <= -GetSoundSpeed()) return PressureOutletStatus::Backflow;
  // Convert the integrated state back to lab-frame momentum, keeping tangential velocity unchanged.
  boundary[0] = state[0];
  for (auto d = 0u; d < nDim; ++d) {
    velocity[d] += (state[2] - normalVelocity) * normal[d];
    boundary[d + 1] = state[0] * velocity[d];
  }
  // Restore kinetic and SST energy to form the total conservative energy.
  boundary[nDim + 1] =
      state[0] * (GetStaticEnergy() + .5 * GeometryToolbox::SquaredNorm(nDim, velocity) + turbulentK);
  return PressureOutletStatus::Valid;
}
