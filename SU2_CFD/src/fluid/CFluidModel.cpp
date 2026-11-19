/*!
 * \file CFluidModel.cpp
 * \brief Source of the fluid model base class containing thermo-physical subroutines.
 * \author S.Vitale, M.Pini, G.Gori, A.Guardone, P.Colonna, T. Economon
 * \version 8.6.0 "Harrier"
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

#include "../../../Common/include/toolboxes/geometry_toolbox.hpp"

#include "../../include/fluid/CConstantConductivity.hpp"
#include "../../include/fluid/CConstantConductivityRANS.hpp"
#include "../../include/fluid/CConstantDiffusivity.hpp"
#include "../../include/fluid/CConstantPrandtl.hpp"
#include "../../include/fluid/CConstantPrandtlRANS.hpp"
#include "../../include/fluid/CConstantSchmidt.hpp"
#include "../../include/fluid/CConstantViscosity.hpp"
#include "../../include/fluid/CFluidScalar.hpp"
#include "../../include/fluid/CFluidCantera.hpp"
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
    case VISCOSITYMODEL::CANTERA:
      /*--- Viscosity is obtained from the Cantera library---*/
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
    case CONDUCTIVITYMODEL::CANTERA:
      /*--- Conductivity is obtained from Cantera library---*/
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
    case DIFFUSIVITYMODEL::CANTERA:
      /*--- Diffusivity is obtained from CANTERA. Transport model option determines the diffusivity model---*/
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
CFluidModel::PressureOutletStatus CFluidModel::ComputePressureOutletState(unsigned short nDim,
                                                                          const su2double* interior,
                                                                          const su2double* normal, su2double gridNormal,
                                                                          su2double targetPressure, su2double* boundary,
                                                                          su2double turbulentK, su2double idealGamma) {
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
  const auto sound2 = GetSoundSpeed2();
  if (!IsPressureOutletStateValid()) return PressureOutletStatus::InvalidState;
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

  PressureOutletState state = {rho, energy, normalVelocity};
  const auto status = idealGamma > 1 ? ComputeIdealGasPressureOutletState(targetPressure, idealGamma, state)
                                     : IntegratePressureOutletIsentrope(targetPressure, state);
  if (status != PressureOutletStatus::Valid) return status;
  // Convert the outlet state back to lab-frame momentum, keeping tangential velocity unchanged.
  boundary[0] = state[0];
  for (auto d = 0u; d < nDim; ++d) {
    velocity[d] += (state[2] - normalVelocity) * normal[d];
    boundary[d + 1] = state[0] * velocity[d];
  }
  // Restore kinetic and SST energy to form the total conservative energy.
  boundary[nDim + 1] = state[0] * (GetStaticEnergy() + .5 * GeometryToolbox::SquaredNorm(nDim, velocity) + turbulentK);
  return PressureOutletStatus::Valid;
}

// Ideal gases have a closed-form solution, so no numerical integration is needed.
CFluidModel::PressureOutletStatus CFluidModel::ComputeIdealGasPressureOutletState(su2double targetPressure,
                                                                                  su2double idealGamma,
                                                                                  PressureOutletState& state) {
  const auto rho = state[0], normalVelocity = state[2];
  const auto pressure = GetPressure();
  const su2double sound = sqrt(GetSoundSpeed2());
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
  if (!IsPressureOutletStateValid()) return PressureOutletStatus::InvalidState;
  state = {density, GetStaticEnergy(), q};
  return PressureOutletStatus::Valid;
}

// Use log pressure so each integration step represents a relative pressure change.
CFluidModel::PressureOutletStatus CFluidModel::IntegratePressureOutletIsentrope(su2double targetPressure,
                                                                                PressureOutletState& state,
                                                                                bool sonicStop) {
  const auto rho = state[0], energy = state[1];
  const auto pressure = GetPressure(), sound2 = GetSoundSpeed2();
  const su2double sound = sqrt(sound2);
  // Start at the interior pressure and choose the direction toward the requested pressure.
  su2double x = log(pressure), end = log(targetPressure);
  const auto direction = end >= x ? 1.0 : -1.0;
  su2double h = direction * .05;
  bool finished = false;
  for (auto iteration = 0u; iteration < 10000; ++iteration) {
    if (fabs(end - x) < 1e-13) {
      // Finish the tiny remaining pressure change while retaining its sensitivity to the inputs.
      PressureOutletState derivative;
      if (!EvaluatePressureOutletIsentrope(x, state, derivative, sound2)) return PressureOutletStatus::InvalidState;
      for (auto k = 0u; k < 3; ++k) state[k] += (end - x) * derivative[k];
      x = end;
      finished = true;
      break;
    }
    // Shorten the last step so it cannot pass the requested pressure.
    if (fabs(h) > fabs(end - x)) h = end - x;
    // Compare one full step with two half steps and retain the finer solution when accepted.
    PressureOutletState coarse, half, fine;
    bool valid = StepPressureOutletIsentrope(x, state, h, coarse, sound2) &&
                 StepPressureOutletIsentrope(x, state, .5 * h, half, sound2) &&
                 StepPressureOutletIsentrope(x + .5 * h, half, .5 * h, fine, sound2);
    su2double error = 1;
    // Scale each error by a typical density, energy or speed so their sizes can be compared.
    if (valid)
      error = fmax(fabs(fine[0] - coarse[0]) / (rho + fabs(fine[0])),
                   fmax(fabs(fine[1] - coarse[1]) / (fabs(energy) + sound2 + fabs(fine[1])),
                        fabs(fine[2] - coarse[2]) / (sound + fabs(fine[2]))));
    // Retry with a smaller step if the EOS rejects a stage or the integration error is too large.
    if (!valid || error > 1e-10) {
      h *= .5;
      if (fabs(h) < 1e-12) return PressureOutletStatus::InvalidState;
      continue;
    }
    const su2double nextSound = sqrt(GetSoundSpeed2());
    // If this step reaches sonic speed, stop expansion at the crossing instead of at the target pressure.
    if (sonicStop && fine[2] >= nextSound) {
      if (!FindPressureOutletSonicState(x, state, h, sound2)) return PressureOutletStatus::InvalidState;
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
  if (!IsPressureOutletStateValid()) return PressureOutletStatus::InvalidState;
  // Check the accumulated pressure error at the target or sonic stop, not at RK predictors.
  if (fabs(GetPressure() / exp(x) - 1) > 1e-6) return PressureOutletStatus::InvalidState;
  if (sonicStop && state[2] <= -GetSoundSpeed()) return PressureOutletStatus::Backflow;
  return PressureOutletStatus::Valid;
}

bool CFluidModel::FindPressureOutletSonicState(su2double& x, PressureOutletState& state, su2double h,
                                               su2double referenceSound2) {
  PressureOutletState half, fine;
  su2double lower = 0, upper = 1;
  // Bisect the fraction of this step, reintegrating from the last accepted state each time.
  for (auto k = 0u; k < 40; ++k) {
    const su2double fraction = .5 * (lower + upper);
    if (!StepPressureOutletIsentrope(x, state, fraction * h * .5, half, referenceSound2) ||
        !StepPressureOutletIsentrope(x + fraction * h * .5, half, fraction * h * .5, fine, referenceSound2))
      return false;
    if (fine[2] >= sqrt(GetSoundSpeed2()))
      upper = fraction;
    else
      lower = fraction;
  }
  // Keep the state and pressure at the sonic crossing.
  x += .5 * (lower + upper) * h;
  state = fine;
  // Bisection finds the sonic point; a Newton correction recovers how it moves with the inputs.
  PressureOutletState derivative;
  if (!EvaluatePressureOutletIsentrope(x, state, derivative, referenceSound2)) return false;
  const su2double sonicResidual = state[2] - sqrt(GetSoundSpeed2());
  // Sample sound speed on either side to estimate how q-a changes along the pressure path.
  constexpr double delta = 1e-5;
  PressureOutletState shifted = state, unused;
  for (auto k = 0u; k < 3; ++k) shifted[k] += delta * derivative[k];
  if (!EvaluatePressureOutletIsentrope(x + delta, shifted, unused, referenceSound2)) return false;
  const double soundPlus = SU2_TYPE::GetValue(GetSoundSpeed());
  shifted = state;
  for (auto k = 0u; k < 3; ++k) shifted[k] -= delta * derivative[k];
  if (!EvaluatePressureOutletIsentrope(x - delta, shifted, unused, referenceSound2)) return false;
  const double soundMinus = SU2_TYPE::GetValue(GetSoundSpeed());
  // Only the residual needs automatic differentiation here; the slope's derivative has negligible weight.
  const double slope = SU2_TYPE::GetValue(derivative[2]) - (soundPlus - soundMinus) / (2 * delta);
  if (!std::isfinite(slope) || fabs(slope) < 1e-12) return false;
  // Apply the correction to both the state and its stopping pressure.
  const su2double correction = -sonicResidual / slope;
  for (auto k = 0u; k < 3; ++k) state[k] += correction * derivative[k];
  x += correction;
  return true;
}

// Classical fourth-order Runge-Kutta; each stage must remain inside the EOS validity range.
bool CFluidModel::StepPressureOutletIsentrope(su2double x, const PressureOutletState& y, su2double h,
                                              PressureOutletState& out, su2double referenceSound2) {
  PressureOutletState k1, k2, k3, k4, trial;
  // Sample the slope at the start and at two estimates of the midpoint.
  if (!EvaluatePressureOutletIsentrope(x, y, k1, referenceSound2)) return false;
  for (auto k = 0u; k < 3; ++k) trial[k] = y[k] + .5 * h * k1[k];
  if (!EvaluatePressureOutletIsentrope(x + .5 * h, trial, k2, referenceSound2)) return false;
  for (auto k = 0u; k < 3; ++k) trial[k] = y[k] + .5 * h * k2[k];
  if (!EvaluatePressureOutletIsentrope(x + .5 * h, trial, k3, referenceSound2)) return false;
  // Sample the endpoint, then combine all four slopes to advance the state.
  for (auto k = 0u; k < 3; ++k) trial[k] = y[k] + h * k3[k];
  if (!EvaluatePressureOutletIsentrope(x + h, trial, k4, referenceSound2)) return false;
  for (auto k = 0u; k < 3; ++k) out[k] = y[k] + h * (k1[k] + 2 * k2[k] + 2 * k3[k] + k4[k]) / 6;
  // Validate the combined state and leave the EOS at that endpoint.
  return EvaluatePressureOutletIsentrope(x + h, out, trial, referenceSound2);
}

// Return how density, thermal energy and relative normal velocity change with log pressure.
bool CFluidModel::EvaluatePressureOutletIsentrope(su2double logPressure, const PressureOutletState& y,
                                                  PressureOutletState& derivative, su2double referenceSound2) {
  if (!(y[0] > 0) || !std::isfinite(SU2_TYPE::GetValue(y[0] + y[1] + y[2]))) return false;
  const su2double p = exp(logPressure);
  SetTDState_rhoe(y[0], y[1]);
  // Intermediate stages can have pressure error, but the EOS must not silently change their density or energy.
  if (!IsPressureOutletStateValid() || fabs(GetDensity() / y[0] - 1) > 1e-6 ||
      fabs(GetStaticEnergy() - y[1]) > 1e-6 * (referenceSound2 + fabs(y[1])))
    return false;
  const auto a2 = GetSoundSpeed2();
  // Constant-entropy derivatives: d(rho)/dp=1/a^2, de/dp=p/(rho^2*a^2), dq/dp=-1/(rho*a).
  derivative = {p / a2, p * p / (y[0] * y[0] * a2), -p / (y[0] * sqrt(a2))};
  return true;
}

bool CFluidModel::IsPressureOutletStateValid() const {
  const auto p = SU2_TYPE::GetValue(GetPressure()), t = SU2_TYPE::GetValue(GetTemperature());
  const auto a2 = SU2_TYPE::GetValue(GetSoundSpeed2());
  return p > 0 && t > 0 && a2 > 0 && std::isfinite(p + t + a2);
}

/*!
 * \brief Build the inflow state of a static-pressure outlet whose backflow comes from a reservoir at rest.
 *
 * For subsonic inflow (-a < q < 0, q relative to the grid) only the q+a wave leaves the domain. The boundary
 * state therefore takes the reservoir total enthalpy and entropy and the given flow direction, and meets the
 * interior outgoing acoustic relation dq = -dp/(rho*a), integrated along the interior isentrope. This is the
 * total-condition inlet closure, but with the interior invariant evaluated on the interior entropy as in
 * ComputePressureOutletState. Flow from a reservoir cannot exceed the sonic speed at the boundary, so a state that
 * would be supersonic, or interior flow that already enters through every characteristic, takes the sonic state
 * of the reservoir. If the interior invariant cannot drive flow in even at the reservoir total pressure, the
 * boundary takes the reservoir state at rest.
 *
 * The pressure is found by a safeguarded Newton search on passive values. One Newton correction with active
 * values then gives the derivative of the root with respect to the inputs (implicit function theorem).
 */
CFluidModel::PressureOutletStatus CFluidModel::ComputeReservoirInflowState(
    unsigned short nDim, const su2double* interior, const su2double* normal, su2double gridNormal,
    su2double totalPressure, su2double totalTemperature, const su2double* flowDirection, su2double* boundary,
    su2double turbulentK, su2double idealGamma, const su2double* tangentialVelocity, su2double restPressure) {
  for (auto v = 0u; v < nDim + 2; ++v)
    if (!std::isfinite(SU2_TYPE::GetValue(interior[v]))) return PressureOutletStatus::InvalidState;
  const auto rho = interior[0];
  if (!(rho > 0) || !(totalPressure > 0 || restPressure > 0) || !(totalTemperature > 0) ||
      !std::isfinite(SU2_TYPE::GetValue(totalPressure + totalTemperature + restPressure)))
    return PressureOutletStatus::InvalidState;

  // Tangential velocity carried into the inflow (its normal component is removed); the inflow is then normal.
  su2double tangential[3] = {};
  if (tangentialVelocity != nullptr) {
    for (auto d = 0u; d < nDim; ++d) tangential[d] = tangentialVelocity[d];
    const su2double normalPart = GeometryToolbox::DotProduct(nDim, tangential, normal);
    for (auto d = 0u; d < nDim; ++d) tangential[d] -= normalPart * normal[d];
    if (flowDirection != nullptr && GeometryToolbox::Norm(nDim, flowDirection) > 0)
      return PressureOutletStatus::InvalidState;
  }
  // Guard the square root at zero so that automatic differentiation does not produce 0 * infinity.
  const su2double tangentialSpeed2 = GeometryToolbox::SquaredNorm(nDim, tangential);
  const su2double tangentialSpeed = tangentialSpeed2 > 0 ? su2double(sqrt(tangentialSpeed2)) : su2double(0.0);

  // A zero direction means inflow along the inward normal.
  su2double direction[3] = {};
  if (flowDirection != nullptr)
    for (auto d = 0u; d < nDim; ++d) direction[d] = flowDirection[d];
  su2double length = GeometryToolbox::Norm(nDim, direction);
  if (!(length > 0)) {
    for (auto d = 0u; d < nDim; ++d) direction[d] = -normal[d];
    length = 1;
  }
  for (auto d = 0u; d < nDim; ++d) direction[d] /= length;

  ReservoirInflow data;
  data.idealGamma = idealGamma > 1 ? idealGamma : su2double(0.0);
  data.gridNormal = gridNormal;
  data.directionNormal = GeometryToolbox::DotProduct(nDim, direction, normal);
  // A direction tangent to the boundary or leaving the domain cannot carry inflow.
  if (!(data.directionNormal < 0)) return PressureOutletStatus::InvalidState;

  su2double velocity[3] = {};
  for (auto d = 0u; d < nDim; ++d) velocity[d] = interior[d + 1] / rho;
  data.density = rho;
  data.energy = interior[nDim + 1] / rho - .5 * GeometryToolbox::SquaredNorm(nDim, velocity) - turbulentK;
  SetTDState_rhoe(data.density, data.energy);
  if (!IsPressureOutletStateValid()) return PressureOutletStatus::InvalidState;
  data.pressure = GetPressure();
  data.sound = sqrt(GetSoundSpeed2());
  data.normalVelocity = GeometryToolbox::DotProduct(nDim, velocity, normal) - gridNormal;

  // Total state of the reservoir. A zero total pressure means that the reservoir has the rest pressure where the
  // inflow has no normal velocity, i.e. at the tangential speed (a plenum at rest when there is no tangential
  // velocity).
  if (totalPressure > 0 ? !SetReservoirTotalState(totalPressure, totalTemperature, data)
                        : !FindReservoirTotalPressure(restPressure, totalTemperature, tangentialSpeed, data))
    return PressureOutletStatus::InvalidState;
  if (tangentialSpeed > 0) {
    // The normal motion sees the reservoir expanded to the tangential speed: same entropy, and the total enthalpy
    // reduced by the tangential kinetic energy, which is added back to the boundary state at the end.
    su2double pressure, density, energy, speed, sound;
    if (!FindReservoirSpeedPressure(data, tangentialSpeed, pressure) ||
        !EvaluateReservoirState(data, pressure, density, energy, speed, sound))
      return PressureOutletStatus::InvalidState;
    data.totalPressure = pressure;
    data.totalDensity = density;
    data.totalEnergy = energy;
    data.totalEnthalpy = energy + pressure / density;
  }

  // Search on passive copies so that only the final correction is recorded by automatic differentiation.
  const auto passive = PassiveCopy(data);

  su2double sonicPressure;
  if (!FindReservoirSonicPressure(passive, sonicPressure)) return PressureOutletStatus::InvalidState;

  enum class Kind { Root, Sonic, Rest } kind = Kind::Root;
  su2double lower = log(sonicPressure), upper = log(passive.totalPressure), x = 0;
  su2double residualLower = 0, residualUpper = 0, speed, slope, density, energy, sound;
  if (passive.normalVelocity <= -passive.sound) {
    // The outgoing invariant also enters the domain: no interior information remains.
    kind = Kind::Sonic;
  } else {
    if (!EvaluateInteriorInvariantSpeed(passive, passive.totalPressure, residualUpper, slope))
      return PressureOutletStatus::InvalidState;
    if (residualUpper <= 0) {
      kind = Kind::Rest;
    } else {
      if (!EvaluateInteriorInvariantSpeed(passive, sonicPressure, residualLower, slope) ||
          !EvaluateReservoirState(passive, sonicPressure, density, energy, speed, sound))
        return PressureOutletStatus::InvalidState;
      residualLower -= speed;
      if (residualLower >= 0) kind = Kind::Sonic;
    }
  }

  if (kind == Kind::Root) {
    // The interior invariant speed grows and the reservoir speed falls with pressure, so the residual is monotone.
    // Newton steps in log pressure are safeguarded by bisection (rtsafe): bisect when a step leaves the bracket or
    // does not halve the previous one, so the bracket shrinks even if the residual is not perfectly smooth.
    x = lower - residualLower * (upper - lower) / (residualUpper - residualLower);
    su2double step = upper - lower, previousStep = step;
    bool converged = false;
    for (auto iteration = 0u; iteration < 200 && !converged; ++iteration) {
      su2double invariantSpeed;
      const su2double p = exp(x);
      if (!EvaluateInteriorInvariantSpeed(passive, p, invariantSpeed, slope) ||
          !EvaluateReservoirState(passive, p, density, energy, speed, sound))
        return PressureOutletStatus::InvalidState;
      const su2double residual = invariantSpeed - speed;
      if (fabs(residual) < 1e-12 * passive.sound) {
        converged = true;
        break;
      }
      if (residual < 0)
        lower = x;
      else
        upper = x;
      const su2double derivative = speed > 0 ? su2double(slope + p / (density * speed)) : su2double(0.0);
      const bool bisect = !(derivative > 0) ||
                          ((x - upper) * derivative - residual) * ((x - lower) * derivative - residual) > 0 ||
                          fabs(2 * residual) > fabs(previousStep * derivative);
      previousStep = step;
      step = bisect ? su2double(x - .5 * (lower + upper)) : su2double(residual / derivative);
      x -= step;
      converged = fabs(step) < 1e-13 || upper - lower < 1e-14;
    }
    if (!converged) return PressureOutletStatus::InvalidState;

    // Repeat the last evaluation with active inputs and apply one Newton correction to carry the derivatives.
    su2double invariantSpeed;
    const su2double p = SU2_TYPE::GetValue(exp(x));
    if (!EvaluateInteriorInvariantSpeed(data, p, invariantSpeed, slope) ||
        !EvaluateReservoirState(data, p, density, energy, speed, sound) || !(speed > 0))
      return PressureOutletStatus::InvalidState;
    const double derivative = SU2_TYPE::GetValue(slope + p / (density * speed));
    const su2double pressure = p * exp(-(invariantSpeed - speed) / derivative);
    if (!EvaluateReservoirState(data, pressure, density, energy, speed, sound))
      return PressureOutletStatus::InvalidState;
  } else if (kind == Kind::Sonic) {
    // The sonic state depends on the reservoir only: recompute its pressure with active reservoir data.
    if (data.idealGamma > 1) {
      FindReservoirSonicPressure(data, sonicPressure);
    } else {
      // One Newton step on V^2 - a^2 = 0 in log pressure, with a passive central-difference slope.
      su2double excess[2];
      for (auto side = 0u; side < 2; ++side) {
        if (!EvaluateReservoirState(passive, sonicPressure * exp(side ? 1e-6 : -1e-6), density, energy, speed, sound))
          return PressureOutletStatus::InvalidState;
        excess[side] = speed * speed - sound * sound;
      }
      const double slope = SU2_TYPE::GetValue((excess[1] - excess[0]) / 2e-6);
      if (!EvaluateReservoirState(data, sonicPressure, density, energy, speed, sound) || !(slope < 0))
        return PressureOutletStatus::InvalidState;
      sonicPressure = sonicPressure * exp(-(speed * speed - sound * sound) / slope);
    }
    if (!EvaluateReservoirState(data, sonicPressure, density, energy, speed, sound))
      return PressureOutletStatus::InvalidState;
  } else {
    density = data.totalDensity;
    energy = data.totalEnergy;
    speed = 0;
  }

  SetTDState_rhoe(density, energy);
  if (!IsPressureOutletStateValid()) return PressureOutletStatus::InvalidState;
  boundary[0] = density;
  for (auto d = 0u; d < nDim; ++d) boundary[d + 1] = density * (speed * direction[d] + tangential[d]);
  boundary[nDim + 1] = density * (energy + .5 * (speed * speed + tangentialSpeed2) + turbulentK);
  return PressureOutletStatus::Valid;
}

CFluidModel::ReservoirInflow CFluidModel::PassiveCopy(const ReservoirInflow& data) {
  ReservoirInflow passive = data;
  for (auto* value : {&passive.totalPressure, &passive.totalEnthalpy, &passive.totalEnergy, &passive.totalDensity,
                      &passive.idealGamma, &passive.density, &passive.energy, &passive.normalVelocity,
                      &passive.pressure, &passive.sound, &passive.directionNormal, &passive.gridNormal})
    *value = SU2_TYPE::GetValue(*value);
  return passive;
}

bool CFluidModel::SetReservoirTotalState(su2double totalPressure, su2double totalTemperature, ReservoirInflow& data) {
  SetTDState_PT(totalPressure, totalTemperature);
  if (!IsPressureOutletStateValid()) return false;
  // Real-gas (p,T) inversions stop at a loose tolerance from a history-dependent guess. Polish the density with
  // Newton steps at fixed temperature, from a passive guess, so the reservoir state is reproducible and its
  // derivatives come only from the converged relation p(rho, T_t) = p_t. The slope is a passive central
  // difference because not every EOS provides dT/drho at constant energy.
  su2double totalDensity = SU2_TYPE::GetValue(GetDensity());
  for (auto iteration = 0u; iteration < 4; ++iteration) {
    const double rho = SU2_TYPE::GetValue(totalDensity), step = 1e-6 * rho;
    const double temperature = SU2_TYPE::GetValue(totalTemperature);
    SetTDState_rhoT(rho + step, temperature);
    const double pressurePlus = SU2_TYPE::GetValue(GetPressure());
    SetTDState_rhoT(rho - step, temperature);
    const double slope = (pressurePlus - SU2_TYPE::GetValue(GetPressure())) / (2 * step);
    SetTDState_rhoT(totalDensity, totalTemperature);
    if (!(slope > 0) || !IsPressureOutletStateValid()) return false;
    totalDensity -= (GetPressure() - totalPressure) / slope;
  }
  SetTDState_rhoT(totalDensity, totalTemperature);
  if (!IsPressureOutletStateValid()) return false;
  data.totalPressure = GetPressure();
  data.totalDensity = GetDensity();
  data.totalEnergy = GetStaticEnergy();
  data.totalEnthalpy = data.totalEnergy + data.totalPressure / data.totalDensity;
  return true;
}

bool CFluidModel::FindReservoirSpeedPressure(const ReservoirInflow& data, su2double speed, su2double& pressure) {
  if (!(speed > 0)) {
    pressure = data.totalPressure;
    return true;
  }
  if (data.idealGamma > 1) {
    const su2double gamma = data.idealGamma, enthalpy = data.totalEnthalpy - .5 * speed * speed;
    if (!(enthalpy > 0)) return false;
    pressure = data.totalPressure * pow(enthalpy / data.totalEnthalpy, gamma / (gamma - 1));
    return true;
  }
  // Newton in log pressure on V(p)^2 = speed^2 along the reservoir isentrope, with d(V^2)/d(log p) = -2 p / rho.
  const auto passive = PassiveCopy(data);
  const su2double target = SU2_TYPE::GetValue(speed * speed);
  su2double density, energy, velocity, sound;
  su2double x = log(passive.totalPressure) - .5 * target * passive.totalDensity / passive.totalPressure;
  bool converged = false;
  su2double previousStep = 1;
  for (auto iteration = 0u; iteration < 50 && !converged; ++iteration) {
    const su2double p = exp(x);
    if (!EvaluateReservoirState(passive, p, density, energy, velocity, sound)) return false;
    const su2double step = (velocity * velocity - target) * density / (2 * p);
    x += fmax(su2double(-.5), fmin(su2double(.5), step));
    // Stop at round-off, or where the steps stop shrinking at the accuracy of the EOS integration.
    converged = fabs(step) < 1e-13 || (fabs(step) < 1e-8 && fabs(step) > .5 * previousStep);
    previousStep = fabs(step);
  }
  if (!converged) return false;
  // One Newton correction with active inputs carries the derivatives of the root.
  const su2double p = SU2_TYPE::GetValue(exp(x));
  if (!EvaluateReservoirState(data, p, density, energy, velocity, sound)) return false;
  const double slope = -2 * SU2_TYPE::GetValue(p / density);
  pressure = p * exp(-(velocity * velocity - speed * speed) / slope);
  return true;
}

bool CFluidModel::FindReservoirTotalPressure(su2double restPressure, su2double totalTemperature, su2double speed,
                                             ReservoirInflow& data) {
  if (!SetReservoirTotalState(restPressure, totalTemperature, data)) return false;
  if (!(speed > 0)) return true;
  if (data.idealGamma > 1) {
    // Ideal gas: the total enthalpy depends on the total temperature only.
    const su2double gamma = data.idealGamma, ratio = 1 - .5 * speed * speed / data.totalEnthalpy;
    if (!(ratio > 0)) return false;
    return SetReservoirTotalState(restPressure * pow(ratio, -gamma / (gamma - 1)), totalTemperature, data);
  }
  // Secant iterations in log total pressure on log p(speed) - log restPressure, with passive values.
  const double target = log(SU2_TYPE::GetValue(restPressure));
  const su2double passiveSpeed = SU2_TYPE::GetValue(speed), passiveTemperature = SU2_TYPE::GetValue(totalTemperature);
  auto residual = [&](double y, double& value) {
    ReservoirInflow trial = data;
    su2double pressure;
    if (!SetReservoirTotalState(exp(y), passiveTemperature, trial) ||
        !FindReservoirSpeedPressure(PassiveCopy(trial), passiveSpeed, pressure))
      return false;
    value = log(SU2_TYPE::GetValue(pressure)) - target;
    return true;
  };
  double y0 = target, y1, r0, r1;
  if (!residual(y0, r0)) return false;
  y1 = y0 - r0;  // The residual grows about one-to-one with the log total pressure.
  if (!residual(y1, r1)) return false;
  for (auto iteration = 0u; iteration < 50 && fabs(y1 - y0) > 1e-13 && r1 != r0; ++iteration) {
    const double y2 = y1 - r1 * (y1 - y0) / (r1 - r0);
    y0 = y1;
    r0 = r1;
    y1 = y2;
    if (!residual(y1, r1)) return false;
  }
  double rPlus;
  if (fabs(r1) > 1e-9 || !residual(y1 + 1e-5, rPlus)) return false;
  const double slope = (rPlus - r1) / 1e-5;
  if (!(slope > 0)) return false;
  // One Newton correction with active inputs carries the derivatives of the root.
  su2double pressure;
  if (!SetReservoirTotalState(exp(y1), totalTemperature, data) || !FindReservoirSpeedPressure(data, speed, pressure))
    return false;
  const su2double y = y1 - (log(pressure) - log(restPressure)) / slope;
  return SetReservoirTotalState(exp(y), totalTemperature, data);
}

bool CFluidModel::EvaluateReservoirState(const ReservoirInflow& data, su2double pressure, su2double& density,
                                         su2double& energy, su2double& speed, su2double& sound) {
  if (!(pressure > 0) || !std::isfinite(SU2_TYPE::GetValue(pressure))) return false;
  su2double enthalpy;
  if (data.idealGamma > 1) {
    // Isentropic ideal gas: h = h_t (p/p_t)^((gamma-1)/gamma), rho = rho_t (p/p_t)^(1/gamma).
    const su2double gamma = data.idealGamma, ratio = pressure / data.totalPressure;
    enthalpy = data.totalEnthalpy * pow(ratio, (gamma - 1) / gamma);
    density = data.totalDensity * pow(ratio, 1 / gamma);
    energy = enthalpy / gamma;
    sound = sqrt((gamma - 1) * enthalpy);
  } else {
    // Follow the reservoir isentrope with the outlet integrator; its velocity component is not used.
    PressureOutletState state = {data.totalDensity, data.totalEnergy, 0};
    SetTDState_rhoe(state[0], state[1]);
    if (IntegratePressureOutletIsentrope(pressure, state, false) != PressureOutletStatus::Valid) return false;
    density = state[0];
    energy = state[1];
    enthalpy = energy + GetPressure() / density;
    sound = sqrt(GetSoundSpeed2());
  }
  // Energy conservation from the reservoir at rest; at the total pressure the speed is zero.
  const su2double speed2 = 2 * (data.totalEnthalpy - enthalpy);
  speed = speed2 > 0 ? su2double(sqrt(speed2)) : su2double(0.0);
  return density > 0 && std::isfinite(SU2_TYPE::GetValue(density + energy + sound));
}

bool CFluidModel::EvaluateInteriorInvariantSpeed(const ReservoirInflow& data, su2double pressure, su2double& speed,
                                                 su2double& slope) {
  su2double density, sound, normalVelocity;
  if (data.idealGamma > 1) {
    // Ideal gas: q + 2a/(gamma-1) is constant along the interior isentrope.
    const su2double gamma = data.idealGamma, ratio = pressure / data.pressure;
    density = data.density * pow(ratio, 1 / gamma);
    sound = data.sound * pow(ratio, (gamma - 1) / (2 * gamma));
    normalVelocity = data.normalVelocity + 2 * (data.sound - sound) / (gamma - 1);
  } else {
    // Without the sonic stop the invariant is a smooth, monotone function of pressure for the root search.
    PressureOutletState state = {data.density, data.energy, data.normalVelocity};
    SetTDState_rhoe(data.density, data.energy);
    if (IntegratePressureOutletIsentrope(pressure, state, false) != PressureOutletStatus::Valid) return false;
    density = state[0];
    sound = sqrt(GetSoundSpeed2());
    normalVelocity = state[2];
  }
  // Lab-frame speed along the inflow direction whose normal component matches the invariant.
  speed = (normalVelocity + data.gridNormal) / data.directionNormal;
  // dq/d(log p) = -p/(rho*a), divided by the (negative) normal component of the direction.
  slope = -pressure / (density * sound * data.directionNormal);
  return true;
}

bool CFluidModel::FindReservoirSonicPressure(const ReservoirInflow& data, su2double& pressure) {
  if (data.idealGamma > 1) {
    const su2double gamma = data.idealGamma;
    pressure = data.totalPressure * pow(2 / (gamma + 1), gamma / (gamma - 1));
    return true;
  }
  // Bracket and bisect V^2 - a^2, which falls from positive at low pressure to -a_t^2 at the total pressure.
  su2double density, energy, speed, sound;
  auto excess = [&](su2double x, su2double& value) {
    if (!EvaluateReservoirState(data, exp(x), density, energy, speed, sound)) return false;
    value = speed * speed - sound * sound;
    return true;
  };
  su2double upper = log(data.totalPressure), lower = upper, value = -1;
  for (auto step = 0u; step < 60 && value <= 0; ++step) {
    upper = lower;
    lower -= .5;
    if (!excess(lower, value)) return false;
  }
  if (value <= 0) return false;
  for (auto iteration = 0u; iteration < 100 && upper - lower > 1e-14; ++iteration) {
    const su2double middle = .5 * (lower + upper);
    if (!excess(middle, value)) return false;
    if (value > 0)
      lower = middle;
    else
      upper = middle;
  }
  pressure = exp(.5 * (lower + upper));
  return true;
}

/*!
 * \brief Fallback for backflow that enters through every characteristic when no reservoir data is given.
 *
 * The interior outgoing acoustic relation is followed to the imposed pressure as for subsonic backflow
 * (ComputePressureOutletState without its check on the interior normal speed). If the result is subsonic it is
 * used; otherwise the boundary keeps the imposed pressure and the interior entropy and tangential velocity, and the
 * inward normal velocity relative to the grid is limited to the local sound speed. Both limits of subsonic backflow
 * (interior speed reaching -a, boundary speed reaching -a) are therefore continuous.
 */
CFluidModel::PressureOutletStatus CFluidModel::ComputeSonicBackflowState(unsigned short nDim, const su2double* interior,
                                                                         const su2double* normal, su2double gridNormal,
                                                                         su2double targetPressure, su2double* boundary,
                                                                         su2double turbulentK, su2double idealGamma) {
  for (auto v = 0u; v < nDim + 2; ++v)
    if (!std::isfinite(SU2_TYPE::GetValue(interior[v]))) return PressureOutletStatus::InvalidState;
  const auto rho = interior[0];
  if (!(rho > 0) || !(targetPressure > 0) || !std::isfinite(SU2_TYPE::GetValue(targetPressure)))
    return PressureOutletStatus::InvalidState;
  su2double velocity[3] = {};
  for (auto d = 0u; d < nDim; ++d) velocity[d] = interior[d + 1] / rho;
  const su2double energy = interior[nDim + 1] / rho - .5 * GeometryToolbox::SquaredNorm(nDim, velocity) - turbulentK;
  SetTDState_rhoe(rho, energy);
  if (!IsPressureOutletStateValid()) return PressureOutletStatus::InvalidState;
  const su2double normalVelocity = GeometryToolbox::DotProduct(nDim, velocity, normal) - gridNormal;
  const su2double interiorPressure = GetPressure();

  PressureOutletState state = {rho, energy, normalVelocity};
  const auto status = idealGamma > 1 ? ComputeIdealGasPressureOutletState(targetPressure, idealGamma, state)
                                     : IntegratePressureOutletIsentrope(targetPressure, state);
  if (status == PressureOutletStatus::InvalidState) return status;
  if (status == PressureOutletStatus::Backflow) {
    // Supersonic inflow: thermodynamic state at the imposed pressure on the interior isentrope, sonic inflow.
    if (idealGamma > 1)
      SetTDState_Prho(targetPressure, rho * pow(targetPressure / interiorPressure, 1 / idealGamma));
    else
      SetTDState_rhoe(state[0], state[1]);  // The integration ended at the imposed pressure.
    if (!IsPressureOutletStateValid() || fabs(GetPressure() / targetPressure - 1) > 1e-6)
      return PressureOutletStatus::InvalidState;
    state = {GetDensity(), GetStaticEnergy(), su2double(-sqrt(GetSoundSpeed2()))};
  }
  for (auto d = 0u; d < nDim; ++d) velocity[d] += (state[2] - normalVelocity) * normal[d];
  boundary[0] = state[0];
  for (auto d = 0u; d < nDim; ++d) boundary[d + 1] = state[0] * velocity[d];
  boundary[nDim + 1] = state[0] * (state[1] + .5 * GeometryToolbox::SquaredNorm(nDim, velocity) + turbulentK);
  return PressureOutletStatus::Valid;
}
