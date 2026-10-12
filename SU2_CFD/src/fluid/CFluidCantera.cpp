/*!
 * \file CFluidCantera.cpp
 * \brief Defines the multicomponent incompressible Ideal Gas model for reacting flows.
 * \author T. Economon, Cristopher Morales Ubal
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

#include "../../include/fluid/CFluidCantera.hpp"

#include <algorithm>
#include <numeric>
#include <cmath>

#ifdef HAVE_CANTERA
#include <cantera/core.h>
#include <cantera/kinetics/Reaction.h>

using namespace Cantera;

namespace {
std::shared_ptr<Solution> LoadSolution(const CConfig* config) {
  try {
    return newSolution(config->GetChemical_MechanismFile(), config->GetPhase_Name(), config->GetTransport_Model());
  } catch (const CanteraError& error) {
    SU2_MPI::Error(string("Cantera error during initialization: ") + error.what(), CURRENT_FUNCTION);
    return nullptr;
  }
}

std::vector<double> GetMolarMasses(Solution& solution) {
  std::vector<double> molarMasses(solution.thermo()->nSpecies());
  solution.thermo()->getMolecularWeights(molarMasses.data());
  return molarMasses;
}

/*!
 * \brief Newton-Raphson on temperature for value(T) = target, kept inside [minT, maxT].
 * \note Without a positive initial temperature, the start is one Newton step from the reference temperature.
 * \return Whether the iteration converged.
 */
template <class SetState, class Value, class Slope>
bool SolveTemperature(double target, double& temperature, double minT, double maxT, const SetState& setState,
                      const Value& value, const Slope& slope) {
  /*--- Temperature tolerance in [K], high accuracy needed for restarts. ---*/
  constexpr double tolerance = 1e-5;
  constexpr int maxIterations = 20;

  if (!(temperature > 0.0)) {
    setState(STD_REF_TEMP);
    temperature = STD_REF_TEMP + (target - value()) / slope();
  }
  temperature = std::clamp(temperature, minT, maxT);

  for (int iter = 0; iter < maxIterations; iter++) {
    setState(temperature);
    const double delta = (target - value()) / slope();
    if (!std::isfinite(delta)) return false;

    const double next = std::clamp(temperature + delta, minT, maxT);
    const bool inRange = (next == temperature + delta);
    temperature = next;
    if (!inRange) return false;
    if (std::abs(delta) <= tolerance) return true;
  }
  return false;
}
}  // namespace

CFluidCantera::CFluidCantera(su2double value_pressure_operating, const CConfig* config)
    : CFluidModel(),
      n_species_mixture(config->GetnSpecies() + 1),
      Pressure_Thermodynamic(value_pressure_operating),
      Prandtl_Turb_Number(config->GetPrandtl_Turb()),
      Schmidt_Turb_Number(config->GetSchmidt_Number_Turbulent()),
      Transport_Model(config->GetTransport_Model()),
      Chemical_MechanismFile(config->GetChemical_MechanismFile()),
      Phase_Name(config->GetPhase_Name()),
      Combustion(config->GetCanteraCombustion()),
      Chemistry_Min_Temperature(config->GetCantera_DC_Min_Temp()),
      Correction_Velocity(config->GetCantera_Correction_Velocity()),
      Turbulent_Conductivity(config->GetKind_Regime() == ENUM_REGIME::INCOMPRESSIBLE),
      sol(LoadSolution(config)),
      molarMasses(GetMolarMasses(*sol)) {
  try {
    const auto& thermo = *sol->thermo();
    const size_t nSpeciesMechanism = thermo.nSpecies();

    massFractions.assign(nSpeciesMechanism, 0.0);
    massDiffusivity.assign(nSpeciesMechanism, 0.0);
    netProductionRates.assign(nSpeciesMechanism, 0.0);
    destructionRates.assign(nSpeciesMechanism, 0.0);
    enthalpiesSpecies.assign(nSpeciesMechanism, 0.0);
    specificHeatSpecies.assign(nSpeciesMechanism, 0.0);
    /*--- Moderate extrapolation of the thermodynamic polynomials is allowed beyond the range of the data. ---*/
    minTemperature = 0.5 * thermo.minTemp();
    maxTemperature = 2.0 * thermo.maxTemp();

    speciesIndices.resize(n_species_mixture);
    if (config->GetnCantera_SpeciesNames() == 0) {
      if (n_species_mixture != nSpeciesMechanism) {
        SU2_MPI::Error("Without CANTERA_SPECIES_NAMES all " + to_string(nSpeciesMechanism) + " species of phase '" +
                           Phase_Name + "' in " + Chemical_MechanismFile + " are transported, which requires " +
                           to_string(nSpeciesMechanism - 1) + " entries in SPECIES_INIT.",
                       CURRENT_FUNCTION);
      }
      std::iota(speciesIndices.begin(), speciesIndices.end(), size_t(0));
    }
    for (size_t iVar = 0; iVar < config->GetnCantera_SpeciesNames(); iVar++) {
      const string name = config->GetChemical_GasComposition(iVar);
      const size_t index = thermo.speciesIndex(name);
      if (index == npos) {
        SU2_MPI::Error("Species '" + name + "' of CANTERA_SPECIES_NAMES is not part of phase '" + Phase_Name +
                           "' in " + Chemical_MechanismFile + ".",
                       CURRENT_FUNCTION);
      }
      if (std::find(speciesIndices.begin(), speciesIndices.begin() + iVar, index) != speciesIndices.begin() + iVar) {
        SU2_MPI::Error("Species '" + name + "' appears more than once in CANTERA_SPECIES_NAMES.", CURRENT_FUNCTION);
      }
      speciesIndices[iVar] = index;
    }

    if (SU2_MPI::GetRank() == MASTER_NODE && n_species_mixture < nSpeciesMechanism) {
      std::vector<bool> transported(nSpeciesMechanism, false);
      for (const auto k : speciesIndices) transported[k] = true;
      string excluded;
      for (size_t k = 0; k < nSpeciesMechanism; k++) {
        if (!transported[k]) excluded += (excluded.empty() ? "" : ", ") + thermo.speciesName(k);
      }
      cout << "WARNING: " << nSpeciesMechanism - n_species_mixture << " species of phase '" << Phase_Name
           << "' are not in CANTERA_SPECIES_NAMES and are not transported (zero mass fraction, no source terms): "
           << excluded << "." << endl;
    }

    chemicalSourceTerm.assign(n_species_mixture, 0.0);
    chemicalSourceJacobian.assign(n_species_mixture, 0.0);
    enthalpyFormation.assign(n_species_mixture, 0.0);
    if (Combustion) SetEnthalpyFormation();
  } catch (const CanteraError& error) {
    SU2_MPI::Error(string("Cantera error during initialization: ") + error.what(), CURRENT_FUNCTION);
  }
}

void CFluidCantera::SetEnthalpyFormation() {
  sol->thermo()->setState_TP(STD_REF_TEMP, SU2_TYPE::GetValue(Pressure_Thermodynamic));
  sol->thermo()->getEnthalpy_RT_ref(enthalpiesSpecies.data());
  for (size_t iVar = 0; iVar < n_species_mixture; iVar++) {
    const size_t k = speciesIndices[iVar];
    enthalpyFormation[iVar] = GasConstant * STD_REF_TEMP * enthalpiesSpecies[k] / molarMasses[k];
  }
}

void CFluidCantera::ComputeChemicalSourceTerm() {
  if (Temperature < Chemistry_Min_Temperature) {
    std::fill(chemicalSourceTerm.begin(), chemicalSourceTerm.end(), 0.0);
    std::fill(chemicalSourceJacobian.begin(), chemicalSourceJacobian.end(), 0.0);
    Heat_Release = 0.0;
    return;
  }
  try {
    sol->kinetics()->getNetProductionRates(netProductionRates.data());
    sol->kinetics()->getDestructionRates(destructionRates.data());
  } catch (const CanteraError& error) {
    SU2_MPI::Error(string("Cantera error in the kinetics evaluation: ") + error.what(), CURRENT_FUNCTION);
  }
  Heat_Release = 0.0;
  for (size_t iVar = 0; iVar < n_species_mixture; iVar++) {
    const size_t k = speciesIndices[iVar];
    chemicalSourceTerm[iVar] = molarMasses[k] * netProductionRates[k];
    /*--- The destruction rate is taken proportional to the mass fraction; the derivative is w.r.t. rho*Y. ---*/
    chemicalSourceJacobian[iVar] =
        -molarMasses[k] * destructionRates[k] / (Density * std::max(massFractions[k], 1e-10));
    Heat_Release -= enthalpyFormation[iVar] * chemicalSourceTerm[iVar];
  }
}

void CFluidCantera::GetEnthalpyDiffusivity(su2double* enthalpy_diffusions) const {
  UpdateDiffusivity();
  sol->thermo()->getEnthalpy_RT_ref(enthalpiesSpecies.data());
  const su2double RT = GasConstant * Temperature;
  const size_t last = speciesIndices[n_species_mixture - 1];
  /*--- With the correction velocity every species enthalpy is taken relative to the mixture enthalpy. ---*/
  const su2double h_ref = Correction_Velocity ? Enthalpy : su2double(0.0);
  const su2double h_last = RT * enthalpiesSpecies[last] / molarMasses[last] - h_ref;
  for (size_t iVar = 0; iVar < n_species_mixture - 1; iVar++) {
    const size_t k = speciesIndices[iVar];
    const su2double h_k = RT * enthalpiesSpecies[k] / molarMasses[k] - h_ref;
    enthalpy_diffusions[iVar] = Density * (h_k * massDiffusivity[k] - h_last * massDiffusivity[last]) +
                                Mu_Turb * (h_k - h_last) / Schmidt_Turb_Number;
  }
}

void CFluidCantera::GetGradEnthalpyDiffusivity(su2double* grad_enthalpy_diffusions) const {
  UpdateDiffusivity();
  sol->thermo()->getCp_R_ref(specificHeatSpecies.data());
  const size_t last = speciesIndices[n_species_mixture - 1];
  const su2double cp_ref = Correction_Velocity ? Cp : su2double(0.0);
  const su2double cp_last = GasConstant * specificHeatSpecies[last] / molarMasses[last] - cp_ref;
  for (size_t iVar = 0; iVar < n_species_mixture - 1; iVar++) {
    const size_t k = speciesIndices[iVar];
    const su2double cp_k = GasConstant * specificHeatSpecies[k] / molarMasses[k] - cp_ref;
    grad_enthalpy_diffusions[iVar] = Density * (cp_k * massDiffusivity[k] - cp_last * massDiffusivity[last]) +
                                     Mu_Turb * (cp_k - cp_last) / Schmidt_Turb_Number;
  }
}

void CFluidCantera::SetMassFractions(const su2double* val_scalars) {
  double scalarsSum = 0.0;
  std::fill(massFractions.begin(), massFractions.end(), 0.0);
  for (size_t iScalar = 0; iScalar < n_species_mixture - 1; iScalar++) {
    const double value = SU2_TYPE::GetValue(val_scalars[iScalar]);
    massFractions[speciesIndices[iScalar]] = value;
    scalarsSum += value;
  }
  massFractions[speciesIndices[n_species_mixture - 1]] = 1.0 - scalarsSum;
}

void CFluidCantera::SetThermoState(double val_temperature) {
  sol->thermo()->setState_TP(val_temperature, SU2_TYPE::GetValue(Pressure_Thermodynamic));
  ReadThermoState();
}

void CFluidCantera::ReadThermoState() {
  const auto& thermo = *sol->thermo();
  Temperature = thermo.temperature();
  Density = thermo.density();
  Pressure = thermo.pressure();
  Enthalpy = thermo.enthalpy_mass();
  StaticEnergy = thermo.intEnergy_mass();
  Entropy = thermo.entropy_mass();
  Cp = thermo.cp_mass();
  Cv = thermo.cv_mass();

  /*--- Ideal gas mixture at frozen composition: P = rho*R*T and e = e(T). ---*/
  const su2double gasConstant = Pressure / (Density * Temperature);
  dPdrho_e = gasConstant * Temperature;
  dPde_rho = Density * gasConstant / Cv;
  dTdrho_e = 0.0;
  dTde_rho = 1.0 / Cv;
  SoundSpeed2 = Cp / Cv * gasConstant * Temperature;

  transportValid = false;
  diffusivityValid = false;
}

void CFluidCantera::SetComposition(const su2double* val_scalars) {
  SetMassFractions(val_scalars);
  try {
    sol->thermo()->setMassFractions(massFractions.data());
  } catch (const CanteraError& error) {
    SU2_MPI::Error(string("Cantera error in SetComposition: ") + error.what(), CURRENT_FUNCTION);
  }
  transportValid = false;
  diffusivityValid = false;
}

void CFluidCantera::SetTDState_rhoe(su2double rho, su2double e) {
  double temperature = temperatureGuess;
  temperatureGuess = 0.0;
  bool converged = false;

  try {
    auto& thermo = *sol->thermo();
    const double density = SU2_TYPE::GetValue(rho);
    converged = SolveTemperature(
        SU2_TYPE::GetValue(e), temperature, minTemperature, maxTemperature,
        [&](double T) { thermo.setState_TD(T, density); }, [&] { return thermo.intEnergy_mass(); },
        [&] { return thermo.cv_mass(); });
    thermo.setState_TD(temperature, density);
    ReadThermoState();
  } catch (const CanteraError& error) {
    SU2_MPI::Error(string("Cantera error in SetTDState_rhoe: ") + error.what(), CURRENT_FUNCTION);
  }
  temperatureIterationFailed = !converged;
}

void CFluidCantera::SetTDState_PT(su2double P, su2double T) {
  try {
    sol->thermo()->setState_TP(SU2_TYPE::GetValue(T), SU2_TYPE::GetValue(P));
    ReadThermoState();
  } catch (const CanteraError& error) {
    SU2_MPI::Error(string("Cantera error in SetTDState_PT: ") + error.what(), CURRENT_FUNCTION);
  }
  temperatureIterationFailed = false;
}

void CFluidCantera::SetTDState_Prho(su2double P, su2double rho) {
  try {
    auto& thermo = *sol->thermo();
    const double density = SU2_TYPE::GetValue(rho);
    const double gasConstant = GasConstant / thermo.meanMolecularWeight();
    thermo.setState_TD(SU2_TYPE::GetValue(P) / (density * gasConstant), density);
    ReadThermoState();
  } catch (const CanteraError& error) {
    SU2_MPI::Error(string("Cantera error in SetTDState_Prho: ") + error.what(), CURRENT_FUNCTION);
  }
  temperatureIterationFailed = false;
}

void CFluidCantera::SetTDState_rhoT(su2double rho, su2double T) {
  try {
    sol->thermo()->setState_TD(SU2_TYPE::GetValue(T), SU2_TYPE::GetValue(rho));
    ReadThermoState();
  } catch (const CanteraError& error) {
    SU2_MPI::Error(string("Cantera error in SetTDState_rhoT: ") + error.what(), CURRENT_FUNCTION);
  }
  temperatureIterationFailed = false;
}

void CFluidCantera::SetTDState_hs(su2double h, su2double s) {
  double temperature = 0.0;
  bool converged = false;

  try {
    auto& thermo = *sol->thermo();
    const double pressureRef = OneAtm;
    converged = SolveTemperature(
        SU2_TYPE::GetValue(h), temperature, minTemperature, maxTemperature,
        [&](double T) { thermo.setState_TP(T, pressureRef); }, [&] { return thermo.enthalpy_mass(); },
        [&] { return thermo.cp_mass(); });
    thermo.setState_TP(temperature, pressureRef);

    /*--- s(T, P) = s(T, P_ref) - R ln(P / P_ref). ---*/
    const double gasConstant = GasConstant / thermo.meanMolecularWeight();
    const double pressure = pressureRef * std::exp((thermo.entropy_mass() - SU2_TYPE::GetValue(s)) / gasConstant);
    thermo.setState_TP(temperature, pressure);
    ReadThermoState();
  } catch (const CanteraError& error) {
    SU2_MPI::Error(string("Cantera error in SetTDState_hs: ") + error.what(), CURRENT_FUNCTION);
  }
  temperatureIterationFailed = !converged;
}

void CFluidCantera::SetTDState_Ps(su2double P, su2double s) {
  try {
    sol->thermo()->setState_SP(SU2_TYPE::GetValue(s), SU2_TYPE::GetValue(P));
    ReadThermoState();
  } catch (const CanteraError& error) {
    SU2_MPI::Error(string("Cantera error in SetTDState_Ps: ") + error.what(), CURRENT_FUNCTION);
  }
  temperatureIterationFailed = false;
}

void CFluidCantera::UpdateTransport() const {
  if (transportValid) return;
  try {
    laminarViscosity = sol->transport()->viscosity();
    laminarConductivity = sol->transport()->thermalConductivity();
  } catch (const CanteraError& error) {
    SU2_MPI::Error(string("Cantera error in the transport evaluation: ") + error.what(), CURRENT_FUNCTION);
  }
  transportValid = true;
}

void CFluidCantera::UpdateDiffusivity() const {
  if (diffusivityValid) return;
  try {
    sol->transport()->getMixDiffCoeffsMass(massDiffusivity.data());
  } catch (const CanteraError& error) {
    SU2_MPI::Error(string("Cantera error in the diffusivity evaluation: ") + error.what(), CURRENT_FUNCTION);
  }
  diffusivityValid = true;
}

void CFluidCantera::SetTDState_T(const su2double val_temperature, const su2double* val_scalars) {
  SetMassFractions(val_scalars);
  try {
    sol->thermo()->setMassFractions(massFractions.data());
    SetThermoState(SU2_TYPE::GetValue(val_temperature));
  } catch (const CanteraError& error) {
    SU2_MPI::Error(string("Cantera error in SetTDState_T: ") + error.what(), CURRENT_FUNCTION);
  }
  temperatureIterationFailed = false;
}

void CFluidCantera::SetTDState_h(const su2double val_enthalpy, const su2double* val_scalars) {
  SetMassFractions(val_scalars);
  double temperature = temperatureGuess;
  temperatureGuess = 0.0;
  bool converged = false;

  try {
    auto& thermo = *sol->thermo();
    thermo.setMassFractions(massFractions.data());
    const double pressure = SU2_TYPE::GetValue(Pressure_Thermodynamic);
    converged = SolveTemperature(
        SU2_TYPE::GetValue(val_enthalpy), temperature, minTemperature, maxTemperature,
        [&](double T) { thermo.setState_TP(T, pressure); }, [&] { return thermo.enthalpy_mass(); },
        [&] { return thermo.cp_mass(); });
    SetThermoState(temperature);
  } catch (const CanteraError& error) {
    SU2_MPI::Error(string("Cantera error in SetTDState_h: ") + error.what(), CURRENT_FUNCTION);
  }
  temperatureIterationFailed = !converged;
}

#else
CFluidCantera::CFluidCantera(su2double value_pressure_operating, const CConfig* config) {
  SU2_MPI::Error(
      "FLUID_CANTERA requires SU2 to be compiled with Cantera (-Denable-cantera=true). "
      "Automatic differentiation builds do not support Cantera.",
      CURRENT_FUNCTION);
}
#endif
