/*!
 * \file CFluidCantera.cpp
 * \brief Defines the multicomponent incompressible Ideal Gas model for reacting flows.
 * \author T. Economon, Cristopher Morales Ubal
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

#include "../../include/fluid/CFluidCantera.hpp"

#include <algorithm>
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
      Combustion(config->GetCombustion()),
      Chemistry_Min_Temperature(config->GetCantera_DC_Min_Temp()),
      Correction_Velocity(config->GetCantera_Correction_Velocity()),
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
    for (size_t iVar = 0; iVar < n_species_mixture; iVar++) {
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

    chemicalSourceTerm.assign(n_species_mixture, 0.0);
    chemicalSourceJacobian.assign(n_species_mixture, 0.0);
    enthalpyFormation.assign(n_species_mixture, 0.0);
    if (Combustion) SetEnthalpyFormation(config);
  } catch (const CanteraError& error) {
    SU2_MPI::Error(string("Cantera error during initialization: ") + error.what(), CURRENT_FUNCTION);
  }
}

void CFluidCantera::SetEnthalpyFormation(const CConfig* config) {
  SetMassFractions(config->GetSpecies_Init());
  sol->thermo()->setMassFractions(massFractions.data());
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
  auto& thermo = *sol->thermo();
  thermo.setState_TP(val_temperature, SU2_TYPE::GetValue(Pressure_Thermodynamic));
  Temperature = val_temperature;
  Density = thermo.density();
  Enthalpy = thermo.enthalpy_mass();
  Cp = thermo.cp_mass();
  Cv = thermo.cv_mass();
  transportValid = false;
  diffusivityValid = false;
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
  /*--- Temperature tolerance in [K], high accuracy needed for restarts. ---*/
  constexpr double tolerance = 1e-5;
  constexpr int maxIterations = 20;

  SetMassFractions(val_scalars);
  const double enthalpy = SU2_TYPE::GetValue(val_enthalpy);
  double temperature = temperatureGuess;
  temperatureGuess = 0.0;
  bool converged = false;

  try {
    sol->thermo()->setMassFractions(massFractions.data());

    /*--- Without a guess, start from a single Newton step at the reference temperature. ---*/
    if (!(temperature > 0.0)) {
      sol->thermo()->setState_TP(STD_REF_TEMP, SU2_TYPE::GetValue(Pressure_Thermodynamic));
      temperature = STD_REF_TEMP + (enthalpy - sol->thermo()->enthalpy_mass()) / sol->thermo()->cp_mass();
    }
    temperature = std::clamp(temperature, minTemperature, maxTemperature);

    /*--- Newton-Raphson on temperature, kept inside the temperature bounds. ---*/
    for (int iter = 0; iter < maxIterations; iter++) {
      sol->thermo()->setState_TP(temperature, SU2_TYPE::GetValue(Pressure_Thermodynamic));
      const double delta = (enthalpy - sol->thermo()->enthalpy_mass()) / sol->thermo()->cp_mass();
      if (!std::isfinite(delta)) break;

      const double next = std::clamp(temperature + delta, minTemperature, maxTemperature);
      const bool inRange = (next == temperature + delta);
      temperature = next;
      if (!inRange) break;
      if (std::abs(delta) <= tolerance) {
        converged = true;
        break;
      }
    }
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
