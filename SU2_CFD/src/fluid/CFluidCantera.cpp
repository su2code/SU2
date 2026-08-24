/*!
 * \file CFluidCantera.cpp
 * \brief Defines the multicomponent incompressible Ideal Gas model for reacting flows.
 * \author T. Economon, Cristopher Morales Ubal
 * \version 8.4.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2024, SU2 Contributors (cf. AUTHORS.md)
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

#include <cmath>

// #include <numeric>
#ifdef USE_CANTERA
#include <cantera/core.h>
#include <cantera/kinetics/Reaction.h>

using namespace Cantera;

CFluidCantera::CFluidCantera(su2double value_pressure_operating, const CConfig* config)
    : CFluidModel(),
      n_species_mixture(config->GetnSpecies() + 1),
      Pressure_Thermodynamic(value_pressure_operating),
      GasConstant_Ref(config->GetGas_Constant_Ref()),
      Prandtl_Turb_Number(config->GetPrandtl_Turb()),
      Schmidt_Turb_Number(config->GetSchmidt_Number_Turbulent()),
      Transport_Model(config->GetTransport_Model()),
      Chemical_MechanismFile(config->GetChemical_MechanismFile()),
      Phase_Name(config->GetPhase_Name()),
      Combustion(config->GetCombustion()) {
  if (n_species_mixture > ARRAYSIZE) {
    SU2_MPI::Error("Too many species, increase ARRAYSIZE", CURRENT_FUNCTION);
  }
  sol = std::shared_ptr<Cantera::Solution>(newSolution(Chemical_MechanismFile, Phase_Name, Transport_Model));
  sol->thermo()->getMolecularWeights(&molarMasses[0]);
  for (int iVar = 0; iVar < n_species_mixture; iVar++) {
    speciesIndices[iVar] = sol->thermo()->speciesIndex(config->GetChemical_GasComposition(iVar));
  }
  enthalpiesSpecies.resize(sol->thermo()->nSpecies());
  specificHeatSpecies.resize(sol->thermo()->nSpecies());
  netProductionRates.resize(sol->thermo()->nSpecies());
  if (Combustion) SetEnthalpyFormation(config);
}

void CFluidCantera::SetEnthalpyFormation(const CConfig* config) {
  SetMassFractions(config->GetSpecies_Init());
  sol->thermo()->setMassFractions(massFractions.data());
  sol->thermo()->setState_TP(STD_REF_TEMP, Pressure_Thermodynamic);
  sol->thermo()->getEnthalpy_RT_ref(enthalpiesSpecies.data());
  for (int iVar = 0; iVar < n_species_mixture; iVar++) {
    enthalpyFormation[iVar] =
        GasConstant * STD_REF_TEMP * enthalpiesSpecies[speciesIndices[iVar]] / molarMasses[speciesIndices[iVar]];
  }
}

void CFluidCantera::ComputeChemicalSourceTerm() {
  sol->kinetics()->getNetProductionRates(netProductionRates.data());
  for (int iVar = 0; iVar < n_species_mixture - 1.0; iVar++) {
    chemicalSourceTerm[iVar] = molarMasses[speciesIndices[iVar]] * netProductionRates[speciesIndices[iVar]];
  }
}

void CFluidCantera::ComputeHeatRelease() {
  sol->kinetics()->getNetProductionRates(netProductionRates.data());
  Heat_Release = 0.0;
  for (int iVar = 0; iVar < n_species_mixture; iVar++) {
    Heat_Release +=
        -1.0 * enthalpyFormation[iVar] * molarMasses[speciesIndices[iVar]] * netProductionRates[speciesIndices[iVar]];
  }
}

void CFluidCantera::GetEnthalpyDiffusivity(su2double* enthalpy_diffusions) const {
  sol->thermo()->getEnthalpy_RT_ref(enthalpiesSpecies.data());
  for (int iVar = 0; iVar < n_species_mixture - 1; iVar++) {
    enthalpy_diffusions[iVar] =
        Density * GasConstant * Temperature *
        ((enthalpiesSpecies[speciesIndices[iVar]] * massDiffusivity[speciesIndices[iVar]] / molarMasses[speciesIndices[iVar]]) -
         (enthalpiesSpecies[speciesIndices[n_species_mixture - 1]] * massDiffusivity[speciesIndices[n_species_mixture - 1]] /
          molarMasses[speciesIndices[n_species_mixture - 1]]));
    enthalpy_diffusions[iVar] +=
        Mu_Turb * GasConstant * Temperature *
        ((enthalpiesSpecies[speciesIndices[iVar]] / molarMasses[speciesIndices[iVar]]) -
         (enthalpiesSpecies[speciesIndices[n_species_mixture - 1]] /
          molarMasses[speciesIndices[n_species_mixture - 1]])) / Schmidt_Turb_Number;
  }
}

void CFluidCantera::GetGradEnthalpyDiffusivity(su2double* grad_enthalpy_diffusions) const {
  sol->thermo()->getCp_R_ref(specificHeatSpecies.data());
  for (int iVar = 0; iVar < n_species_mixture - 1; iVar++) {
    grad_enthalpy_diffusions[iVar] =
        Density * GasConstant *
        ((specificHeatSpecies[speciesIndices[iVar]] * massDiffusivity[speciesIndices[iVar]] / molarMasses[speciesIndices[iVar]]) -
         (specificHeatSpecies[speciesIndices[n_species_mixture - 1]] * massDiffusivity[speciesIndices[n_species_mixture - 1]] /
          molarMasses[speciesIndices[n_species_mixture - 1]]));
    grad_enthalpy_diffusions[iVar] += Mu_Turb * GasConstant *
                                      ((specificHeatSpecies[speciesIndices[iVar]] / molarMasses[speciesIndices[iVar]]) -
                                       (specificHeatSpecies[speciesIndices[n_species_mixture - 1]] /
                                        molarMasses[speciesIndices[n_species_mixture - 1]])) /
                                      Schmidt_Turb_Number;
  }
}

void CFluidCantera::SetMassFractions(const su2double* val_scalars) {
  su2double val_scalars_sum{0.0};
  massFractions.fill(0.0);
  for (int i_scalar = 0; i_scalar < n_species_mixture - 1; i_scalar++) {
    massFractions[speciesIndices[i_scalar]] = val_scalars[i_scalar];
    val_scalars_sum += val_scalars[i_scalar];
  }
  massFractions[speciesIndices[n_species_mixture - 1]] = 1.0 - val_scalars_sum;
}

void CFluidCantera::SetTDState_T(const su2double val_temperature, const su2double* val_scalars) {
  Temperature = val_temperature;
  SetMassFractions(val_scalars);
  
  sol->thermo()->setMassFractions(massFractions.data());
  sol->thermo()->setState_TP(Temperature, Pressure_Thermodynamic);
  Density = sol->thermo()->density();
  Enthalpy = sol->thermo()->enthalpy_mass();
  Cp = sol->thermo()->cp_mass();
  Cv = sol->thermo()->cv_mass();
  Mu = sol->transport()->viscosity();
  Kt = sol->transport()->thermalConductivity();

  sol->transport()->getMixDiffCoeffsMass(massDiffusivity.data());
  if (Combustion) ComputeHeatRelease();
}

void CFluidCantera::SetTDState_h(const su2double val_enthalpy, const su2double* val_scalars) {
  Enthalpy = val_enthalpy;
  /*--- convergence criterion for temperature in [K], high accuracy needed for restarts. ---*/
  su2double toll = 1e-5;
  su2double temp_iter = 300.0;
  su2double delta_temp_iter = 1e10;
  su2double delta_enthalpy_iter;
  const int counter_limit = 20;

  int counter = 0;

  /*--- Set mass fractions. ---*/
  SetMassFractions(val_scalars);
  sol->thermo()->setMassFractions(massFractions.data());

  /*--- Computing temperature given enthalpy and species mass fractions using Newton-Raphson. ---*/
  while ((abs(delta_temp_iter) > toll) && (counter++ < counter_limit)) {
    /*--- Set thermodynamic state based on the current value of temperature. ---*/
    sol->thermo()->setState_TP(temp_iter, Pressure_Thermodynamic);

    su2double Enthalpy_iter = sol->thermo()->enthalpy_mass();
    su2double Cp_iter = sol->thermo()->cp_mass();

    delta_enthalpy_iter = Enthalpy - Enthalpy_iter;

    delta_temp_iter = delta_enthalpy_iter / Cp_iter;

    temp_iter += delta_temp_iter;
    if (temp_iter < 0.0) {
      cout << "Warning: Negative temperature has been found during Newton-Raphson" << endl;
      break;
    }
  }
  Temperature = temp_iter;
  if (counter == counter_limit) {
    cout << "Warning Newton-Raphson exceed number of max iteration in temperature computation" << endl;
  }
  SetTDState_T(Temperature, val_scalars);
}

#else
CFluidCantera::CFluidCantera(su2double value_pressure_operating, const CConfig* config) {
  SU2_MPI::Error("SU2 was not compiled with Cantera(-Denable-cantera=true)", CURRENT_FUNCTION);
}
#endif