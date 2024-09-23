/*!
 * \file CFluidCantera.cpp
 * \brief Defines the multicomponent incompressible Ideal Gas model for reacting flows.
 * \author T. Economon, Cristopher Morales Ubal
 * \version 8.0.1 "Harrier"
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

#include <numeric>
#ifdef USE_CANTERA
#include "../../Common/include/basic_types/ad_structure.hpp"

#include "/home/cristopher/codes/cantera/include/cantera/core.h"
#include <fstream>
#include <iostream>

using namespace Cantera;
using namespace SU2_TYPE;
#endif

CFluidCantera::CFluidCantera(su2double value_pressure_operating, const CConfig* config)
    : CFluidModel(),
      n_species_mixture(config->GetnSpecies() + 1),
      Pressure_Thermodynamic(value_pressure_operating),
      GasConstant_Ref(config->GetGas_Constant_Ref()),
      Prandtl_Number(config->GetPrandtl_Turb()),
      Transport_Model(config->GetTransport_Model()),
      Chemical_MechanismFile(config->GetChemical_MechanismFile()),
      Phase_Name(config->GetPhase_Name()) {
  if (n_species_mixture > ARRAYSIZE) {
    SU2_MPI::Error("Too many species, increase ARRAYSIZE", CURRENT_FUNCTION);
  }

  #ifdef USE_CANTERA
  for (int iVar = 0; iVar < n_species_mixture; iVar++) { 
    gasComposition[iVar]=config->GetChemical_GasComposition(iVar);
  }
  #endif

  SetMassDiffusivityModel(config);
}

void CFluidCantera::SetMassDiffusivityModel(const CConfig* config) {
  for (int iVar = 0; iVar < n_species_mixture; iVar++) {
    MassDiffusivityPointers[iVar] = MakeMassDiffusivityModel(config, iVar);
  }
}

void CFluidCantera::ComputeMassDiffusivity() {
  for (int iVar = 0; iVar < n_species_mixture; iVar++) {
    MassDiffusivityPointers[iVar]->SetDiffusivity(Density, Mu, Cp, Kt);
    massDiffusivity[iVar] = MassDiffusivityPointers[iVar]->GetDiffusivity();
  }
}

#ifdef USE_CANTERA
string CFluidCantera::DictionaryChemicalComposition(const su2double* val_scalars) {
  su2double val_scalars_sum{0.0};
  chemical_composition="";
  for (int i_scalar = 0; i_scalar < n_species_mixture - 1; i_scalar++) {
    chemical_composition.append(gasComposition[i_scalar] + ":" + to_string(val_scalars[i_scalar])+", ");
    val_scalars_sum += val_scalars[i_scalar];
  }
  chemical_composition.append(gasComposition[n_species_mixture - 1] + ":" + to_string(1.0 - val_scalars_sum));

  return chemical_composition;
}

void CFluidCantera::SetTDState_T(const su2double val_temperature, const su2double* val_scalars) {
  auto sol = newSolution(Chemical_MechanismFile, Phase_Name, Transport_Model);
  auto gas = sol->thermo();
  DictionaryChemicalComposition(val_scalars);
  Temperature = val_temperature;
  gas->setState_TPY(GetValue(Temperature), GetValue(Pressure_Thermodynamic), chemical_composition);
  Density = gas->density();
  Cp = gas->cp_mass();
  Cv = gas->cv_mass();
  Mu = sol->transport()->viscosity();
  Kt = sol->transport()->thermalConductivity();

  ComputeMassDiffusivity();
}
#endif