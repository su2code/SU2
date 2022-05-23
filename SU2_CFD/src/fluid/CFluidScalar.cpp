/*!
 * \file CFluidScalar.hpp
 * \brief Defines the incompressible Ideal Gas model.
 * \author T. Economon, Mark Heimgartner, Cristopher Morales Ubal 
 * \version 7.3.1 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2022, SU2 Contributors (cf. AUTHORS.md)
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

#include <cmath>
#include <math.h>
#include <vector>
#include <numeric>
#include <iostream>

#include "../../include/fluid/CFluidScalar.hpp"

CFluidScalar::CFluidScalar(su2double val_Cp, su2double val_gas_constant, const su2double value_pressure_operating, CConfig *config) : CFluidModel() {
  
  n_scalars = config->GetnSpecies();
  n_species_mixture = n_scalars + 1;

  specificHeat.resize(n_species_mixture);
  molarMasses.resize(n_species_mixture);
  massFractions.resize(n_species_mixture);
  moleFractions.resize(n_species_mixture);
  laminarViscosity.resize(n_species_mixture);
  laminarthermalConductivity.resize(n_species_mixture);

  for (int iVar = 0; iVar < n_species_mixture; iVar++) {
    molarMasses[iVar] = config->GetMolecular_Weight(iVar);
    specificHeat[iVar] = 1004.703; 
  }

  wilke = false;
  davidson = true;

  Pressure_Thermodynamic = value_pressure_operating;
  Gas_Constant = val_gas_constant; 
  Gamma = 1.0;
  Cp = val_Cp;
  Cv = Cp;
  SetLaminarViscosityModel(config);
  SetThermalConductivityModel(config);
}

void CFluidScalar::SetTDState_T(su2double val_temperature, const su2double* val_scalars){
  const su2double MeanMolecularWeight = ComputeMeanMolecularWeight(molarMasses, val_scalars);
  Temperature = val_temperature;
  Density = Pressure_Thermodynamic / ((Temperature * UNIVERSAL_GAS_CONSTANT*MeanMolecularWeight)); 
}
