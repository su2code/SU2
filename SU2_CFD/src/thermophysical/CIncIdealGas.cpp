/*!
 * CIncIdealGas.cpp
 * \brief Source of the incompressible Ideal Gas model.
 * \author T. Economon
 * \version 7.0.4 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2020, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../include/fluid_model.hpp"

CIncIdealGas::CIncIdealGas() : CFluidModel() {
  Pressure = 0.0;
  Gamma = 0.0;
  Gas_Constant = 0.0;
  Cp = 0.0;
  Cv = 0.0;
}

CIncIdealGas::CIncIdealGas(su2double val_Cp, su2double val_gas_constant, su2double val_operating_pressure)
    : CFluidModel() {
  /*--- In the incompressible ideal gas model, the thermodynamic pressure
    is decoupled from the governing equations and held constant. The
    density is therefore only a function of temperature variations. ---*/

  Gas_Constant = val_gas_constant;
  Pressure = val_operating_pressure;
  Gamma = 1.0;
  Cp = val_Cp;
  Cv = Cp;
}

CIncIdealGas::~CIncIdealGas(void) {}

void CIncIdealGas::SetTDState_T(su2double val_temperature) {
  /*--- The EoS only depends upon temperature. ---*/

  Temperature = val_temperature;
  Density = Pressure / (Temperature * Gas_Constant);
}
