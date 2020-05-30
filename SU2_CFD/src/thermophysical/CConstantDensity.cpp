/*!
 * CConstantDensity.cpp
 * \brief Source of the incompressible constant density model.
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

CConstantDensity::CConstantDensity() : CFluidModel() {
  Density = 0.0;
  Cp = 0.0;
  Cv = 0.0;
}

CConstantDensity::CConstantDensity(su2double val_Density, su2double val_Cp) : CFluidModel() {
  Density = val_Density;
  Cp = val_Cp;
  Cv = val_Cp;
}

CConstantDensity::~CConstantDensity(void) {}

void CConstantDensity::SetTDState_T(su2double val_Temperature) {
  /*--- Density is constant and thermodynamic pressure is
   not required for incompressible, constant density flows,
   but the energy equation can still be computed as a
   decoupled equation. Hence, we update the value.
   Note Cp = Cv (gamma = 1). ---*/

  Temperature = val_Temperature;
}
