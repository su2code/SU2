/*!
 * fluid_model_inc.cpp
 * \brief Source of the incompressible fluid models.
 * \author T. Economon
 * \version 5.0.0 "Raven"
 *
 * SU2 Original Developers: Dr. Francisco D. Palacios.
 *                          Dr. Thomas D. Economon.
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *                 Prof. Edwin van der Weide's group at the University of Twente.
 *                 Prof. Vincent Terrapon's group at the University of Liege.
 *
 * Copyright (C) 2012-2017 SU2, the open-source CFD code.
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

#include "../include/fluid_model.hpp"

CConstantDensity::CConstantDensity() : CFluidModel() {
  Density = 0.0;
  Cp      = 0.0;
  Cv      = 0.0;
}

CConstantDensity::CConstantDensity(su2double val_Density, su2double val_Cp, su2double val_Cv) : CFluidModel() {
  Density = val_Density;
  Cp      = val_Cp;
  Cv      = val_Cv;
}

CConstantDensity::~CConstantDensity(void) { }

void CConstantDensity::SetTDState_T (su2double val_Temperature) {
  
  /*--- Density is constant and thermodynamic pressure is
    not required for incompressible, constant density flows,
    but the energy equation can still be computed as a 
    decoupled equation. Hence, we update the value. ---*/

  Temperature = val_Temperature;

}

CIncIdealGas::CIncIdealGas() : CFluidModel() {
  Pressure        = 0.0;
  Gamma           = 0.0;
  Gamma_Minus_One = 0.0;
  Gas_Constant    = 0.0;
  Cp              = 0.0;
  Cv              = 0.0;
}

CIncIdealGas::CIncIdealGas(su2double val_Gamma, su2double val_Gas_Constant, su2double val_Pressure) : CFluidModel() {

  /*--- In the incompressible ideal gas model, the thermodynamic pressure
    is decoupled from the governing equations and held constant. The 
    density is therefore only a function of temperature variations. ---*/

  Gamma        = val_Gamma;
  Gas_Constant = val_Gas_Constant;
  Pressure     = val_Pressure;

  /*--- Derived gas quantities from input. ---*/

  Gamma_Minus_One = Gamma - 1.0;
  Cp              = Gas_Constant*Gamma/Gamma_Minus_One;
  Cv              = Cp/Gamma;

}

CIncIdealGas::~CIncIdealGas(void) { }

void CIncIdealGas::SetTDState_T(su2double val_Temperature) {

 /*--- The EoS only depends upon temperature. ---*/

  Temperature  = val_Temperature;
  Density      = Pressure/(Temperature*Gas_Constant);
  StaticEnergy = Temperature*Gas_Constant/Gamma_Minus_One;
  SoundSpeed2  = Gamma*Pressure/Density;
  Entropy      = (1.0/Gamma_Minus_One*log(Temperature) + log(1.0/Density))*Gas_Constant;
  dPdrho_e     = Gamma_Minus_One*StaticEnergy;
  dPde_rho     = Gamma_Minus_One*Density;
  dTdrho_e     = 0.0;
  dTde_rho     = Gamma_Minus_One/Gas_Constant;

}
