/*!
 * \file CIdealGas.cpp
 * \brief Source of the ideal gas model.
 * \author S. Vitale, G. Gori, M. Pini, A. Guardone, P. Colonna
 * \version 8.0.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2023, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../include/fluid/CIdealGas.hpp"

CIdealGas::CIdealGas(su2double gamma, su2double R, bool CompEntropy) : CFluidModel() {
  Gamma = gamma;
  Gamma_Minus_One = Gamma - 1.0;
  Gas_Constant = R;
  Cp = Gamma / Gamma_Minus_One * Gas_Constant;
  Cv = Cp - R;

  ComputeEntropy = CompEntropy;
}

void CIdealGas::SetTDState_rhoe(su2double rho, su2double e) {
  Density = rho;
  StaticEnergy = e;
  Pressure = Gamma_Minus_One * Density * StaticEnergy;
  Temperature = Gamma_Minus_One * StaticEnergy / Gas_Constant;
  SoundSpeed2 = Gamma * Pressure / Density;
  dPdrho_e = Gamma_Minus_One * StaticEnergy;
  dPde_rho = Gamma_Minus_One * Density;
  dTdrho_e = 0.0;
  dTde_rho = Gamma_Minus_One / Gas_Constant;

  if (ComputeEntropy) Entropy = (1.0 / Gamma_Minus_One * log(Temperature) + log(1.0 / Density)) * Gas_Constant;
}

void CIdealGas::SetTDState_PT(su2double P, su2double T) {
  su2double e = T * Gas_Constant / Gamma_Minus_One;
  su2double rho = P / (T * Gas_Constant);
  SetTDState_rhoe(rho, e);
}

void CIdealGas::SetTDState_Prho(su2double P, su2double rho) {
  su2double e = P / (Gamma_Minus_One * rho);
  SetTDState_rhoe(rho, e);
}

void CIdealGas::SetEnergy_Prho(su2double P, su2double rho) { StaticEnergy = P / (rho * Gamma_Minus_One); }

void CIdealGas::SetTDState_hs(su2double h, su2double s) {
  su2double T = h * Gamma_Minus_One / Gas_Constant / Gamma;
  su2double e = h / Gamma;
  su2double v = exp(-1 / Gamma_Minus_One * log(T) + s / Gas_Constant);

  SetTDState_rhoe(1 / v, e);
}

void CIdealGas::SetTDState_Ps(su2double P, su2double s) {
  su2double T = exp(Gamma_Minus_One / Gamma * (s / Gas_Constant + log(P) - log(Gas_Constant)));
  su2double rho = P / (T * Gas_Constant);

  SetTDState_Prho(P, rho);
}

void CIdealGas::SetTDState_rhoT(su2double rho, su2double T) {
  su2double e = T * Gas_Constant / Gamma_Minus_One;
  SetTDState_rhoe(rho, e);
}

void CIdealGas::ComputeDerivativeNRBC_Prho(su2double P, su2double rho) {
  su2double dPdT_rho, dPdrho_T, dPds_rho;

  SetTDState_Prho(P, rho);

  dPdT_rho = Gas_Constant * rho;
  dPdrho_T = Gas_Constant * Temperature;

  dhdrho_P = -dPdrho_e / dPde_rho - P / rho / rho;
  dhdP_rho = 1.0 / dPde_rho + 1.0 / rho;
  dPds_rho = rho * rho * (SoundSpeed2 - dPdrho_T) / dPdT_rho;
  dsdP_rho = 1.0 / dPds_rho;
  dsdrho_P = -SoundSpeed2 / dPds_rho;
}
