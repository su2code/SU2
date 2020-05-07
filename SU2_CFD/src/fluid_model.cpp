/*!
 * fluid_model.cpp
 * \brief Source of the main thermo-physical subroutines of the SU2 solvers.
 * \author S.Vitale, M.Pini, G.Gori, A.Guardone, P.Colonna
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


#include "../include/fluid_model.hpp"

CFluidModel::CFluidModel(void) {

  /*--- Attributes initialization ---*/

  StaticEnergy = 0.0;
  Entropy = 0.0;
  Density = 0.0;
  Pressure = 0.0;
  SoundSpeed2 = 0.0;
  Temperature = 0.0;
  dPdrho_e = 0.0;
  dPde_rho = 0.0;
  dTdrho_e = 0.0;
  dTde_rho = 0.0;
  Cp       = 0.0;
  Cv       = 0.0;
  Mu       = 0.0;
  Mu_Turb  = 0.0;

  LaminarViscosity = NULL;
  ThermalConductivity = NULL;

}

CFluidModel::~CFluidModel(void) {
  if (LaminarViscosity!= NULL) delete LaminarViscosity;
  if (ThermalConductivity!= NULL) delete ThermalConductivity;
}

void CFluidModel::SetLaminarViscosityModel (CConfig *config) {
  
  switch (config->GetKind_ViscosityModel()) {
    case CONSTANT_VISCOSITY:
      LaminarViscosity = new CConstantViscosity(config->GetMu_ConstantND());
      break;
    case SUTHERLAND:
      LaminarViscosity = new CSutherland(config->GetMu_RefND(), config->GetMu_Temperature_RefND(), config->GetMu_SND());
      break;
    case POLYNOMIAL_VISCOSITY:
      LaminarViscosity = new CPolynomialViscosity(config->GetnPolyCoeffs(), config->GetMu_PolyCoeffND());
      break;
    default:
      SU2_MPI::Error("Viscosity model not available.", CURRENT_FUNCTION);
      break;
  }
  
}

void CFluidModel::SetThermalConductivityModel (CConfig *config) {
  
  switch (config->GetKind_ConductivityModel()) {
    case CONSTANT_CONDUCTIVITY:
      if (config->GetKind_ConductivityModel_Turb() == CONSTANT_PRANDTL_TURB) {
        ThermalConductivity = new CConstantConductivityRANS(config->GetKt_ConstantND(), config->GetPrandtl_Turb());
      } else {
        ThermalConductivity = new CConstantConductivity(config->GetKt_ConstantND());
      }
      break;
    case CONSTANT_PRANDTL:
      if (config->GetKind_ConductivityModel_Turb() == CONSTANT_PRANDTL_TURB) {
        ThermalConductivity = new CConstantPrandtlRANS(config->GetPrandtl_Lam(), config->GetPrandtl_Turb());
      } else {
        ThermalConductivity = new CConstantPrandtl(config->GetPrandtl_Lam());
      }
      break;
    case POLYNOMIAL_CONDUCTIVITY:
      if (config->GetKind_ConductivityModel_Turb() == CONSTANT_PRANDTL_TURB) {
        ThermalConductivity = new CPolynomialConductivityRANS(config->GetnPolyCoeffs(), config->GetKt_PolyCoeffND(), config->GetPrandtl_Turb());
      } else {
        ThermalConductivity = new CPolynomialConductivity(config->GetnPolyCoeffs(), config->GetKt_PolyCoeffND());
      }
      break;
    default:
      SU2_MPI::Error("Conductivity model not available.", CURRENT_FUNCTION);
      break;
  }
  
}

