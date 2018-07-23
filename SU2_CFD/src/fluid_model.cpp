/*!
 * fluid_model.cpp
 * \brief Source of the main thermo-physical subroutines of the SU2 solvers.
 * \author S.Vitale, M.Pini, G.Gori, A.Guardone, P.Colonna
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
  case LUT_VISCOSITY:
	LaminarViscosity = new CLookUpTable_Viscosity(config, config->GetBoolDimensionalLUTViscosity());
	break;

    }
}

void CFluidModel::SetThermalConductivityModel (CConfig *config) {
  
  switch (config->GetKind_ConductivityModel()) {
  case CONSTANT_CONDUCTIVITY:
    ThermalConductivity = new CConstantConductivity(config->GetKt_ConstantND());
    break;
  case CONSTANT_PRANDTL:
    ThermalConductivity = new CConstantPrandtl(config->GetPrandtl_Lam());
    break;
  case LUT_CONDUCTIVITY:
		ThermalConductivity = new CLookUpTable_Conductivity(config);
  break;

  }
  
}

