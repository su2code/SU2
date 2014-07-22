/*!
 * fluid_model.cpp
 * \brief Source of the main thermo-physical subroutines of the SU2 solvers.
 * \author S.Vitale, M.Pini, G.Gori, A.Guardone, P.Colonna
 * \version 3.2.0 "eagle"
 *
 * SU2, Copyright (C) 2012-2014 Aerospace Design Laboratory (ADL).
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

	DynamicViscosity = NULL;
	ThermalConductivity = NULL;

}

CFluidModel::~CFluidModel(void) {

  }

void CFluidModel::SetViscosityModel (CConfig *config){
	switch (config->GetKind_ViscosityModel()) {

	case CONSTANT_VISCOSITY:
		DynamicViscosity = new CConstantViscosity(config->GetMu_ConstantND());
		break;
	case SUTHERLAND:
		DynamicViscosity = new CSutherland(config->GetMu_RefND(), config->GetMu_Temperature_RefND(), config->GetMu_SND());
		break;

	}
}

//void CFluidModel::SetThermalConductivityModel (CConfig *config){
//	switch (config->GetKind_ThermalConductivityModel()) {
//
//	case CONSTANT_THERMALCONDUCTIVITY:
//		ThermalConductivity = new CConstantThermalConductivity(double kt_const);
//		break;
//	case CONSTANT_PRANDTL:
//		ThermalConductivity = new CConstantPrandtl();
//		break;
//
//	}
//}

