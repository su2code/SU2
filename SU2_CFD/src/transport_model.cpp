/*!
 * transport_model.cpp
 * \brief Source of the main transport properties subroutines of the SU2 solvers.
 * \author: S.Vitale, M.Pini, G.Gori, A.Guardone, P.Colonna
 * \version 3.2.4 "eagle"
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

#include "../include/transport_model.hpp"


/* ------------------------------------------------- */
/* ----------- Dynamic Viscosity Models ------------ */
/* ------------------------------------------------- */

CViscosityModel::CViscosityModel(void) {

  /*--- Attributes initialization ---*/

	Mu = 0.0;
	dmudrho_T = 0.0;
	dmudT_rho = 0.0;

}

CViscosityModel::~CViscosityModel(void) { }


CConstantViscosity::CConstantViscosity(void) : CViscosityModel() { }

CConstantViscosity::CConstantViscosity(double mu_const) : CViscosityModel() {

  /*--- Attributes initialization ---*/

	Mu = mu_const;

}

CConstantViscosity::~CConstantViscosity(void) { }




CSutherland::CSutherland(void) : CViscosityModel() {
	Mu_ref = 0.0;
	T_ref = 0.0;
	S = 0.0;

}

CSutherland::CSutherland(double mu_ref, double t_ref, double s) : CViscosityModel() {

	Mu_ref = mu_ref;
	T_ref = t_ref;
	S = s;
}

CSutherland::~CSutherland(void) { }


void CSutherland::SetViscosity(double T, double rho) {

	Mu = Mu_ref*pow((T/T_ref),(3.0/2.0))*((T_ref + S)/(T + S));
	dmudrho_T = 0.0;
	dmudT_rho = 0.0;

}


/* ------------------------------------------------- */
/* ---------- Thermal Conductivity Models ---------- */
/* ------------------------------------------------- */

CThermalConductivityModel::CThermalConductivityModel(void) {

  /*--- Attributes initialization ---*/

	Kt = 0.0;
	dktdrho_T = 0.0;
	dktdT_rho = 0.0;

}

CThermalConductivityModel::~CThermalConductivityModel(void) { }


CConstantThermalConductivity::CConstantThermalConductivity(void) : CThermalConductivityModel() { }

CConstantThermalConductivity::CConstantThermalConductivity(double kt_const) : CThermalConductivityModel() {

  /*--- Attributes initialization ---*/

	Kt = kt_const;

}

CConstantThermalConductivity::~CConstantThermalConductivity(void) { }


CConstantPrandtl::CConstantPrandtl(void) : CThermalConductivityModel() { }

CConstantPrandtl::CConstantPrandtl(double pr_const) : CThermalConductivityModel() {

  /*--- Attributes initialization ---*/

	Pr_const = pr_const;

}

void CConstantPrandtl::SetThermalConductivity(double par1, double par2) {

	double Cp = par1;
	double Mu = par2;

	Kt = Mu*Cp/Pr_const;
	dktdrho_T = 0.0;
	dktdT_rho = 0.0;

}

CConstantPrandtl::~CConstantPrandtl(void) { }

