/*!
 * transport_model.cpp
 * \brief Source of the main transport properties subroutines of the SU2 solvers.
 * \author S. Vitale, M. Pini, G. Gori, A. Guardone, P. Colonna
 * \version 3.2.8.1 "eagle"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (fpalacios@stanford.edu).
 *                      Dr. Thomas D. Economon (economon@stanford.edu).
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
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
	dmudrho_T = 0.0;
	dmudT_rho = 0.0;

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

}

void CSutherland::SetDerViscosity(double T, double rho) {

	dmudrho_T = 0.0;
	dmudT_rho = Mu_ref*( (3.0/2.0)*pow( (T/T_ref),(1.0/2.0) )*( (T_ref + S)/(T + S) )
			    -pow( (T/T_ref),(3.0/2.0) )*(T_ref + S)/(T + S)/(T + S) );

}

/* ------------------------------------------------- */
/* ---------- Thermal Conductivity Models ---------- */
/* ------------------------------------------------- */

CConductivityModel::CConductivityModel(void) {

  /*--- Attributes initialization ---*/

	Kt = 0.0;
	dktdrho_T = 0.0;
	dktdT_rho = 0.0;

}

CConductivityModel::~CConductivityModel(void) { }


CConstantConductivity::CConstantConductivity(void) : CConductivityModel() { }

CConstantConductivity::CConstantConductivity(double kt_const) : CConductivityModel() {

  /*--- Attributes initialization ---*/

	Kt = kt_const;
	dktdrho_T = 0.0;
	dktdT_rho = 0.0;

}

CConstantConductivity::~CConstantConductivity(void) { }


CConstantPrandtl::CConstantPrandtl(void) : CConductivityModel() { }

CConstantPrandtl::CConstantPrandtl(double pr_const) : CConductivityModel() {

  /*--- Attributes initialization ---*/

	Pr_const = pr_const;

}

void CConstantPrandtl::SetConductivity(double T, double rho, double mu, double cp) {

	Kt = mu*cp/Pr_const;

}

void CConstantPrandtl::SetDerConductivity(double T, double rho, double dmudrho_T, double dmudT_rho, double cp) {

	dktdrho_T = dmudrho_T*cp/Pr_const;
	dktdT_rho = dmudT_rho*cp/Pr_const;

}

CConstantPrandtl::~CConstantPrandtl(void) { }

