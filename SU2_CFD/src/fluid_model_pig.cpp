/*!
 * fluid_model_pig.cpp
 * \brief Source of the ideal gas model.
 * \author S. Vitale, G. Gori, M. Pini, A. Guardone, P. Colonna
 * \version 4.0.0 "Cardinal"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
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

#include "../include/fluid_model.hpp"

CIdealGas::CIdealGas() : CFluidModel() {

	Gamma = 0.0;
	Gamma_Minus_One = 0.0;
	Gas_Constant = 0.0;
	Cp = 0.0;
}


CIdealGas::CIdealGas(double gamma, double R ) : CFluidModel() {
	Gamma = gamma;
	Gamma_Minus_One = Gamma - 1.0;
	Gas_Constant = R;
	Cp = Gamma/Gamma_Minus_One*Gas_Constant;
}


CIdealGas::~CIdealGas(void) {

}

void CIdealGas::SetTDState_rhoe (double rho, double e ) {
  
	Density = rho;
	StaticEnergy = e;
	Pressure = Gamma_Minus_One*Density*StaticEnergy;
	Temperature = Gamma_Minus_One*StaticEnergy/Gas_Constant;
	SoundSpeed2 = Gamma*Pressure/Density;
	Entropy = (1.0/Gamma_Minus_One*log(Temperature) + log(1.0/Density))*Gas_Constant;
	dPdrho_e = Gamma_Minus_One*StaticEnergy;
	dPde_rho = Gamma_Minus_One*Density;
	dTdrho_e = 0.0;
	dTde_rho = Gamma_Minus_One/Gas_Constant;

}

void CIdealGas::SetTDState_PT (double P, double T ) {
	double e = T*Gas_Constant/Gamma_Minus_One;
	double rho = P/(T*Gas_Constant);
	SetTDState_rhoe(rho, e);

}

void CIdealGas::SetTDState_Prho (double P, double rho ) {
	double e = P/(Gamma_Minus_One*rho);
	SetTDState_rhoe(rho, e);

}

void CIdealGas::SetEnergy_Prho (double P, double rho ) {
	StaticEnergy = P/(rho*Gamma_Minus_One);

}

void CIdealGas::SetTDState_hs (double h, double s ) {

	double T = h*Gamma_Minus_One/Gas_Constant/Gamma;
	double e = h/Gamma;
	double v = exp(-1/Gamma_Minus_One*log(T) + s/Gas_Constant);

	SetTDState_rhoe(1/v, e);

}

void CIdealGas::SetTDState_rhoT (double rho, double T ) {

	double e = T*Gas_Constant/Gamma_Minus_One;
	SetTDState_rhoe(rho, e);

}









