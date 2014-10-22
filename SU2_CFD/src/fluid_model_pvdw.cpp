/*!
 * fluid_model_pvdw.cpp
 * \brief Source of the Polytropic Van der Waals model.
 * \author: S.Vitale, G.Gori, M.Pini, A.Guardone, P.Colonna
 * \version 3.2.3 "eagle"
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

CVanDerWaalsGas::CVanDerWaalsGas() : CIdealGas() {
	a = 0.0;
	b = 0.0;
}


CVanDerWaalsGas::CVanDerWaalsGas(double gamma, double R, double Pstar, double Tstar): CIdealGas(gamma, R) {
    a = 27.0/64.0*Gas_Constant*Gas_Constant*Tstar*Tstar/Pstar;
	b = 1.0/8.0*Gas_Constant*Tstar/Pstar;


}


CVanDerWaalsGas::~CVanDerWaalsGas(void) { }


void CVanDerWaalsGas::SetTDState_rhoe (double rho, double e ) {
	Density = rho;
	StaticEnergy = e;

	Pressure = Gamma_Minus_One*Density/(1.0-Density*b)*(StaticEnergy + Density*a) - a*Density*Density;
    Temperature = (Pressure+Density*Density*a)*((1-Density*b)/(Density*Gas_Constant));
    Entropy = Gas_Constant *( log(Temperature)/Gamma_Minus_One + log(1/Density - b));


    dPde_rho = Density*Gamma_Minus_One/(1.0 - Density*b);
    dPdrho_e = Gamma_Minus_One/(1.0 - Density*b)*((StaticEnergy + 2*Density*a) + Density*b*(StaticEnergy + Density*a)/(1.0 - Density*b)) - 2*Density*a;
    dTdrho_e = Gamma_Minus_One/Gas_Constant*a;
    dTde_rho = Gamma_Minus_One/Gas_Constant;

    SoundSpeed2 = dPdrho_e + Pressure/(Density*Density)*dPde_rho;

}


void CVanDerWaalsGas::SetTDState_PT (double P, double T ) {
	double toll= 1e-9;
	double A, B, Z, DZ, F, F1;
	A= a*P/(T*Gas_Constant)/(T*Gas_Constant);
	B= b*P/(T*Gas_Constant);

    Z= max(B, 1.01);
	DZ= 1.0;
	do{
		F = Z*Z*Z - Z*Z*(B+1.0) + Z*A - A*B;
		F1 = 3*Z*Z - 2*Z*(B+1.0) + A;
		DZ = F/F1;
		Z-= DZ;
	}while(DZ>toll);
	Density = P/(Z*Gas_Constant*T);

    double e = T*Gas_Constant/Gamma_Minus_One - a*Density;
	SetTDState_rhoe(Density, e);

}

void CVanDerWaalsGas::SetTDState_Prho (double P, double rho ) {

	SetEnergy_Prho(P, rho);

	SetTDState_rhoe(rho, StaticEnergy);

}

void CVanDerWaalsGas::SetTDState_hs (double h, double s ){

    double v, T, dv, f, f1;
    double toll = 1e-9;

    T = h*Gamma_Minus_One/Gas_Constant/Gamma;
    v = exp(-1/Gamma_Minus_One*log(T) + s/Gas_Constant);
	do{
		f=  log(v-b) - s/Gas_Constant + log(T)/Gamma_Minus_One;
		f1= 1/(v-b);
		dv= f/f1;
		v-= dv;
		T= (h+ 2*a/v)/Gas_Constant/(1/Gamma_Minus_One+ v/(v-b));
	}while(abs(dv) > toll);

	Density = 1/v;
	Temperature = T;
    Pressure = Gas_Constant*Temperature*Density / (1 - Density*b) - a*Density*Density;

    SetTDState_Prho(Pressure, Density);

}

void CVanDerWaalsGas::SetEnergy_Prho (double P, double rho ) {
	double T = (P+rho*rho*a)*(1-rho*b)/(rho*Gas_Constant);
	StaticEnergy = T*Gas_Constant/Gamma_Minus_One - rho*a;

}

void CVanDerWaalsGas::SetTDState_rhoT (double rho, double T) {

	double e = T*Gas_Constant/Gamma_Minus_One - a*rho;
	SetTDState_rhoe(rho, e);

}






