/*!
 * fluid_model_pig.cpp
 * \brief Source of the ideal gas model.
 * \author S. Vitale, G. Gori, M. Pini, A. Guardone, P. Colonna
 * \version 5.0.0 "Raven"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
 *                      Dr. Thomas D. Economon (economon@stanford.edu).
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

CIdealGas::CIdealGas() : CFluidModel() {

  Gamma = 0.0;
  Gamma_Minus_One = 0.0;
  Gas_Constant = 0.0;
  Cp = 0.0;
}


CIdealGas::CIdealGas(su2double gamma, su2double R ) : CFluidModel() {
  Gamma = gamma;
  Gamma_Minus_One = Gamma - 1.0;
  Gas_Constant = R;
  Cp = Gamma/Gamma_Minus_One*Gas_Constant;
}


CIdealGas::~CIdealGas(void) {

}


void CIdealGas::SetTDState_rhoe (su2double rho, su2double e ) {
  
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



void CIdealGas::SetTDState_PT (su2double P, su2double T ) {
  su2double e = T*Gas_Constant/(Gamma_Minus_One);
  su2double rho = P/(T*Gas_Constant);
  SetTDState_rhoe(rho, e);
}




void CIdealGas::SetTDState_Prho (su2double P, su2double rho ) {
  su2double e = P/(Gamma_Minus_One*rho);
  SetTDState_rhoe(rho, e);

}



void CIdealGas::SetEnergy_Prho (su2double P, su2double rho ) {
  StaticEnergy = P/(rho*Gamma_Minus_One);

}



void CIdealGas::SetTDState_hs (su2double h, su2double s ) {

  su2double T = h*Gamma_Minus_One/Gas_Constant/Gamma;
  su2double e = h/Gamma;
  su2double v = exp(-1/Gamma_Minus_One*log(T) + s/Gas_Constant);

  SetTDState_rhoe(1/v, e);

}




void CIdealGas::SetTDState_Ps (su2double P, su2double s ) {

  su2double T   = exp(Gamma_Minus_One/Gamma* (s/Gas_Constant +log(P) -log(Gas_Constant)) );
  su2double rho = P/(T*Gas_Constant);

  SetTDState_Prho(P, rho);

}



void CIdealGas::SetTDState_rhoT (su2double rho, su2double T ) {

  su2double e = T*Gas_Constant/Gamma_Minus_One;
  SetTDState_rhoe(rho, e);

}


void CIdealGas_Generic::SetTDState_rhoT (su2double rho, su2double T ) {

  HeatCapacity->Set_Cv0 (T);
  Cv = HeatCapacity->Get_Cv0 ();

  su2double e = T*Cv;
  SetTDState_rhoe(rho, e);


}


void CIdealGas_Generic::SetGamma_Trho (su2double Temperature, su2double Density) {

  Cp = Cv + Gas_Constant;
  Gamma = Cp/Cv;

}




CIdealGas_Generic::CIdealGas_Generic() : CIdealGas() {

  Gamma = 0.0;
  Gamma_Minus_One = 0.0;
  Gas_Constant = 0.0;
  Cp = 0.0;
}


CIdealGas_Generic::CIdealGas_Generic(su2double gamma, su2double R ) : CIdealGas(gamma, R ) {
  Gamma = gamma;
  Gamma_Minus_One = Gamma - 1.0;
  Gas_Constant = R;
  Cp = Gamma/Gamma_Minus_One*Gas_Constant;
}


CIdealGas_Generic::~CIdealGas_Generic(void) {

}



void CIdealGas_Generic::SetTDState_rhoe (su2double rho, su2double e ) {

  su2double Temperature_new, error;
  su2double toll = 1e-3;
  unsigned int ITMAX = 1000, count_T = 0;

  Density = rho;
  StaticEnergy = e;

  Temperature = Gamma_Minus_One/Gas_Constant*StaticEnergy;

  do{
  	Temperature_new = Temperature;
  	HeatCapacity->Set_Cv0 (Temperature_new);
  	Cv = HeatCapacity->Get_Cv0 ();
  	Temperature = Temperature_new * 0.5 + 0.5 * StaticEnergy/ Cv;
  	error = abs(Temperature_new-Temperature)/Temperature_new;

  	count_T++;

  	}
  while(error > toll && count_T<ITMAX);

  if (count_T==ITMAX) {
//	 cout <<"Too many iterations for T in e-rho call" << endl;
  }

  SetGamma_Trho (Temperature, rho);
  Gamma_Minus_One = Gamma - 1;

  Pressure = Gas_Constant*Density*Temperature;
  SoundSpeed2 = Gamma*Pressure/Density;
  Entropy = Cv*log(Temperature) + log(1.0/Density)*Gas_Constant;
  dPdrho_e = Gamma_Minus_One*StaticEnergy;
  dPde_rho = Gamma_Minus_One*Density;
  dTdrho_e = 0.0;
  dTde_rho = Gamma_Minus_One/Gas_Constant;

}



void CIdealGas_Generic::SetTDState_PT (su2double P, su2double T ) {
  HeatCapacity->Set_Cv0 (T);
  Cv = HeatCapacity->Get_Cv0 ();
  su2double e = T*Cv;
  su2double rho = P/(T*Gas_Constant);
  SetTDState_rhoe(rho, e);

}

void CIdealGas_Generic::SetTDState_Prho (su2double P, su2double rho ) {
  su2double T = P/(Gas_Constant*rho);
  SetTDState_PT(P, T);

}


void CIdealGas_Generic::SetEnergy_Prho (su2double P, su2double rho ) {
	su2double T = P/(Gas_Constant*rho);
	HeatCapacity->Set_Cv0 (T);
	Cv = HeatCapacity->Get_Cv0 ();
	StaticEnergy = T*Cv;

}

void CIdealGas_Generic::SetTDState_hs (su2double h, su2double s ) {

	  su2double error, T;
	  su2double toll = 1e-3;
	  unsigned int ITMAX = 1000, count_T = 0;

	  su2double T_new = h*Gamma_Minus_One/Gas_Constant/Gamma;

	  do{
	  	T = T_new;
	  	HeatCapacity->Set_Cv0 (T);
	  	Cv = HeatCapacity->Get_Cv0 ();
	  	Cp = Cv + Gas_Constant;
	  	T_new = h/ Cp;
	  	error = abs(T_new-T)/T_new;
	  	count_T++;
	  	}

	  while(error > toll && count_T<ITMAX);


      su2double e = T*Cv;
      su2double  v = exp((s - Cv*log(T)) / Gas_Constant) ;

      SetTDState_rhoe(1/v, e);

}

void CIdealGas_Generic::SetTDState_Ps (su2double P, su2double s ) {

	su2double T_new   = exp(Gamma_Minus_One/Gamma* (s/Gas_Constant +log(P) -log(Gas_Constant)) );
	su2double T, v, error, toll = 1e-3;
	unsigned int ITMAX = 1000, count_T = 0;

	do{
	   	  T = T_new;
	   	  HeatCapacity->Set_Cv0 (T);
		  Cv = HeatCapacity->Get_Cv0 ();

		  v =  exp((s - Cv*log(T))/Gas_Constant);
		  T_new = P*v / Gas_Constant;

		  error = abs(T - T_new)/T;
		  count_T++;
	}while(error >toll && count_T < ITMAX);

    SetTDState_Prho(P, 1/v);
}
