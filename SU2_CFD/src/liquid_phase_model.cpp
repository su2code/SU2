/*! * liquid_phase_model.cpp
 * \brief Source of the main transport properties subroutines of the SU2 solvers.
 * \author S. Vitale, M. Pini, G. Gori, A. Guardone, P. Colonna
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

#include "../include/liquid_phase_model.hpp"


CLiquidModel::CLiquidModel() {
}

CLiquidModel::CLiquidModel(CConfig *config) {

  /*--- Attributes initialization ---*/

  rho_l  = 0.0;
  h_l    = 0.0;
  dGibbs = 0.0;
  Tsat   = 0.0;
  Psat   = 0.0;
  sigma  = 0.0;


  Gas_Constant = config->GetGas_Constant();
  Tstar = config->GetTemperature_Critical();

}

CLiquidModel::~CLiquidModel(void) { }

void CLiquidModel::Set_LiquidProp(su2double P, su2double T, su2double rho, su2double h_v, su2double Rcritical, su2double Rdroplet, su2double mom3) {

	// guess for R critical, a loop is required to evaluate the right properties
	//Rc = 1e-12;

	//SetRadius(Two_Phase_Var);

	SetTsat(P);
	SetPsat(T);

	if (Tsat > T) {

		SetSurfaceTension(T, Rdroplet);

		SetTLiquid(T, Rcritical, Rdroplet);
		SetLiquidDensity();

		SetLiquidEnthalpy(h_v);

		SetRCritical(P, T);

		//Rc = min(Rc, Rdroplet);

		SetDensity_Mixture(rho, mom3);
	} else {

		sigma = 0;
		T_l = T;
		rho_l = rho;
		Rc = 0;
		h_l = h_v;
		rho_m = rho;
		Rdroplet = 0;

	}

}

void CLiquidModel::SetRCritical(su2double P, su2double T) {

    if (Psat < P) {
		dGibbs = T * Gas_Constant * log(P/Psat);
		Rc = 2.0*sigma / (rho_l * dGibbs);
    }

/*    if (Tsat > T) {
		dGibbs = P * (Tsat - T) / Tsat;
		Rc = 2.0*sigma / (rho_l * dGibbs);
    }
*/
    else { Rc = 0.0;}

}

void CLiquidModel::SetDensity_Mixture(su2double rho, su2double mom3) {

	    if (mom3 == 0) {
	    	rho_m = rho;
	    } else {
			y = mom3*(rho_l - rho);
			y = y + 0.75 * rho / 3.1415;
			y = mom3*rho_l / y;

			rho_m = y/ rho_l + (1.0 - y)/ rho;
			rho_m = 1.0/ rho_m;
	    }

}



CWater::CWater(CConfig *config) : CLiquidModel(config) {

	coeff_saturation  = new su2double [10];
	coeff_latent_heat = new su2double [4];

    coeff_saturation[0]  =  0.11670521452767e4;
    coeff_saturation[1]  = -0.72421316703206e6;
    coeff_saturation[2]  = -0.17073846940092e2;
    coeff_saturation[3]  =  0.12020824702470e5;
    coeff_saturation[4]  = -0.32325550322333e7;
    coeff_saturation[5]  =  0.14915108613530e2;
    coeff_saturation[6]  = -0.48232657361591e4;
    coeff_saturation[7]  =  0.40511340542057e6;
    coeff_saturation[8]  = -0.23855557567849e0;
    coeff_saturation[9]  =  0.65017534844798e3;

    coeff_latent_heat[0] =  3.8788e6;
    coeff_latent_heat[1] = -5.9196e6;
    coeff_latent_heat[2] =  8.8253e6;
    coeff_latent_heat[3] = -5.9584e6;

    Asat = 0;
    Bsat = 0;
    Csat = 0;
    Dsat = 0;
    Esat = 0;
    Fsat = 0;
    Gsat = 0;
    Hsat = 0;

    Gas_Constant = config->GetGas_Constant();
    Tstar = config->GetTemperature_Critical();

}



CWater::~CWater(void) {

	delete [] coeff_saturation;
	delete [] coeff_latent_heat;

}


void CWater::SetTsat(su2double P) {

    Hsat = pow((P/1E6),0.25);
    Esat =  Hsat*Hsat + coeff_saturation[2] * Hsat + coeff_saturation[5];
    Fsat =  coeff_saturation[0]*Hsat*Hsat + coeff_saturation[3]*Hsat + coeff_saturation[6];
    Gsat =  coeff_saturation[1]*Hsat*Hsat + coeff_saturation[4]*Hsat + coeff_saturation[7];
    Dsat =  2.0 * Gsat /(-Fsat -sqrt(Fsat*Fsat - 4.0 * Esat * Gsat));

    Tsat = pow((coeff_saturation[9] + Dsat), 2) - 4.0 * (coeff_saturation[8] + coeff_saturation[9] * Dsat);
    Tsat = -sqrt(Tsat) + coeff_saturation[9] + Dsat;

    Tsat =  Tsat * 0.5;

}

void CWater::SetPsat(su2double T) {


	su2double  T_limited = min(T, Tstar);

	Hsat = T_limited + coeff_saturation[8] / (T_limited - coeff_saturation[9] ) ;
	Asat = Hsat*Hsat + coeff_saturation[0]*Hsat + coeff_saturation[1];
	Bsat = coeff_saturation[2]*Hsat*Hsat + coeff_saturation[3]*Hsat + coeff_saturation[4];
	Csat = coeff_saturation[5]*Hsat*Hsat + coeff_saturation[6]*Hsat + coeff_saturation[7];

	Psat    = -Bsat + sqrt(Bsat*Bsat- 4.0*Asat*Csat);
	Psat    = pow((2.0*Csat/Psat), 4);

	Psat    = Psat *1E6;

}

void CWater::SetSurfaceTension(su2double T, su2double Rdroplet) {

    su2double  T_limited = min(T, Tstar);

	sigma =  1.0 - T/Tstar;
	sigma =  235.8e-3 * (0.375 + 0.625*(1.0-sigma)) * pow(sigma,1.256);

    /*if (T < 373.15) sigma = 0.122 -0.17e-3*T;
    else sigma = 0.138- 0.212e-3 * T;


	sigma = sigma*0.87;
    */


	//if (Rdroplet > 0) sigma = sigma / (1 + 2e-10/Rdroplet);


}

void CWater::SetTLiquid(su2double T, su2double Rcritical, su2double Rdroplet) {

		if (Rdroplet!=0.0) T_l   = max(T, Tsat  - (Tsat - T)*Rcritical/Rdroplet);
		else      T_l = T;

}

void CWater::SetLiquidDensity() {

	    su2double  T_limited = T_l;

		rho_l = 928.08 + 464.63*T_limited/Tstar - 568.46*(T_limited/Tstar)*(T_limited/Tstar)
				- 255.17*(T_limited/Tstar)*(T_limited/Tstar)*(T_limited/Tstar);

}

void CWater::SetLiquidEnthalpy(su2double h_v) {

	    //enthalpy
		h_l = coeff_latent_heat[0] + coeff_latent_heat[1]*(T_l/Tstar)  ;
		h_l = h_l + coeff_latent_heat[2]*pow((T_l/Tstar), 2) +
					coeff_latent_heat[3]*pow((T_l/Tstar), 3);
		h_l = h_v - h_l;
}


CCO2::CCO2(CConfig *config) : CLiquidModel(config) {

	Ei = new su2double [8];

    Ei[0] =  0.9287838;
    Ei[1] = -0.5910856e-2;
    Ei[2] = -0.5047727e-4;
    Ei[3] = -0.9781836e-6;
    Ei[4] = -0.2475923e-7;
    Ei[5] = -0.5433246e-9;
    Ei[6] = -0.6854221e-11;
    Ei[7] = -0.3503978e-13;

    Gas_Constant = config->GetGas_Constant();
    Tstar = config->GetTemperature_Critical();

}



CCO2::~CCO2(void) {

	delete [] Ei;

}


void CCO2::SetTsat(su2double P) {

	su2double P_limited;

	P_limited = P/1e5;
	Tsat = -6E-06 * pow(P_limited, 4)  + 0.0013 * pow(P_limited,3)
	       -0.0981 * pow(P_limited,2) + 4.2985 * P_limited + 198.22;

}

void CCO2::SetPsat(su2double T) {

	su2double  T_limited = T;

	Psat = 4E-05 * pow(T_limited,3) -0.0213 * pow(T_limited,2) + 4.1933 * T_limited - 284.83;
	Psat = Psat * 1e5;

}

void CCO2::SetSurfaceTension(su2double T, su2double Rdroplet) {

    su2double  T_limited = min(T, Tstar);

    sigma = 1.0 - 1.0 * T_limited/Tstar;
    sigma = 84.72 * pow(sigma, 1.281);
    sigma = sigma * 1.0e-3;
}

void CCO2::SetTLiquid(su2double T, su2double Rcritical, su2double Rdroplet) {

		if (Rdroplet!=0.0) T_l   = max(T, Tsat  - (Tsat - T)*Rcritical/Rdroplet);
		else      T_l = T;

}

void CCO2::SetLiquidDensity() {

	    su2double  T_limited = T_l;

        rho_l = 0.0;

        for ( i=0; i< 8; i++) {
             rho_l += 1000.0 * Ei[i] * pow(T_limited-273.15,i);
        }

}

void CCO2::SetLiquidEnthalpy(su2double h_v) {

	su2double T_limited = T_l;

	h_l = -8E-09 * pow(T_limited, 6) + 1E-05 * pow(T_limited, 5) -0.0078 * pow(T_limited, 4)
	      + 2.6673 * pow(T_limited, 3) - 508.95 *pow(T_limited, 2) + 51679 * T_limited - 2E+06;

    h_l = h_l * 1e3 + 252943.0559;   // H ref evaluate at 1bar , 300 K  from Fluidprop, stanmix


}

CR22::CR22(CConfig *config) : CLiquidModel(config) {

	Ei = new su2double [8];

    Ei[0] =  0.1284436e1;
    Ei[1] = -0.3440780e-2;
    Ei[2] = -0.7891991e-5;
    Ei[3] = -0.4743679e-7;
    Ei[4] = -0.6503566e-9;
    Ei[5] = +0.5864434e-11;
    Ei[6] = +0.1317490e-13;
    Ei[7] = -0.1788847e-14;

    Gas_Constant = config->GetGas_Constant();
    Tstar = config->GetTemperature_Critical();

}



CR22::~CR22(void) {

	delete [] Ei;

}


void CR22::SetTsat(su2double P) {

	su2double P_limited;

	P_limited = P/1e5;

	Tsat = 4E-06 * pow(P_limited, 5) - 0.0006 * pow(P_limited, 4)  + 0.0346 * pow(P_limited,3)
	       -0.9326 * pow(P_limited,2) + 14.359 * P_limited + 217.93;


}

void CR22::SetPsat(su2double T) {


	su2double  T_limited = T;

	Psat = 1E-05 * pow(T_limited,3) - 0.0078 * pow(T_limited,2) + 1.7371 * T_limited - 132.32;
	Psat = Psat * 1e5;

}

void CR22::SetSurfaceTension(su2double T, su2double Rdroplet) {

    su2double  T_limited = min(T, Tstar);

    sigma = 1.0 - 1.0 * T_limited/Tstar;
    sigma = 69.93 * (1.0 - 0.154 *pow(sigma, 0.87)) * pow(sigma, 1.285);
    sigma = sigma * 1.0e-3;

}

void CR22::SetTLiquid(su2double T, su2double Rcritical, su2double Rdroplet) {

		if (Rdroplet!=0.0) T_l   = max(T, Tsat  - (Tsat - T)*Rcritical/Rdroplet);
		else      T_l = T;

}

void CR22::SetLiquidDensity() {

	    su2double  T_limited = T_l;

        rho_l = 0.0;

        for ( i=0; i< 8; i++) {
             rho_l += 1000.0 * Ei[i] * pow(T_limited-273.15,i);
        }

}

void CR22::SetLiquidEnthalpy(su2double h_v) {

	su2double T_limited = T_l;

	h_l = -2E-10 * pow(T_limited, 6) + 4E-07 * pow(T_limited, 5) - 0.0003 * pow(T_limited, 4)
	      + 0.1101 * pow(T_limited, 3) - 23.296 *pow(T_limited, 2) + 2614.9 * T_limited - 121670;

    h_l = h_l * 1e3 + 194147.1632;   // H ref evaluate at 1bar , 300 K  from Fluidprop, stanmix


}






CR12::CR12(CConfig *config) : CLiquidModel(config) {

	Ei = new su2double [8];

    Ei[0] =  0.1395549e1;
    Ei[1] = -0.3225727e-2;
    Ei[2] = -0.5998138e-5;
    Ei[3] = -0.1307566e-7;
    Ei[4] = -0.1071878e-8;
    Ei[5] = -0.9172160e-11;
    Ei[6] = +0.3799672e-12;
    Ei[7] = -0.2864822e-14;

    Gas_Constant = config->GetGas_Constant();
    Tstar = config->GetTemperature_Critical();

}



CR12::~CR12(void) {

	delete [] Ei;

}


void CR12::SetTsat(su2double P) {

	su2double P_limited;

	P_limited = P/1e5;

	Tsat = -3E-05 * pow(P_limited, 4)  + 0.0041 * pow(P_limited,3)
	- 0.1972 * pow(P_limited,2) + 5.6176 * P_limited + 192.93;


}

void CR12::SetPsat(su2double T) {


	su2double  T_limited = T;

	Psat = 3E-05 * pow(T_limited,3) - 0.0187 * pow(T_limited,2) + 3.5625 * T_limited - 232.73;
	Psat = Psat * 1e5;

}

void CR12::SetSurfaceTension(su2double T, su2double Rdroplet) {


    su2double  T_limited = min(T, Tstar);

    sigma = 1.0 - 1.0 * T_limited/Tstar;
    sigma = 61.23 * (1.0 - 0.094 *pow(sigma, 0.584)) * pow(sigma, 1.285);
    sigma = sigma * 1.0e-3;

}

void CR12::SetTLiquid(su2double T, su2double Rcritical, su2double Rdroplet) {

		if (Rdroplet!=0.0) T_l   = max(T, Tsat  - (Tsat - T)*Rcritical/Rdroplet);
		else      T_l = T;

}

void CR12::SetLiquidDensity() {

	    su2double  T_limited = T_l;

        rho_l = 0.0;

        for ( i=0; i< 8; i++) {
             rho_l += 1000.0 * Ei[i] * pow(T_limited-273.15,i);
        }

}

void CR12::SetLiquidEnthalpy(su2double h_v) {


	su2double T_limited = T_l;

	h_l = -1E-06 * pow(T_limited, 4) + 0.0009 * pow(T_limited, 3) - 0.3005 * pow(T_limited, 2)
	+ 47.331 * T_limited - 2947.7;

    h_l = h_l * 1e3 + 151401.6891;   // H ref evaluate at 1bar , 300 K  from Fluidprop, stanmix


}

