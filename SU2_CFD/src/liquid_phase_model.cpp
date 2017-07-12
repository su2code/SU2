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
  Fluid = config->GetKind_Liquid_Model();

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

		if (Fluid == WATER)
			SetRCritical(P, T);
		else
//			SetRCritical(P, T);
			SetRCritical(h_v-h_l, P, T);

		//Rc = min(Rc, Rdroplet);

		SetDensity_Mixture(rho, mom3);
	} else {

		sigma = 0;
		T_l = T;
		rho_l = rho;
		Rc = 0;
		h_l = h_v;
		rho_m = rho;

	}


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

    Ptriple = 611;
    Ttriple = 273.15;

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

/*	if (P < Ptriple) {
		cout << "Warning: Pressure lower than triple point" << endl;
		cout << "Saturation conditions are not reliable" << endl;
	}
*/


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

void CWater::SetRCritical(su2double P, su2double T) {

    if (Psat < P) {

		dGibbs = T * Gas_Constant * log(P/Psat);
		Rc = 2.0*sigma / (rho_l * dGibbs);
    }
    else { Rc = 0.0;}

 /*   if (Tsat > T) {
		dGibbs = P * (Tsat - T) / Tsat;
		Rc = 2.0*sigma / (rho_l * dGibbs);
   }
*/


}

void CWater::SetSurfaceTension(su2double T, su2double Rdroplet) {

    su2double  T_limited = min(T, Tstar);

	sigma =  1.0 - T/Tstar;
	sigma =  235.8e-3 * (0.375 + 0.625*(1.0-sigma)) * pow(sigma,1.256);



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
	ac = new su2double [4];
	tc = new su2double [4];
	as = new su2double [3];
    ts = new su2double [3];

    Ei[0] =  0.9287838;
    Ei[1] = -0.5910856e-2;
    Ei[2] = -0.5047727e-4;
    Ei[3] = -0.9781836e-6;
    Ei[4] = -0.2475923e-7;
    Ei[5] = -0.5433246e-9;
    Ei[6] = -0.6854221e-11;
    Ei[7] = -0.3503978e-13;

    ac[0] = -7.0602087;
    ac[1] =  1.9391218;
    ac[2] = -1.6463597;
    ac[3] = -3.2995634;

    tc[0] =  1.0;
    tc[1] =  1.5;
    tc[2] =  2.0;
    tc[3] =  4.0;

    as[0] = -14.740846;
    as[1] =  2.4327015;
    as[2] = -5.3061778;

    ts[0] =  1.0;
    ts[1] =  1.9;
    ts[2] =  2.9;

    Gas_Constant = config->GetGas_Constant();
    Tstar = config->GetTemperature_Critical();
    Pstar = config->GetPressure_Critical();

    MolMass   = config->GetMolecular_Mass() * 1.66054e-27;

    Ptriple = 5.12e5;
    Ttriple = 216.57;

}



CCO2::~CCO2(void) {

	delete [] Ei;
    delete [] ac;
    delete [] tc;
    delete [] as;
    delete [] ts;
}


void CCO2::SetTsat(su2double P) {

	su2double P_limited;

	P_limited = P/Pstar;

	if (P > Pstar)
		Tsat = 0; //outside two phase region
	else {

		Tsat = -6.2185E-01 * pow(P_limited, 4)  + 1.6817 * pow(P_limited,3)
		- 1.7566 * pow(P_limited,2) + 1.0431 * P_limited + 6.5160E-01;

		Tsat = Tsat*Tstar;

		if (P < Ptriple) {
			cout << "Warning: Pressure lower than triple point" << endl;
			cout << "Saturation conditions are not reliable" << endl;
		}
	}
}

void CCO2::SetPsat(su2double T) {

	su2double  T_limited;

	if (T > Tstar)
		Psat = 0; //outside two phase region

	else if (T < Tstar && T > Ttriple) {

		T_limited = T/Tstar;
		Psat = 0.0;
		for (i=0 ; i<4 ; i++) Psat += ac[i] * pow((1- T_limited), tc[i]);

		Psat = Psat / T_limited;
		Psat = exp(Psat) * Pstar;

	} else {

		T_limited = T/Ttriple;

		Psat = 0.0;
		for (i=0 ;i<3; i++) Psat += as[i] * pow((1- T_limited), ts[i]);

		Psat = Psat / T_limited;
		Psat = exp(Psat) * Ptriple;

	}

}

void CCO2::SetRCritical(su2double dh, su2double P, su2double T) {


	su2double dG, dGsat;

	if (Psat < P) {
		dGibbs = T * Gas_Constant * log(P/Psat);
		Rc = sigma / (rho_l * dGibbs);
		Rc = Rc * (1 + sqrt(1+ 2*3.28e-10 * dGibbs * rho_l/sigma) ) -2*3.28e-10;
		Rc = max(Rc, 0.0);

	} else
		Rc = 0;

/*	if (Psat < P) {
		dGibbs = T * Gas_Constant * log(P/Psat);
		Rc = sigma / (rho_l * dGibbs);
		Rc = Rc * (1 + sqrt( pow(MolMass/36/3.1415/rho_l, 1.0/3.0)*dGibbs*rho_l/sigma ));
		Rc = max(Rc, 0.0);

	} else
		Rc = 0;
*/


/*    if (Tsat > T) {
		dGibbs = (dh) * (Tsat - T) / Tsat;
		Rc = 2.0*sigma / (rho_l * dGibbs);
   } else
	    Rc = 0;
*/
    su2double a = 204.3918, b = 6.06e-4;

//    if (Tsat > T) {
/*		dG = -Gas_Constant*T*log(P*(1/rho - b)/Gas_Constant/T) + a/2/sqrt(2)*log( (1/rho+(1-sqrt(2))*b) / (1/rho+(1+sqrt(2))*b));
		dG = dG + P/rho - Gas_Constant*T;
		dG = dG + T * Gas_Constant * log(P);

		dGsat = - Gas_Constant*T*log(Psat*(1/rho - b)/Gas_Constant/T) - a/2/sqrt(2)*log( (1/rho+(1-sqrt(2))*b)/(1/rho+(1+sqrt(2))*b));
		dGsat = dGsat + Psat/rho - Gas_Constant*T;
		dGsat = dGsat + T * Gas_Constant * log(Psat);
*/
/*    	dG = a/b * log(sqrt(2)-1)/2 + P/rho - Gas_Constant*T;
    	dG = dG + T * Gas_Constant * log(P);

    	dGsat = a/b * log(sqrt(2)-1)/2 + Psat/rho - Gas_Constant*T;
    	dG = dG + T * Gas_Constant * log(Psat);

		dGibbs = dG - dGsat;

		Rc = abs(2.0*sigma / (rho_l * dGibbs));
*/
//   } else
//	    Rc = 0;

}

void CCO2::SetRCritical(su2double P, su2double T) {

	if (Psat < P) {
		dGibbs = T * Gas_Constant * log(P/Psat);
		Rc = sigma / (rho_l * dGibbs)*2;
	} else
		Rc = 0;

}

void CCO2::SetSurfaceTension(su2double T, su2double Rdroplet) {

    su2double  T_limited = min(T, Tstar);

    sigma = 1.0 - T_limited/Tstar;
    sigma = 84.72 * pow(sigma, 1.281);
    sigma = sigma * 1.0e-3;

    //if (Rdroplet > 0) sigma = sigma / (1 + 2e-10/Rdroplet);
}

void CCO2::SetTLiquid(su2double T, su2double Rcritical, su2double Rdroplet) {

		if (Rdroplet!=0.0) T_l   = max(T, Tsat  - (Tsat - T)*Rcritical/Rdroplet);
		else T_l = Tsat;

		if (T_l > Tsat) T_l = Tsat;
}

void CCO2::SetLiquidDensity() {

	    su2double  T_limited = T_l;

        rho_l = 0.0;

        for ( i=0; i< 8; i++) {
             rho_l = rho_l + 1000.0 * Ei[i] * pow(T_limited-273.15, i);

        }

}

void CCO2::SetLiquidEnthalpy(su2double h_v) {

	su2double T_limited = T_l;

//	if (T_l > Tstar) {
//		h_l = 0;
//	} else  {

//	    h_l = 8.9641E-06 * pow(T_limited, 4) - 9.0172E-03 * pow(T_limited, 3) + 3.3989 *pow(T_limited, 2)
//	      - 5.6672E+02 * T_limited + 3.5023E+04;
//	    h_l = h_l * 1e3;   // H ref evaluate at 1bar , 300 K  from Fluidprop, stanmix


	    h_l = -5.1612E-07 * pow(T_limited, 5) + 6.5716E-04 *pow(T_limited, 4)
	          - 3.3400E-01 * pow(T_limited, 3) + 8.4682E+01 *pow(T_limited, 2)
	          - 1.0710E+04 * T_limited + 5.4085E+05;

	    h_l = h_v - h_l * 1000;

//		if (T_l < Ttriple) {
//
//		    h_l = 8.9641E-06 * pow(Ttriple, 4) - 9.0172E-03 * pow(Ttriple, 3) + 3.3989 *pow(Ttriple, 2)
//		      - 5.6672E+02 * Ttriple + 3.5023E+04;
//		    h_l = h_l * 1e3;
//			cout << "Warning: Temperature lower than triple point" << endl;
//			cout << "Liquid enthalpy is not reliable" << endl;
//		}
//
//	}

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

