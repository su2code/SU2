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


CLiquidModel::CLiquidModel(void) {

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

void CLiquidModel::Set_LiquidProp(su2double P, su2double T, su2double rho, su2double h_v, su2double *Two_Phase_Var) {

	SetRadius(Two_Phase_Var);

	SetTsat(P);
	SetPsat(T);

	SetSurfaceTension(T);

    SetTLiquid(T);
    SetLiquidDensity();

    SetLiquidEnthalpy(h_v);

    SetRCritical(P, T);

    SetDensity_Mixture(rho, Two_Phase_Var);
}

void CLiquidModel::SetRCritical(su2double P, su2double T) {

    if (Tsat > T) {
		dGibbs = fabs(T * Gas_Constant * log(P/Psat));

		Rc = 2.0*sigma / (rho_l * dGibbs);
    }
    else { Rc = 0.0  	;    }

}

void CLiquidModel::SetDensity_Mixture(su2double rho, su2double *Two_Phase_Var) {

	    y = Two_Phase_Var[3]*(rho_l - rho);
	    y = y + 0.75 * rho / 3.14;
	    y = y*rho_l / y;

	    rho_m = y/ rho_l + (1.0 - y)/ rho;
	    rho_m = 1.0/ rho_m;
}



void CLiquidModel::SetRadius(su2double *Two_Phase_Var) {

	if (Two_Phase_Var[0] == 0) R = 0;
	else 	R = Two_Phase_Var[1]/Two_Phase_Var[0];

}




CWater::CWater() : CLiquidModel() {

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

}



CWater::~CWater(void) {

	delete [] coeff_saturation;
	delete [] coeff_latent_heat;

}


void CWater::SetTsat(su2double P) {


    Hsat = pow((P/1E6),0.25);
    Esat =  Hsat*Hsat + coeff_saturation[2] * Hsat + coeff_saturation[5];
    Fsat =  coeff_saturation[0]*Hsat*Hsat + coeff_saturation[4]*Hsat + coeff_saturation[6];
    Gsat =  coeff_saturation[1]*Hsat*Hsat + coeff_saturation[4]*Hsat + coeff_saturation[7];
    Dsat =  2.0 * Gsat /(-Fsat -sqrt(Fsat*Fsat - 4.0 * Esat * Gsat));

    Tsat = pow((coeff_saturation[9] + Dsat), 2) - 4.0 * (coeff_saturation[8] + coeff_saturation[9] * Dsat);
    Tsat = -sqrt(Tsat) + coeff_saturation[9] + Dsat;

    Tsat =  Tsat * 0.5;
}

void CWater::SetPsat(su2double T) {

	Hsat = T + coeff_saturation[8] / (T - coeff_saturation[9] ) ;
	Asat = Hsat*Hsat + coeff_saturation[0]*Hsat + coeff_saturation[1];
	Bsat = coeff_saturation[2]*Hsat*Hsat + coeff_saturation[3]*Hsat + coeff_saturation[4];
	Csat = coeff_saturation[5]*Hsat*Hsat + coeff_saturation[6]*Hsat + coeff_saturation[7];

	Psat    = -Bsat + sqrt(Bsat*Bsat- 4.0*Asat*Csat);
	Psat    = pow((2.0*Csat/Psat), 4);

	Psat    = Psat *1E6;
}

void CWater::SetSurfaceTension(su2double T) {


	sigma =  1.0 - 1.0 * T/Tstar;
	sigma =  235.8e-3 * (0.375 + 0.625*(1.0-sigma)) * pow(sigma,1.256);
}

void CWater::SetTLiquid(su2double T) {

    	//Guess Rcritical, a loop is required to evaluate the right properties
		Rc = 1e-15;

		if (R!=0) T_l   = Tsat  - (Tsat - T)*Rc/R;
		else      T_l = T;

}

void CWater::SetLiquidDensity() {

		rho_l = 928.08 + 464.63*T_l/Tstar - 568.46*(T_l/Tstar)*(T_l/Tstar)
				- 255.17*(T_l/Tstar)*(T_l/Tstar)*(T_l/Tstar);

}

void CWater::SetLiquidEnthalpy(su2double h_v) {

	    //enthalpy
		h_l = coeff_latent_heat[0] + coeff_latent_heat[1]*(T_l/Tstar)  ;
		h_l = h_l + coeff_latent_heat[2]*pow((T_l/Tstar), 2) +
					coeff_latent_heat[3]*pow((T_l/Tstar), 3);
		h_l = h_v - h_l;
}

