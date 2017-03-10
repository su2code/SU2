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

}

CLiquidModel::~CLiquidModel(void) { }



CWater::CWater(void) : CLiquidModel() {

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

    rho_l  = 0.0;
    h_l    = 0.0;
    dGibbs = 0.0;
    Tsat   = 0.0;
    Psat   = 0.0;
    sigma  = 0.0;

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


void CWater::SetLiquidDensity_PT(su2double P, su2double T, su2double Tstar) {
	rho_l = 928.08 + 464.63*T/Tstar - 568.46*(T/Tstar)*(T/Tstar) &
	                 - 255.17*pow((T/Tstar),3);

}

/*!
* \brief return liquid enthalpy value.
*/
void CWater::SetLiquidEnthalpy_PT(su2double P, su2double T, su2double h_v, su2double Tstar) {

    h_l = coeff_latent_heat[0] + coeff_latent_heat[1]*(T/Tstar)  ;
    h_l = h_l + coeff_latent_heat[2]*pow((T/Tstar), 2) +
	            coeff_latent_heat[3]*pow((T/Tstar), 3);
    h_l = h_v - h_l;

}

/*!
* \brief return surface tension value.
*/
void CWater::SetSurfaceTension_T(su2double T, su2double Tstar) {

    sigma =  1.0 - 1.0 * T/Tstar;
    sigma =  235.8e-3 * (0.375 + 0.625*(1.0-sigma)) * pow(sigma,1.256);

}

/*!
* \brief return saturation temperature value at correspondent pressure.
*/
void CWater::SetTsat_P(su2double P) {
    Hsat = pow((P/1E6),0.25);
    Esat =  Hsat*Hsat + coeff_saturation[2] * Hsat + coeff_saturation[5];
    Fsat =  coeff_saturation[0]*Hsat*Hsat + coeff_saturation[4]*Hsat + coeff_saturation[6];
    Gsat =  coeff_saturation[1]*Hsat*Hsat + coeff_saturation[4]*Hsat + coeff_saturation[7];
    Dsat =  2.0 * Gsat /(-Fsat -sqrt(Fsat*Fsat - 4.0 * Esat * Gsat));

    Tsat = pow((coeff_saturation[9] + Dsat), 2) - 4.0  (coeff_saturation[8] + coeff_saturation[9] * Dsat);
    Tsat = -sqrt(Tsat) + coeff_saturation[9] + Dsat;

    Tsat =  Tsat * 0.5;

}

/*!
* \brief return free gibbs energy variation and Psat at same temperature.
*/
void CWater::SetdGibbs_PT(su2double P, su2double T, su2double Gas_Constant) {
    Hsat = T + coeff_saturation[8] / (T - coeff_saturation[9] ) ;
    Asat = Hsat*Hsat + coeff_saturation[0]*Hsat + coeff_saturation[1];
    Bsat = coeff_saturation[2]*Hsat*Hsat + coeff_saturation[3]*Hsat + coeff_saturation[4];
    Csat = coeff_saturation[5]*Hsat*Hsat + coeff_saturation[6]*Hsat + coeff_saturation[7];

    Psat    = -Bsat + sqrt(Bsat*Bsat- 4.0*Asat*Csat);
    Psat    = pow((2.0*Csat/Psat), 4);

    Psat    = Psat *1E6;

	dGibbs = abs(T * Gas_Constant * log(P/Psat));

}


