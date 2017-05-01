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

#include "../include/heat_capacity.hpp"


CHeatCapacity::CHeatCapacity() {
}

CHeatCapacity::CHeatCapacity(CConfig *config) {

	Gas_Constant = config->GetGas_Constant();
	coeff_Cp0  = new su2double [5];

}

CHeatCapacity::~CHeatCapacity(void) {

	delete [] coeff_Cp0;
}



void CHeatCapacity::Set_Cv0(su2double T) {

	unsigned short i;
    su2double Cp0;

    Cp0 = 0;

    for (i=0;i<5;i++)
       Cp0 += coeff_Cp0[i] * pow(T, i);

    Cv0 = Cp0 - Gas_Constant;

}


CSteam::CSteam() : CHeatCapacity() {

    coeff_Cp0[0] =  34.51/0.018      ;
    coeff_Cp0[1] = - 0.01/0.018      ;
    coeff_Cp0[2] =   4.686e-05/0.018 ;
    coeff_Cp0[3] = - 3.781e-08/0.018 ;
    coeff_Cp0[3] =   1.173e-11/0.018 ;
}



CSteam::~CSteam(void) {

}

CCo2::CCo2() : CHeatCapacity() {

    coeff_Cp0[0] =  17.46/0.044      ;
    coeff_Cp0[1] =   0.09/0.044      ;
    coeff_Cp0[2] = - 1.008e-04/0.044 ;
    coeff_Cp0[3] =   6.539e-08/0.044 ;
    coeff_Cp0[3] = - 1.845e-11/0.044 ;
}



CCo2::~CCo2(void) {

}





