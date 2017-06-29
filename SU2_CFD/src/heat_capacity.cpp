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

	Gas_Constant   = config->GetGas_Constant();
	Gas_ConstantND = config->GetGas_ConstantND();
	Gamma          = config->GetGamma();
	Constant_Gamma = config->Get_ConstantGamma();
	Tref           = config->GetTemperature_Ref();

	coeff_Cp0      = config->GetCoeff_HeatCapacity();

	T_Total        = config->GetRiemann_Var2(config->GetMarker_All_TagBound(1));

}

CHeatCapacity::~CHeatCapacity(void) {

}


CHeatCapacity_Dimensional::CHeatCapacity_Dimensional() : CHeatCapacity() {

}

CHeatCapacity_Dimensional::CHeatCapacity_Dimensional(CConfig *config) : CHeatCapacity(config) {

}

CHeatCapacity_Dimensional::~CHeatCapacity_Dimensional(void) {

}

void CHeatCapacity_Dimensional::Set_Cv0(su2double T) {

	unsigned short i;
    su2double Cp0, T_limited;

    if (coeff_Cp0[0] == 0.0) {
    	cout << "Warning: no value set for HEAT_CAPACITY_MODEL" << endl;
        getchar();
    }

    if (Constant_Gamma) {
    	Cp0 = Gas_Constant * Gamma / (Gamma-1);

    } else {

    	Cp0 = 0;

        for (i=0;i<5;i++)
           Cp0 += coeff_Cp0[i] * pow(T, i);
    }

    Cv0 = Cp0 - Gas_Constant;


}


CHeatCapacity_Dimensionless::CHeatCapacity_Dimensionless() : CHeatCapacity() {

}

CHeatCapacity_Dimensionless::CHeatCapacity_Dimensionless(CConfig *config) : CHeatCapacity(config) {

}

CHeatCapacity_Dimensionless::~CHeatCapacity_Dimensionless(void) {

}

void CHeatCapacity_Dimensionless::Set_Cv0(su2double T) {

	unsigned short i;
    su2double Cp0, T_limited;

    if (coeff_Cp0[0] == 0.0) {
    	cout << "Warning: no value set for HEAT_CAPACITY_MODEL" << endl;
        getchar();
    }


    if (Constant_Gamma) {
    	Cp0 = Gas_Constant * Gamma / (Gamma-1);

    } else {
        Cp0 = 0;

        T_limited = T*Tref;

        for (i=0;i<5;i++) {
           Cp0 += coeff_Cp0[i] * pow(T_limited, i);
        }
    }

    Cp0 = Cp0 / Gas_Constant * Gas_ConstantND;
    Cv0 = Cp0 - Gas_ConstantND;

}





