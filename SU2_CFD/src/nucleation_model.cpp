/*! * nucleation_model.cpp
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

#include "../include/nucleation_model.hpp"


CNucleationModel::CNucleationModel(void) {


}

CNucleationModel::~CNucleationModel(void) { }



CClassicalTheory::CClassicalTheory(CConfig *config, CGeometry *geometry) : CNucleationModel() {

	Gamma = config->GetGamma();
	Gas_Constant = config -> GetGas_ConstantND();
	nDim = geometry->GetnDim();

	Boltzmann = 1.38064852/config->GetBoltzmann_Ref();
	MolMass   = 18.0 * 1.66054e-27 / config ->GetMass_Ref();

}

CClassicalTheory::~CClassicalTheory(void) {


}


void CClassicalTheory::SetNucleationRate (su2double V_i, su2double V_Liquid) {

	Theta = (V_i[nDim + 3] - V_Liquid[2])/Gas_Constant/V_i[0];
	Theta = (Theta-0.5)*Theta;
	Theta = Theta * 2 * (Gamma - 1)/(Gamma + 1);

	if (V_i[nDim + 1] < V_Liquid[3])
		Rc = 0;
	else
	    Rc = 2*V_Liquid[5] / V_Liquid[1] / V_Liquid[6] ;


	J     = -4.0*3.14*V_Liquid[5]*Rc*Rc / 3 / V_i[0] / Boltzmann;

	J = exp(J) * V_i[nDim+2]/V_Liquid[1]*sqrt(2*V_Liquid[5] / 3.14 / MolMass);

	J = J / (1+Theta);

}


void CClassicalTheory::SetGrowthRate (su2double V_i, su2double V_Liquid) {

	if (V_i[0] > V_Liquid[4])
		G = 0;
	else {

		if (V_i[nDim + 1] < V_Liquid[3])
			Rc = 0;
		else
		    Rc = 2*V_Liquid[5] / V_Liquid[1] / V_Liquid[6] ;

		G = V_i[nDim+7] * (V_Liquid[4] - V_i[0]);
		G = G * (1-Rc/(V_Liquid[7] + 1e-30));
		G = G / V_Liquid[1] / (V_i[nDim + 3] - V_Liquid[2]);

		Ni = Gas_Constant * V_Liquid[4] / (V_i[nDim + 3] - V_Liquid[2]);
		Ni = Ni * 0.25 * (Gamma + 1) / (Gamma -1);
		Ni = 0.5 - Ni;
		Ni = Ni * Gas_Constant * V_Liquid[4] / (V_i[nDim+3] - V_Liquid[2]);

		Lambda = 1.5* V_i[nDim+5] * sqrt(Gas_Constant * V_i[0] / V_i[nDim + 1]);
		Pr = Gas_Constant * Gamma / (Gamma - 1) * V_i[nDim+5] / V_i[nDim+7];

		Lambda = V_Liquid[7] + 1.89 * (1-Ni) * Lambda / Pr;
		G   = G / Lambda;
	}

}


