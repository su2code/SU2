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

CNucleationModel::CNucleationModel() {

}

CNucleationModel::CNucleationModel(CConfig *config) {

}

CNucleationModel::~CNucleationModel(void) { }



CClassicalTheory::CClassicalTheory(CConfig *config) : CNucleationModel() {

	Gas_Constant = config -> GetGas_ConstantND();

	Boltzmann = 1.38064852 * 1.0e-23 /config->GetBoltzmann_Ref();
	MolMass   = config->GetMolecular_Mass() * 1.66054e-27 / config ->GetMass_Ref();

	Fluid     = config->GetKind_Liquid_Model();
}

CClassicalTheory::~CClassicalTheory(void) {


}


su2double CClassicalTheory::SetNucleationRate (su2double P, su2double T, su2double rho,
		                      su2double h, su2double k, su2double mu, su2double Gam, su2double Cp, su2double *V_l) {

	su2double J1, J2, J3, J0;
	// V_l = T, rho, h, Psat, Tsat, sigma, Rc, R, rhom, (G)

	Theta = (h - V_l[2])/Gas_Constant/T - 0.5;

	Theta = Theta * (h - V_l[2])/Gas_Constant/T;

	Theta = Theta * 2 * (Gam - 1)/(Gam + 1);

	if (Theta < 0 ) {
		cout << "WARNING: Non-isothermal correction < 0, isothermal J model imposed" << endl;
		Theta = 0.0;
	}

	if (Fluid == WATER) {

			J  = exp(-4.0*3.1415*V_l[5]*V_l[6]*V_l[6] / 3 / T / Boltzmann);
			J = J * rho/V_l[1]*sqrt(2*V_l[5] / 3.1415 / pow(MolMass,3) );
			J = J / (1+Theta);


	} else {

			J0 = log(1.0/(1.0+Theta));

			J1 = -4.0*3.1415*V_l[5]*V_l[6]*V_l[6] / 3 / T / Boltzmann;

			J2 = sqrt(2*V_l[5] / 3.1415 / pow(MolMass,3) );
			J2 = log(J2);

			J3 = log(rho/V_l[1]);

			J = J0 + J1 + J2 + J3;

			J = exp(J);

	}

	if (J < 0 || J != J) {

		cout << "Error in Nucleation rate evaluation! J = " << J << " " << J0 << " " << J1 << " " << J2 << " " << J3 << endl;
		cout << "Last values in function: "<< endl;
		cout << "Theta " << " " << Theta << endl;
		cout << "delta Enthalpy " << " " << h - V_l[2] << endl;
		cout << "Surface tension"  << " " << V_l[5] << endl;

		cout << h << " "<<  V_l[2] << " " << T <<" "<<  V_l[4] << " " << V_l[0]<<" "<< V_l[1] <<" "<< rho << " " << V_l[7] <<" "<< V_l[6] <<endl;
		getchar();
		//getchar();//exit(EXIT_FAILURE);

	}

/*	cout << h << " "<<  V_l[2] << " " << T <<" "<<  V_l[4] << " " << V_l[0]<<" "<< rho << " " << V_l[1] <<" "<< V_l[7] <<" "
			<< V_l[6] << " " << J << " " << V_l[5] << " " << Theta << endl;
	getchar();
*/

	return J;



}

su2double CClassicalTheory::SetGrowthRate (su2double P, su2double T, su2double rho,
		                      su2double h, su2double k, su2double mu, su2double Gam, su2double Cp, su2double *V_l) {

	su2double Delta, Alpha;

	// V_l = T, rho, h, Psat, Tsat, sigma, Rc, R, rho_m, (G)
    if (Fluid == WATER ) {

			G = k * (V_l[4] - T);

			if (V_l[7] > V_l[6] && V_l[7] !=0) G = G * (1-V_l[6]/V_l[7]);

			G = G / V_l[1] / (h - V_l[2]);

			Ni = Gas_Constant * V_l[4] / (h - V_l[2]);
			Ni = Ni * 0.25 * (Gam + 1) / (Gam -1);
			Ni = 0.5 - Ni;
			Ni = Ni * Gas_Constant * V_l[4] / (h - V_l[2]);

			Lambda = 1.5* mu * sqrt(Gas_Constant * T) / P;

			Pr = Gas_Constant * Gam / (Gam - 1) * mu / k;
			Lambda = V_l[7] + 1.89 * (1-Ni) * Lambda / Pr;

			G = G/Lambda;

			if (G < 0 ||  G!= G) {
				cout << "Error in Growth rate evaluation!" << endl;
				cout << "Last values in function: "<< endl;
				cout << "Density(vapor) " << " " << rho << endl;
				cout << "Radius(vapor) "  << " " << V_l[7] << endl;
				cout << "Critical Radius(vapor) "  << " " << V_l[6] << endl;
				cout << "Enthalpy(liquid) " << " " << V_l[2] << endl;

				exit(EXIT_FAILURE);
			}

    } else {

			Delta  = 32.0 / 7.5 / 3.1415 / rho;
			Delta  = Delta / 2.0/ sqrt(2.0 * Boltzmann * T / MolMass /3.1415);
			Alpha  = V_l[7] + 2.0 * sqrt(2.0 * 3.1415) * (k/(Cp + Cp/Gam)) * Delta;
			Alpha  = k/ Alpha;

			G = Alpha * fabs(V_l[0] - T)/ V_l[1];

			//if (h - V_l[2] > 0)
			G = G /(h - V_l[2]);

			if (Fluid == CO2) {G/=1;}

/*			else {
				cout << "Warning: fluidprop error on liquid enthalpy" << endl;
				cout << "Evaluation will continue anyway" << endl;
			}
*/
			if (G < 0 || G != G) {

				cout << "Error in Growth rate evaluation!" << endl;
				cout << "Last values in function: "<< endl;
				cout << "Density(vapor) " << " " << rho << endl;
				cout << "Radius(vapor) "  << " " << V_l[7] << endl;
				cout << "Temperature(liquid) " << " " << V_l[0] << endl;
				cout << "Enthalpy(liquid) " << " " << V_l[2] << endl;

                //break;
				cout << V_l[1] << " " << 2.0 * Boltzmann * T / MolMass /3.1417 << " " << V_l[7] << " " << Cp << " " << Gam <<
				" " << k << " " << Delta << " " << " " << G << " " << h << " " << V_l[2] << " " << V_l[0] - T<< endl;
				getchar();//exit(EXIT_FAILURE);
			}

//			cout << "Grate" << G << endl;
    }

	return G;
}



