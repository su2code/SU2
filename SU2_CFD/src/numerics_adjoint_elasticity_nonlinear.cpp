/*!
 * \file numerics_adjoint_elasticity_nonlinear.cpp
 * \brief This file contains the routines for setting the tangent matrix and residual of an adjoint FEM structural problem.
 * \author R. Sanchez
 * \version 4.1.0 "Cardinal"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
 *                      Dr. Thomas D. Economon (economon@stanford.edu).
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *
 * Copyright (C) 2012-2016 SU2, the open-source CFD code.
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

#include "../include/numerics_structure.hpp"
#include <limits>


CFEM_NeoHookean_Comp_Adj::CFEM_NeoHookean_Comp_Adj(unsigned short val_nDim, unsigned short val_nVar,
                                   CConfig *config) : CFEM_NonlinearElasticity(val_nDim, val_nVar, config) {

	Mu_E = 1 / (2.0*(1.0 + Nu));
	Lambda_E = Nu /((1.0+Nu)*(1.0-2.0*Nu));

}

CFEM_NeoHookean_Comp_Adj::~CFEM_NeoHookean_Comp_Adj(void) {

}

void CFEM_NeoHookean_Comp_Adj::Compute_Plane_Stress_Term(CElement *element, CConfig *config) {

	su2double j_red = 1.0;
	su2double fx = 0.0, fpx = 1.0;
	su2double xkm1 = 1.0, xk = 1.0;
	su2double cte = 0.0;

	unsigned short iNR, nNR;
	su2double NRTOL;

	// Maximum number of iterations and tolerance (relative)
	nNR = 10;
	NRTOL = 1E-25;

	// j_red: reduced jacobian, for the 2x2 submatrix of F
	j_red = F_Mat[0][0] * F_Mat[1][1] - F_Mat[1][0] * F_Mat[0][1];
	// cte: constant term in the NR method
	cte = Lambda*log(j_red) - Mu;

	// f(f33)  = mu*f33^2 + lambda*ln(f33) + (lambda*ln(j_red)-mu) = 0
	// f'(f33) = 2*mu*f33 + lambda/f33

	for (iNR = 0; iNR < nNR; iNR++){
		fx  = Mu*pow(xk,2.0) + Lambda*log(xk) + cte;
		fpx = 2*Mu*xk + (Lambda / xk);
		xkm1 = xk - fx / fpx;
		if (((xkm1 - xk) / xk) < NRTOL) break;
		xk = xkm1;
	}

	f33 = xkm1;

}

void CFEM_NeoHookean_Comp_Adj::Compute_Constitutive_Matrix(CElement *element, CConfig *config) {

	su2double Mu_p = 0.0, Lambda_p = 0.0;

	/*--- This can be done in a better way ---*/
	if (J_F != 0.0){
		Mu_p = (Mu_E - Lambda_E*log(J_F))/J_F;
		Lambda_p = Lambda_E/J_F;
	}

	/*--- Assuming plane strain ---*/

	if (nDim == 2){
	    D_Mat[0][0] = Lambda_p + 2.0 * Mu_p;	D_Mat[0][1] = Lambda_p;					D_Mat[0][2] = 0.0;
	    D_Mat[1][0] = Lambda_p;					D_Mat[1][1] = Lambda_p + 2.0 * Mu_p;	D_Mat[1][2] = 0.0;
	    D_Mat[2][0] = 0.0;						D_Mat[2][1] = 0.0;						D_Mat[2][2] = Mu_p;
	}
	else if (nDim == 3){
	    D_Mat[0][0] = Lambda_p + 2.0 * Mu_p;	D_Mat[0][1] = Lambda_p;					D_Mat[0][2] = Lambda_p;				D_Mat[0][3] = 0.0;	D_Mat[0][4] = 0.0;	D_Mat[0][5] = 0.0;
	    D_Mat[1][0] = Lambda_p;					D_Mat[1][1] = Lambda_p + 2.0 * Mu_p;	D_Mat[1][2] = Lambda_p;				D_Mat[1][3] = 0.0;	D_Mat[1][4] = 0.0;	D_Mat[1][5] = 0.0;
	    D_Mat[2][0] = Lambda_p;					D_Mat[2][1] = Lambda_p;					D_Mat[2][2] = Lambda_p + 2.0 * Mu_p;	D_Mat[2][3] = 0.0;	D_Mat[2][4] = 0.0;	D_Mat[2][5] = 0.0;
	    D_Mat[3][0] = 0.0;						D_Mat[3][1] = 0.0;						D_Mat[3][2] = 0.0;					D_Mat[3][3] = Mu_p;	D_Mat[3][4] = 0.0;	D_Mat[3][5] = 0.0;
	    D_Mat[4][0] = 0.0;						D_Mat[4][1] = 0.0;						D_Mat[4][2] = 0.0;					D_Mat[4][3] = 0.0;	D_Mat[4][4] = Mu_p;	D_Mat[4][5] = 0.0;
	    D_Mat[5][0] = 0.0;						D_Mat[5][1] = 0.0;						D_Mat[5][2] = 0.0;					D_Mat[5][3] = 0.0;	D_Mat[5][4] = 0.0;	D_Mat[5][5] = Mu_p;
	}

}

void CFEM_NeoHookean_Comp_Adj::Compute_Stress_Tensor(CElement *element, CConfig *config) {

	unsigned short iVar,jVar;
	su2double Mu_J = 0.0, Lambda_J = 0.0;
	su2double dij = 0.0;

	/*--- This can be done in a better way ---*/
	if (J_F != 0.0){
		Mu_J = Mu_E/J_F;
		Lambda_J = Lambda_E/J_F;
	}

	for (iVar = 0; iVar < 3; iVar++){
		for (jVar = 0; jVar < 3; jVar++){
			if (iVar == jVar) dij = 1.0;
			else if (iVar != jVar) dij = 0.0;
			Stress_Tensor[iVar][jVar] = Mu_J * (b_Mat[iVar][jVar] - dij) + Lambda_J * log(J_F) * dij;
		}
	}

}




CFEM_DielectricElastomer_Adj::CFEM_DielectricElastomer_Adj(unsigned short val_nDim, unsigned short val_nVar,
                                   CConfig *config) : CFEM_NonlinearElasticity(val_nDim, val_nVar, config) {

	su2double Electric_Field_Mod;
	su2double *Electric_Field_Dir = config->Get_Electric_Field_Dir();
	su2double ref_Efield_mod;
	unsigned short iVar, iDim;

	ke_DE = config->GetDE_Modulus();

	nElectric_Field = config->GetnElectric_Field();
	nDim_Electric_Field = config->GetnDim_Electric_Field();

	if (nDim != nDim_Electric_Field) cout << "DIMENSIONS DON'T AGREE (Fix this)" << endl;

	/*--- Initialize pointer for the electric field ---*/
	Electric_Field_Ref_Unit = new su2double* [nElectric_Field];
	for (iVar = 0; iVar < nElectric_Field; iVar++) {
		Electric_Field_Ref_Unit[iVar] = new su2double[nDim_Electric_Field];
	}

	ref_Efield_mod = 0.0;
	/*--- Normalize the electric field vector ---*/
	for (iDim = 0; iDim < nDim_Electric_Field; iDim++) {
		ref_Efield_mod += Electric_Field_Dir[iDim]*Electric_Field_Dir[iDim];
	}
	ref_Efield_mod = sqrt(ref_Efield_mod);

	if (ref_Efield_mod == 0){
		cout << "The electric field direction is incorrectly defined!!!!!" << endl;
		exit(EXIT_FAILURE);
	}

	/*--- Assign values to the auxiliary Electric_Field structure ---*/
	for (iVar = 0; iVar < nElectric_Field; iVar++) {
		for (iDim = 0; iDim < nDim_Electric_Field; iDim++) {
			Electric_Field_Ref_Unit[iVar][iDim] = Electric_Field_Dir[iDim]/ref_Efield_mod;
		}
	}

	/*--- Auxiliary vector for computing the electric field in the current configuration ---*/
	Electric_Field_Curr = new su2double[nDim_Electric_Field];
	for (iDim = 0; iDim < nDim_Electric_Field; iDim++) {
		Electric_Field_Curr[iDim] = 0.0;
	}

	/*--- Auxiliary vector for hosting the electric field modulus in the reference configuration ---*/
	E_mod = new su2double[nElectric_Field];
	for (iVar = 0; iVar < nElectric_Field; iVar++) {
		E_mod[iVar] = config->Get_Electric_Field_Mod(iVar);
	}


}

CFEM_DielectricElastomer_Adj::~CFEM_DielectricElastomer_Adj(void) {

	unsigned short iVar;

	for (iVar = 0; iVar < nElectric_Field; iVar++){
		delete [] Electric_Field_Ref_Unit[iVar];
	}

	delete [] Electric_Field_Ref_Unit;
	delete [] Electric_Field_Curr;

}

void CFEM_DielectricElastomer_Adj::Compute_Plane_Stress_Term(CElement *element, CConfig *config) {

}

void CFEM_DielectricElastomer_Adj::Compute_Constitutive_Matrix(CElement *element, CConfig *config) {

	/*--- This reduces performance by now, but it is temporal ---*/

	if (nDim == 2){
	    D_Mat[0][0] = 0.0;	D_Mat[0][1] = 0.0;	D_Mat[0][2] = 0.0;
	    D_Mat[1][0] = 0.0;	D_Mat[1][1] = 0.0;	D_Mat[1][2] = 0.0;
	    D_Mat[2][0] = 0.0;	D_Mat[2][1] = 0.0;	D_Mat[2][2] = 0.0;
	}
	else if (nDim == 3){
	    D_Mat[0][0] = 0.0;	D_Mat[0][1] = 0.0;	D_Mat[0][2] = 0.0;	D_Mat[0][3] = 0.0;	D_Mat[0][4] = 0.0;	D_Mat[0][5] = 0.0;
	    D_Mat[1][0] = 0.0;	D_Mat[1][1] = 0.0;	D_Mat[1][2] = 0.0;	D_Mat[1][3] = 0.0;	D_Mat[1][4] = 0.0;	D_Mat[1][5] = 0.0;
	    D_Mat[2][0] = 0.0;	D_Mat[2][1] = 0.0;	D_Mat[2][2] = 0.0;	D_Mat[2][3] = 0.0;	D_Mat[2][4] = 0.0;	D_Mat[2][5] = 0.0;
	    D_Mat[3][0] = 0.0;	D_Mat[3][1] = 0.0;	D_Mat[3][2] = 0.0;	D_Mat[3][3] = 0.0;	D_Mat[3][4] = 0.0;	D_Mat[3][5] = 0.0;
	    D_Mat[4][0] = 0.0;	D_Mat[4][1] = 0.0;	D_Mat[4][2] = 0.0;	D_Mat[4][3] = 0.0;	D_Mat[4][4] = 0.0;	D_Mat[4][5] = 0.0;
	    D_Mat[5][0] = 0.0;	D_Mat[5][1] = 0.0;	D_Mat[5][2] = 0.0;	D_Mat[5][3] = 0.0;	D_Mat[5][4] = 0.0;	D_Mat[5][5] = 0.0;
	}


}

void CFEM_DielectricElastomer_Adj::Compute_Stress_Tensor(CElement *element, CConfig *config) {

	unsigned short iVar, iDim, jDim;
	su2double mod_Curr, mod_Ref;

	su2double E00 = 0.0, E11 = 0.0, E22 = 0.0;
	su2double E01 = 0.0, E02 = 0.0, E12 = 0.0;

//	Compute_Eigenproblem(element, config);

	for (iVar = 0; iVar < nElectric_Field; iVar++){
		for (iDim = 0; iDim < nDim; iDim++){
			Electric_Field_Curr[iDim] = 0.0;
			for (jDim = 0; jDim < nDim; jDim++){
				Electric_Field_Curr[iDim] += F_Mat[iDim][jDim] * Electric_Field_Ref_Unit[iVar][jDim];
			}
		}

		mod_Curr = sqrt(pow(Electric_Field_Curr[0],2)+pow(Electric_Field_Curr[1],2));
		mod_Ref = sqrt(pow(Electric_Field_Ref_Unit[iVar][0],2)+pow(Electric_Field_Ref_Unit[iVar][1],2));
//		cout << endl;
//		cout << "E_Ref(" << iVar << ")  = (" << Electric_Field_Ref_Unit[iVar][0] << "," << Electric_Field_Ref_Unit[iVar][1] << ").  |E_Ref|  = " << mod_Ref << "." << endl;
//		cout << "E_Curr(" << iVar << ") = (" << Electric_Field_Curr[0] << "," << Electric_Field_Curr[1] << "). |E_Curr| = " << mod_Curr << "." << endl;
	}

	mod_Curr = 0.0;
	/*--- Normalize the electric field vector ---*/
	for (iDim = 0; iDim < nDim_Electric_Field; iDim++) {
		mod_Curr += Electric_Field_Curr[iDim]*Electric_Field_Curr[iDim];
	}
	mod_Curr = sqrt(mod_Curr);

	/*--- v = Electric_Field_Curr[iDim] ---*/

	E00 = (2*(pow(Electric_Field_Curr[0],2))-pow(mod_Curr,2))*E_mod[0];
	E11 = (2*(pow(Electric_Field_Curr[1],2))-pow(mod_Curr,2))*E_mod[0];
	E01 = 2*Electric_Field_Curr[0]*Electric_Field_Curr[0]*E_mod[0];
	E22 = -1*pow(mod_Curr,2)*E_mod[0];

	if (nDim == 3){
		E02 = 2*Electric_Field_Curr[0]*Electric_Field_Curr[2]*E_mod[0];
		E12 = 2*Electric_Field_Curr[1]*Electric_Field_Curr[2]*E_mod[0];
		E22+=2*(pow(Electric_Field_Curr[1],2))*E_mod[0];
	}

	Stress_Tensor[0][0] = ke_DE*E00;	Stress_Tensor[0][1] = ke_DE*E01;	Stress_Tensor[0][2] = ke_DE*E02;
	Stress_Tensor[1][0] = ke_DE*E01;	Stress_Tensor[1][1] = ke_DE*E11;	Stress_Tensor[1][2] = ke_DE*E12;
	Stress_Tensor[2][0] = ke_DE*E02;	Stress_Tensor[2][1] = ke_DE*E12;	Stress_Tensor[2][2] = ke_DE*E22;

//	cout << endl << "SXX:" << endl;
//	cout << Stress_Tensor[0][0] << " " << Stress_Tensor[0][1] << " " << Stress_Tensor[0][2] << endl;
//	cout << Stress_Tensor[1][0] << " " << Stress_Tensor[1][1] << " " << Stress_Tensor[1][2] << endl;
//	cout << Stress_Tensor[2][0] << " " << Stress_Tensor[2][1] << " " << Stress_Tensor[2][2] << endl;

}
