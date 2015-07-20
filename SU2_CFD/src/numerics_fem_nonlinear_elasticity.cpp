/*!
 * \file numerics_fem_nonlinear_elasticity.cpp
 * \brief This file contains the routines for setting the tangent matrix and residual of a FEM nonlinear elastic structural problem.
 * \author R. Sanchez
 * \version 4.0.0 "Cardinal"
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

CFEM_NonlinearElasticity::CFEM_NonlinearElasticity(unsigned short val_nDim, unsigned short val_nVar,
                                   CConfig *config) : CFEM_Elasticity(val_nDim, val_nVar, config) {

	unsigned short i;

	F_Mat = new double *[3];
	b_Mat = new double *[3];
	Stress_Tensor = new double *[3];
	for (i = 0; i < 3; i++){
		F_Mat[i] = new double [3];
		b_Mat[i] = new double [3];
		Stress_Tensor[i] = new double [3];
	}

	if (nDim == 2){
		currentCoord = new double* [4];	/*--- As of now, 4 is the maximum number of nodes for 2D problems ---*/
		for (i = 0; i < 4; i++) currentCoord[i] = new double[nDim];
	}
	else if (nDim == 3){
		currentCoord = new double* [8];	/*--- As of now, 4 is the maximum number of nodes for 3D problems ---*/
		for (i = 0; i < 8; i++) currentCoord[i] = new double[nDim];
	}

	J_F = 0.0;


}

CFEM_NonlinearElasticity::~CFEM_NonlinearElasticity(void) {

}


void CFEM_NonlinearElasticity::Compute_Tangent_Matrix(CElement *element){

	unsigned short i, j, k;
	unsigned short iGauss, nGauss;
	unsigned short iNode, jNode, nNode;
	unsigned short iDim;
	unsigned short bDim;

	double Weight, Jac_X;

	double AuxMatrix[6][3];

	/*--- Initialize auxiliary matrices ---*/

	if (nDim == 2) bDim = 3;
	else if (nDim == 3) bDim = 6;

	for (i = 0; i < bDim; i++){
		for (j = 0; j < nDim; j++){
			Ba_Mat[i][j] = 0.0;
			Bb_Mat[i][j] = 0.0;
		}
	}

	for (i = 0; i < 6; i++){
		for (j = 0; j < 3; j++){
			AuxMatrix[i][j] = 0.0;
			AuxMatrix[i][j] = 0.0;
		}
	}

	element->clearElement(); 			/*--- Restarts the element: avoids adding over previous results in other elements --*/
	element->ComputeGrad_Linear();
	nNode = element->GetnNodes();
	nGauss = element->GetnGaussPoints();

	for (iGauss = 0; iGauss < nGauss; iGauss++){

		Weight = element->GetWeight(iGauss);
		Jac_X = element->GetJ_X(iGauss);

		/*--- Initialize the deformation gradient for each Gauss Point ---*/

		for (i = 0; i < nDim; i++){
			for (j = 0; j < nDim; j++){
				F_Mat[i][j] = 0.0;
				b_Mat[i][j] = 0.0;
			}
		}

		/*--- Retrieve the values of the gradients of the shape functions for each node ---*/
		/*--- This avoids repeated operations ---*/

		for (iNode = 0; iNode < nNode; iNode++){

			for (iDim = 0; iDim < nDim; iDim++){
				GradNi_Mat[iNode][iDim] = element->GetGradNi_X(iNode,iGauss,iDim);
				currentCoord[iNode][iDim] = element->GetCurr_Coord(iNode, iDim);
			}

			/*--- Compute the deformation gradient ---*/

			for (i = 0; i < nDim; i++){
				for (j = 0; j < nDim; j++){
					F_Mat[i][j] += currentCoord[iNode][i]*GradNi_Mat[iNode][j];
				}
			}

			/*--- This implies plane strain --> Consider the possible implementation for plane stress --*/
			if (nDim == 2){
				F_Mat[3][3] = 1.0;
			}

		}

		/*--- Determinant of F --> Jacobian of the transformation ---*/

		J_F = 	F_Mat[0][0]*F_Mat[1][1]*F_Mat[2][2]+
				F_Mat[0][1]*F_Mat[1][2]*F_Mat[2][0]+
				F_Mat[0][2]*F_Mat[1][0]*F_Mat[2][1]-
				F_Mat[0][2]*F_Mat[1][1]*F_Mat[2][0]-
				F_Mat[1][2]*F_Mat[2][1]*F_Mat[0][0]-
				F_Mat[2][2]*F_Mat[0][1]*F_Mat[1][0];

		/*--- Compute the left Cauchy deformation tensor ---*/

		for (i = 0; i < 3; i++){
			for (j = 0; j < 3; j++){
				for (k = 0; k < 3; k++){
					b_Mat[i][j] += F_Mat[i][k]*F_Mat[j][k];
				}
			}
		}


		for (iNode = 0; iNode < nNode; iNode++){

			if (nDim == 2){
				Ba_Mat[0][0] = GradNi_Mat[iNode][0];
				Ba_Mat[1][1] = GradNi_Mat[iNode][1];
				Ba_Mat[2][0] = GradNi_Mat[iNode][1];
				Ba_Mat[2][1] = GradNi_Mat[iNode][0];
			}
			else if (nDim ==3){
				Ba_Mat[0][0] = GradNi_Mat[iNode][0];
				Ba_Mat[1][1] = GradNi_Mat[iNode][1];
				Ba_Mat[2][2] = GradNi_Mat[iNode][2];
				Ba_Mat[3][0] = GradNi_Mat[iNode][1];
				Ba_Mat[3][1] = GradNi_Mat[iNode][0];
				Ba_Mat[4][0] = GradNi_Mat[iNode][2];
				Ba_Mat[4][2] = GradNi_Mat[iNode][0];
				Ba_Mat[5][1] = GradNi_Mat[iNode][2];
				Ba_Mat[5][2] = GradNi_Mat[iNode][1];					;
			}


			/*--- Compute the constitutive matrix ---*/

			Compute_Constitutive_Matrix();

		    /*--- Compute the BT.D Matrix ---*/

			for (i = 0; i < nDim; i++){
				for (j = 0; j < bDim; j++){
					AuxMatrix[i][j] = 0.0;
					for (k = 0; k < bDim; k++){
						AuxMatrix[i][j] += Ba_Mat[k][i]*D_Mat[k][j];
					}
				}
			}

			/*--- Assumming symmetry ---*/
			for (jNode = iNode; jNode < nNode; jNode++){
				if (nDim == 2){
					Bb_Mat[0][0] = GradNi_Mat[jNode][0];
					Bb_Mat[1][1] = GradNi_Mat[jNode][1];
					Bb_Mat[2][0] = GradNi_Mat[jNode][1];
					Bb_Mat[2][1] = GradNi_Mat[jNode][0];
				}
				else if (nDim ==3){
					Bb_Mat[0][0] = GradNi_Mat[iNode][0];
					Bb_Mat[1][1] = GradNi_Mat[iNode][1];
					Bb_Mat[2][2] = GradNi_Mat[iNode][2];
					Bb_Mat[3][0] = GradNi_Mat[iNode][1];
					Bb_Mat[3][1] = GradNi_Mat[iNode][0];
					Bb_Mat[4][0] = GradNi_Mat[iNode][2];
					Bb_Mat[4][2] = GradNi_Mat[iNode][0];
					Bb_Mat[5][1] = GradNi_Mat[iNode][2];
					Bb_Mat[5][2] = GradNi_Mat[iNode][1];
				}

				for (i = 0; i < nDim; i++){
					for (j = 0; j < nDim; j++){
						KAux_ab[i][j] = 0.0;
						for (k = 0; k < bDim; k++){
							KAux_ab[i][j] += Weight * AuxMatrix[i][k] * Bb_Mat[k][j] * Jac_X;
						}
					}
				}

				element->Add_Kab(KAux_ab,iNode, jNode);
				/*--- Symmetric terms --*/
				if (iNode != jNode){
					element->Add_Kab_T(KAux_ab, jNode, iNode);
				}

			}

		}

	}

}

CFEM_NeoHookean_Comp::CFEM_NeoHookean_Comp(unsigned short val_nDim, unsigned short val_nVar,
                                   CConfig *config) : CFEM_NonlinearElasticity(val_nDim, val_nVar, config) {


}

CFEM_NeoHookean_Comp::~CFEM_NeoHookean_Comp(void) {

}

void CFEM_NeoHookean_Comp::Compute_Constitutive_Matrix(void) {

	unsigned short i,j;
	double Mu_p, Lambda_p;
	double dij;

	/*--- This can be done in a better way ---*/
	if (J_F != 0.0){
		Mu_p = (Mu - Lambda*log(J_F))/J_F;
		Lambda_p = Lambda/J_F;
	}

    D_Mat[0][0] = Lambda_p + 2 * Mu_p;	D_Mat[0][1] = Lambda_p;				D_Mat[0][2] = Lambda_p;				D_Mat[0][3] = 0.0;	D_Mat[0][4] = 0.0;	D_Mat[0][5] = 0.0;
    D_Mat[1][0] = Lambda_p;				D_Mat[1][1] = Lambda_p + 2 * Mu_p;	D_Mat[1][2] = Lambda_p;				D_Mat[1][3] = 0.0;	D_Mat[1][4] = 0.0;	D_Mat[1][5] = 0.0;
    D_Mat[2][0] = Lambda_p;				D_Mat[2][1] = Lambda_p;				D_Mat[2][2] = Lambda_p + 2 * Mu_p;	D_Mat[2][3] = 0.0;	D_Mat[2][4] = 0.0;	D_Mat[2][5] = 0.0;
    D_Mat[3][0] = 0.0;					D_Mat[3][1] = 0.0;					D_Mat[3][2] = 0.0;					D_Mat[3][3] = Mu_p;	D_Mat[3][4] = 0.0;	D_Mat[3][5] = 0.0;
    D_Mat[4][0] = 0.0;					D_Mat[4][1] = 0.0;					D_Mat[4][2] = 0.0;					D_Mat[4][3] = 0.0;	D_Mat[4][4] = Mu_p;	D_Mat[4][5] = 0.0;
    D_Mat[5][0] = 0.0;					D_Mat[5][1] = 0.0;					D_Mat[5][2] = 0.0;					D_Mat[5][3] = 0.0;	D_Mat[5][4] = 0.0;	D_Mat[5][5] = Mu_p;



}

void CFEM_NeoHookean_Comp::Compute_Stress_Tensor(void) {

	unsigned short i,j;
	double Mu_J, Lambda_J;
	double dij;

	/*--- This can be done in a better way ---*/
	if (J_F != 0.0){
		Mu_J = Mu/J_F;
		Lambda_J = Lambda/J_F;
	}

	for (i = 0; i < 3; i++){
		for (j = 0; j < 3; j++){
			if (i == j) dij = 1.0;
			else if (i != j) dij = 0.0;
			Stress_Tensor[i][j] = Mu_J * (b_Mat[i][j] - dij) + Lambda_J * log(J_F) * dij;
		}
	}


}



