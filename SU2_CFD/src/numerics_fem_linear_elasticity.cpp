/*!
 * \file numerics_fem_elasticity.cpp
 * \brief This file contains the routines for setting the FEM elastic structural problem.
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

CFEM_LinearElasticity::CFEM_LinearElasticity(unsigned short val_nDim, unsigned short val_nVar,
                                   CConfig *config) : CFEM_Elasticity(val_nDim, val_nVar, config) {


	/*--- If it is linear elasticity, D is constant along the calculations ---*/

	Compute_Constitutive_Matrix();

}

CFEM_LinearElasticity::~CFEM_LinearElasticity(void) {

}

void CFEM_LinearElasticity::Compute_Tangent_Matrix(CElement *element){

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

		/*--- Retrieve the values of the gradients of the shape functions for each node ---*/
		/*--- This avoids repeated operations ---*/
		for (iNode = 0; iNode < nNode; iNode++){
			for (iDim = 0; iDim < nDim; iDim++){
				GradNi_Ref_Mat[iNode][iDim] = element->GetGradNi_X(iNode,iGauss,iDim);
			}
		}

		for (iNode = 0; iNode < nNode; iNode++){

			if (nDim == 2){
				Ba_Mat[0][0] = GradNi_Ref_Mat[iNode][0];
				Ba_Mat[1][1] = GradNi_Ref_Mat[iNode][1];
				Ba_Mat[2][0] = GradNi_Ref_Mat[iNode][1];
				Ba_Mat[2][1] = GradNi_Ref_Mat[iNode][0];
			}
			else if (nDim ==3){
				Ba_Mat[0][0] = GradNi_Ref_Mat[iNode][0];
				Ba_Mat[1][1] = GradNi_Ref_Mat[iNode][1];
				Ba_Mat[2][2] = GradNi_Ref_Mat[iNode][2];
				Ba_Mat[3][0] = GradNi_Ref_Mat[iNode][1];
				Ba_Mat[3][1] = GradNi_Ref_Mat[iNode][0];
				Ba_Mat[4][0] = GradNi_Ref_Mat[iNode][2];
				Ba_Mat[4][2] = GradNi_Ref_Mat[iNode][0];
				Ba_Mat[5][1] = GradNi_Ref_Mat[iNode][2];
				Ba_Mat[5][2] = GradNi_Ref_Mat[iNode][1];					;
			}

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
					Bb_Mat[0][0] = GradNi_Ref_Mat[jNode][0];
					Bb_Mat[1][1] = GradNi_Ref_Mat[jNode][1];
					Bb_Mat[2][0] = GradNi_Ref_Mat[jNode][1];
					Bb_Mat[2][1] = GradNi_Ref_Mat[jNode][0];
				}
				else if (nDim ==3){
					Bb_Mat[0][0] = GradNi_Ref_Mat[iNode][0];
					Bb_Mat[1][1] = GradNi_Ref_Mat[iNode][1];
					Bb_Mat[2][2] = GradNi_Ref_Mat[iNode][2];
					Bb_Mat[3][0] = GradNi_Ref_Mat[iNode][1];
					Bb_Mat[3][1] = GradNi_Ref_Mat[iNode][0];
					Bb_Mat[4][0] = GradNi_Ref_Mat[iNode][2];
					Bb_Mat[4][2] = GradNi_Ref_Mat[iNode][0];
					Bb_Mat[5][1] = GradNi_Ref_Mat[iNode][2];
					Bb_Mat[5][2] = GradNi_Ref_Mat[iNode][1];
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


void CFEM_LinearElasticity::Compute_Constitutive_Matrix(void){

	/*--- Assuming plane strain ---*/

	if (nDim == 2){
		D_Mat[0][0] = Lambda + 2.0*Mu;	D_Mat[0][1] = Lambda;            D_Mat[0][2] = 0.0;
		D_Mat[1][0] = Lambda;           D_Mat[1][1] = Lambda + 2.0*Mu;   D_Mat[1][2] = 0.0;
		D_Mat[2][0] = 0.0;              D_Mat[2][1] = 0.0;               D_Mat[2][2] = Mu;
	}
	else if (nDim == 3){

	    D_Mat[0][0] = Lambda + 2.0*Mu;	D_Mat[0][1] = Lambda;			D_Mat[0][2] = Lambda;			D_Mat[0][3] = 0.0;	D_Mat[0][4] = 0.0;	D_Mat[0][5] = 0.0;
	    D_Mat[1][0] = Lambda;			D_Mat[1][1] = Lambda + 2.0*Mu;	D_Mat[1][2] = Lambda;			D_Mat[1][3] = 0.0;	D_Mat[1][4] = 0.0;	D_Mat[1][5] = 0.0;
	    D_Mat[2][0] = Lambda;			D_Mat[2][1] = Lambda;			D_Mat[2][2] = Lambda + 2.0*Mu;	D_Mat[2][3] = 0.0;	D_Mat[2][4] = 0.0;	D_Mat[2][5] = 0.0;
	    D_Mat[3][0] = 0.0;				D_Mat[3][1] = 0.0;				D_Mat[3][2] = 0.0;				D_Mat[3][3] = Mu;	D_Mat[3][4] = 0.0;	D_Mat[3][5] = 0.0;
	    D_Mat[4][0] = 0.0;				D_Mat[4][1] = 0.0;				D_Mat[4][2] = 0.0;				D_Mat[4][3] = 0.0;	D_Mat[4][4] = Mu;	D_Mat[4][5] = 0.0;
	    D_Mat[5][0] = 0.0;				D_Mat[5][1] = 0.0;				D_Mat[5][2] = 0.0;				D_Mat[5][3] = 0.0;	D_Mat[5][4] = 0.0;	D_Mat[5][5] = Mu;

	}

}
