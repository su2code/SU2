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

	KAux_t_a = new double [nDim];

	KAux_P_ab = new double* [nDim];
	for (i = 0; i < nDim; i++) {
		KAux_P_ab[i] = new double[nDim];
	}

	if (nDim == 2){
		currentCoord = new double* [4];	/*--- As of now, 4 is the maximum number of nodes for 2D problems ---*/
		for (i = 0; i < 4; i++) currentCoord[i] = new double[nDim];
	}
	else if (nDim == 3){
		currentCoord = new double* [8];	/*--- As of now, 8 is the maximum number of nodes for 3D problems ---*/
		for (i = 0; i < 8; i++) currentCoord[i] = new double[nDim];
	}

	J_F = 0.0;


}

CFEM_NonlinearElasticity::~CFEM_NonlinearElasticity(void) {

	unsigned short iVar, jVar;

	for (iVar = 0; iVar < 3; iVar++){
		delete [] F_Mat[iVar];
		delete [] b_Mat[iVar];
		delete [] Stress_Tensor[iVar];
	}

	for (iVar = 0; iVar < nDim; iVar++){
		delete [] KAux_P_ab[iVar];
	}

	if (nDim == 2){
		for (iVar = 0; iVar < 4; iVar++){
			delete [] currentCoord[iVar];
		}
	}
	else if (nDim == 3){
		for (iVar = 0; iVar < 8; iVar++){
			delete [] currentCoord[iVar];
		}
	}

	delete [] F_Mat;
	delete [] b_Mat;
	delete [] Stress_Tensor;
	delete [] KAux_t_a;
	delete [] KAux_P_ab;
	delete [] currentCoord;

}


void CFEM_NonlinearElasticity::Compute_Tangent_Matrix(CElement *element){

	unsigned short i, j, k;
	unsigned short iGauss, nGauss;
	unsigned short iNode, jNode, nNode;
	unsigned short iDim, bDim;

	double Ks_Aux_ab;

	double Weight, Jac_X, Jac_x;

	double AuxMatrixKc[6][3];
	double AuxMatrixKs[3];

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
			AuxMatrixKc[i][j] = 0.0;
		}
	}

	for (i = 0; i < 3; i++){
		AuxMatrixKs[i] = 0.0;
	}

	element->clearElement(); 			/*--- Restarts the element: avoids adding over previous results in other elements --*/
	element->ComputeGrad_NonLinear();

	nNode = element->GetnNodes();
	nGauss = element->GetnGaussPoints();

	/*--- Full integration of the constitutive and stress term ---*/

	for (iGauss = 0; iGauss < nGauss; iGauss++){

		Weight = element->GetWeight(iGauss);
		Jac_X = element->GetJ_X(iGauss);
		Jac_x = element->GetJ_x(iGauss);

		/*--- Initialize the deformation gradient for each Gauss Point ---*/

		for (i = 0; i < 3; i++){
			for (j = 0; j < 3; j++){
				F_Mat[i][j] = 0.0;
				b_Mat[i][j] = 0.0;
			}
		}

		/*--- Retrieve the values of the gradients of the shape functions for each node ---*/
		/*--- This avoids repeated operations ---*/

		for (iNode = 0; iNode < nNode; iNode++){

			for (iDim = 0; iDim < nDim; iDim++){
				GradNi_Ref_Mat[iNode][iDim] = element->GetGradNi_X(iNode,iGauss,iDim);
				GradNi_Curr_Mat[iNode][iDim] = element->GetGradNi_x(iNode,iGauss,iDim);
				currentCoord[iNode][iDim] = element->GetCurr_Coord(iNode, iDim);
			}

			/*--- Compute the deformation gradient ---*/

			for (i = 0; i < nDim; i++){
				for (j = 0; j < nDim; j++){
					F_Mat[i][j] += currentCoord[iNode][i]*GradNi_Ref_Mat[iNode][j];
				}
			}

			/*--- This implies plane strain --> Consider the possible implementation for plane stress --*/
			if (nDim == 2){
				F_Mat[2][2] = 1.0;
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

		/*--- Compute the constitutive matrix ---*/

		Compute_Constitutive_Matrix(element);
		Compute_Stress_Tensor(element);


		for (iNode = 0; iNode < nNode; iNode++){

			/*--------------------------------------------------------------------------------*/
			/*---------------------------- NODAL STRESS TERM ---------------------------------*/
			/*--------------------------------------------------------------------------------*/
		    /*--- Compute the nodal stress term for each gaussian point and for each node, ---*/
		    /*--- and add it to the element structure to be retrieved from the solver      ---*/

			for (i = 0; i < nDim; i++){
				KAux_t_a[i] = 0.0;
				for (j = 0; j < nDim; j++){
					KAux_t_a[i] += Weight * Stress_Tensor[i][j] * GradNi_Curr_Mat[iNode][j] * Jac_x;
				}
			}

			element->Add_Kt_a(KAux_t_a, iNode);

			/*--------------------------------------------------------------------------------*/
			/*----------------------- CONSTITUTIVE AND STRESS TERM ---------------------------*/
			/*--------------------------------------------------------------------------------*/

			if (nDim == 2){
				Ba_Mat[0][0] = GradNi_Curr_Mat[iNode][0];
				Ba_Mat[1][1] = GradNi_Curr_Mat[iNode][1];
				Ba_Mat[2][0] = GradNi_Curr_Mat[iNode][1];
				Ba_Mat[2][1] = GradNi_Curr_Mat[iNode][0];
			}
			else if (nDim ==3){
				Ba_Mat[0][0] = GradNi_Curr_Mat[iNode][0];
				Ba_Mat[1][1] = GradNi_Curr_Mat[iNode][1];
				Ba_Mat[2][2] = GradNi_Curr_Mat[iNode][2];
				Ba_Mat[3][0] = GradNi_Curr_Mat[iNode][1];
				Ba_Mat[3][1] = GradNi_Curr_Mat[iNode][0];
				Ba_Mat[4][0] = GradNi_Curr_Mat[iNode][2];
				Ba_Mat[4][2] = GradNi_Curr_Mat[iNode][0];
				Ba_Mat[5][1] = GradNi_Curr_Mat[iNode][2];
				Ba_Mat[5][2] = GradNi_Curr_Mat[iNode][1];
			}

		    /*--- Compute the BT.D Matrix ---*/

			for (i = 0; i < nDim; i++){
				for (j = 0; j < bDim; j++){
					AuxMatrixKc[i][j] = 0.0;
					for (k = 0; k < bDim; k++){
						AuxMatrixKc[i][j] += Ba_Mat[k][i]*D_Mat[k][j];
					}
				}
			}

		    /*--- Compute the BT.D Matrix ---*/

			for (i = 0; i < nDim; i++){
				AuxMatrixKs[i] = 0.0;
				for (j = 0; j < nDim; j++){
					AuxMatrixKs[i] += GradNi_Curr_Mat[iNode][j]*Stress_Tensor[j][i]; // DOUBLE CHECK
				}
			}

			/*--- Assumming symmetry ---*/
			for (jNode = iNode; jNode < nNode; jNode++){
				if (nDim == 2){
					Bb_Mat[0][0] = GradNi_Curr_Mat[jNode][0];
					Bb_Mat[1][1] = GradNi_Curr_Mat[jNode][1];
					Bb_Mat[2][0] = GradNi_Curr_Mat[jNode][1];
					Bb_Mat[2][1] = GradNi_Curr_Mat[jNode][0];
				}
				else if (nDim ==3){
					Bb_Mat[0][0] = GradNi_Curr_Mat[jNode][0];
					Bb_Mat[1][1] = GradNi_Curr_Mat[jNode][1];
					Bb_Mat[2][2] = GradNi_Curr_Mat[jNode][2];
					Bb_Mat[3][0] = GradNi_Curr_Mat[jNode][1];
					Bb_Mat[3][1] = GradNi_Curr_Mat[jNode][0];
					Bb_Mat[4][0] = GradNi_Curr_Mat[jNode][2];
					Bb_Mat[4][2] = GradNi_Curr_Mat[jNode][0];
					Bb_Mat[5][1] = GradNi_Curr_Mat[jNode][2];
					Bb_Mat[5][2] = GradNi_Curr_Mat[jNode][1];
				}

				/*--- KAux_ab is the term for the constitutive part of the tangent matrix ---*/
				for (i = 0; i < nDim; i++){
					for (j = 0; j < nDim; j++){
						KAux_ab[i][j] = 0.0;
						for (k = 0; k < bDim; k++){
							KAux_ab[i][j] += Weight * AuxMatrixKc[i][k] * Bb_Mat[k][j] * Jac_x;
						}
					}
				}

				/*--- Ks_Aux_ab is the term for the constitutive part of the tangent matrix ---*/
				Ks_Aux_ab = 0.0;
				for (i = 0; i < nDim; i++){
					Ks_Aux_ab += Weight * AuxMatrixKs[i] * GradNi_Curr_Mat[jNode][i] * Jac_x;
				}

				element->Add_Kab(KAux_ab,iNode, jNode);
				element->Add_Ks_ab(Ks_Aux_ab,iNode, jNode);
				/*--- Symmetric terms --*/
				if (iNode != jNode){
					element->Add_Kab_T(KAux_ab, jNode, iNode);
					element->Add_Ks_ab(Ks_Aux_ab,jNode, iNode);
				}

			}

		}

	}

}

void CFEM_NonlinearElasticity::Compute_MeanDilatation_Term(CElement *element){

	unsigned short i, j, k;
	unsigned short iGauss, nGauss;
	unsigned short iNode, jNode, nNode;
	double Weight, Jac_X, Jac_x;
	unsigned short iDim ;

	double GradNi_Mat_Term;
	double Vol_current, Vol_reference;
	double Avg_kappa;
	double el_Pressure;


	/*--- Under integration of the pressure term, if the calculations assume incompressibility or near incompressibility ---*/

	element->ComputeGrad_Pressure(); // Check if we can take this out!

	/*--- nGauss is here the number of Gaussian Points for the pressure term ---*/
	nGauss = element->GetnGaussPointsP();
	nNode = element->GetnNodes();

	/*--- Initialize the Gradient auxiliary Matrix ---*/
	for (iNode = 0; iNode < nNode; iNode++){
		for (iDim = 0; iDim < nDim; iDim++){
			GradNi_Curr_Mat[iNode][iDim] = 0.0;
		}
	}

	Vol_current = 0.0;
	Vol_reference = 0.0;

	/*--------------------------------------------------------------------------------*/
	/*-------------------------- INCOMPRESSIBLE TERM ---------------------------------*/
	/*--------------------------------------------------------------------------------*/

	for (iGauss = 0; iGauss < nGauss; iGauss++){

		Weight = element->GetWeight_P(iGauss);
		Jac_X = element->GetJ_X_P(iGauss);
		Jac_x = element->GetJ_x_P(iGauss);

		/*--- Retrieve the values of the gradients of the shape functions for each node ---*/
		/*--- This avoids repeated operations ---*/

		/*--- We compute the average gradient ---*/
		for (iNode = 0; iNode < nNode; iNode++){
			for (iDim = 0; iDim < nDim; iDim++){
				GradNi_Mat_Term = element->GetGradNi_x_P(iNode,iGauss,iDim);
				GradNi_Curr_Mat[iNode][iDim] += Weight * GradNi_Mat_Term * Jac_x;
			}
		}

		Vol_reference += Weight * Jac_X;
		Vol_current += Weight * Jac_x;

	}

	if ((Vol_current != 0.0) && (Vol_reference != 0.0)) {

		/*--- It is necessary to divide over the current volume to obtain the averaged gradients ---*/
		/*--- TODO: Check this operation and add exit if the volumes are 0. ---*/
		for (iNode = 0; iNode < nNode; iNode++){
			for (iDim = 0; iDim < nDim; iDim++){
				GradNi_Curr_Mat[iNode][iDim] = GradNi_Curr_Mat[iNode][iDim] / Vol_current;
			}
		}

		Avg_kappa = Kappa * Vol_current / Vol_reference;

		el_Pressure = Kappa * ((Vol_current / Vol_reference) - 1);

		element->SetElement_Pressure(el_Pressure);

//		cout << "ELEMENT PRESSURE: " << el_Pressure << endl;

	}

	for (iNode = 0; iNode < nNode; iNode++){

		for (jNode = 0; jNode < nNode; jNode++){

			/*--- KAux_P_ab is the term for the incompressibility part of the tangent matrix ---*/
			for (i = 0; i < nDim; i++){
				for (j = 0; j < nDim; j++){
					KAux_P_ab[i][j] = Avg_kappa * Vol_current * GradNi_Curr_Mat[iNode][i] * GradNi_Curr_Mat[jNode][j];
				}
			}

			element->Set_Kk_ab(KAux_P_ab,iNode, jNode);

		}

	}

}


void CFEM_NonlinearElasticity::Compute_NodalStress_Term(CElement *element){

	unsigned short i, j, k;
	unsigned short iGauss, nGauss;
	unsigned short iNode, jNode, nNode;
	unsigned short iDim, bDim;

	double Ks_Aux_ab;

	double Weight, Jac_X, Jac_x;

	double AuxMatrixKt[3];

	/*--- Initialize auxiliary matrices ---*/

	for (i = 0; i < 3; i++){
		AuxMatrixKt[i] = 0.0;
	}

	element->clearElement(); 			/*--- Restarts the element: avoids adding over previous results in other elements --*/
	element->ComputeGrad_NonLinear();	/*--- Check if we can take this out... so we don't have to do it twice ---*/

	nNode = element->GetnNodes();
	nGauss = element->GetnGaussPoints();

	/*--- Full integration of the nodal stress ---*/

	for (iGauss = 0; iGauss < nGauss; iGauss++){

		Weight = element->GetWeight(iGauss);
		Jac_X = element->GetJ_X(iGauss);
		Jac_x = element->GetJ_x(iGauss);

		/*--- Initialize the deformation gradient for each Gauss Point ---*/

		for (i = 0; i < 3; i++){
			for (j = 0; j < 3; j++){
				F_Mat[i][j] = 0.0;
				b_Mat[i][j] = 0.0;
			}
		}

		/*--- Retrieve the values of the gradients of the shape functions for each node ---*/
		/*--- This avoids repeated operations ---*/

		for (iNode = 0; iNode < nNode; iNode++){

			for (iDim = 0; iDim < nDim; iDim++){
				GradNi_Ref_Mat[iNode][iDim] = element->GetGradNi_X(iNode,iGauss,iDim);
				GradNi_Curr_Mat[iNode][iDim] = element->GetGradNi_x(iNode,iGauss,iDim);
				currentCoord[iNode][iDim] = element->GetCurr_Coord(iNode, iDim);
			}

			/*--- Compute the deformation gradient ---*/

			for (i = 0; i < nDim; i++){
				for (j = 0; j < nDim; j++){
					F_Mat[i][j] += currentCoord[iNode][i]*GradNi_Ref_Mat[iNode][j];
				}
			}

			/*--- This implies plane strain --> Consider the possible implementation for plane stress --*/
			if (nDim == 2){
				F_Mat[2][2] = 1.0;
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

		/*--- Compute the stress tensor ---*/

		Compute_Stress_Tensor(element);


		for (iNode = 0; iNode < nNode; iNode++){

		    /*--- Compute the nodal stress term for each gaussian point and for each node, ---*/
		    /*--- and add it to the element structure to be retrieved from the solver      ---*/

			for (i = 0; i < nDim; i++){
				KAux_t_a[i] = 0.0;
				for (j = 0; j < nDim; j++){
					KAux_t_a[i] += Weight * Stress_Tensor[i][j] * GradNi_Curr_Mat[iNode][j] * Jac_x;
				}
			}

			element->Add_Kt_a(KAux_t_a, iNode);

		}

	}

}

void CFEM_NonlinearElasticity::Compute_Averaged_NodalStress(CElement *element){

	unsigned short i, j, k;
	unsigned short iGauss, nGauss;
	unsigned short iNode, jNode, nNode;
	unsigned short iDim, bDim;

	double Ks_Aux_ab;

	double Weight, Jac_X, Jac_x;

	element->clearStress();
	element->ComputeGrad_Linear();

	nNode = element->GetnNodes();
	nGauss = element->GetnGaussPoints();

	/*--- Computation of the deformation gradient ---*/

	for (iGauss = 0; iGauss < nGauss; iGauss++){

		Weight = element->GetWeight(iGauss);
		Jac_X = element->GetJ_X(iGauss);
		Jac_x = element->GetJ_x(iGauss);

		/*--- Initialize the deformation gradient for each Gauss Point ---*/

		for (i = 0; i < 3; i++){
			for (j = 0; j < 3; j++){
				F_Mat[i][j] = 0.0;
				b_Mat[i][j] = 0.0;
			}
		}

		/*--- Retrieve the values of the gradients of the shape functions for each node ---*/
		/*--- This avoids repeated operations ---*/

		for (iNode = 0; iNode < nNode; iNode++){

			for (iDim = 0; iDim < nDim; iDim++){
				GradNi_Ref_Mat[iNode][iDim] = element->GetGradNi_X(iNode,iGauss,iDim);
				currentCoord[iNode][iDim] = element->GetCurr_Coord(iNode, iDim);
			}

			/*--- Compute the deformation gradient ---*/

			for (i = 0; i < nDim; i++){
				for (j = 0; j < nDim; j++){
					F_Mat[i][j] += currentCoord[iNode][i]*GradNi_Ref_Mat[iNode][j];
				}
			}

			/*--- This implies plane strain --> Consider the possible implementation for plane stress --*/
			if (nDim == 2){
				F_Mat[2][2] = 1.0;
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

		/*--- Compute the stress tensor ---*/

		Compute_Stress_Tensor(element);

		for (iNode = 0; iNode < nNode; iNode++){
			element->Add_NodalStress(Stress_Tensor[0][0] * element->GetNi_Extrap(iNode, iGauss), iNode, 0);
			element->Add_NodalStress(Stress_Tensor[1][1] * element->GetNi_Extrap(iNode, iGauss), iNode, 1);
			element->Add_NodalStress(Stress_Tensor[0][1] * element->GetNi_Extrap(iNode, iGauss), iNode, 2);
			if (nDim == 3){
				element->Add_NodalStress(Stress_Tensor[2][2] * element->GetNi_Extrap(iNode, iGauss), iNode, 3);
				element->Add_NodalStress(Stress_Tensor[0][2] * element->GetNi_Extrap(iNode, iGauss), iNode, 4);
				element->Add_NodalStress(Stress_Tensor[1][2] * element->GetNi_Extrap(iNode, iGauss), iNode, 5);
			}
		}

	}


}


CFEM_NeoHookean_Comp::CFEM_NeoHookean_Comp(unsigned short val_nDim, unsigned short val_nVar,
                                   CConfig *config) : CFEM_NonlinearElasticity(val_nDim, val_nVar, config) {


}

CFEM_NeoHookean_Comp::~CFEM_NeoHookean_Comp(void) {

}

void CFEM_NeoHookean_Comp::Compute_Constitutive_Matrix(CElement *element) {

	unsigned short i,j;
	double Mu_p, Lambda_p;
	double dij;

	/*--- This can be done in a better way ---*/
	if (J_F != 0.0){
		Mu_p = (Mu - Lambda*log(J_F))/J_F;
		Lambda_p = Lambda/J_F;
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

void CFEM_NeoHookean_Comp::Compute_Stress_Tensor(CElement *element) {

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

CFEM_NeoHookean_Incomp::CFEM_NeoHookean_Incomp(unsigned short val_nDim, unsigned short val_nVar,
                                   CConfig *config) : CFEM_NonlinearElasticity(val_nDim, val_nVar, config) {


}

CFEM_NeoHookean_Incomp::~CFEM_NeoHookean_Incomp(void) {

}

void CFEM_NeoHookean_Incomp::Compute_Constitutive_Matrix(CElement *element) {

	unsigned short i,j;
	double dij, el_P;
	double Ib = 0.0, Jft;

	/*--- First invariant of b -> Ib = tr(b) ---*/
	for (i = 0; i < 3; i++){
		Ib += b_Mat[i][i];
	}

	/*--- Retrieve element pressure ---*/

	el_P = element->GetElement_Pressure();

	/*--- J^(-5/3) ---*/
	Jft = pow(J_F, -1.666666666666667);

	if (nDim == 2){

		/*--- Diagonal terms ---*/
		D_Mat[0][0] = 2.0 * Mu * Jft * ((4.0 / 9.0) * Ib - (2.0 / 3.0) * b_Mat[0][0]) - el_P;
		D_Mat[1][1] = 2.0 * Mu * Jft * ((4.0 / 9.0) * Ib - (2.0 / 3.0) * b_Mat[1][1]) - el_P;

		D_Mat[2][2] = (1.0 / 3.0) * Mu * Jft * Ib - el_P;

		/*--- Q2 off diagonal terms (and symmetric) ---*/

		D_Mat[0][1] = (-2.0 / 3.0) * Mu * Jft * b_Mat[0][1];
		D_Mat[1][0] = (-2.0 / 3.0) * Mu * Jft * b_Mat[0][1];

		D_Mat[0][2] = 0.0;
		D_Mat[2][0] = 0.0;

	}
	else if (nDim == 3){

		/*--- Diagonal terms ---*/
		D_Mat[0][0] = 2.0 * Mu * Jft * ((4.0 / 9.0) * Ib - (2.0 / 3.0) * b_Mat[0][0]) - el_P;
		D_Mat[1][1] = 2.0 * Mu * Jft * ((4.0 / 9.0) * Ib - (2.0 / 3.0) * b_Mat[1][1]) - el_P;
		D_Mat[2][2] = 2.0 * Mu * Jft * ((4.0 / 9.0) * Ib - (2.0 / 3.0) * b_Mat[2][2]) - el_P;

		D_Mat[3][3] = (1.0 / 3.0) * Mu * Jft * Ib - el_P;
		D_Mat[4][4] = (1.0 / 3.0) * Mu * Jft * Ib - el_P;
		D_Mat[5][5] = (1.0 / 3.0) * Mu * Jft * Ib - el_P;

		/*--- Q1 off diagonal terms (and symmetric) ---*/

		D_Mat[0][1] = 2.0 * Mu * Jft * ((1.0 / 9.0) * Ib - (1.0 / 3.0) * b_Mat[0][0] - (1.0 / 3.0) * b_Mat[1][1]) + el_P;
		D_Mat[0][2] = 2.0 * Mu * Jft * ((1.0 / 9.0) * Ib - (1.0 / 3.0) * b_Mat[0][0] - (1.0 / 3.0) * b_Mat[2][2]) + el_P;
		D_Mat[1][2] = 2.0 * Mu * Jft * ((1.0 / 9.0) * Ib - (1.0 / 3.0) * b_Mat[1][1] - (1.0 / 3.0) * b_Mat[2][2]) + el_P;

		D_Mat[1][0] = 2.0 * Mu * Jft * ((1.0 / 9.0) * Ib - (1.0 / 3.0) * b_Mat[0][0] - (1.0 / 3.0) * b_Mat[1][1]) + el_P;
		D_Mat[2][0] = 2.0 * Mu * Jft * ((1.0 / 9.0) * Ib - (1.0 / 3.0) * b_Mat[0][0] - (1.0 / 3.0) * b_Mat[2][2]) + el_P;
		D_Mat[2][1] = 2.0 * Mu * Jft * ((1.0 / 9.0) * Ib - (1.0 / 3.0) * b_Mat[1][1] - (1.0 / 3.0) * b_Mat[2][2]) + el_P;

		/*--- Q2 off diagonal terms (and symmetric) ---*/

		D_Mat[0][3] = (-2.0 / 3.0) * Mu * Jft * b_Mat[0][1];
		D_Mat[1][3] = (-2.0 / 3.0) * Mu * Jft * b_Mat[0][1];
		D_Mat[2][3] = (-2.0 / 3.0) * Mu * Jft * b_Mat[0][1];

		D_Mat[0][4] = (-2.0 / 3.0) * Mu * Jft * b_Mat[0][2];
		D_Mat[1][4] = (-2.0 / 3.0) * Mu * Jft * b_Mat[0][2];
		D_Mat[2][4] = (-2.0 / 3.0) * Mu * Jft * b_Mat[0][2];

		D_Mat[0][5] = (-2.0 / 3.0) * Mu * Jft * b_Mat[1][2];
		D_Mat[1][5] = (-2.0 / 3.0) * Mu * Jft * b_Mat[1][2];
		D_Mat[2][5] = (-2.0 / 3.0) * Mu * Jft * b_Mat[1][2];


		D_Mat[3][0] = (-2.0 / 3.0) * Mu * Jft * b_Mat[0][1];
		D_Mat[3][1] = (-2.0 / 3.0) * Mu * Jft * b_Mat[0][1];
		D_Mat[3][2] = (-2.0 / 3.0) * Mu * Jft * b_Mat[0][1];

		D_Mat[4][0] = (-2.0 / 3.0) * Mu * Jft * b_Mat[0][2];
		D_Mat[4][1] = (-2.0 / 3.0) * Mu * Jft * b_Mat[0][2];
		D_Mat[4][2] = (-2.0 / 3.0) * Mu * Jft * b_Mat[0][2];

		D_Mat[5][0] = (-2.0 / 3.0) * Mu * Jft * b_Mat[1][2];
		D_Mat[5][1] = (-2.0 / 3.0) * Mu * Jft * b_Mat[1][2];
		D_Mat[5][2] = (-2.0 / 3.0) * Mu * Jft * b_Mat[1][2];

		/*--- Q3 off diagonal terms (and symmetric) are all zeros ---*/

		D_Mat[3][4] = 0.0;
		D_Mat[3][5] = 0.0;
		D_Mat[4][5] = 0.0;

		D_Mat[4][3] = 0.0;
		D_Mat[5][3] = 0.0;
		D_Mat[5][4] = 0.0;

	}

}

void CFEM_NeoHookean_Incomp::Compute_Stress_Tensor(CElement *element) {

	unsigned short i,j;
	double dij, el_P;
	double Ib = 0.0, Jft;

	/*--- First invariant of b -> Ib = tr(b) ---*/
	for (i = 0; i < 3; i++){
		Ib += b_Mat[i][i];
	}

	/*--- Retrieve element pressure ---*/

	el_P = element->GetElement_Pressure();

	/*--- J^(-5/3) ---*/
	Jft = pow(J_F, -1.666666666666667);

	for (i = 0; i < 3; i++){
		for (j = 0; j < 3; j++){
			if (i == j) dij = 1.0;
			else if (i != j) dij = 0.0;
			Stress_Tensor[i][j] = Mu *  Jft * (b_Mat[i][j] - ((1.0 / 3.0) * Ib * dij )) + el_P * dij;
		}
	}


}



