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

	unsigned short iVar;

	F_Mat = new su2double *[3];
	b_Mat = new su2double *[3];
	Stress_Tensor = new su2double *[3];
	for (iVar = 0; iVar < 3; iVar++){
		F_Mat[iVar] = new su2double [3];
		b_Mat[iVar] = new su2double [3];
		Stress_Tensor[iVar] = new su2double [3];
	}

	KAux_t_a = new su2double [nDim];

	KAux_P_ab = new su2double* [nDim];
	for (iVar = 0; iVar < nDim; iVar++) {
		KAux_P_ab[iVar] = new su2double[nDim];
	}

	if (nDim == 2){
		currentCoord = new su2double* [4];	/*--- As of now, 4 is the maximum number of nodes for 2D problems ---*/
		for (iVar = 0; iVar < 4; iVar++) currentCoord[iVar] = new su2double[nDim];
	}
	else if (nDim == 3){
		currentCoord = new su2double* [8];	/*--- As of now, 8 is the maximum number of nodes for 3D problems ---*/
		for (iVar = 0; iVar < 8; iVar++) currentCoord[iVar] = new su2double[nDim];
	}

	J_F = 0.0;
	f33 = 1.0;


}

CFEM_NonlinearElasticity::~CFEM_NonlinearElasticity(void) {

	unsigned short iVar;

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


void CFEM_NonlinearElasticity::Compute_Tangent_Matrix(CElement *element, CConfig *config){

	unsigned short iVar, jVar, kVar;
	unsigned short iGauss, nGauss;
	unsigned short iNode, jNode, nNode;
	unsigned short iDim, bDim;

	su2double Ks_Aux_ab;

	su2double Weight, Jac_x;

	su2double AuxMatrixKc[3][6];
	su2double AuxMatrixKs[3];

	/*--- Initialize auxiliary matrices ---*/

	if (nDim == 2) bDim = 3;
	else bDim = 6;

	for (iVar = 0; iVar < bDim; iVar++){
		for (jVar = 0; jVar < nDim; jVar++){
			Ba_Mat[iVar][jVar] = 0.0;
			Bb_Mat[iVar][jVar] = 0.0;
		}
	}

	for (iVar = 0; iVar < 3; iVar++){
		for (jVar = 0; jVar < 6; jVar++){
			AuxMatrixKc[iVar][jVar] = 0.0;
		}
	}

	for (iVar = 0; iVar < 3; iVar++){
		AuxMatrixKs[iVar] = 0.0;
	}

	element->clearElement(); 			/*--- Restarts the element: avoids adding over previous results in other elements --*/
	element->ComputeGrad_NonLinear();

	nNode = element->GetnNodes();
	nGauss = element->GetnGaussPoints();

	/*--- Full integration of the constitutive and stress term ---*/

	for (iGauss = 0; iGauss < nGauss; iGauss++){

		Weight = element->GetWeight(iGauss);
		Jac_x = element->GetJ_x(iGauss);

		/*--- Initialize the deformation gradient for each Gauss Point ---*/

		for (iVar = 0; iVar < 3; iVar++){
			for (jVar = 0; jVar < 3; jVar++){
				F_Mat[iVar][jVar] = 0.0;
				b_Mat[iVar][jVar] = 0.0;
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

			for (iVar = 0; iVar < nDim; iVar++){
				for (jVar = 0; jVar < nDim; jVar++){
					F_Mat[iVar][jVar] += currentCoord[iNode][iVar]*GradNi_Ref_Mat[iNode][jVar];
				}
			}

			/*--- This implies plane strain --> Consider the possible implementation for plane stress --*/
			if (nDim == 2){
				F_Mat[2][2] = 1.0;
			}

		}

		if (nDim == 2) {
			if (plane_stress){
				// Compute the value of the term 33 for the deformation gradient
				Compute_Plane_Stress_Term(element, config);
				F_Mat[2][2] = f33;
			}
			else{
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

		for (iVar = 0; iVar < 3; iVar++){
			for (jVar = 0; jVar < 3; jVar++){
				for (kVar = 0; kVar < 3; kVar++){
					b_Mat[iVar][jVar] += F_Mat[iVar][kVar]*F_Mat[jVar][kVar];
				}
			}
		}

		/*--- Compute the constitutive matrix ---*/

		Compute_Constitutive_Matrix(element, config);
		Compute_Stress_Tensor(element, config);


		for (iNode = 0; iNode < nNode; iNode++){

			/*--------------------------------------------------------------------------------*/
			/*---------------------------- NODAL STRESS TERM ---------------------------------*/
			/*--------------------------------------------------------------------------------*/
		    /*--- Compute the nodal stress term for each gaussian point and for each node, ---*/
		    /*--- and add it to the element structure to be retrieved from the solver      ---*/

			for (iVar = 0; iVar < nDim; iVar++){
				KAux_t_a[iVar] = 0.0;
				for (jVar = 0; jVar < nDim; jVar++){
					KAux_t_a[iVar] += Weight * Stress_Tensor[iVar][jVar] * GradNi_Curr_Mat[iNode][jVar] * Jac_x;
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

			for (iVar = 0; iVar < nDim; iVar++){
				for (jVar = 0; jVar < bDim; jVar++){
					AuxMatrixKc[iVar][jVar] = 0.0;
					for (kVar = 0; kVar < bDim; kVar++){
						AuxMatrixKc[iVar][jVar] += Ba_Mat[kVar][iVar]*D_Mat[kVar][jVar];
					}
				}
			}

		    /*--- Compute the BT.D Matrix ---*/

			for (iVar = 0; iVar < nDim; iVar++){
				AuxMatrixKs[iVar] = 0.0;
				for (jVar = 0; jVar < nDim; jVar++){
					AuxMatrixKs[iVar] += GradNi_Curr_Mat[iNode][jVar]*Stress_Tensor[jVar][iVar]; // DOUBLE CHECK
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
				for (iVar = 0; iVar < nDim; iVar++){
					for (jVar = 0; jVar < nDim; jVar++){
						KAux_ab[iVar][jVar] = 0.0;
						for (kVar = 0; kVar < bDim; kVar++){
							KAux_ab[iVar][jVar] += Weight * AuxMatrixKc[iVar][kVar] * Bb_Mat[kVar][jVar] * Jac_x;
						}
					}
				}

				/*--- Ks_Aux_ab is the term for the constitutive part of the tangent matrix ---*/
				Ks_Aux_ab = 0.0;
				for (iVar = 0; iVar < nDim; iVar++){
					Ks_Aux_ab += Weight * AuxMatrixKs[iVar] * GradNi_Curr_Mat[jNode][iVar] * Jac_x;
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

void CFEM_NonlinearElasticity::Compute_MeanDilatation_Term(CElement *element, CConfig *config){

	unsigned short iVar, jVar;
	unsigned short iGauss, nGauss;
	unsigned short iNode, jNode, nNode;
	su2double Weight, Jac_X, Jac_x;
	unsigned short iDim ;

	su2double GradNi_Mat_Term;
	su2double Vol_current, Vol_reference;
	su2double Avg_kappa;
	su2double el_Pressure;


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

	if ((Vol_current > 0.0) && (Vol_reference > 0.0)) {

		/*--- It is necessary to divide over the current volume to obtain the averaged gradients ---*/
		for (iNode = 0; iNode < nNode; iNode++){
			for (iDim = 0; iDim < nDim; iDim++){
				GradNi_Curr_Mat[iNode][iDim] = GradNi_Curr_Mat[iNode][iDim] / Vol_current;
			}
		}

		Avg_kappa = Kappa * Vol_current / Vol_reference;

		el_Pressure = Kappa * ((Vol_current / Vol_reference) - 1);

		element->SetElement_Pressure(el_Pressure);

	}
	else {
		cout << "Warning: Negative volume computed during FE structural analysis. Exiting..." << endl;
		exit(EXIT_FAILURE);
	}

	for (iNode = 0; iNode < nNode; iNode++){

		for (jNode = 0; jNode < nNode; jNode++){

			/*--- KAux_P_ab is the term for the incompressibility part of the tangent matrix ---*/
			for (iVar = 0; iVar < nDim; iVar++){
				for (jVar = 0; jVar < nDim; jVar++){
					KAux_P_ab[iVar][jVar] = Avg_kappa * Vol_current * GradNi_Curr_Mat[iNode][iVar] * GradNi_Curr_Mat[jNode][jVar];
				}
			}

			element->Set_Kk_ab(KAux_P_ab,iNode, jNode);

		}

	}

}


void CFEM_NonlinearElasticity::Compute_NodalStress_Term(CElement *element, CConfig *config){

	unsigned short iVar, jVar, kVar;
	unsigned short iGauss, nGauss;
	unsigned short iNode, nNode;
	unsigned short iDim;

	su2double Weight, Jac_x;

	element->clearElement(); 			/*--- Restarts the element: avoids adding over previous results in other elements --*/
	element->ComputeGrad_NonLinear();	/*--- Check if we can take this out... so we don't have to do it twice ---*/

	nNode = element->GetnNodes();
	nGauss = element->GetnGaussPoints();

	/*--- Full integration of the nodal stress ---*/

	for (iGauss = 0; iGauss < nGauss; iGauss++){

		Weight = element->GetWeight(iGauss);
		Jac_x = element->GetJ_x(iGauss);

		/*--- Initialize the deformation gradient for each Gauss Point ---*/

		for (iVar = 0; iVar < 3; iVar++){
			for (jVar = 0; jVar < 3; jVar++){
				F_Mat[iVar][jVar] = 0.0;
				b_Mat[iVar][jVar] = 0.0;
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

			for (iVar = 0; iVar < nDim; iVar++){
				for (jVar = 0; jVar < nDim; jVar++){
					F_Mat[iVar][jVar] += currentCoord[iNode][iVar]*GradNi_Ref_Mat[iNode][jVar];
				}
			}

			/*--- This implies plane strain --> Consider the possible implementation for plane stress --*/
			if (nDim == 2){
				F_Mat[2][2] = 1.0;
			}

		}

		if (nDim == 2) {
			if (plane_stress){
				// Compute the value of the term 33 for the deformation gradient
				Compute_Plane_Stress_Term(element, config);
				F_Mat[2][2] = f33;
			}
			else{
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

		for (iVar = 0; iVar < 3; iVar++){
			for (jVar = 0; jVar < 3; jVar++){
				for (kVar = 0; kVar < 3; kVar++){
					b_Mat[iVar][jVar] += F_Mat[iVar][kVar]*F_Mat[jVar][kVar];
				}
			}
		}

		/*--- Compute the stress tensor ---*/

		Compute_Stress_Tensor(element, config);


		for (iNode = 0; iNode < nNode; iNode++){

		    /*--- Compute the nodal stress term for each gaussian point and for each node, ---*/
		    /*--- and add it to the element structure to be retrieved from the solver      ---*/

			for (iVar = 0; iVar < nDim; iVar++){
				KAux_t_a[iVar] = 0.0;
				for (jVar = 0; jVar < nDim; jVar++){
					KAux_t_a[iVar] += Weight * Stress_Tensor[iVar][jVar] * GradNi_Curr_Mat[iNode][jVar] * Jac_x;
				}
			}

			element->Add_Kt_a(KAux_t_a, iNode);

		}

	}

}

void CFEM_NonlinearElasticity::Compute_Tangent_Matrix_DE(CElement *element, CConfig *config){

	unsigned short iVar, jVar, kVar;
	unsigned short iGauss, nGauss;
	unsigned short iNode, jNode, nNode;
	unsigned short iDim, bDim;

	su2double Ks_Aux_ab;

	su2double Weight, Jac_x;

	su2double AuxMatrixKs[3];

	/*--- Initialize auxiliary matrix ---*/

	for (iVar = 0; iVar < 3; iVar++){
		AuxMatrixKs[iVar] = 0.0;
	}

//	cout << "I ENTER HERE" << endl;

	//element->clearElement(); 			/*--- Restarts the element: avoids adding over previous results in other elements --*/
	element->ComputeGrad_Pressure();    /*--- Compute the gradients for the subintegrated element ---*/

	nNode = element->GetnNodes();
	nGauss = element->GetnGaussPoints();

	/*--- Reduced integration of the stress term for the tangent matrix and for the residual ---*/

	for (iGauss = 0; iGauss < nGauss; iGauss++){

//		cout << "And I compute this " << iGaussP+1 << " time. " << endl;

		Weight = element->GetWeight(iGauss);
		Jac_x = element->GetJ_x(iGauss);

		/*--- Initialize the deformation gradient for each Gauss Point ---*/

		for (iVar = 0; iVar < 3; iVar++){
			for (jVar = 0; jVar < 3; jVar++){
				F_Mat[iVar][jVar] = 0.0;
				b_Mat[iVar][jVar] = 0.0;
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

			for (iVar = 0; iVar < nDim; iVar++){
				for (jVar = 0; jVar < nDim; jVar++){
					F_Mat[iVar][jVar] += currentCoord[iNode][iVar]*GradNi_Ref_Mat[iNode][jVar];
				}
			}

			/*--- This implies plane strain --> Consider the possible implementation for plane stress --*/
			if (nDim == 2){
				F_Mat[2][2] = 1.0;
			}

		}

		cout << "Node 1 (x, y): " << currentCoord[0][0] << ", " << currentCoord[0][1] << endl;
		cout << "Node 2 (x, y): " << currentCoord[1][0] << ", " << currentCoord[1][1] << endl;
		cout << "Node 3 (x, y): " << currentCoord[2][0] << ", " << currentCoord[2][1] << endl;
		cout << "Node 4 (x, y): " << currentCoord[3][0] << ", " << currentCoord[3][1] << endl;

		// As a first approach, we assume plane strain

//		if (nDim == 2) {
//			if (plane_stress){
//				// Compute the value of the term 33 for the deformation gradient
//				Compute_Plane_Stress_Term(element);
//				F_Mat[2][2] = f33;
//			}
//			else{
//				F_Mat[2][2] = 1.0;
//			}
//		}

		/*--- Determinant of F --> Jacobian of the transformation ---*/

		J_F = 	F_Mat[0][0]*F_Mat[1][1]*F_Mat[2][2]+
				F_Mat[0][1]*F_Mat[1][2]*F_Mat[2][0]+
				F_Mat[0][2]*F_Mat[1][0]*F_Mat[2][1]-
				F_Mat[0][2]*F_Mat[1][1]*F_Mat[2][0]-
				F_Mat[1][2]*F_Mat[2][1]*F_Mat[0][0]-
				F_Mat[2][2]*F_Mat[0][1]*F_Mat[1][0];

		/*--- Compute the left Cauchy deformation tensor ---*/

		for (iVar = 0; iVar < 3; iVar++){
			for (jVar = 0; jVar < 3; jVar++){
				for (kVar = 0; kVar < 3; kVar++){
					b_Mat[iVar][jVar] += F_Mat[iVar][kVar]*F_Mat[jVar][kVar];
				}
			}
		}

//		cout << b_Mat[0][0] << " " << b_Mat[0][1] << " " << b_Mat[0][2] << endl;
//		cout << b_Mat[1][0] << " " << b_Mat[1][1] << " " << b_Mat[1][2] << endl;
//		cout << b_Mat[2][0] << " " << b_Mat[2][1] << " " << b_Mat[2][2] << endl;

		/*--- Compute the constitutive matrix ---*/

		Compute_Eigenproblem(element, config);

		Compute_Stress_Tensor(element, config);

		for (iNode = 0; iNode < nNode; iNode++){

			/*--------------------------------------------------------------------------------*/
			/*---------------------------- NODAL STRESS TERM ---------------------------------*/
			/*--------------------------------------------------------------------------------*/
		    /*--- Compute the nodal stress term for each gaussian point and for each node, ---*/
		    /*--- and add it to the element structure to be retrieved from the solver      ---*/

			for (iVar = 0; iVar < nDim; iVar++){
				KAux_t_a[iVar] = 0.0;
				for (jVar = 0; jVar < nDim; jVar++){
					KAux_t_a[iVar] += Weight * Stress_Tensor[iVar][jVar] * GradNi_Curr_Mat[iNode][jVar] * Jac_x;
				}
			}

			/*--------------------------------------------------------------------------------*/
			/*-------------------------------  STRESS TERM  ----------------------------------*/
			/*--------------------------------------------------------------------------------*/

		    /*--- Compute the BT.D Matrix ---*/

			for (iVar = 0; iVar < nDim; iVar++){
				AuxMatrixKs[iVar] = 0.0;
				for (jVar = 0; jVar < nDim; jVar++){
					AuxMatrixKs[iVar] += GradNi_Curr_Mat[iNode][jVar]*Stress_Tensor[jVar][iVar]; // DOUBLE CHECK
				}
			}

			/*--- Assumming symmetry ---*/
			for (jNode = iNode; jNode < nNode; jNode++){

				/*--- Ks_Aux_ab is the term for the constitutive part of the tangent matrix ---*/
				Ks_Aux_ab = 0.0;
				for (iVar = 0; iVar < nDim; iVar++){
					Ks_Aux_ab += Weight * AuxMatrixKs[iVar] * GradNi_Curr_Mat[jNode][iVar] * Jac_x;
				}

			}

		}

	}

}

void CFEM_NonlinearElasticity::Compute_Eigenproblem(CElement *element, CConfig *config){

	su2double l1, l2, J1, J2;
	su2double v12_1, v12_2;
	su2double v21_1, v21_2;

	unsigned short iVar, jVar, kVar;
	double C_Mat[3][3];

	// Define 2 unit vectors E1 and E2 in the reference configuration (1,0) and (0,1)
	// The vectors E1_def and E2_def are going to be the deformed of E1 and E2

	su2double E1[2] = {1.0,0.0}, E2[2] = {0.0,1.0};
	su2double E1_def[2] = {0.0,0.0}, E2_def[2] = {0.0,0.0};

	E1_def[0] = F_Mat[0][0]*E1[0]+F_Mat[0][1]*E1[1];
	E1_def[1] = F_Mat[1][0]*E1[0]+F_Mat[1][1]*E1[1];

	E2_def[0] = F_Mat[0][0]*E2[0]+F_Mat[0][1]*E2[1];
	E2_def[1] = F_Mat[1][0]*E2[0]+F_Mat[1][1]*E2[1];

	cout << "Vector (1,0) projects into (" << E1_def[0] << "," <<  E1_def[1] << ") and vector (0,1) projects into (" << E2_def[0] << "," <<  E2_def[1] << ")." <<endl;


	/*--- Compute the right Cauchy deformation tensor ---*/

//	for (iVar = 0; iVar < 3; iVar++){
//		for (jVar = 0; jVar < 3; jVar++){
//			C_Mat[iVar][jVar] = 0.0;
//			for (kVar = 0; kVar < 3; kVar++){
//				C_Mat[iVar][jVar] += F_Mat[kVar][iVar] * F_Mat[kVar][jVar];
//			}
//		}
//	}



//	cout << "----------------------INSIDE THE EIGENPROBLEM ROUTINE-------------------------" << endl;
//
//	cout << C_Mat[0][0] << " " << C_Mat[0][1] << " " << C_Mat[0][2] << endl;
//	cout << C_Mat[1][0] << " " << C_Mat[1][1] << " " << C_Mat[1][2] << endl;
//	cout << C_Mat[2][0] << " " << C_Mat[2][1] << " " << C_Mat[2][2] << endl;
////
//	cout << "------------------------------------------------------------------------------" << endl;

//	if (nDim == 2){
//
//		l1 = 0.5*(b_Mat[0][0]+b_Mat[1][1]) + 0.5*sqrt(pow((b_Mat[0][0] + b_Mat[1][1]),2) - ( 4* (b_Mat[0][0]*b_Mat[1][1] - b_Mat[0][1]*b_Mat[1][0])));
//		l2 = 0.5*(b_Mat[0][0]+b_Mat[1][1]) - 0.5*sqrt(pow((b_Mat[0][0] + b_Mat[1][1]),2) - ( 4* (b_Mat[0][0]*b_Mat[1][1] - b_Mat[0][1]*b_Mat[1][0])));
//
//		J1 = 0.5*(C_Mat[0][0]+C_Mat[1][1]) + 0.5*sqrt(pow((C_Mat[0][0] + C_Mat[1][1]),2) - ( 4* (C_Mat[0][0]*C_Mat[1][1] - C_Mat[0][1]*C_Mat[1][0])));
//		J2 = 0.5*(C_Mat[0][0]+C_Mat[1][1]) - 0.5*sqrt(pow((C_Mat[0][0] + C_Mat[1][1]),2) - ( 4* (C_Mat[0][0]*C_Mat[1][1] - C_Mat[0][1]*C_Mat[1][0])));
//
//		v12_1 = - 1 * b_Mat[1][0] / (b_Mat[1][1]-l1);
//		v12_2 = - 1 * (b_Mat[0][0]-l1) / b_Mat[1][0];
//
//		v21_1 = - 1 * b_Mat[1][0] / (b_Mat[0][0]-l2);
//		v21_2 = - 1 * (b_Mat[1][1]-l2) / b_Mat[1][0];
//
//		cout << "Eigenvector b 1: (1," << v12_1 << ") or (1," << v12_2 << ")" << endl;
//		cout << "Eigenvector b 2: (" << v21_1 << ",1) or (" << v21_2 << ",1)" << endl;
//
//		cout << "LAMBDA_1^2=" << l1 << " and LAMBDA_2^2=" << l2 << ". " << "LAMBDA_1=" << sqrt(l1) << " and LAMBDA_2=" << sqrt(l2) << endl;
//
//		v12_1 = - 1 * C_Mat[1][0] / (C_Mat[1][1]-J1);
//		v12_2 = - 1 * (C_Mat[0][0]-J1) / C_Mat[1][0];
//
//		v21_1 = - 1 * C_Mat[1][0] / (C_Mat[0][0]-J2);
//		v21_2 = - 1 * (C_Mat[1][1]-J2) / C_Mat[1][0];
//
//		cout << "Eigenvector C 1: (1," << v12_1 << ") or (1," << v12_2 << ")" << endl;
//		cout << "Eigenvector C 2: (" << v21_1 << ",1) or (" << v21_2 << ",1)" << endl;
//
//		cout << "LAMBDA_1^2_C=" << J1 << " and LAMBDA_2^2_C=" << J2 << ". " << "LAMBDA_1_C=" << sqrt(J1) << " and LAMBDA_2_C=" << sqrt(J2) << endl;
//
//
////		cout << b_Mat[1][0] << endl;
////		cout << "Check L1=" << (b_Mat[0][0]-l1)*(b_Mat[1][1]-l1)-(b_Mat[0][1])*(b_Mat[1][0]) << endl;
////		cout << "Check L2=" << (b_Mat[0][0]-l2)*(b_Mat[1][1]-l2)-(b_Mat[0][1])*(b_Mat[1][0]) << endl;
//	}



}

void CFEM_NonlinearElasticity::Compute_NodalStress_Term_DE(CElement *element, CConfig *config){

//	unsigned short iVar, jVar, kVar;
//	unsigned short iGauss, nGauss;
//	unsigned short iNode, nNode;
//	unsigned short iDim;
//
//	su2double Weight, Jac_x;
//
//	element->clearElement(); 			/*--- Restarts the element: avoids adding over previous results in other elements --*/
//	element->ComputeGrad_NonLinear();	/*--- Check if we can take this out... so we don't have to do it twice ---*/
//
//	nNode = element->GetnNodes();
//	nGauss = element->GetnGaussPoints();
//
//	/*--- Full integration of the nodal stress ---*/
//
//	for (iGauss = 0; iGauss < nGauss; iGauss++){
//
//		Weight = element->GetWeight(iGauss);
//		Jac_x = element->GetJ_x(iGauss);
//
//		/*--- Initialize the deformation gradient for each Gauss Point ---*/
//
//		for (iVar = 0; iVar < 3; iVar++){
//			for (jVar = 0; jVar < 3; jVar++){
//				F_Mat[iVar][jVar] = 0.0;
//				b_Mat[iVar][jVar] = 0.0;
//			}
//		}
//
//		/*--- Retrieve the values of the gradients of the shape functions for each node ---*/
//		/*--- This avoids repeated operations ---*/
//
//		for (iNode = 0; iNode < nNode; iNode++){
//
//			for (iDim = 0; iDim < nDim; iDim++){
//				GradNi_Ref_Mat[iNode][iDim] = element->GetGradNi_X(iNode,iGauss,iDim);
//				GradNi_Curr_Mat[iNode][iDim] = element->GetGradNi_x(iNode,iGauss,iDim);
//				currentCoord[iNode][iDim] = element->GetCurr_Coord(iNode, iDim);
//			}
//
//			/*--- Compute the deformation gradient ---*/
//
//			for (iVar = 0; iVar < nDim; iVar++){
//				for (jVar = 0; jVar < nDim; jVar++){
//					F_Mat[iVar][jVar] += currentCoord[iNode][iVar]*GradNi_Ref_Mat[iNode][jVar];
//				}
//			}
//
//			/*--- This implies plane strain --> Consider the possible implementation for plane stress --*/
//			if (nDim == 2){
//				F_Mat[2][2] = 1.0;
//			}
//
//		}
//
//		if (nDim == 2) {
//			if (plane_stress){
//				// Compute the value of the term 33 for the deformation gradient
//				Compute_Plane_Stress_Term(element);
//				F_Mat[2][2] = f33;
//			}
//			else{
//				F_Mat[2][2] = 1.0;
//			}
//		}
//
//		/*--- Determinant of F --> Jacobian of the transformation ---*/
//
//		J_F = 	F_Mat[0][0]*F_Mat[1][1]*F_Mat[2][2]+
//				F_Mat[0][1]*F_Mat[1][2]*F_Mat[2][0]+
//				F_Mat[0][2]*F_Mat[1][0]*F_Mat[2][1]-
//				F_Mat[0][2]*F_Mat[1][1]*F_Mat[2][0]-
//				F_Mat[1][2]*F_Mat[2][1]*F_Mat[0][0]-
//				F_Mat[2][2]*F_Mat[0][1]*F_Mat[1][0];
//
//		/*--- Compute the left Cauchy deformation tensor ---*/
//
//		for (iVar = 0; iVar < 3; iVar++){
//			for (jVar = 0; jVar < 3; jVar++){
//				for (kVar = 0; kVar < 3; kVar++){
//					b_Mat[iVar][jVar] += F_Mat[iVar][kVar]*F_Mat[jVar][kVar];
//				}
//			}
//		}
//
//		/*--- Compute the stress tensor ---*/
//
//		Compute_Stress_Tensor(element);
//
//
//		for (iNode = 0; iNode < nNode; iNode++){
//
//		    /*--- Compute the nodal stress term for each gaussian point and for each node, ---*/
//		    /*--- and add it to the element structure to be retrieved from the solver      ---*/
//
//			for (iVar = 0; iVar < nDim; iVar++){
//				KAux_t_a[iVar] = 0.0;
//				for (jVar = 0; jVar < nDim; jVar++){
//					KAux_t_a[iVar] += Weight * Stress_Tensor[iVar][jVar] * GradNi_Curr_Mat[iNode][jVar] * Jac_x;
//				}
//			}
//
//			element->Add_Kt_a(KAux_t_a, iNode);
//
//		}
//
//	}

}

void CFEM_NonlinearElasticity::Compute_Averaged_NodalStress(CElement *element, CConfig *config){

	unsigned short iVar, jVar, kVar;
	unsigned short iGauss, nGauss;
	unsigned short iDim, iNode, nNode;

	su2double Weight, Jac_x;

	element->clearStress();
	element->clearElement(); 			/*--- Restarts the element: avoids adding over previous results in other elements --*/
	element->ComputeGrad_NonLinear();

	nNode = element->GetnNodes();
	nGauss = element->GetnGaussPoints();

	/*--- Computation of the deformation gradient ---*/

	for (iGauss = 0; iGauss < nGauss; iGauss++){

		Weight = element->GetWeight(iGauss);
		Jac_x = element->GetJ_x(iGauss);

		/*--- Initialize the deformation gradient for each Gauss Point ---*/

		for (iVar = 0; iVar < 3; iVar++){
			for (jVar = 0; jVar < 3; jVar++){
				F_Mat[iVar][jVar] = 0.0;
				b_Mat[iVar][jVar] = 0.0;
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

			for (iVar = 0; iVar < nDim; iVar++){
				for (jVar = 0; jVar < nDim; jVar++){
					F_Mat[iVar][jVar] += currentCoord[iNode][iVar]*GradNi_Ref_Mat[iNode][jVar];
				}
			}

			/*--- This implies plane strain --> Consider the possible implementation for plane stress --*/
			if (nDim == 2){
				F_Mat[2][2] = 1.0;
			}

		}

		if (nDim == 2) {
			if (plane_stress){
				// Compute the value of the term 33 for the deformation gradient
				Compute_Plane_Stress_Term(element, config);
				F_Mat[2][2] = f33;
			}
			else{
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

		for (iVar = 0; iVar < 3; iVar++){
			for (jVar = 0; jVar < 3; jVar++){
				for (kVar = 0; kVar < 3; kVar++){
					b_Mat[iVar][jVar] += F_Mat[iVar][kVar]*F_Mat[jVar][kVar];
				}
			}
		}

		/*--- Compute the stress tensor ---*/

		Compute_Stress_Tensor(element, config);

		for (iNode = 0; iNode < nNode; iNode++){


		    /*--- Compute the nodal stress term for each gaussian point and for each node, ---*/
		    /*--- and add it to the element structure to be retrieved from the solver      ---*/

			for (iVar = 0; iVar < nDim; iVar++){
				KAux_t_a[iVar] = 0.0;
				for (jVar = 0; jVar < nDim; jVar++){
					KAux_t_a[iVar] += Weight * Stress_Tensor[iVar][jVar] * GradNi_Curr_Mat[iNode][jVar] * Jac_x;
				}
			}

			element->Add_Kt_a(KAux_t_a, iNode);

		    /*--- Compute the average nodal stresses for each node ---*/

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

void CFEM_NeoHookean_Comp::Compute_Plane_Stress_Term(CElement *element, CConfig *config) {

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

void CFEM_NeoHookean_Comp::Compute_Constitutive_Matrix(CElement *element, CConfig *config) {

	su2double Mu_p = 0.0, Lambda_p = 0.0;

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

void CFEM_NeoHookean_Comp::Compute_Stress_Tensor(CElement *element, CConfig *config) {

	unsigned short iVar,jVar;
	su2double Mu_J = 0.0, Lambda_J = 0.0;
	su2double dij = 0.0;

	/*--- This can be done in a better way ---*/
	if (J_F != 0.0){
		Mu_J = Mu/J_F;
		Lambda_J = Lambda/J_F;
	}

	for (iVar = 0; iVar < 3; iVar++){
		for (jVar = 0; jVar < 3; jVar++){
			if (iVar == jVar) dij = 1.0;
			else if (iVar != jVar) dij = 0.0;
			Stress_Tensor[iVar][jVar] = Mu_J * (b_Mat[iVar][jVar] - dij) + Lambda_J * log(J_F) * dij;
		}
	}

}

CFEM_NeoHookean_Incomp::CFEM_NeoHookean_Incomp(unsigned short val_nDim, unsigned short val_nVar,
                                   CConfig *config) : CFEM_NonlinearElasticity(val_nDim, val_nVar, config) {


}

CFEM_NeoHookean_Incomp::~CFEM_NeoHookean_Incomp(void) {

}

void CFEM_NeoHookean_Incomp::Compute_Plane_Stress_Term(CElement *element, CConfig *config) {

}

void CFEM_NeoHookean_Incomp::Compute_Constitutive_Matrix(CElement *element, CConfig *config) {

	unsigned short iVar;
	su2double el_P;
	su2double Ib = 0.0, Jft;

	/*--- First invariant of b -> Ib = tr(b) ---*/
	for (iVar = 0; iVar < 3; iVar++){
		Ib += b_Mat[iVar][iVar];
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

void CFEM_NeoHookean_Incomp::Compute_Stress_Tensor(CElement *element, CConfig *config) {

	unsigned short iDim,jDim;
	su2double dij, el_P;
	su2double Ib = 0.0, Jft;

	/*--- First invariant of b -> Ib = tr(b) ---*/
	for (iDim = 0; iDim < 3; iDim++){
		Ib += b_Mat[iDim][iDim];
	}

	/*--- Retrieve element pressure ---*/

	el_P = element->GetElement_Pressure();

	/*--- J^(-5/3) ---*/
	Jft = pow(J_F, -1.666666666666667);

	for (iDim = 0; iDim < 3; iDim++){
		for (jDim = 0; jDim < 3; jDim++){
			if (iDim == jDim) dij = 1.0;
			else if (iDim != jDim) dij = 0.0;
			Stress_Tensor[iDim][jDim] = Mu *  Jft * (b_Mat[iDim][jDim] - ((1.0 / 3.0) * Ib * dij )) + el_P * dij;
		}
	}


}

CFEM_DielectricElastomer::CFEM_DielectricElastomer(unsigned short val_nDim, unsigned short val_nVar,
                                   CConfig *config) : CFEM_NonlinearElasticity(val_nDim, val_nVar, config) {

	su2double Electric_Field_Mod;
	su2double *Electric_Field_Dir = config->Get_Electric_Field_Dir();
	unsigned short iVar, iDim;

	ke_DE = config->GetDE_Modulus();

	nElectric_Field = config->GetnElectric_Field();
	nDim_Electric_Field = config->GetnDim_Electric_Field();

	if (nDim != nDim_Electric_Field) cout << "DIMENSIONS DON'T AGREE (Fix this)" << endl;

	/*--- Initialize pointer for the electric field ---*/
	Electric_Field_Ref = new su2double* [nElectric_Field];
	for (iVar = 0; iVar < nElectric_Field; iVar++) {
		Electric_Field_Ref[iVar] = new su2double[nDim_Electric_Field];
	}

	/*--- Assign values to the auxiliary Electric_Field structure ---*/
	for (iVar = 0; iVar < nElectric_Field; iVar++) {
		for (iDim = 0; iDim < nDim_Electric_Field; iDim++) {
			Electric_Field_Ref[iVar][iDim] = Electric_Field_Dir[iDim]*config->Get_Electric_Field_Mod(iVar);
		}
	}

	/*--- Auxiliary vector for computing the electric field in the current configuration ---*/
	Electric_Field_Curr = new su2double[nDim_Electric_Field];
	for (iDim = 0; iDim < nDim_Electric_Field; iDim++) {
		Electric_Field_Curr[iDim] = 0.0;
	}


}

CFEM_DielectricElastomer::~CFEM_DielectricElastomer(void) {

	unsigned short iVar;

	for (iVar = 0; iVar < nElectric_Field; iVar++){
		delete [] Electric_Field_Ref[iVar];
	}

	delete [] Electric_Field_Ref;
	delete [] Electric_Field_Curr;

}

void CFEM_DielectricElastomer::Compute_Plane_Stress_Term(CElement *element, CConfig *config) {

}

void CFEM_DielectricElastomer::Compute_Constitutive_Matrix(CElement *element, CConfig *config) {

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

void CFEM_DielectricElastomer::Compute_Stress_Tensor(CElement *element, CConfig *config) {

	unsigned short iVar, iDim, jDim;
	su2double mod_Curr, mod_Ref;

	su2double E0 = 0.0, E1 = 0.0, E2 = 0.0;
	su2double E0_2 = 0.0, E1_2 = 0.0, E2_2 = 0.0;
	su2double E_2 = 0.0;

//	Compute_Eigenproblem(element, config);

	for (iVar = 0; iVar < nElectric_Field; iVar++){
		for (iDim = 0; iDim < nDim; iDim++){
			Electric_Field_Curr[iDim] = 0.0;
			for (jDim = 0; jDim < nDim; jDim++){
				Electric_Field_Curr[iDim] += F_Mat[iDim][jDim] * Electric_Field_Ref[iVar][jDim];
			}
		}

//		mod_Curr = sqrt(pow(Electric_Field_Curr[0],2)+pow(Electric_Field_Curr[1],2));
//		mod_Ref = sqrt(pow(Electric_Field_Ref[iVar][0],2)+pow(Electric_Field_Ref[iVar][1],2));
//		cout << endl;
//		cout << "E_Ref(" << iVar << ")  = (" << Electric_Field_Ref[iVar][0] << "," << Electric_Field_Ref[iVar][1] << ").  |E_Ref|  = " << mod_Ref << "." << endl;
//		cout << "E_Curr(" << iVar << ") = (" << Electric_Field_Curr[0] << "," << Electric_Field_Curr[1] << "). |E_Curr| = " << mod_Curr << "." << endl;
	}

	E0 = Electric_Field_Curr[0];	E0_2 = pow(E0,2);
	E1 = Electric_Field_Curr[1];	E1_2 = pow(E1,2);
	if (nDim == 3) {E2 = Electric_Field_Curr[2];	E2_2 = pow(E2,2);}


//	cout << "E0: " << E0 << " E1: " << E1 << " E0_2: " << E0_2 << " E1_2: " << E1_2 << endl;
//	E_2 = E0_2+E1_2+E2_2;

	Stress_Tensor[0][0] = ke_DE*(E0_2-0.5*E_2);	Stress_Tensor[0][1] = ke_DE*E0*E1;			Stress_Tensor[0][2] = ke_DE*E0*E2;
	Stress_Tensor[1][0] = ke_DE*E1*E0;			Stress_Tensor[1][1] = ke_DE*(E1_2-0.5*E_2);	Stress_Tensor[1][2] = ke_DE*E1*E2;
	Stress_Tensor[2][0] = ke_DE*E2*E0;			Stress_Tensor[2][1] = ke_DE*E2*E1;			Stress_Tensor[2][2] = ke_DE*(E2_2-0.5*E_2);

//	cout << endl << "SXX:" << endl;
//	cout << Stress_Tensor[0][0] << " " << Stress_Tensor[0][1] << " " << Stress_Tensor[0][2] << endl;
//	cout << Stress_Tensor[1][0] << " " << Stress_Tensor[1][1] << " " << Stress_Tensor[1][2] << endl;
//	cout << Stress_Tensor[2][0] << " " << Stress_Tensor[2][1] << " " << Stress_Tensor[2][2] << endl;

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
		cout << endl;
		cout << "E_Ref(" << iVar << ")  = (" << Electric_Field_Ref_Unit[iVar][0] << "," << Electric_Field_Ref_Unit[iVar][1] << ").  |E_Ref|  = " << mod_Ref << "." << endl;
		cout << "E_Curr(" << iVar << ") = (" << Electric_Field_Curr[0] << "," << Electric_Field_Curr[1] << "). |E_Curr| = " << mod_Curr << "." << endl;
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

	cout << endl << "SXX:" << endl;
	cout << Stress_Tensor[0][0] << " " << Stress_Tensor[0][1] << " " << Stress_Tensor[0][2] << endl;
	cout << Stress_Tensor[1][0] << " " << Stress_Tensor[1][1] << " " << Stress_Tensor[1][2] << endl;
	cout << Stress_Tensor[2][0] << " " << Stress_Tensor[2][1] << " " << Stress_Tensor[2][2] << endl;

}

