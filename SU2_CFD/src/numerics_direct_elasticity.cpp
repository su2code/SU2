/*!
 * \file numerics_direct_elasticity.cpp
 * \brief This file contains the routines for setting the tangent matrix and residual of a FEM linear elastic structural problem.
 * \author R. Sanchez
 * \version 4.2.0 "Cardinal"
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

CFEM_Elasticity::CFEM_Elasticity(unsigned short val_nDim, unsigned short val_nVar,
                                   CConfig *config) : CNumerics(val_nDim, val_nVar, config) {

	bool body_forces = config->GetDeadLoad();	// Body forces (dead loads).

	E = config->GetElasticyMod();
	Nu = config->GetPoissonRatio();
	Rho_s = config->GetMaterialDensity();
	Mu = E / (2.0*(1.0 + Nu));
	Lambda = Nu*E/((1.0+Nu)*(1.0-2.0*Nu));
	Kappa = Lambda + (2/3)*Mu;
	//Kappa = config->GetBulk_Modulus_Struct();

	// Auxiliary vector for body forces (dead load)
	if (body_forces) FAux_Dead_Load = new su2double [nDim]; else FAux_Dead_Load = NULL;

	plane_stress = (config->GetElas2D_Formulation() == PLANE_STRESS);

	unsigned short iVar;

	KAux_ab = new su2double* [nDim];
	for (iVar = 0; iVar < nDim; iVar++) {
		KAux_ab[iVar] = new su2double[nDim];
	}


	if (nDim == 2){
		Ba_Mat = new su2double* [3];
		Bb_Mat = new su2double* [3];
		D_Mat  = new su2double* [3];
		Ni_Vec  = new su2double [4];			/*--- As of now, 4 is the maximum number of nodes for 2D problems ---*/
		GradNi_Ref_Mat = new su2double* [4];	/*--- As of now, 4 is the maximum number of nodes for 2D problems ---*/
		GradNi_Curr_Mat = new su2double* [4];	/*--- As of now, 4 is the maximum number of nodes for 2D problems ---*/
		for (iVar = 0; iVar < 3; iVar++) {
			Ba_Mat[iVar]  		= new su2double[nDim];
			Bb_Mat[iVar]  		= new su2double[nDim];
			D_Mat[iVar]	   	= new su2double[3];
		}
		for (iVar = 0; iVar < 4; iVar++) {
			GradNi_Ref_Mat[iVar] 	= new su2double[nDim];
			GradNi_Curr_Mat[iVar] 	= new su2double[nDim];
		}
	}
	else if (nDim == 3){
		Ba_Mat = new su2double* [6];
		Bb_Mat = new su2double* [6];
		D_Mat  = new su2double* [6];
		Ni_Vec  = new su2double [8];			/*--- As of now, 8 is the maximum number of nodes for 3D problems ---*/
		GradNi_Ref_Mat = new su2double* [8];	/*--- As of now, 8 is the maximum number of nodes for 3D problems ---*/
		GradNi_Curr_Mat = new su2double* [8];	/*--- As of now, 8 is the maximum number of nodes for 3D problems ---*/
		for (iVar = 0; iVar < 6; iVar++) {
			Ba_Mat[iVar]  		= new su2double[nDim];
			Bb_Mat[iVar]  		= new su2double[nDim];
			D_Mat[iVar]      	= new su2double[6];
		}
		for (iVar = 0; iVar < 8; iVar++) {
			GradNi_Ref_Mat[iVar] 	= new su2double[nDim];
			GradNi_Curr_Mat[iVar] 	= new su2double[nDim];
		}
	}
}

CFEM_Elasticity::~CFEM_Elasticity(void) {

	unsigned short iVar;

	for (iVar = 0; iVar < nDim; iVar++){
		delete [] KAux_ab[iVar];
	}

	if (nDim == 2){
		for (iVar = 0; iVar < 3; iVar++){
			delete [] Ba_Mat[iVar];
			delete [] Bb_Mat[iVar];
			delete [] D_Mat[iVar];
		}
		for (iVar = 0; iVar < 4; iVar++){
			delete [] GradNi_Ref_Mat[iVar];
			delete [] GradNi_Curr_Mat[iVar];
		}
	}
	else if (nDim == 3){
		for (iVar = 0; iVar < 6; iVar++){
			delete [] Ba_Mat[iVar];
			delete [] Bb_Mat[iVar];
			delete [] D_Mat[iVar];
		}
		for (iVar = 0; iVar < 8; iVar++){
			delete [] GradNi_Ref_Mat[iVar];
			delete [] GradNi_Curr_Mat[iVar];
		}
	}

	delete [] KAux_ab;
	delete [] Ba_Mat;
	delete [] Bb_Mat;
	delete [] D_Mat;
	delete [] GradNi_Ref_Mat;
	delete [] GradNi_Curr_Mat;

	if (FAux_Dead_Load 		!= NULL) delete [] FAux_Dead_Load;

}

void CFEM_Elasticity::Compute_Mass_Matrix(CElement *element, CConfig *config){

	unsigned short iGauss, nGauss;
	unsigned short iNode, jNode, nNode;

	su2double Weight, Jac_X;

	su2double val_Mab;

	element->clearElement(); 			/*--- Restarts the element: avoids adding over previous results in other elements --*/
	element->ComputeGrad_Linear();		/*--- Need to compute the gradients to obtain the Jacobian ---*/

	nNode = element->GetnNodes();
	nGauss = element->GetnGaussPoints();

	for (iGauss = 0; iGauss < nGauss; iGauss++){

		Weight = element->GetWeight(iGauss);
		Jac_X = element->GetJ_X(iGauss);			/*--- The mass matrix is computed in the reference configuration ---*/

		/*--- Retrieve the values of the shape functions for each node ---*/
		/*--- This avoids repeated operations ---*/
		for (iNode = 0; iNode < nNode; iNode++){
			Ni_Vec[iNode] = element->GetNi(iNode,iGauss);
		}

		for (iNode = 0; iNode < nNode; iNode++){

			/*--- Assumming symmetry ---*/
			for (jNode = iNode; jNode < nNode; jNode++){

				val_Mab = Weight * Ni_Vec[iNode] * Ni_Vec[jNode] * Jac_X * Rho_s;

				element->Add_Mab(val_Mab,iNode, jNode);
				/*--- Symmetric terms --*/
				if (iNode != jNode){
					element->Add_Mab(val_Mab, jNode, iNode);
				}

			}

		}

	}

}

void CFEM_Elasticity::Compute_Dead_Load(CElement *element, CConfig *config){

	unsigned short iGauss, nGauss;
	unsigned short iNode, iDim, nNode;

	su2double Weight, Jac_X;

	/* -- Gravity directionality:
	 * -- For 2D problems, we assume the direction for gravity is -y
	 * -- For 3D problems, we assume the direction for gravity is -z
	 */
	su2double g_force[3] = {0.0,0.0,0.0};

	if (nDim == 2) g_force[1] = -1*STANDART_GRAVITY;
	else if (nDim == 3) g_force[2] = -1*STANDART_GRAVITY;

	element->clearElement(); 			/*--- Restarts the element: avoids adding over previous results in other elements and sets initial values to 0--*/
	element->ComputeGrad_Linear();		/*--- Need to compute the gradients to obtain the Jacobian ---*/

	nNode = element->GetnNodes();
	nGauss = element->GetnGaussPoints();

	for (iGauss = 0; iGauss < nGauss; iGauss++){

		Weight = element->GetWeight(iGauss);
		Jac_X = element->GetJ_X(iGauss);			/*--- The dead load is computed in the reference configuration ---*/

		/*--- Retrieve the values of the shape functions for each node ---*/
		/*--- This avoids repeated operations ---*/
		for (iNode = 0; iNode < nNode; iNode++){
			Ni_Vec[iNode] = element->GetNi(iNode,iGauss);
		}

		for (iNode = 0; iNode < nNode; iNode++){

			for (iDim = 0; iDim < nDim; iDim++){
				FAux_Dead_Load[iDim] = Weight * Ni_Vec[iNode] * Jac_X * Rho_s * g_force[iDim];
			}

			element->Add_FDL_a(FAux_Dead_Load,iNode);

		}

	}

}

