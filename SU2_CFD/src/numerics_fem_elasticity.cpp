/*!
 * \file numerics_fem_linear_elasticity.cpp
 * \brief This file contains the routines for setting the tangent matrix and residual of a FEM linear elastic structural problem.
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

CFEM_Elasticity::CFEM_Elasticity(unsigned short val_nDim, unsigned short val_nVar,
                                   CConfig *config) : CNumerics(val_nDim, val_nVar, config) {

	E = config->GetElasticyMod();
	Nu = config->GetPoissonRatio();
	Rho_s = config->GetMaterialDensity();
	Mu = E / (2.0*(1.0 + Nu));
	Lambda = Nu*E/((1.0+Nu)*(1.0-2.0*Nu));
	Kappa = config->GetBulk_Modulus_Struct();

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

	unsigned short iVar, jVar;

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

}

void CFEM_Elasticity::Compute_Mass_Matrix(CElement *element){

	unsigned short iGauss, nGauss;
	unsigned short iNode, jNode, nNode;
	unsigned short iDim;
	unsigned short bDim;

	su2double Weight, Jac_X;

	su2double val_Mab;

	element->clearElement(); 			/*--- Restarts the element: avoids adding over previous results in other elements --*/
	element->ComputeGrad_Linear();		/*--- Need to compute the gradients to obtain the Jacobian: TODO: this may be improved (another method only to compute J_X) ---*/

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

