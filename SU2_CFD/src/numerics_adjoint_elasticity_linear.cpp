/*!
 * \file numerics_adjoint_elasticity_linear_linear.cpp
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

CFEM_LinearElasticity_Adj::CFEM_LinearElasticity_Adj(unsigned short val_nDim, unsigned short val_nVar,
                                   CConfig *config) : CFEM_LinearElasticity(val_nDim, val_nVar, config) {

	/*--- If it is linear elasticity, dD/dE is constant along the calculations ---*/

	Compute_Constitutive_Matrix();

}

CFEM_LinearElasticity_Adj::~CFEM_LinearElasticity_Adj(void) {

}

void CFEM_LinearElasticity_Adj::Compute_Constitutive_Matrix(void){

     /*--- Compute the dD/dE Matrix (for plane stress and 2-D)---*/

	if (nDim == 2){
		if (plane_stress){

			/*--- We enable plane stress cases ---*/

			D_Mat[0][0] = 1.0/(1.0-Nu*Nu);	  		D_Mat[0][1] = Nu/(1.0-Nu*Nu);  	 D_Mat[0][2] = 0.0;
			D_Mat[1][0] = Nu/(1.0-Nu*Nu);    		D_Mat[1][1] = 1.0/(1.0-Nu*Nu);   D_Mat[1][2] = 0.0;
			D_Mat[2][0] = 0.0;               	    D_Mat[2][1] = 0.0;               D_Mat[2][2] = (1.0-Nu)/(2.0*(1.0-Nu*Nu));
		}
		else{
			/*--- Assuming plane strain as a general case ---*/

			D_Mat[0][0] = (1.0-Nu)/((1.0+Nu)*(1.0-2.0*Nu));	D_Mat[0][1] = Nu/((1.0+Nu)*(1.0-2.0*Nu));       	D_Mat[0][2] = 0.0;
			D_Mat[1][0] = Nu/((1.0+Nu)*(1.0-2.0*Nu));       D_Mat[1][1] = (1.0-Nu)/((1.0+Nu)*(1.0-2.0*Nu));  	D_Mat[1][2] = 0.0;
			D_Mat[2][0] = 0.0;              				D_Mat[2][1] = 0.0;               					D_Mat[2][2] = 1.0/(2.0*(1.0+Nu));
		}

	}
	else if (nDim == 3){

		D_Mat[0][0] = (1.0-Nu)/((1.0+Nu)*(1.0-2.0*Nu));	D_Mat[0][1] = Nu/((1.0+Nu)*(1.0-2.0*Nu));		D_Mat[0][2] = Nu/((1.0+Nu)*(1.0-2.0*Nu));
		D_Mat[1][0] = Nu/((1.0+Nu)*(1.0-2.0*Nu));		D_Mat[1][1] = (1.0-Nu)/((1.0+Nu)*(1.0-2.0*Nu));	D_Mat[1][2] = Nu/((1.0+Nu)*(1.0-2.0*Nu));
		D_Mat[2][0] = Nu/((1.0+Nu)*(1.0-2.0*Nu));		D_Mat[2][1] = Nu/((1.0+Nu)*(1.0-2.0*Nu));		D_Mat[2][2] = (1.0-Nu)/((1.0+Nu)*(1.0-2.0*Nu));

		D_Mat[3][3] = 1.0/(2.0*(1.0+Nu));				D_Mat[4][4] = 1.0/(2.0*(1.0+Nu));				D_Mat[5][5] = 1.0/(2.0*(1.0+Nu));

		D_Mat[0][3] = 0.0;	D_Mat[0][4] = 0.0;	D_Mat[0][5] = 0.0;
		D_Mat[1][3] = 0.0;	D_Mat[1][4] = 0.0;	D_Mat[1][5] = 0.0;
		D_Mat[2][3] = 0.0;	D_Mat[2][4] = 0.0;	D_Mat[2][5] = 0.0;
		D_Mat[3][0] = 0.0;	D_Mat[3][1] = 0.0;	D_Mat[3][2] = 0.0;
		D_Mat[4][0] = 0.0;	D_Mat[4][1] = 0.0;	D_Mat[4][2] = 0.0;
		D_Mat[5][0] = 0.0;	D_Mat[5][1] = 0.0;	D_Mat[5][2] = 0.0;

		D_Mat[3][4] = 0.0;	D_Mat[3][5] = 0.0;
		D_Mat[4][3] = 0.0;	D_Mat[4][5] = 0.0;
		D_Mat[5][3] = 0.0;	D_Mat[5][4] = 0.0;

	}

}
