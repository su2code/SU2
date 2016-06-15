/*!
 * \file numerics_direct_poisson.cpp
 * \brief This file contains all the convective term discretization.
 * \author F. Palacios
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

CGalerkin_Flow::CGalerkin_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;
}

CGalerkin_Flow::~CGalerkin_Flow(void) { }

void CGalerkin_Flow::ComputeResidual(su2double **val_stiffmatrix_elem, CConfig *config) {
  
	su2double a[4], b[4], c[4], d[4], Area, B_Matrix[4][4];
	unsigned short iVar, jVar;
  
	if (nDim == 2) {
		for (unsigned short iDim = 0; iDim < nDim; iDim++) {
			a[iDim] = Coord_0[iDim]-Coord_2[iDim];
			b[iDim] = Coord_1[iDim]-Coord_2[iDim];
		}
    
		Area = 0.5*fabs(a[0]*b[1]-a[1]*b[0]);	/* Norm of the normal component of area, area = 1/2*cross(a, b) */
    
		a[0] = 0.5 * (Coord_1[0]*Coord_2[1]-Coord_2[0]*Coord_1[1]) / Area;
		a[1] = 0.5 * (Coord_2[0]*Coord_0[1]-Coord_0[0]*Coord_2[1]) / Area;
		a[2] = 0.5 * (Coord_0[0]*Coord_1[1]-Coord_1[0]*Coord_0[1]) / Area;
    
		b[0] = 0.5 * (Coord_1[1]-Coord_2[1]) / Area;
		b[1] = 0.5 * (Coord_2[1]-Coord_0[1]) / Area;
		b[2] = 0.5 * (Coord_0[1]-Coord_1[1]) / Area;
    
		c[0] = 0.5 * (Coord_2[0]-Coord_1[0]) / Area;
		c[1] = 0.5 * (Coord_0[0]-Coord_2[0]) / Area;
		c[2] = 0.5 * (Coord_1[0]-Coord_0[0]) / Area;
    
		/* Compute the stiffness matrix, K & multiply it by the Area */
    
		val_stiffmatrix_elem[0][0] = Area * (b[0]*b[0]+c[0]*c[0]);
		val_stiffmatrix_elem[0][1] = Area * (b[0]*b[1]+c[0]*c[1]);
		val_stiffmatrix_elem[0][2] = Area * (b[0]*b[2]+c[0]*c[2]);
		val_stiffmatrix_elem[1][0] = Area * (b[0]*b[1]+c[0]*c[1]);
		val_stiffmatrix_elem[1][1] = Area * (b[1]*b[1]+c[1]*c[1]);
		val_stiffmatrix_elem[1][2] = Area * (b[1]*b[2]+c[1]*c[2]);
		val_stiffmatrix_elem[2][0] = Area * (b[0]*b[2]+c[0]*c[2]);
		val_stiffmatrix_elem[2][1] = Area * (b[1]*b[2]+c[1]*c[2]);
		val_stiffmatrix_elem[2][2] = Area * (b[2]*b[2]+c[2]*c[2]);
	}
  
	if (nDim == 3) {
		su2double Volume = 0.0;
		Volume -= Determinant_3x3(Coord_1[0], Coord_1[1], Coord_1[2], Coord_2[0], Coord_2[1], Coord_2[2], Coord_3[0], Coord_3[1], Coord_3[2]);
		Volume += Determinant_3x3(Coord_0[0], Coord_0[1], Coord_0[2], Coord_2[0], Coord_2[1], Coord_2[2], Coord_3[0], Coord_3[1], Coord_3[2]);
		Volume -= Determinant_3x3(Coord_0[0], Coord_0[1], Coord_0[2], Coord_1[0], Coord_1[1], Coord_1[2], Coord_3[0], Coord_3[1], Coord_3[2]);
		Volume += Determinant_3x3(Coord_0[0], Coord_0[1], Coord_0[2], Coord_1[0], Coord_1[1], Coord_1[2], Coord_2[0], Coord_2[1], Coord_2[2]);
		Volume = fabs(Volume / 6.0);
    
		a[0] = Determinant_3x3(Coord_1[0], Coord_1[1], Coord_1[2], Coord_2[0], Coord_2[1], Coord_2[2], Coord_3[0], Coord_3[1], Coord_3[2])/(6.0*Volume);
		b[0] = -Determinant_3x3(1.0, Coord_1[1], Coord_1[2],1.0, Coord_2[1], Coord_2[2],1.0, Coord_3[1], Coord_3[2])/(6.0*Volume);
		c[0] = -Determinant_3x3(Coord_1[0],1.0, Coord_1[2], Coord_2[0],1.0, Coord_2[2], Coord_3[0],1.0, Coord_3[2])/(6.0*Volume);
		d[0] = -Determinant_3x3(Coord_1[0], Coord_1[1],1.0, Coord_2[0], Coord_2[1],1.0, Coord_3[0], Coord_3[1],1.0)/(6.0*Volume);
    
		a[1] = -Determinant_3x3(Coord_2[0], Coord_2[1], Coord_2[2], Coord_3[0], Coord_3[1], Coord_3[2], Coord_0[0], Coord_0[1], Coord_0[2])/(6.0*Volume);
		b[1] = Determinant_3x3(1.0, Coord_2[1], Coord_2[2],1.0, Coord_3[1], Coord_3[2],1.0, Coord_0[1], Coord_0[2])/(6.0*Volume);
		c[1] = Determinant_3x3(Coord_2[0],1.0, Coord_2[2], Coord_3[0],1.0, Coord_3[2], Coord_0[0],1.0, Coord_0[2])/(6.0*Volume);
		d[1] = Determinant_3x3(Coord_2[0], Coord_2[1],1.0, Coord_3[0], Coord_3[1],1.0, Coord_0[0], Coord_0[1],1.0)/(6.0*Volume);
    
		a[2] = Determinant_3x3(Coord_3[0], Coord_3[1], Coord_3[2], Coord_0[0], Coord_0[1], Coord_0[2], Coord_1[0], Coord_1[1], Coord_1[2])/(6.0*Volume);
		b[2] = -Determinant_3x3(1.0, Coord_3[1], Coord_3[2],1.0, Coord_0[1], Coord_0[2],1.0, Coord_1[1], Coord_1[2])/(6.0*Volume);
		c[2] = -Determinant_3x3(Coord_3[0],1.0, Coord_3[2], Coord_0[0],1.0, Coord_0[2], Coord_1[0],1.0, Coord_1[2])/(6.0*Volume);
		d[2] = -Determinant_3x3(Coord_3[0], Coord_3[1],1.0, Coord_0[0], Coord_0[1],1.0, Coord_1[0], Coord_1[1],1.0)/(6.0*Volume);
    
		a[3] = -Determinant_3x3(Coord_0[0], Coord_0[1], Coord_0[2], Coord_1[0], Coord_1[1], Coord_1[2], Coord_2[0], Coord_2[1], Coord_2[2])/(6.0*Volume);
		b[3] = Determinant_3x3(1.0, Coord_0[1], Coord_0[2],1.0, Coord_1[1], Coord_1[2],1.0, Coord_2[1], Coord_2[2])/(6.0*Volume);
		c[3] = Determinant_3x3(Coord_0[0],1.0, Coord_0[2], Coord_1[0],1.0, Coord_1[2], Coord_2[0],1.0, Coord_2[2])/(6.0*Volume);
		d[3] = Determinant_3x3(Coord_0[0], Coord_0[1],1.0, Coord_1[0], Coord_1[1],1.0, Coord_2[0], Coord_2[1],1.0)/(6.0*Volume);
    
		/*--- Compute the B Matrix = grad N_j, dot grad N_i  ---*/
		B_Matrix[0][0] = b[0]*b[0] + c[0]*c[0] + d[0]*d[0];
		B_Matrix[0][1] = b[0]*b[1] + c[0]*c[1] + d[0]*d[1];
		B_Matrix[0][2] = b[0]*b[2] + c[0]*c[2] + d[0]*d[2];
		B_Matrix[0][3] = b[0]*b[3] + c[0]*c[3] + d[0]*d[3];
    
		B_Matrix[1][0] = b[1]*b[0] + c[1]*c[0] + d[1]*d[0];
		B_Matrix[1][1] = b[1]*b[1] + c[1]*c[1] + d[1]*d[1];
		B_Matrix[1][2] = b[1]*b[2] + c[1]*c[2] + d[1]*d[2];
		B_Matrix[1][3] = b[1]*b[3] + c[1]*c[3] + d[1]*d[3];
    
		B_Matrix[2][0] = b[2]*b[0] + c[2]*c[0] + d[2]*d[0];
		B_Matrix[2][1] = b[2]*b[1] + c[2]*c[1] + d[2]*d[1];
		B_Matrix[2][2] = b[2]*b[2] + c[2]*c[2] + d[2]*d[2];
		B_Matrix[2][3] = b[2]*b[3] + c[2]*c[3] + d[2]*d[3];
    
		B_Matrix[3][0] = b[3]*b[0] + c[3]*c[0] + d[3]*d[0];
		B_Matrix[3][1] = b[3]*b[1] + c[3]*c[1] + d[3]*d[1];
		B_Matrix[3][2] = b[3]*b[2] + c[3]*c[2] + d[3]*d[2];
		B_Matrix[3][3] = b[3]*b[3] + c[3]*c[3] + d[3]*d[3];
    
    
		/*--- Compute the BT.D.B Matrix (stiffness matrix) ---*/
		for (iVar = 0; iVar < 4; iVar++)
			for (jVar = 0; jVar < 4; jVar++)
				val_stiffmatrix_elem[iVar][jVar] = Volume * B_Matrix[iVar][jVar];
    
	}
}
