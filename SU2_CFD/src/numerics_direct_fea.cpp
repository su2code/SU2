/*!
 * \file numerics_direct_fea.cpp
 * \brief This file contains all the convective term discretization.
 * \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
 * \version 2.0.6
 *
 * Stanford University Unstructured (SU2) Code
 * Copyright (C) 2012 Aerospace Design Laboratory
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "../include/numerics_structure.hpp"
#include <limits>

CGalerkin_FEA::CGalerkin_FEA(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
	E = config->GetElasticyMod();
	Nu = config->GetPoissonRatio();
	Mu = E / (2.0*(1.0 + Nu));
	Lambda = Nu*E/((1.0+Nu)*(1.0-2.0*Nu));
}

CGalerkin_FEA::~CGalerkin_FEA(void) { }

void CGalerkin_FEA::ComputeResidual(double **val_stiffmatrix_elem, CConfig *config) {
  
	double a[4], b[4], c[4], d[4], Area, B_Matrix[6][12], BT_Matrix[12][6], D_Matrix[6][6], Aux_Matrix[12][6];
	unsigned short iDim, iVar, jVar, kVar;
  
	if (nDim == 2) {
    
		for (iDim = 0; iDim < nDim; iDim++) {
			a[iDim] = Coord_0[iDim]-Coord_2[iDim];
			b[iDim] = Coord_1[iDim]-Coord_2[iDim];
		}
    
		Area = 0.5*fabs(a[0]*b[1]-a[1]*b[0]);
    
		a[0] = 0.5 * (Coord_1[0]*Coord_2[1]-Coord_2[0]*Coord_1[1]) / Area;
		a[1] = 0.5 * (Coord_2[0]*Coord_0[1]-Coord_0[0]*Coord_2[1]) / Area;
		a[2] = 0.5 * (Coord_0[0]*Coord_1[1]-Coord_1[0]*Coord_0[1]) / Area;
    
		b[0] = 0.5 * (Coord_1[1]-Coord_2[1]) / Area;
		b[1] = 0.5 * (Coord_2[1]-Coord_0[1]) / Area;
		b[2] = 0.5 * (Coord_0[1]-Coord_1[1]) / Area;
    
		c[0] = 0.5 * (Coord_2[0]-Coord_1[0]) / Area;
		c[1] = 0.5 * (Coord_0[0]-Coord_2[0]) / Area;
		c[2] = 0.5 * (Coord_1[0]-Coord_0[0]) / Area;
    
		/*		 cout << a[0]+b[0]*Coord_0[0]+c[0]*Coord_0[1] << endl;
		 cout << a[0]+b[0]*Coord_1[0]+c[0]*Coord_1[1] << endl;
		 cout << a[0]+b[0]*Coord_2[0]+c[0]*Coord_2[1] << endl;
     
		 cout << a[1]+b[1]*Coord_0[0]+c[1]*Coord_0[1] << endl;
		 cout << a[1]+b[1]*Coord_1[0]+c[1]*Coord_1[1] << endl;
		 cout << a[1]+b[1]*Coord_2[0]+c[1]*Coord_2[1] << endl;
     
		 cout << a[2]+b[2]*Coord_0[0]+c[2]*Coord_0[1]+d[2] << endl;
		 cout << a[2]+b[2]*Coord_1[0]+c[2]*Coord_1[1]+d[2] << endl;
		 cout << a[2]+b[2]*Coord_2[0]+c[2]*Coord_2[1]+d[2] << endl;
     
		 cin.get();
		 */
    
		/*--- Compute the B Matrix ---*/
		B_Matrix[0][0] = b[0];	B_Matrix[0][1] = 0.0;		B_Matrix[0][2] = b[1];	B_Matrix[0][3] = 0.0;		B_Matrix[0][4] = b[2];	B_Matrix[0][5] = 0.0;
		B_Matrix[1][0] = 0.0;		B_Matrix[1][1] = c[0];	B_Matrix[1][2] = 0.0;		B_Matrix[1][3] = c[1];	B_Matrix[1][4] = 0.0;		B_Matrix[1][5] = c[2];
		B_Matrix[2][0] = c[0];	B_Matrix[2][1] = b[0];	B_Matrix[2][2] = c[1];	B_Matrix[2][3] = b[1];	B_Matrix[2][4] = c[2];	B_Matrix[2][5] = b[2];
    
		for (iVar = 0; iVar < 3; iVar++)
			for (jVar = 0; jVar < 6; jVar++)
				BT_Matrix[jVar][iVar] = B_Matrix[iVar][jVar];
    
		/*--- Compute the D Matrix (for plane strain and 3-D)---*/
		D_Matrix[0][0] = Lambda + 2.0*Mu;		D_Matrix[0][1] = Lambda;						D_Matrix[0][2] = 0.0;
		D_Matrix[1][0] = Lambda;						D_Matrix[1][1] = Lambda + 2.0*Mu;		D_Matrix[1][2] = 0.0;
		D_Matrix[2][0] = 0.0;								D_Matrix[2][1] = 0.0;								D_Matrix[2][2] = Mu;
    
		/*		double coeff = E/(1.0-Nu*Nu);
     D_Matrix[0][0] = 1.0;		D_Matrix[0][1] = Nu;		D_Matrix[0][2] = 0.0;
     D_Matrix[1][0] = Nu;		D_Matrix[1][1] = 1.0;		D_Matrix[1][2] = 0.0;
     D_Matrix[2][0] = 0.0;		D_Matrix[2][1] = 0.0;		D_Matrix[2][2] = (1.0-Nu)/2.0;	*/
    
		/*--- Compute the BT.D Matrix ---*/
		for (iVar = 0; iVar < 6; iVar++) {
			for (jVar = 0; jVar < 3; jVar++) {
				Aux_Matrix[iVar][jVar] = 0.0;
				for (kVar = 0; kVar < 3; kVar++)
					Aux_Matrix[iVar][jVar] += BT_Matrix[iVar][kVar]*D_Matrix[kVar][jVar];
			}
		}
    
		/*--- Compute the BT.D.B Matrix (stiffness matrix) ---*/
		for (iVar = 0; iVar < 6; iVar++) {
			for (jVar = 0; jVar < 6; jVar++) {
				val_stiffmatrix_elem[iVar][jVar] = 0.0;
				for (kVar = 0; kVar < 3; kVar++)
					val_stiffmatrix_elem[iVar][jVar] += Area * Aux_Matrix[iVar][kVar]*B_Matrix[kVar][jVar];
			}
		}
    
	}
  
	if (nDim == 3) {
    
		double Volume = 0.0;
		Volume -= Determinant_3x3(Coord_1[0],Coord_1[1],Coord_1[2],Coord_2[0],Coord_2[1],Coord_2[2],Coord_3[0],Coord_3[1],Coord_3[2]);
		Volume += Determinant_3x3(Coord_0[0],Coord_0[1],Coord_0[2],Coord_2[0],Coord_2[1],Coord_2[2],Coord_3[0],Coord_3[1],Coord_3[2]);
		Volume -= Determinant_3x3(Coord_0[0],Coord_0[1],Coord_0[2],Coord_1[0],Coord_1[1],Coord_1[2],Coord_3[0],Coord_3[1],Coord_3[2]);
		Volume += Determinant_3x3(Coord_0[0],Coord_0[1],Coord_0[2],Coord_1[0],Coord_1[1],Coord_1[2],Coord_2[0],Coord_2[1],Coord_2[2]);
		Volume = fabs(Volume / 6.0);
    
		a[0] = Determinant_3x3(Coord_1[0],Coord_1[1],Coord_1[2],Coord_2[0],Coord_2[1],Coord_2[2],Coord_3[0],Coord_3[1],Coord_3[2])/(6.0*Volume);
		b[0] = -Determinant_3x3(1.0,Coord_1[1],Coord_1[2],1.0,Coord_2[1],Coord_2[2],1.0,Coord_3[1],Coord_3[2])/(6.0*Volume);
		c[0] = -Determinant_3x3(Coord_1[0],1.0,Coord_1[2],Coord_2[0],1.0,Coord_2[2],Coord_3[0],1.0,Coord_3[2])/(6.0*Volume);
		d[0] = -Determinant_3x3(Coord_1[0],Coord_1[1],1.0,Coord_2[0],Coord_2[1],1.0,Coord_3[0],Coord_3[1],1.0)/(6.0*Volume);
    
		a[1] = -Determinant_3x3(Coord_2[0],Coord_2[1],Coord_2[2],Coord_3[0],Coord_3[1],Coord_3[2],Coord_0[0],Coord_0[1],Coord_0[2])/(6.0*Volume);
		b[1] = Determinant_3x3(1.0,Coord_2[1],Coord_2[2],1.0,Coord_3[1],Coord_3[2],1.0,Coord_0[1],Coord_0[2])/(6.0*Volume);
		c[1] = Determinant_3x3(Coord_2[0],1.0,Coord_2[2],Coord_3[0],1.0,Coord_3[2],Coord_0[0],1.0,Coord_0[2])/(6.0*Volume);
		d[1] = Determinant_3x3(Coord_2[0],Coord_2[1],1.0,Coord_3[0],Coord_3[1],1.0,Coord_0[0],Coord_0[1],1.0)/(6.0*Volume);
    
		a[2] = Determinant_3x3(Coord_3[0],Coord_3[1],Coord_3[2],Coord_0[0],Coord_0[1],Coord_0[2],Coord_1[0],Coord_1[1],Coord_1[2])/(6.0*Volume);
		b[2] = -Determinant_3x3(1.0,Coord_3[1],Coord_3[2],1.0,Coord_0[1],Coord_0[2],1.0,Coord_1[1],Coord_1[2])/(6.0*Volume);
		c[2] = -Determinant_3x3(Coord_3[0],1.0,Coord_3[2],Coord_0[0],1.0,Coord_0[2],Coord_1[0],1.0,Coord_1[2])/(6.0*Volume);
		d[2] = -Determinant_3x3(Coord_3[0],Coord_3[1],1.0,Coord_0[0],Coord_0[1],1.0,Coord_1[0],Coord_1[1],1.0)/(6.0*Volume);
    
		a[3] = -Determinant_3x3(Coord_0[0],Coord_0[1],Coord_0[2],Coord_1[0],Coord_1[1],Coord_1[2],Coord_2[0],Coord_2[1],Coord_2[2])/(6.0*Volume);
		b[3] = Determinant_3x3(1.0,Coord_0[1],Coord_0[2],1.0,Coord_1[1],Coord_1[2],1.0,Coord_2[1],Coord_2[2])/(6.0*Volume);
		c[3] = Determinant_3x3(Coord_0[0],1.0,Coord_0[2],Coord_1[0],1.0,Coord_1[2],Coord_2[0],1.0,Coord_2[2])/(6.0*Volume);
		d[3] = Determinant_3x3(Coord_0[0],Coord_0[1],1.0,Coord_1[0],Coord_1[1],1.0,Coord_2[0],Coord_2[1],1.0)/(6.0*Volume);
    
		/*
     if (Volume < 0.0) {
     cout << Volume << endl;
     
     cout << a[0]+b[0]*Coord_0[0]+c[0]*Coord_0[1]+d[0]*Coord_0[2] << endl;
     cout << a[0]+b[0]*Coord_1[0]+c[0]*Coord_1[1]+d[0]*Coord_1[2] << endl;
     cout << a[0]+b[0]*Coord_2[0]+c[0]*Coord_2[1]+d[0]*Coord_2[2] << endl;
     cout << a[0]+b[0]*Coord_3[0]+c[0]*Coord_3[1]+d[0]*Coord_3[2] << endl;
     
     cout << a[1]+b[1]*Coord_0[0]+c[1]*Coord_0[1]+d[1]*Coord_0[2] << endl;
     cout << a[1]+b[1]*Coord_1[0]+c[1]*Coord_1[1]+d[1]*Coord_1[2] << endl;
     cout << a[1]+b[1]*Coord_2[0]+c[1]*Coord_2[1]+d[1]*Coord_2[2] << endl;
     cout << a[1]+b[1]*Coord_3[0]+c[1]*Coord_3[1]+d[1]*Coord_3[2] << endl;
     
     cout << a[2]+b[2]*Coord_0[0]+c[2]*Coord_0[1]+d[2]*Coord_0[2] << endl;
     cout << a[2]+b[2]*Coord_1[0]+c[2]*Coord_1[1]+d[2]*Coord_1[2] << endl;
     cout << a[2]+b[2]*Coord_2[0]+c[2]*Coord_2[1]+d[2]*Coord_2[2] << endl;
     cout << a[2]+b[2]*Coord_3[0]+c[2]*Coord_3[1]+d[2]*Coord_3[2] << endl;
     
     cout << a[3]+b[3]*Coord_0[0]+c[3]*Coord_0[1]+d[3]*Coord_0[2] << endl;
     cout << a[3]+b[3]*Coord_1[0]+c[3]*Coord_1[1]+d[3]*Coord_1[2] << endl;
     cout << a[3]+b[3]*Coord_2[0]+c[3]*Coord_2[1]+d[3]*Coord_2[2] << endl;
     cout << a[3]+b[3]*Coord_3[0]+c[3]*Coord_3[1]+d[3]*Coord_3[2] << endl;
     
     cin.get();
     }*/
    
		/*--- Compute the B Matrix ---*/
		B_Matrix[0][0] = b[0];	B_Matrix[0][1] = 0.0;		B_Matrix[0][2] = 0.0;		B_Matrix[0][3] = b[1];	B_Matrix[0][4] = 0.0;		B_Matrix[0][5] = 0.0;		B_Matrix[0][6] = b[2];	B_Matrix[0][7] = 0.0;		B_Matrix[0][8] = 0.0;		B_Matrix[0][9] = b[3];	B_Matrix[0][10] = 0.0;	B_Matrix[0][11] = 0.0;
		B_Matrix[1][0] = 0.0;		B_Matrix[1][1] = c[0];	B_Matrix[1][2] = 0.0;		B_Matrix[1][3] = 0.0;		B_Matrix[1][4] = c[1];	B_Matrix[1][5] = 0.0;		B_Matrix[1][6] = 0.0;		B_Matrix[1][7] = c[2];	B_Matrix[1][8] = 0.0;		B_Matrix[1][9] = 0.0;		B_Matrix[1][10] = c[3];	B_Matrix[1][11] = 0.0;
		B_Matrix[2][0] = 0.0;		B_Matrix[2][1] = 0.0;		B_Matrix[2][2] = d[0];	B_Matrix[2][3] = 0.0;		B_Matrix[2][4] = 0.0;		B_Matrix[2][5] = d[1];	B_Matrix[2][6] = 0.0;		B_Matrix[2][7] = 0.0;		B_Matrix[2][8] = d[2];	B_Matrix[2][9] = 0.0;		B_Matrix[2][10] = 0.0;	B_Matrix[2][11] = d[3];
		B_Matrix[3][0] = c[0];	B_Matrix[3][1] = b[0];	B_Matrix[3][2] = 0.0;		B_Matrix[3][3] = c[1];	B_Matrix[3][4] = b[1];	B_Matrix[3][5] = 0.0;		B_Matrix[3][6] = c[2];	B_Matrix[3][7] = b[2];	B_Matrix[3][8] = 0.0;		B_Matrix[3][9] = c[3];	B_Matrix[3][10] = b[3];	B_Matrix[3][11] = 0.0;
		B_Matrix[4][0] = 0.0;		B_Matrix[4][1] = d[0];	B_Matrix[4][2] = c[0];	B_Matrix[4][3] = 0.0;		B_Matrix[4][4] = d[1];	B_Matrix[4][5] = c[1];	B_Matrix[4][6] = 0.0;		B_Matrix[4][7] = d[2];	B_Matrix[4][8] = c[2];	B_Matrix[4][9] = 0.0;		B_Matrix[4][10] = d[3];	B_Matrix[4][11] = c[3];
		B_Matrix[5][0] = d[0];	B_Matrix[5][1] = 0.0;		B_Matrix[5][2] = b[0];	B_Matrix[5][3] = d[1];	B_Matrix[5][4] = 0.0;		B_Matrix[5][5] = b[1];	B_Matrix[5][6] = d[2];	B_Matrix[5][7] = 0.0;		B_Matrix[5][8] = b[2];	B_Matrix[5][9] = d[3];	B_Matrix[5][10] = 0.0;	B_Matrix[5][11] = b[3];
    
		for (iVar = 0; iVar < 6; iVar++)
			for (jVar = 0; jVar < 12; jVar++)
				BT_Matrix[jVar][iVar] = B_Matrix[iVar][jVar];
    
		/*--- Compute the D Matrix (for plane strain and 3-D)---*/
		D_Matrix[0][0] = Lambda + 2.0*Mu;	D_Matrix[0][1] = Lambda;					D_Matrix[0][2] = Lambda;					D_Matrix[0][3] = 0.0;	D_Matrix[0][4] = 0.0;	D_Matrix[0][5] = 0.0;
		D_Matrix[1][0] = Lambda;					D_Matrix[1][1] = Lambda + 2.0*Mu;	D_Matrix[1][2] = Lambda;					D_Matrix[1][3] = 0.0;	D_Matrix[1][4] = 0.0;	D_Matrix[1][5] = 0.0;
		D_Matrix[2][0] = Lambda;					D_Matrix[2][1] = Lambda;					D_Matrix[2][2] = Lambda + 2.0*Mu;	D_Matrix[2][3] = 0.0;	D_Matrix[2][4] = 0.0;	D_Matrix[2][5] = 0.0;
		D_Matrix[3][0] = 0.0;							D_Matrix[3][1] = 0.0;							D_Matrix[3][2] = 0.0;							D_Matrix[3][3] = Mu;	D_Matrix[3][4] = 0.0;	D_Matrix[3][5] = 0.0;
		D_Matrix[4][0] = 0.0;							D_Matrix[4][1] = 0.0;							D_Matrix[4][2] = 0.0;							D_Matrix[4][3] = 0.0;	D_Matrix[4][4] = Mu;	D_Matrix[4][5] = 0.0;
		D_Matrix[5][0] = 0.0;							D_Matrix[5][1] = 0.0;							D_Matrix[5][2] = 0.0;							D_Matrix[5][3] = 0.0;	D_Matrix[5][4] = 0.0;	D_Matrix[5][5] = Mu;
    
		/*--- Compute the BT.D Matrix ---*/
		for (iVar = 0; iVar < 12; iVar++) {
			for (jVar = 0; jVar < 6; jVar++) {
				Aux_Matrix[iVar][jVar] = 0.0;
				for (kVar = 0; kVar < 6; kVar++)
					Aux_Matrix[iVar][jVar] += BT_Matrix[iVar][kVar]*D_Matrix[kVar][jVar];
			}
		}
    
		/*--- Compute the BT.D.B Matrix (stiffness matrix) ---*/
		for (iVar = 0; iVar < 12; iVar++) {
			for (jVar = 0; jVar < 12; jVar++) {
				val_stiffmatrix_elem[iVar][jVar] = 0.0;
				for (kVar = 0; kVar < 6; kVar++)
					val_stiffmatrix_elem[iVar][jVar] += Volume * Aux_Matrix[iVar][kVar]*B_Matrix[kVar][jVar];
			}
		}
	}
}