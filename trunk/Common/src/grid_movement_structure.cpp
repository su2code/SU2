/*!
 * \file grid_movement_structure.cpp
 * \brief Subroutines for doing the grid movement using different strategies.
 * \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
 * \version 2.0.
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

#ifndef NO_OMP
#include <omp.h>
#endif

#include <ctime>
#include "../include/grid_movement_structure.hpp"

using namespace std;

CGridMovement::CGridMovement(void) { }

CGridMovement::~CGridMovement(void) { }

CVolumetricMovement::CVolumetricMovement(CGeometry *geometry) : CGridMovement() {
	unsigned short iVar, jVar;
	unsigned long iElem, iElem_aux;
	
	C_tor = 1E-6;
	C_lin = 1.0;		
	niter = 999999999;
	nDim = geometry->GetnDim();
	
	/*--- This method only use triangles, it is necessary to divide everything into triangles ---*/
	nElem = 0;
	for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
		if (geometry->elem[iElem]->GetVTK_Type() == TRIANGLE) nElem = nElem + 1;
		if (geometry->elem[iElem]->GetVTK_Type() == RECTANGLE) nElem = nElem + 4;
		if (geometry->elem[iElem]->GetVTK_Type() == TETRAHEDRON) nElem = nElem + 4; 
	}
	
	nodes = new unsigned long* [nElem];
	for (iElem = 0; iElem < nElem; iElem++)
		nodes[iElem] = new unsigned long [3];
	
	iElem_aux = 0;
	for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
		if (geometry->elem[iElem]->GetVTK_Type() == TRIANGLE) {
			nodes[iElem_aux][0] = geometry->elem[iElem]->GetNode(0);
			nodes[iElem_aux][1] = geometry->elem[iElem]->GetNode(1);
			nodes[iElem_aux][2] = geometry->elem[iElem]->GetNode(2);
			iElem_aux++;
		}
		if (geometry->elem[iElem]->GetVTK_Type() == RECTANGLE) {
			nodes[iElem_aux][0] = geometry->elem[iElem]->GetNode(0);
			nodes[iElem_aux][1] = geometry->elem[iElem]->GetNode(1);
			nodes[iElem_aux][2] = geometry->elem[iElem]->GetNode(2);
			iElem_aux++;
			
			nodes[iElem_aux][0] = geometry->elem[iElem]->GetNode(0);
			nodes[iElem_aux][1] = geometry->elem[iElem]->GetNode(2);
			nodes[iElem_aux][2] = geometry->elem[iElem]->GetNode(3);
			iElem_aux++;

			nodes[iElem_aux][0] = geometry->elem[iElem]->GetNode(1);
			nodes[iElem_aux][1] = geometry->elem[iElem]->GetNode(2);
			nodes[iElem_aux][2] = geometry->elem[iElem]->GetNode(3);
			iElem_aux++;
			
			nodes[iElem_aux][0] = geometry->elem[iElem]->GetNode(3);
			nodes[iElem_aux][1] = geometry->elem[iElem]->GetNode(0);
			nodes[iElem_aux][2] = geometry->elem[iElem]->GetNode(1);
			iElem_aux++;
		}
		if (geometry->elem[iElem]->GetVTK_Type() == TETRAHEDRON) {
			nodes[iElem_aux][0] = geometry->elem[iElem]->GetNode(0);
			nodes[iElem_aux][1] = geometry->elem[iElem]->GetNode(1);
			nodes[iElem_aux][2] = geometry->elem[iElem]->GetNode(2);
			iElem_aux++;
			
			nodes[iElem_aux][0] = geometry->elem[iElem]->GetNode(3);
			nodes[iElem_aux][1] = geometry->elem[iElem]->GetNode(2);
			nodes[iElem_aux][2] = geometry->elem[iElem]->GetNode(1);
			iElem_aux++;
			
			nodes[iElem_aux][0] = geometry->elem[iElem]->GetNode(0);
			nodes[iElem_aux][1] = geometry->elem[iElem]->GetNode(2);
			nodes[iElem_aux][2] = geometry->elem[iElem]->GetNode(3);
			iElem_aux++;
			
			nodes[iElem_aux][0] = geometry->elem[iElem]->GetNode(0);
			nodes[iElem_aux][1] = geometry->elem[iElem]->GetNode(3);
			nodes[iElem_aux][2] = geometry->elem[iElem]->GetNode(1);
			iElem_aux++;
		}
	}
	
	kijk = new double** [3*nDim];
	for (iVar = 0; iVar < 3*nDim; iVar++) {
		kijk[iVar] = new double* [3*nDim];
		for (jVar = 0; jVar < 3*nDim; jVar++)
			kijk[iVar][jVar] = new double [nElem]; 
	}
	
	oldgcenter = new double* [nDim];
	for (iVar = 0; iVar < nDim; iVar++)
		oldgcenter[iVar] = new double [nElem];
	
	x = new double [nDim*geometry->GetnPoint()];
	Initial_Boundary = new double [nDim*geometry->GetnPoint()+1];
	diagk = new double [nDim*geometry->GetnPoint()+1];
	p = new double [nDim*geometry->GetnPoint()+1];
	r = new double [nDim*geometry->GetnPoint()+1];
	z = new double [nDim*geometry->GetnPoint()+1];
}

CVolumetricMovement::~CVolumetricMovement(void) {
unsigned short iVar, jVar;

	for (iVar = 0; iVar < 6; iVar++) {
		for (jVar = 0; jVar < 6; jVar++)
			delete [] kijk[iVar][jVar]; 
		delete [] kijk[iVar];
	}
	delete [] kijk;

	for (iVar = 0; iVar < 2; iVar++)
		delete [] oldgcenter[iVar];
	delete [] oldgcenter;

	delete [] x;
	delete [] Initial_Boundary;
	delete [] diagk;
	delete [] p;
	delete [] r;
	delete [] z;
}

void CVolumetricMovement::Set2DMatrix_Structure(CGeometry *geometry) {
unsigned short iVar, jVar, kVar, iDim;
unsigned long iPoint_0, iPoint_1, iPoint_2, iPoint, iElem;	
double ar, rx12, ry12, l12, a12, b12, a21, b21, rx13, ry13, l13, a13, b13, a31, b31,
	   rx23, ry23, l23, a23, b23, a32, b32, cosgam, singam;


/*--- 1st step : Calculation of the stiffnesses, based on the former grid  ---*/
	for (iVar = 0; iVar < 6; iVar++)
		for (jVar = 0; jVar < 6; jVar++)
			for (iElem = 0; iElem < nElem; iElem++)
				kijk[iVar][jVar][iElem] = 0;

	for (iPoint = 0; iPoint < 2*geometry->GetnPoint(); iPoint++) {
		x[iPoint] = 0; 
		diagk[iPoint] = 0;
	}

	for (iElem = 0; iElem < nElem; iElem ++) {

		iPoint_0 = nodes[iElem][0];
		iPoint_1 = nodes[iElem][1];
		iPoint_2 = nodes[iElem][2];

		for (iDim = 0; iDim < 2; iDim++) {
			vec_a[iDim] = geometry->node[iPoint_2]->GetCoord(iDim)
				-geometry->node[iPoint_0]->GetCoord(iDim);
			vec_b[iDim] = geometry->node[iPoint_1]->GetCoord(iDim)
				-geometry->node[iPoint_0]->GetCoord(iDim);
		}

		Dim_Normal[0] = (vec_a[1]*vec_b[2]-vec_a[2]*vec_b[1])*0.5;
		Dim_Normal[1] = -(vec_a[0]*vec_b[2]-vec_a[2]*vec_b[0])*0.5;
		Dim_Normal[2] = (vec_a[0]*vec_b[1]-vec_a[1]*vec_b[0])*0.5;

		ar = sqrt(Dim_Normal[0]*Dim_Normal[0]+Dim_Normal[1]*Dim_Normal[1]+Dim_Normal[2]*Dim_Normal[2]);
		
		rx12 = geometry->node[iPoint_1]->GetCoord(0) - geometry->node[iPoint_0]->GetCoord(0);
		ry12 = geometry->node[iPoint_1]->GetCoord(1) - geometry->node[iPoint_0]->GetCoord(1);
		l12 = sqrt(rx12*rx12 + ry12*ry12);
		a12 = rx12/(l12*l12);
		b12 = ry12/(l12*l12);
		a21 = -a12;
		b21 = -b12;

		rx13 = geometry->node[iPoint_2]->GetCoord(0) - geometry->node[iPoint_0]->GetCoord(0);
		ry13 = geometry->node[iPoint_2]->GetCoord(1) - geometry->node[iPoint_0]->GetCoord(1);
		l13 = sqrt(rx13*rx13 + ry13*ry13);
		a13 = rx13/(l13*l13); 
		b13 = ry13/(l13*l13); 
		a31 = -a13; 
		b31 = -b13; 

		rx23 = geometry->node[iPoint_2]->GetCoord(0) - geometry->node[iPoint_1]->GetCoord(0);
		ry23 = geometry->node[iPoint_2]->GetCoord(1) - geometry->node[iPoint_1]->GetCoord(1);
		l23 = sqrt(rx23*rx23 + ry23*ry23);
		a23 = rx23/(l23*l23);
		b23 = ry23/(l23*l23);
		a32 = -a23; 
		b32 = -b23;

		/*--- Torsional stiffness matrix calculation ---*/
		C_mat[0][0] = (l12*l12)*(l13*l13)/(4.0*ar*ar);
		C_mat[0][1] = 0.0;
		C_mat[0][2] = 0.0;
		C_mat[1][0] = 0.0;
		C_mat[1][1] = (l12*l12)*(l23*l23)/(4.0*ar*ar);
		C_mat[1][2] = 0.0;
		C_mat[2][0] = 0.0;
		C_mat[2][1] = 0.0;
		C_mat[2][2] = (l13*l13)*(l23*l23)/(4.0*ar*ar);

		/*--- Torsional kinematic matrix calculation ---*/
		R_mat[0][0] = b13-b12; 
		R_mat[0][1] = a12-a13; 
		R_mat[0][2] = b12;
		R_mat[0][3] = -a12; 
		R_mat[0][4] = -b13; 
		R_mat[0][5] = a13;
		R_mat[1][0] = -b21; 
		R_mat[1][1] =  a21; 
		R_mat[1][2] = b21-b23;
		R_mat[1][3] = a23-a21; 
		R_mat[1][4] = b23; 
		R_mat[1][5] = -a23;
		R_mat[2][0] = b31; 
		R_mat[2][1] = -a31; 
		R_mat[2][2] = -b32;
		R_mat[2][3] = a32; 
		R_mat[2][4] = b32-b31; 
		R_mat[2][5] = a31-a32;

		/*--- Force-reduced torsional matrix calculation ---*/
		for (iVar = 0; iVar < 3; iVar ++)
			for (jVar = 0; jVar < 6; jVar ++)
				Rt_mat[jVar][iVar] = R_mat[iVar][jVar];

		for(iVar = 0; iVar < 3; iVar++)
			for(jVar = 0; jVar < 6; jVar++) {
				Aux_mat[iVar][jVar] = 0.0;
				for(kVar = 0; kVar < 3; kVar++)
					Aux_mat[iVar][jVar] += C_mat[iVar][kVar] * R_mat[kVar][jVar]; 
			}

		for(iVar = 0; iVar < 6; iVar++)
			for(jVar = 0; jVar < 6; jVar++) {
				Ktor_mat[iVar][jVar] = 0.0;
				for(kVar = 0; kVar < 3; kVar++)
					Ktor_mat[iVar][jVar] += Rt_mat[iVar][kVar] * Aux_mat[kVar][jVar];
			}

		for(iVar = 0; iVar < 6; iVar++)
			for(jVar = 0; jVar < 6; jVar++) {
				Ktor_mat[iVar][jVar] = C_tor * Ktor_mat[iVar][jVar];
			}

		/*--- Linear stiffness matrix calculation ---*/
		cosgam = rx12/l12;
		singam = ry12/l12;
		Klin_mat[0][0] = cosgam*cosgam;
		Klin_mat[0][1] = singam*cosgam;
		Klin_mat[1][0] = Klin_mat[0][1];
		Klin_mat[1][1] = singam*singam;

		for(iVar = 0; iVar < 2; iVar++)
			for(jVar = 0; jVar < 2; jVar++) {
				Klin_mat[iVar][jVar+2] = -Klin_mat[iVar][jVar];
				Klin_mat[iVar+2][jVar] = Klin_mat[iVar][jVar+2];
				Klin_mat[iVar+2][jVar+2] = Klin_mat[iVar][jVar];
			}

		for(iVar = 0; iVar < 4; iVar++)
			for(jVar = 0; jVar < 4; jVar++)
				Ktor_mat[iVar][jVar] = Ktor_mat[iVar][jVar] + C_lin/(2.0*l12)*Klin_mat[iVar][jVar];


		cosgam = rx23/l23; 
		singam = ry23/l23;
		Klin_mat[0][0] = cosgam*cosgam;
		Klin_mat[0][1] = singam*cosgam;
		Klin_mat[1][0] = Klin_mat[0][1];
		Klin_mat[1][1] = singam*singam;

		for(iVar = 0; iVar < 2; iVar++)
			for(jVar = 0; jVar < 2; jVar++) {
				Klin_mat[iVar][jVar+2] = -Klin_mat[iVar][jVar];
				Klin_mat[iVar+2][jVar] = Klin_mat[iVar][jVar+2];
				Klin_mat[iVar+2][jVar+2] = Klin_mat[iVar][jVar];
			}

		for(iVar = 0; iVar < 4; iVar++)
			for(jVar = 0; jVar < 4; jVar++)
				Ktor_mat[iVar+2][jVar+2] = Ktor_mat[iVar+2][jVar+2] + C_lin/(2.0*l23)*Klin_mat[iVar][jVar];

		cosgam = rx13/l13; 
		singam = ry13/l13;
		Klin_mat[0][0] = cosgam*cosgam; 
		Klin_mat[0][1] = singam*cosgam;
		Klin_mat[1][0] = Klin_mat[0][1]; 
		Klin_mat[1][1] = singam*singam;

		for(iVar = 0; iVar < 2; iVar++)
			for(jVar = 0; jVar < 2; jVar++) {
				Klin_mat[iVar][jVar+2] = -Klin_mat[iVar][jVar];
				Klin_mat[iVar+2][jVar] = Klin_mat[iVar][jVar+2];
				Klin_mat[iVar+2][jVar+2] = Klin_mat[iVar][jVar];
			}

		for(iVar = 0; iVar < 2; iVar++)
			for(jVar = 0; jVar < 2; jVar++) {
				Ktor_mat[iVar][jVar] = Ktor_mat[iVar][jVar] + C_lin/(2.0*l13)*Klin_mat[iVar][jVar];
				Ktor_mat[iVar][jVar+4] = Ktor_mat[iVar][jVar+4] + C_lin/(2.0*l13)*Klin_mat[iVar][jVar+2];
				Ktor_mat[iVar+4][jVar] = Ktor_mat[iVar+4][jVar] + C_lin/(2.0*l13)*Klin_mat[iVar+2][jVar];
				Ktor_mat[iVar+4][jVar+4] = Ktor_mat[iVar+4][jVar+4] + C_lin/(2.0*l13)*Klin_mat[iVar+2][jVar+2];
			}

		/*--- Storage ---*/
		for(iVar = 0; iVar < 6; iVar++)
			for(jVar = 0; jVar < 6; jVar++)
				kijk[iVar][jVar][iElem] = Ktor_mat[iVar][jVar];

		/*--- Calculation of the diagonal part of the "assembled" matrix ---*/
		diagk[2*iPoint_0] = diagk[2*iPoint_0] + kijk[0][0][iElem];
		diagk[2*iPoint_0+1] = diagk[2*iPoint_0+1] + kijk[1][1][iElem];
		diagk[2*iPoint_1] = diagk[2*iPoint_1] + kijk[2][2][iElem];
		diagk[2*iPoint_1+1] = diagk[2*iPoint_1+1] + kijk[3][3][iElem];
		diagk[2*iPoint_2] = diagk[2*iPoint_2] + kijk[4][4][iElem];
		diagk[2*iPoint_2+1] = diagk[2*iPoint_2+1] + kijk[5][5][iElem];
	}
}

void CVolumetricMovement::Set3DMatrix_Structure(CGeometry *geometry) {
	unsigned short iVar, jVar, iDim;
	unsigned long iPoint_0, iPoint_1, iPoint_2, iPoint, iElem;	
	double ar, rx12, ry12, rz12, l12, rx13, ry13, rz13, l13, rx23, ry23, rz23, l23, l, m, n;

	/*--- 1st step : Calculation of the stiffnesses, based on the former grid ---*/
	for (iVar = 0; iVar < 9; iVar++)
		for (jVar = 0; jVar < 9; jVar++)
			for (iElem = 0; iElem < nElem; iElem++)
				kijk[iVar][jVar][iElem] = 0;
	
	for (iPoint = 0; iPoint < 3*geometry->GetnPoint(); iPoint++) {
		x[iPoint] = 0; diagk[iPoint] = 0; }
	
	for (iElem = 0; iElem < nElem; iElem ++) {
		
		iPoint_0 = nodes[iElem][0]; iPoint_1 = nodes[iElem][1]; iPoint_2 = nodes[iElem][2];
		
		for (iDim = 0; iDim < 3; iDim++) {
			vec_a[iDim] = geometry->node[iPoint_2]->GetCoord(iDim) - geometry->node[iPoint_0]->GetCoord(iDim);
			vec_b[iDim] = geometry->node[iPoint_1]->GetCoord(iDim) - geometry->node[iPoint_0]->GetCoord(iDim);
		}
		
		Dim_Normal[0] = (vec_a[1]*vec_b[2]-vec_a[2]*vec_b[1])*0.5;
		Dim_Normal[1] = -(vec_a[0]*vec_b[2]-vec_a[2]*vec_b[0])*0.5;
		Dim_Normal[2] = (vec_a[0]*vec_b[1]-vec_a[1]*vec_b[0])*0.5;
		
		ar = sqrt(Dim_Normal[0]*Dim_Normal[0]+Dim_Normal[1]*Dim_Normal[1]+Dim_Normal[2]*Dim_Normal[2]);
		
		rx12 = geometry->node[iPoint_1]->GetCoord(0) - geometry->node[iPoint_0]->GetCoord(0);
		ry12 = geometry->node[iPoint_1]->GetCoord(1) - geometry->node[iPoint_0]->GetCoord(1);
		rz12 = geometry->node[iPoint_1]->GetCoord(2) - geometry->node[iPoint_0]->GetCoord(2);
		l12 = sqrt(rx12*rx12 + ry12*ry12 + rz12*rz12);
		
		rx13 = geometry->node[iPoint_2]->GetCoord(0) - geometry->node[iPoint_0]->GetCoord(0);
		ry13 = geometry->node[iPoint_2]->GetCoord(1) - geometry->node[iPoint_0]->GetCoord(1);
		rz13 = geometry->node[iPoint_2]->GetCoord(2) - geometry->node[iPoint_0]->GetCoord(2);
		l13 = sqrt(rx13*rx13 + ry13*ry13 + rz13*rz13);

		
		rx23 = geometry->node[iPoint_2]->GetCoord(0) - geometry->node[iPoint_1]->GetCoord(0);
		ry23 = geometry->node[iPoint_2]->GetCoord(1) - geometry->node[iPoint_1]->GetCoord(1);
		rz23 = geometry->node[iPoint_2]->GetCoord(2) - geometry->node[iPoint_1]->GetCoord(2);
		l23 = sqrt(rx23*rx23 + ry23*ry23 + rz23*rz23);
		
		for(iVar = 0; iVar < 9; iVar++)
			for(jVar = 0; jVar < 9; jVar++)
				Ktor_mat[iVar][jVar] = 0.0;
		
		/*--- Linear stiffness matrix calculation
		 |	l^2		lm		ln		-l^2	-lm		-ln		|
		 |	lm		m^2		mn		-lm		-m^2	-mn		|
		 |	ln		mn		n^2		-ln		-mn		-n^2	|
		 |	-l^2	-lm		-ln		l^2		lm		ln		|
		 |	-lm		-m^2	-mn		lm		m^2		mn		|
		 |	-ln		-mn		-n^2	ln		mn		n^2		|
		 
		 l_e = sqrt ((x_2-x_1)^2+(y_2-y_1)^2+(z_2-z_1)^2)
		 l = cos(theta) = (x_2-x_1) / l_e
		 m = cos(phi) = (y_2-y_1) / l_e
		 n = cos(psi) = (z_2-z_1) / l_e ---*/
		
		l = rx12/l12; m = ry12/l12; n = rz12/l12;
		Klin_mat[0][0] = l*l;
		Klin_mat[0][1] = l*m;
		Klin_mat[0][2] = l*n;
		Klin_mat[1][0] = Klin_mat[0][1];
		Klin_mat[1][1] = m*m;
		Klin_mat[1][2] = m*n;
		Klin_mat[2][0] = Klin_mat[0][2];
		Klin_mat[2][1] = Klin_mat[1][2];
		Klin_mat[2][2] = n*n;
		
		for(iVar = 0; iVar < 3; iVar++)
			for(jVar = 0; jVar < 3; jVar++) {
				Klin_mat[iVar][jVar+3] = -Klin_mat[iVar][jVar];
				Klin_mat[iVar+3][jVar] = Klin_mat[iVar][jVar+3];
				Klin_mat[iVar+3][jVar+3] = Klin_mat[iVar][jVar];
			}
		
		for(iVar = 0; iVar < 6; iVar++)
			for(jVar = 0; jVar < 6; jVar++)
				Ktor_mat[iVar][jVar] = 1.0/(2.0*l12)*Klin_mat[iVar][jVar];
		
		
		l = rx23/l23; m = ry23/l23; n = rz23/l23;
		Klin_mat[0][0] = l*l;
		Klin_mat[0][1] = l*m;
		Klin_mat[0][2] = l*n;
		Klin_mat[1][0] = Klin_mat[0][1];
		Klin_mat[1][1] = m*m;
		Klin_mat[1][2] = m*n;
		Klin_mat[2][0] = Klin_mat[0][2];
		Klin_mat[2][1] = Klin_mat[1][2];
		Klin_mat[2][2] = n*n;
		
		for(iVar = 0; iVar < 3; iVar++)
			for(jVar = 0; jVar < 3; jVar++) {
				Klin_mat[iVar][jVar+3] = -Klin_mat[iVar][jVar];
				Klin_mat[iVar+3][jVar] = Klin_mat[iVar][jVar+3];
				Klin_mat[iVar+3][jVar+3] = Klin_mat[iVar][jVar];
			}
		
		for(iVar = 0; iVar < 6; iVar++)
			for(jVar = 0; jVar < 6; jVar++)
				Ktor_mat[iVar+3][jVar+3] = C_lin/(2.0*l23)*Klin_mat[iVar][jVar];
		
		l = rx13/l13; m = ry13/l13; n = rz13/l13;
		Klin_mat[0][0] = l*l;
		Klin_mat[0][1] = l*m;
		Klin_mat[0][2] = l*n;
		Klin_mat[1][0] = Klin_mat[0][1];
		Klin_mat[1][1] = m*m;
		Klin_mat[1][2] = m*n;
		Klin_mat[2][0] = Klin_mat[0][2];
		Klin_mat[2][1] = Klin_mat[1][2];
		Klin_mat[2][2] = n*n;
		
		for(iVar = 0; iVar < 3; iVar++)
			for(jVar = 0; jVar < 3; jVar++) {
				Klin_mat[iVar][jVar+3] = -Klin_mat[iVar][jVar];
				Klin_mat[iVar+3][jVar] = Klin_mat[iVar][jVar+3];
				Klin_mat[iVar+3][jVar+3] = Klin_mat[iVar][jVar];
			}
		
		for(iVar = 0; iVar < 3; iVar++)
			for(jVar = 0; jVar < 3; jVar++) {
				Ktor_mat[iVar][jVar] = Ktor_mat[iVar][jVar] + C_lin/(2.0*l13)*Klin_mat[iVar][jVar];
				Ktor_mat[iVar][jVar+6] = Ktor_mat[iVar][jVar+6] + C_lin/(2.0*l13)*Klin_mat[iVar][jVar+3];
				Ktor_mat[iVar+6][jVar] = Ktor_mat[iVar+6][jVar] + C_lin/(2.0*l13)*Klin_mat[iVar+3][jVar];
				Ktor_mat[iVar+6][jVar+6] = Ktor_mat[iVar+6][jVar+6] + C_lin/(2.0*l13)*Klin_mat[iVar+3][jVar+3];
			}
		
		/*--- Storage ---*/
		for(iVar = 0; iVar < 9; iVar++)
			for(jVar = 0; jVar < 9; jVar++)
				kijk[iVar][jVar][iElem] = Ktor_mat[iVar][jVar];
		
		/*--- Calculation of the diagonal part of the "assembled" matrix ---*/
		diagk[3*iPoint_0] = diagk[3*iPoint_0] + kijk[0][0][iElem];
		diagk[3*iPoint_0+1] = diagk[3*iPoint_0+1] + kijk[1][1][iElem];
		diagk[3*iPoint_0+2] = diagk[3*iPoint_0+2] + kijk[2][2][iElem];
		diagk[3*iPoint_1] = diagk[3*iPoint_1] + kijk[3][3][iElem];
		diagk[3*iPoint_1+1] = diagk[3*iPoint_1+1] + kijk[4][4][iElem];
		diagk[3*iPoint_1+2] = diagk[3*iPoint_1+2] + kijk[5][5][iElem];
		diagk[3*iPoint_2] = diagk[3*iPoint_2] + kijk[6][6][iElem];
		diagk[3*iPoint_2+1] = diagk[3*iPoint_2+1] + kijk[7][7][iElem];
		diagk[3*iPoint_2+2] = diagk[3*iPoint_2+2] + kijk[8][8][iElem];
	}
}

void CVolumetricMovement::SetBoundary(CGeometry *geometry, CConfig *config) {
	unsigned long iVertex, Point;
	unsigned short iMarker, iDim;
	double *VarCoord;
	unsigned short nDim = geometry->GetnDim();
	
	for (iMarker=0; iMarker < config->GetnMarker_All(); iMarker++)
		if (config->GetMarker_All_Moving(iMarker) == YES)
			for(iVertex = 0; iVertex<geometry->nVertex[iMarker]; iVertex++) {
				Point = geometry->vertex[iMarker][iVertex]->GetNode();
				VarCoord = geometry->vertex[iMarker][iVertex]->GetVarCoord();
				for (iDim = 0; iDim < nDim; iDim++)
					x[nDim*Point + iDim] = VarCoord[iDim];
			}
}

void CVolumetricMovement::SetInitial_Boundary(CGeometry *geometry, CConfig *config) {
	unsigned short iMarker, Boundary;
	unsigned long iVertex, Point;

	for (iMarker=0; iMarker < config->GetnMarker_All(); iMarker++) {
		Boundary = config->GetMarker_All_Boundary(iMarker);
		for(iVertex = 0; iVertex<geometry->nVertex[iMarker]; iVertex++) {
			Point = geometry->vertex[iMarker][iVertex]->GetNode();
			double *Coord = geometry->node[Point]->GetCoord();
			Initial_Boundary[2*Point] = Coord[0]; 
			Initial_Boundary[2*Point+1] = Coord[1]; 
		}
	}
}	

void CVolumetricMovement::SetBoundary_Smooth(CGeometry *geometry, CConfig *config) {
	unsigned short iMarker, Boundary, iDim;
	unsigned long iVertex, Point;
	double smooth_coeff = 15.0;
	unsigned short nSmooth = 1000;
	unsigned long jVertex;
	tol = config->GetGridDef_Error();
	double **Coordinate_Old = NULL, **Coordinate_Sum = NULL, **Coordinate = NULL;

	for (iMarker=0; iMarker < config->GetnMarker_All(); iMarker++) {
		Boundary = config->GetMarker_All_Boundary(iMarker);
		switch(Boundary) {
		case (EULER_WALL): case (NO_SLIP_WALL):
			Coordinate_Old = new double* [geometry->nVertex[iMarker]];
			Coordinate_Sum = new double* [geometry->nVertex[iMarker]];
			Coordinate = new double* [geometry->nVertex[iMarker]];
			for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
				Coordinate_Old[iVertex] = new double [geometry->GetnDim()];
				Coordinate_Sum[iVertex] = new double [geometry->GetnDim()];
				Coordinate[iVertex] = new double [geometry->GetnDim()];
			}

			for(iVertex = 0; iVertex<geometry->nVertex[iMarker]; iVertex++) {
				unsigned long iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
				for(iDim = 0; iDim < geometry->GetnDim(); iDim++) {
					Coordinate_Old[iVertex][iDim] = geometry->node[iPoint]->GetCoord(iDim);
					Coordinate[iVertex][iDim] = geometry->node[iPoint]->GetCoord(iDim);
				}
			}

			for (unsigned short iSmooth = 0; iSmooth < nSmooth; iSmooth++) {
			double max_error = 0;

				for(iVertex = 0; iVertex<geometry->nVertex[iMarker]; iVertex++)
					for(iDim = 0; iDim < geometry->GetnDim(); iDim++)
						Coordinate_Sum[iVertex][iDim] = 0.0;

				for(iVertex = 0; iVertex<geometry->nVertex[iMarker]; iVertex++) {
					unsigned long iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
					for (unsigned short iNeighbor = 0; iNeighbor < geometry->node[iPoint]->GetnPoint(); iNeighbor++) {
						unsigned long jPoint = geometry->node[iPoint]->GetPoint(iNeighbor);
						if (geometry->node[jPoint]->GetBoundary()) {
							for(jVertex = 0; jVertex<geometry->nVertex[iMarker]; jVertex++) {
								unsigned long kPoint = geometry->vertex[iMarker][jVertex]->GetNode();
								if (kPoint==jPoint) break;
							}
							for(iDim = 0; iDim < geometry->GetnDim(); iDim++) 
								Coordinate_Sum[iVertex][iDim] += Coordinate[jVertex][iDim];
						}
					}
				}


				for(unsigned long iVertex = 0; iVertex<geometry->nVertex[iMarker]; iVertex++) {
					const double nneigh = 2.0; //Para 2D

					double error_x = 
						Coordinate[iVertex][0] - (Coordinate_Old[iVertex][0] 
						+ smooth_coeff*Coordinate_Sum[iVertex][0])/(1.0 + smooth_coeff*nneigh);

					double error_y = 
						Coordinate[iVertex][1] - (Coordinate_Old[iVertex][1] 
						+ smooth_coeff*Coordinate_Sum[iVertex][1])/(1.0 + smooth_coeff*nneigh);

					double total_error = sqrt (error_x*error_x+error_y*error_y);

					if (max_error < total_error) max_error = total_error;

					for(iDim = 0; iDim < geometry->GetnDim(); iDim++)
						Coordinate[iVertex][iDim] =(Coordinate_Old[iVertex][iDim] 
							+ smooth_coeff*Coordinate_Sum[iVertex][iDim])/(1.0 + smooth_coeff*nneigh);
				}

				if (max_error < tol) break;

			}

			break;
		}
	}


	for (iMarker=0; iMarker < config->GetnMarker_All(); iMarker++) {
		Boundary = config->GetMarker_All_Boundary(iMarker);
		switch(Boundary) {
		case (EULER_WALL): case (NO_SLIP_WALL):
			for(iVertex = 0; iVertex<geometry->nVertex[iMarker]; iVertex++) {
				Point = geometry->vertex[iMarker][iVertex]->GetNode();
				x[2*Point] = Coordinate[iVertex][0]-geometry->node[Point]->GetCoord(0); 
				x[2*Point+1] = Coordinate[iVertex][1]-geometry->node[Point]->GetCoord(1);
			}
			for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
				delete [] Coordinate_Old[iVertex];
				delete [] Coordinate_Sum[iVertex];
				delete [] Coordinate[iVertex];
			}
			delete [] Coordinate_Old;
			delete [] Coordinate_Sum;
			delete [] Coordinate;
			break;

		}	
	}

}

void CVolumetricMovement::SetSolution_Smoothing(CGeometry *geometry, CConfig *config) {
	unsigned long val_nSmooth = 1;
	double val_smooth_coeff = 0.25;
	
	/*--- Perform a Jacobi approximation to an implicit residual smoothing ---*/
	double **coord_old, **coord_sum, **coord;
	unsigned short nDim = 2;

	coord_old = new double* [geometry->GetnPoint()];
	coord_sum = new double* [geometry->GetnPoint()];
	coord = new double* [geometry->GetnPoint()];	
	for (unsigned long iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
		coord_old[iPoint] = new double [nDim];
		coord_sum[iPoint] = new double [nDim];
		coord[iPoint] = new double [nDim];		
	}

	unsigned short iDim;
	unsigned long iPoint, iEdge;
	unsigned long iVertex;
	unsigned short iMarker;	
	
	/*--- Copy the initial grid ---*/
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
		for (iDim = 0; iDim < nDim; iDim++) {
			coord_old[iPoint][iDim] = geometry->node[iPoint]->GetCoord(iDim);
			coord[iPoint][iDim] = geometry->node[iPoint]->GetCoord(iDim);
		}
	}
	
	/*--- Copy the deformed boundary ---*/
	for (iMarker=0; iMarker < config->GetnMarker_All(); iMarker++)
		if (config->GetMarker_All_Moving(iMarker) == YES)
			for(iVertex = 0; iVertex<geometry->nVertex[iMarker]; iVertex++) {
				unsigned long Point = geometry->vertex[iMarker][iVertex]->GetNode();
				for (iDim = 0; iDim < nDim; iDim++) {
					coord_old[Point][iDim] = geometry->node[Point]->GetCoord(iDim);
					coord[Point][iDim] = geometry->node[Point]->GetCoord(iDim);
				}
			}


	
	/*--- Jacobi iterations ---*/
	for (unsigned short iSmooth = 0; iSmooth < val_nSmooth; iSmooth++) {
		
		for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
			for (iDim = 0; iDim < nDim; iDim++)
				coord_sum[iPoint][iDim]= 0.0;
		
		/*--- Loop over Interior edges ---*/
		for(iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {	
			const unsigned long Point_0 = geometry->edge[iEdge]->GetNode(0);			
			const unsigned long Point_1 = geometry->edge[iEdge]->GetNode(1);
			
			/*--- Accumulate nearest neighbor residual to Res_sum for each variable ---*/
			for (iDim = 0; iDim < nDim; iDim++) {
				coord_sum[Point_0][iDim]= coord_sum[Point_0][iDim] + coord[Point_1][iDim];
				coord_sum[Point_1][iDim]= coord_sum[Point_1][iDim] + coord[Point_0][iDim];
			}
		}
		
		/*--- Loop over all mesh points (Update Residuals with averaged sum) ---*/
		for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
			const unsigned short nneigh = geometry->node[iPoint]->GetnPoint();
			for (iDim = 0; iDim < nDim; iDim++) {
				coord[iPoint][iDim] =(coord_old[iPoint][iDim] + 
								 val_smooth_coeff*coord_sum[iPoint][iDim])
				/(1.0 + val_smooth_coeff*double(nneigh));
			}
		}
		
		/*--- Copy the deformed boundary ---*/
		for (iMarker=0; iMarker < config->GetnMarker_All(); iMarker++)
			for(iVertex = 0; iVertex<geometry->nVertex[iMarker]; iVertex++) {
				unsigned long Point = geometry->vertex[iMarker][iVertex]->GetNode();
				for (iDim = 0; iDim < nDim; iDim++)
					coord[Point][iDim] = coord_old[Point][iDim];
			}		
		
	}
	
	for (unsigned long iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
		for (unsigned short iDim = 0; iDim < geometry->GetnDim(); iDim++) {
			geometry->node[iPoint]->SetCoord(iDim, coord[iPoint][iDim]);
		}
	
}



void CVolumetricMovement::SetSolution(CGeometry *geometry, CConfig *config) {
	
	unsigned long iPoint, Point, iVertex;
	unsigned short nDim = geometry->GetnDim();
	unsigned long iter;
	unsigned short iMarker;
	
	tol = config->GetGridDef_Error();
	
	Setkijk_times(geometry,x,r);
	
	for (iPoint = 0; iPoint < nDim*geometry->GetnPoint(); iPoint++) {
		r[iPoint] = -r[iPoint];
		z[iPoint] = r[iPoint] / diagk[iPoint];
	}
	
	for (iter = 0; iter < niter; iter++) {
		
		for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
			for(iVertex = 0; iVertex<geometry->nVertex[iMarker]; iVertex++) {
				Point = geometry->vertex[iMarker][iVertex]->GetNode();
				for(iPoint = nDim*Point; iPoint <= nDim*Point+1; iPoint++) {
					r[iPoint] = 0.0; 
					z[iPoint] = 0.0;
				}
			}	
		
		bknum=0.0;
		for (iPoint = 0; iPoint < nDim*geometry->GetnPoint(); iPoint++) {
			bknum = bknum + z[iPoint] * r[iPoint];
		}
		
		if (iter == 0) {
			for (iPoint = 0; iPoint < nDim*geometry->GetnPoint(); iPoint++) {
				p[iPoint] = z[iPoint]; 
			}
		}
		else {
			bk = bknum / bkden;
			for (iPoint = 0; iPoint < nDim*geometry->GetnPoint(); iPoint++)
				p[iPoint] = bk * p[iPoint] + z[iPoint]; 
		}
		
		bkden = bknum;
		
		Setkijk_times(geometry,p,z);
		
		akden=0.0;
		for (iPoint = 0; iPoint < nDim*geometry->GetnPoint(); iPoint++)
			akden = akden + z[iPoint]*p[iPoint];
		
		ak = bknum/akden;
		for (iPoint = 0; iPoint < nDim*geometry->GetnPoint(); iPoint++) {
			x[iPoint] = x[iPoint] + ak*p[iPoint];
			r[iPoint] = r[iPoint] - ak*z[iPoint];
			z[iPoint] = r[iPoint] / diagk[iPoint];
		}
		
		err = 0.0;
		for (iPoint = 0; iPoint < nDim*geometry->GetnPoint(); iPoint++)
			err = err + r[iPoint]*r[iPoint];
		
		err = sqrt(err);
		
		if (err < tol) {
			cout <<"Total number of iterations: "<< iter <<", with a final error of " << err <<"."<<endl;
			return;
		}
	}
}

void CVolumetricMovement::Setkijk_times(CGeometry *geometry, double *vect, double *res) {
	
	unsigned long iPoint_0, iPoint_1, iPoint_2, iPoint, iElem;
	unsigned short nDim = geometry->GetnDim();
	
	for (iPoint = 0; iPoint < nDim*geometry->GetnPoint(); iPoint++)
		res[iPoint] = 0.0;
	
	if (nDim == 2) {
		for (iElem = 0; iElem < nElem; iElem ++) {
			iPoint_0 = nodes[iElem][0];
			iPoint_1 = nodes[iElem][1];
			iPoint_2 = nodes[iElem][2];
			
			res[2*iPoint_0] = res[2*iPoint_0] +
			kijk[0][0][iElem]*vect[2*iPoint_0] +  kijk[0][1][iElem]*vect[2*iPoint_0+1] + 
			kijk[0][2][iElem]*vect[2*iPoint_1] +  kijk[0][3][iElem]*vect[2*iPoint_1+1] + 
			kijk[0][4][iElem]*vect[2*iPoint_2] +  kijk[0][5][iElem]*vect[2*iPoint_2+1];
			res[2*iPoint_0+1] = res[2*iPoint_0+1] +
			kijk[1][0][iElem]*vect[2*iPoint_0] +  kijk[1][1][iElem]*vect[2*iPoint_0+1] + 
			kijk[1][2][iElem]*vect[2*iPoint_1] +  kijk[1][3][iElem]*vect[2*iPoint_1+1] + 
			kijk[1][4][iElem]*vect[2*iPoint_2] +  kijk[1][5][iElem]*vect[2*iPoint_2+1];
			res[2*iPoint_1] = res[2*iPoint_1] +
			kijk[2][0][iElem]*vect[2*iPoint_0] +  kijk[2][1][iElem]*vect[2*iPoint_0+1] + 
			kijk[2][2][iElem]*vect[2*iPoint_1] +  kijk[2][3][iElem]*vect[2*iPoint_1+1] + 
			kijk[2][4][iElem]*vect[2*iPoint_2] +  kijk[2][5][iElem]*vect[2*iPoint_2+1];
			res[2*iPoint_1+1] = res[2*iPoint_1+1] +
			kijk[3][0][iElem]*vect[2*iPoint_0] +  kijk[3][1][iElem]*vect[2*iPoint_0+1] + 
			kijk[3][2][iElem]*vect[2*iPoint_1] +  kijk[3][3][iElem]*vect[2*iPoint_1+1] + 
			kijk[3][4][iElem]*vect[2*iPoint_2] +  kijk[3][5][iElem]*vect[2*iPoint_2+1];
			res[2*iPoint_2] = res[2*iPoint_2] +
			kijk[4][0][iElem]*vect[2*iPoint_0] +   kijk[4][1][iElem]*vect[2*iPoint_0+1] + 
			kijk[4][2][iElem]*vect[2*iPoint_1] +   kijk[4][3][iElem]*vect[2*iPoint_1+1] + 
			kijk[4][4][iElem]*vect[2*iPoint_2] +   kijk[4][5][iElem]*vect[2*iPoint_2+1];
			res[2*iPoint_2+1] = res[2*iPoint_2+1] +
			kijk[5][0][iElem]*vect[2*iPoint_0] +   kijk[5][1][iElem]*vect[2*iPoint_0+1] + 
			kijk[5][2][iElem]*vect[2*iPoint_1] +   kijk[5][3][iElem]*vect[2*iPoint_1+1] + 
			kijk[5][4][iElem]*vect[2*iPoint_2] +   kijk[5][5][iElem]*vect[2*iPoint_2+1];
		}
	}
	else {
		for (iElem = 0; iElem < nElem; iElem ++) {
			iPoint_0 = nodes[iElem][0];
			iPoint_1 = nodes[iElem][1];
			iPoint_2 = nodes[iElem][2];
			
			res[3*iPoint_0] = res[3*iPoint_0] +
			kijk[0][0][iElem]*vect[3*iPoint_0] +  kijk[0][1][iElem]*vect[3*iPoint_0+1] +  kijk[0][2][iElem]*vect[3*iPoint_0+2] + 
			kijk[0][3][iElem]*vect[3*iPoint_1] +  kijk[0][4][iElem]*vect[3*iPoint_1+1] +  kijk[0][5][iElem]*vect[3*iPoint_0+2] + 
			kijk[0][6][iElem]*vect[3*iPoint_2] +  kijk[0][7][iElem]*vect[3*iPoint_2+1] +  kijk[0][8][iElem]*vect[3*iPoint_0+2];
			res[3*iPoint_0+1] = res[3*iPoint_0+1] +
			kijk[1][0][iElem]*vect[3*iPoint_0] +  kijk[1][1][iElem]*vect[3*iPoint_0+1] +  kijk[1][2][iElem]*vect[3*iPoint_0+2] + 
			kijk[1][3][iElem]*vect[3*iPoint_1] +  kijk[1][4][iElem]*vect[3*iPoint_1+1] +  kijk[1][5][iElem]*vect[3*iPoint_0+2] + 
			kijk[1][6][iElem]*vect[3*iPoint_2] +  kijk[1][7][iElem]*vect[3*iPoint_2+1] +  kijk[1][8][iElem]*vect[3*iPoint_0+2];
			res[3*iPoint_0+2] = res[3*iPoint_0+2] +
			kijk[2][0][iElem]*vect[3*iPoint_0] +  kijk[2][1][iElem]*vect[3*iPoint_0+1] +  kijk[2][2][iElem]*vect[3*iPoint_0+2] + 
			kijk[2][3][iElem]*vect[3*iPoint_1] +  kijk[2][4][iElem]*vect[3*iPoint_1+1] +  kijk[2][5][iElem]*vect[3*iPoint_0+2] + 
			kijk[2][6][iElem]*vect[3*iPoint_2] +  kijk[2][7][iElem]*vect[3*iPoint_2+1] +  kijk[2][8][iElem]*vect[3*iPoint_0+2];
			
			
			res[3*iPoint_1] = res[3*iPoint_1] +
			kijk[3][0][iElem]*vect[3*iPoint_0] +  kijk[3][1][iElem]*vect[3*iPoint_0+1] +  kijk[3][2][iElem]*vect[3*iPoint_0+2] + 
			kijk[3][3][iElem]*vect[3*iPoint_1] +  kijk[3][4][iElem]*vect[3*iPoint_1+1] +  kijk[3][5][iElem]*vect[3*iPoint_0+2] + 
			kijk[3][6][iElem]*vect[3*iPoint_2] +  kijk[3][7][iElem]*vect[3*iPoint_2+1] +  kijk[3][8][iElem]*vect[3*iPoint_0+2];
			res[3*iPoint_1+1] = res[3*iPoint_1+1] +
			kijk[4][0][iElem]*vect[3*iPoint_0] +  kijk[4][1][iElem]*vect[3*iPoint_0+1] +  kijk[4][2][iElem]*vect[3*iPoint_0+2] + 
			kijk[4][3][iElem]*vect[3*iPoint_1] +  kijk[4][4][iElem]*vect[3*iPoint_1+1] +  kijk[4][5][iElem]*vect[3*iPoint_0+2] + 
			kijk[4][6][iElem]*vect[3*iPoint_2] +  kijk[4][7][iElem]*vect[3*iPoint_2+1] +  kijk[4][8][iElem]*vect[3*iPoint_0+2];
			res[3*iPoint_1+2] = res[3*iPoint_1+2] +
			kijk[5][0][iElem]*vect[3*iPoint_0] +  kijk[5][1][iElem]*vect[3*iPoint_0+1] +  kijk[5][2][iElem]*vect[3*iPoint_0+2] + 
			kijk[5][3][iElem]*vect[3*iPoint_1] +  kijk[5][4][iElem]*vect[3*iPoint_1+1] +  kijk[5][5][iElem]*vect[3*iPoint_0+2] + 
			kijk[5][6][iElem]*vect[3*iPoint_2] +  kijk[5][7][iElem]*vect[3*iPoint_2+1] +  kijk[5][8][iElem]*vect[3*iPoint_0+2];
			
			
			res[3*iPoint_2] = res[3*iPoint_2] +
			kijk[6][0][iElem]*vect[3*iPoint_0] +  kijk[6][1][iElem]*vect[3*iPoint_0+1] +  kijk[6][2][iElem]*vect[3*iPoint_0+2] + 
			kijk[6][3][iElem]*vect[3*iPoint_1] +  kijk[6][4][iElem]*vect[3*iPoint_1+1] +  kijk[6][5][iElem]*vect[3*iPoint_0+2] + 
			kijk[6][6][iElem]*vect[3*iPoint_2] +  kijk[6][7][iElem]*vect[3*iPoint_2+1] +  kijk[6][8][iElem]*vect[3*iPoint_0+2];
			res[3*iPoint_2+1] = res[3*iPoint_2+1] +
			kijk[7][0][iElem]*vect[3*iPoint_0] +  kijk[7][1][iElem]*vect[3*iPoint_0+1] +  kijk[7][2][iElem]*vect[3*iPoint_0+2] + 
			kijk[7][3][iElem]*vect[3*iPoint_1] +  kijk[7][4][iElem]*vect[3*iPoint_1+1] +  kijk[7][5][iElem]*vect[3*iPoint_0+2] + 
			kijk[7][6][iElem]*vect[3*iPoint_2] +  kijk[7][7][iElem]*vect[3*iPoint_2+1] +  kijk[7][8][iElem]*vect[3*iPoint_0+2];
			res[3*iPoint_2+2] = res[3*iPoint_2+2] +
			kijk[8][0][iElem]*vect[3*iPoint_0] +  kijk[8][1][iElem]*vect[3*iPoint_0+1] +  kijk[8][2][iElem]*vect[3*iPoint_0+2] + 
			kijk[8][3][iElem]*vect[3*iPoint_1] +  kijk[8][4][iElem]*vect[3*iPoint_1+1] +  kijk[8][5][iElem]*vect[3*iPoint_0+2] + 
			kijk[8][6][iElem]*vect[3*iPoint_2] +  kijk[8][7][iElem]*vect[3*iPoint_2+1] +  kijk[8][8][iElem]*vect[3*iPoint_0+2];
			
		}
	}
}

void CVolumetricMovement::UpdateGrid(CGeometry *geometry, CConfig *config) {
	double new_coord;
	
	unsigned long iPoint;
	unsigned short iDim;
	unsigned short nDim = geometry->GetnDim();
	
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
		for (iDim = 0; iDim < nDim; iDim++) {
			new_coord = geometry->node[iPoint]->GetCoord(iDim) + x[nDim*iPoint + iDim];
			geometry->node[iPoint]->SetCoord(iDim,new_coord);
			//			velocity = x[nDim*iPoint + iDim] / config->GetDelta_UnstTimeND();
			//			geometry->node[iPoint]->SetGridVel(iDim,velocity);
		}
	
	//	geometry->SetCG();
	//	geometry->SetControlVolume(config,UPDATE);
	//	geometry->SetBoundControlVolume(config, UPDATE);
}

void CVolumetricMovement::UpdateMultiGrid(CGeometry **geometry, CConfig *config) {
	unsigned long Fine_Point, Coarse_Point;
	unsigned short iDim, iChildren;
	double Area_Parent, Area_Children;
	unsigned short nDim = geometry[0]->GetnDim();
	double *GridVel_fine, *GridVel;
	GridVel = new double[nDim];

	for (unsigned short iMGlevel = 1; iMGlevel <=config->GetMGLevels(); iMGlevel++) {
		geometry[iMGlevel]->SetControlVolume(config,geometry[iMGlevel-1], UPDATE);
		geometry[iMGlevel]->SetBoundControlVolume(config,geometry[iMGlevel-1], UPDATE);
	}

	for (unsigned short iMesh = 0; iMesh < config->GetMGLevels(); iMesh++) {
		for (Coarse_Point = 0; Coarse_Point < geometry[iMesh+1]->GetnPoint(); Coarse_Point++) {
			Area_Parent = geometry[iMesh+1]->node[Coarse_Point]->GetVolume();

			for (iDim = 0; iDim < nDim; iDim++) GridVel[iDim] = 0.0;

			for (iChildren = 0; iChildren < 
				geometry[iMesh+1]->node[Coarse_Point]->GetnChildren_CV(); iChildren++) {
	
				Fine_Point = geometry[iMesh+1]->node[Coarse_Point]->GetChildren_CV(iChildren);
				Area_Children = geometry[iMesh]->node[Fine_Point]->GetVolume();
				GridVel_fine = geometry[iMesh]->node[Fine_Point]->GetGridVel();

				for (iDim = 0; iDim < nDim; iDim++)
					GridVel[iDim] += GridVel_fine[iDim]*Area_Children/Area_Parent;  
			}
			geometry[iMesh+1]->node[Coarse_Point]->SetGridVel(GridVel);
		}
	}
	delete [] GridVel;
}

void CVolumetricMovement::GetBoundary(CGeometry *geometry, CConfig *config, string val_filename) {

	unsigned long Point, iVertex;
	double xCoord, *Normal, Proj_Def;
	unsigned short iMarker, Boundary;

	ofstream Surface_file;
	char *cstr; cstr = new char [val_filename.size()+1];
	strcpy (cstr, val_filename.c_str());
	Surface_file.open(cstr, ios::out);

	Surface_file <<  "\"x_coord\",\"Deformation_x\",\"Deformation_y\",\"Projected_Def\" " << endl;

	for (iMarker=0; iMarker < config->GetnMarker_All(); iMarker++) {
		Boundary = config->GetMarker_All_Boundary(iMarker);
		switch(Boundary) {
		case (EULER_WALL): case (NO_SLIP_WALL):
			for(iVertex = 0; iVertex<geometry->nVertex[iMarker]; iVertex++) {
				Point = geometry->vertex[iMarker][iVertex]->GetNode();
				Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
				Proj_Def = x[2*Point]*Normal[0]+x[2*Point+1]*Normal[1];
				xCoord = geometry->node[Point]->GetCoord(0);
				Surface_file << xCoord << "," << x[2*Point] << "," << x[2*Point+1] << "," << Proj_Def << endl;
			}
			break;
		}	
	}
}


void CVolumetricMovement::Initialize_StiffMatrix_Structure(CGeometry *geometry) {
	unsigned long iPoint, nPoint = geometry->GetnPoint(), nPointDomain = geometry->GetnPointDomain();
  unsigned long *row_ptr, *col_ind, *vneighs, index, nnz;
	unsigned short iNeigh, nNeigh;
	unsigned short nSub_blocks;
  unsigned short *Sub_block_sizes;
	unsigned short nDim = geometry->GetnDim();
  
	row_ptr = new unsigned long [nPoint+1];
	
	/*--- +1 -> to include diagonal element	---*/
	row_ptr[0] = 0;
	for (iPoint = 0; iPoint < nPoint; iPoint++)
		row_ptr[iPoint+1] = row_ptr[iPoint]+(geometry->node[iPoint]->GetnPoint()+1);
	nnz = row_ptr[nPoint]; 
	
	col_ind = new unsigned long [nnz];
	vneighs = new unsigned long [MAX_NEIGHBORS];
	
	/*--- Neighbors to point iPoint to include relation of the point with 
	 itself (matrix diagonal) ---*/
	for (iPoint = 0; iPoint < nPoint; iPoint++) {
		nNeigh = geometry->node[iPoint]->GetnPoint();
		
		for (iNeigh = 0; iNeigh < nNeigh; iNeigh++)
			vneighs[iNeigh] = geometry->node[iPoint]->GetPoint(iNeigh);
		
		vneighs[nNeigh] = iPoint;
		sort(vneighs,vneighs+nNeigh+1);
		index = row_ptr[iPoint];
		for (iNeigh = 0; iNeigh <= nNeigh; iNeigh++) {
			col_ind[index] = vneighs[iNeigh];
			index++;
		}
	}
	
	/*--- solve with preconditioned conjugate gradient method ---*/
	bool preconditioner = true; 
  bool blockDiagonalJacobian = false;
  nSub_blocks = 1;
  Sub_block_sizes = new unsigned short [nSub_blocks];
  Sub_block_sizes[0] = geometry->GetnDim();
	StiffMatrix.SetIndexes(nPoint, nPointDomain, nDim, row_ptr, col_ind, nnz, preconditioner, blockDiagonalJacobian, nSub_blocks, Sub_block_sizes);

	/*--- Set StiffMatrix entries to zero ---*/
	StiffMatrix.SetValZero();

	rhs  = new double [geometry->GetnPoint() * geometry->GetnDim()];
	usol = new double [geometry->GetnPoint() * geometry->GetnDim()];
	
	/*--- Don't delete *row_ptr, *col_ind because they are asigned to the Jacobian structure ---*/
	delete[] vneighs;
}

void CVolumetricMovement::Deallocate_StiffMatrix_Structure(CGeometry *geometry) {
	
	delete[] rhs;
	delete[] usol;
	
}


double CVolumetricMovement::SetSpringMethodContributions_Edges(CGeometry *geometry) {
	unsigned short iDim, jDim, nDim = geometry->GetnDim();
	unsigned long iEdge, Point_0, Point_1;
	double *Coord_0, *Coord_1, *Edge_Vector, *Unit_Vector;
	double Length, kij, **Smatrix, MinLength;

	Edge_Vector = new double [nDim];
	Unit_Vector = new double [nDim];
	MinLength = 1E10;
	
	Smatrix = new double* [nDim];
	for (iDim = 0; iDim < nDim; iDim++)
		Smatrix[iDim] = new double [nDim];

	/*--- Compute contributions of the basic edge spring method ---*/
	for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {

		/*--- Points in edge and coordinates ---*/
		Point_0 = geometry->edge[iEdge]->GetNode(0);
		Point_1 = geometry->edge[iEdge]->GetNode(1);
		Coord_0 = geometry->node[Point_0]->GetCoord();
		Coord_1 = geometry->node[Point_1]->GetCoord();

		/*--- Compute Edge_Vector ---*/
		Length = 0;
		for (iDim = 0; iDim < nDim; iDim++) {
			Edge_Vector[iDim] = Coord_1[iDim] - Coord_0[iDim];
			Length += Edge_Vector[iDim]*Edge_Vector[iDim];
		}
		Length = sqrt(Length);
		MinLength = min(Length, MinLength);
		
		if (Length == 0.0) {
			cout << Point_0 << endl;
			cout << Point_1 << endl;
			cout << Coord_0[0] <<" "<< Coord_0[1] <<" "<< Coord_0[2] << endl;
			cout << Coord_1[0] <<" "<< Coord_1[1] <<" "<< Coord_1[2] << endl;
			cin.get();
		}

		/*--- Compute Unit_Vector ---*/
		for (iDim = 0; iDim < nDim; iDim++)
			Unit_Vector[iDim] = Edge_Vector[iDim]/Length;

		/*--- Compute spring stiffness (kij) and point-to-point matrix ---*/
		kij = 1.0/Length;
		
		for (iDim = 0; iDim < nDim; iDim++)
			for (jDim = 0; jDim < nDim; jDim++)
				Smatrix[iDim][jDim] = kij*Unit_Vector[iDim]*Unit_Vector[jDim];

		/*--- Add and substract contributions to the global matrix ---*/
		StiffMatrix.AddBlock(Point_0, Point_0, Smatrix);
		StiffMatrix.SubtractBlock(Point_0, Point_1, Smatrix);
		StiffMatrix.SubtractBlock(Point_1, Point_0, Smatrix);
		StiffMatrix.AddBlock(Point_1, Point_1, Smatrix);
	}
	
	for (iDim = 0; iDim < nDim; iDim++)
		delete [] Smatrix [iDim];
	delete [] Smatrix;
	delete [] Unit_Vector;
	delete [] Edge_Vector;
	
	return MinLength;
}

void CVolumetricMovement::SetBoundaryDisplacements(CGeometry *geometry, CConfig *config) {

	unsigned short iDim, nDim = geometry->GetnDim(), iMarker, axis = 0;
	unsigned long iPoint, total_index, iVertex;
	double *VarCoord, MeanCoord[3];
	
	/*--- As initialization, set to zero displacements of all the surfaces except the symmetry 
	 plane, and the MPI boundaries (receive). ---*/
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
		if ((config->GetMarker_All_Boundary(iMarker) != SYMMETRY_PLANE) &&
				(!((config->GetMarker_All_Boundary(iMarker) == SEND_RECEIVE) 
					 && (config->GetMarker_All_SendRecv(iMarker) > 0))))
			for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
				iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
				for (iDim = 0; iDim < nDim; iDim++) {
					total_index = iPoint*nDim + iDim;
					rhs[total_index]  = 0.0;
					usol[total_index] = 0.0;
					StiffMatrix.DeleteValsRowi(total_index);
				}
			}
	
	/*--- Set the known displacements and modify the linear system ---*/
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
		if (config->GetMarker_All_Moving(iMarker) == YES)  {
			for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
				iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
				VarCoord = geometry->vertex[iMarker][iVertex]->GetVarCoord();
				for (iDim = 0; iDim < nDim; iDim++) {
					total_index = iPoint*nDim + iDim;
					rhs[total_index]  = VarCoord[iDim];
					usol[total_index] = VarCoord[iDim];
				}
			}
    }

	/*--- Set to zero displacements of the normal component for the Symmetry plane condition ---*/
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
		if (config->GetMarker_All_Boundary(iMarker) == SYMMETRY_PLANE) {

			for (iDim = 0; iDim < nDim; iDim++) MeanCoord[iDim] = 0.0;
			for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
				iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
				VarCoord = geometry->node[iPoint]->GetCoord();
				for (iDim = 0; iDim < nDim; iDim++)
					MeanCoord[iDim] += VarCoord[iDim]*VarCoord[iDim];
			}
			for (iDim = 0; iDim < nDim; iDim++) MeanCoord[iDim] = sqrt(MeanCoord[iDim]);
			
			if ((MeanCoord[0] <= MeanCoord[1]) && (MeanCoord[0] <= MeanCoord[2])) axis = 0;
			if ((MeanCoord[1] <= MeanCoord[0]) && (MeanCoord[1] <= MeanCoord[2])) axis = 1;
			if ((MeanCoord[2] <= MeanCoord[0]) && (MeanCoord[2] <= MeanCoord[1])) axis = 2;
			
			if (nDim == 2) axis = 1;
			
			for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
				iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
				total_index = iPoint*nDim + axis;
				rhs[total_index] = 0.0; 
				usol[total_index] = 0.0;
				StiffMatrix.DeleteValsRowi(total_index);
			}
		}
	}
	
	/*--- Set to zero displacements of the near-field ---*/
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
		if (config->GetMarker_All_Boundary(iMarker) == NEARFIELD_BOUNDARY)
			for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
				iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
				for (iDim = 0; iDim < nDim; iDim++) {
					total_index = iPoint*nDim + iDim;
					rhs[total_index]  = 0.0;
					usol[total_index] = 0.0;
					StiffMatrix.DeleteValsRowi(total_index);
				}
			}
}

void CVolumetricMovement::SetDomainDisplacements(CGeometry *geometry, CConfig *config) {
	unsigned short iDim, nDim = geometry->GetnDim();
	unsigned long iPoint, total_index;
	double *Coord, MinCoordValues[3], MaxCoordValues[3], *Hold_GridFixed_Coord;
	
	Hold_GridFixed_Coord = config->GetHold_GridFixed_Coord();
	
	MinCoordValues[0] = Hold_GridFixed_Coord[0];
	MinCoordValues[1] = Hold_GridFixed_Coord[1];
	MinCoordValues[2] = Hold_GridFixed_Coord[2];
	MaxCoordValues[0] = Hold_GridFixed_Coord[3];
	MaxCoordValues[1] = Hold_GridFixed_Coord[4];
	MaxCoordValues[2] = Hold_GridFixed_Coord[5];

	/*--- Set to zero displacements of all the points that are not going to be moved
	 except the surfaces ---*/
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
		Coord = geometry->node[iPoint]->GetCoord();
		for (iDim = 0; iDim < nDim; iDim++) {
			if ((Coord[iDim] < MinCoordValues[iDim]) || (Coord[iDim] > MaxCoordValues[iDim])) {
				total_index = iPoint*nDim + iDim;
				rhs[total_index]  = 0.0;
				usol[total_index] = 0.0;
				StiffMatrix.DeleteValsRowi(total_index);
			}
		}
	}
	
	//  /*--- Hold some regions of the mesh fixed ---*/
  //	for (unsigned long iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {		
  //		unsigned long total_index = iPoint*nDim;
  //		rhs[total_index] = 0.0; 
  //		usol[total_index] = 0.0;
  //		StiffMatrix.DeleteValsRowi(total_index);	
  //    double *Coord = geometry->node[iPoint]->GetCoord();
  //    double Radius = 0.0;
  //    for (unsigned short iDim = 0; iDim < nDim; iDim++)
  //      Radius += Coord[iDim]*Coord[iDim];
  //    Radius = sqrt(Radius);
  //    
  //		if (Radius > 1.95) {
  //			for (unsigned short iDim = 0; iDim < nDim; iDim++) {
  //				unsigned long total_index = iPoint*nDim + iDim;
  //				rhs[total_index]  = 0.0;
  //				usol[total_index] = 0.0;
  //				StiffMatrix.DeleteValsRowi(total_index);
  //			}
  //		}
  //  }
	
}

void CVolumetricMovement::UpdateSpringGrid(CGeometry *geometry, CConfig *config) {
	unsigned long iPoint, total_index;
	double new_coord;

	unsigned short iDim, nDim = geometry->GetnDim();
	unsigned long nPoint = geometry->GetnPoint();
  
	for (iPoint = 0; iPoint < nPoint; iPoint++)
		for (iDim = 0; iDim < nDim; iDim++) {
			total_index = iPoint*nDim + iDim;
			new_coord = geometry->node[iPoint]->GetCoord(iDim) + usol[total_index];
			if (fabs(new_coord) < EPS*EPS) new_coord = 0.0;
			geometry->node[iPoint]->SetCoord(iDim, new_coord);
		}

}


void CVolumetricMovement::AlgebraicMethod(CGeometry *geometry, CConfig *config, bool UpdateGeo) {
	unsigned long iPoint, jPoint, iVertex;
	unsigned short iMarker, iDim;
	double *iCoord, *jCoord, *DeltaS, *jNorm, jArea, jProjVec, *NuFunc, jLength, Up[3], 
	Low, *NewCoordx, *NewCoordy, *NewCoordz, dist, alpha;
	
	/*--- This parameter limit the area of influence of a perturbation ---*/
	alpha = 100.0;

	/*--- Compute the nu function (The nu function is defined so that the deformation is smoothly 
	 reduced as the distance to the body increases ---*/
	NuFunc = new double [geometry->GetnPoint()]; 
	NewCoordx = new double [geometry->GetnPoint()]; 
	NewCoordy = new double [geometry->GetnPoint()]; 
	NewCoordz = new double [geometry->GetnPoint()]; 
	
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
		iCoord = geometry->node[iPoint]->GetCoord();
		dist = (iCoord[nDim-1]+0.5)*(iCoord[nDim-1]+0.5);
		dist = sqrt(dist);
		if (dist <= 0.2) NuFunc[iPoint] = 1.0;
		if ((dist > 0.2) && (dist < 0.4))
			NuFunc[iPoint] = -5.0*dist + 2.0;
		if (dist >= 0.4) NuFunc[iPoint] = 0.0;		
	}
	
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
		iCoord = geometry->node[iPoint]->GetCoord();
		
		Up[0] = 0.0; Up[1] = 0.0; Up[2] = 0.0; Low = 0.0;
		for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
			if (config->GetMarker_All_Moving(iMarker) == YES)		
				for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
					jPoint = geometry->vertex[iMarker][iVertex]->GetNode();
					jCoord = geometry->node[jPoint]->GetCoord();
					jNorm = geometry->vertex[iMarker][iVertex]->GetNormal();				
					DeltaS = geometry->vertex[iMarker][iVertex]->GetVarCoord();
					
					jArea = 0.0; 
					for (iDim = 0; iDim < geometry->GetnDim(); iDim++)
						jArea += jNorm[iDim]*jNorm[iDim];
					jArea = sqrt(jArea); 
						
					jProjVec = 0.0; jLength = 0.0;
					for (iDim = 0; iDim < geometry->GetnDim(); iDim++) {
						jProjVec += (iCoord[iDim]-jCoord[iDim])*jNorm[iDim]/jArea;
						jLength += (iCoord[iDim]-jCoord[iDim])*(iCoord[iDim]-jCoord[iDim]);
					}
					jLength = sqrt(jLength);		
					
					if (jPoint != iPoint) {
						
						/*--- Compute the numerator of the algeabric method ---*/
						for (iDim = 0; iDim < geometry->GetnDim(); iDim++)
							Up[iDim] += jArea * DeltaS[iDim] * (jProjVec + 1.0) / pow(jLength, alpha);
						
						/*--- Compute the denominator of the algeabric method ---*/
						Low += jArea * (jProjVec + 1.0) / pow(jLength, alpha);
					}
					
				}
						
		/*--- Evaluate the new position of the interior point ---*/
		NewCoordx[iPoint] = geometry->node[iPoint]->GetCoord(0) + NuFunc[iPoint]*Up[0] / Low;
		NewCoordy[iPoint] = geometry->node[iPoint]->GetCoord(1) + NuFunc[iPoint]*Up[1] / Low;
		if (geometry->GetnDim() == 3)
			NewCoordz[iPoint] = geometry->node[iPoint]->GetCoord(2) + NuFunc[iPoint]*Up[2] / Low;
		
	}
	
/*	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
		if (config->GetMarker_All_Moving(iMarker) == YES)		
			for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
				iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
				DeltaS = geometry->vertex[iMarker][iVertex]->GetVarCoord();
				
				NewCoordx[iPoint] = geometry->node[iPoint]->GetCoord(0) + DeltaS[0];
				NewCoordy[iPoint] = geometry->node[iPoint]->GetCoord(1) + DeltaS[1];
				if (geometry->GetnDim() == 3)
					NewCoordz[iPoint] = geometry->node[iPoint]->GetCoord(2) + DeltaS[2];
			}		*/		
	
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
		geometry->node[iPoint]->SetCoord(0, NewCoordx[iPoint]);
		geometry->node[iPoint]->SetCoord(1, NewCoordy[iPoint]);
		if (geometry->GetnDim() == 3)
			geometry->node[iPoint]->SetCoord(2, NewCoordz[iPoint]);
	}
		
	if (UpdateGeo) {
		geometry->SetCG();
		geometry->SetControlVolume(config, UPDATE);
		geometry->SetBoundControlVolume(config, UPDATE);
	}
	
	delete [] NuFunc; 
	delete [] NewCoordx; 
	delete [] NewCoordy; 
	delete [] NewCoordz; 

}

void CVolumetricMovement::SpringMethod(CGeometry *geometry, CConfig *config, bool UpdateGeo) {
	double MinLength, NumError;
	unsigned long iPoint, total_index;
	unsigned short iDim;
	
	Initialize_StiffMatrix_Structure(geometry);
	
	StiffMatrix.SetValZero();
	
	MinLength = SetSpringMethodContributions_Edges(geometry);
	
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++)
		for (iDim = 0; iDim < geometry->GetnDim(); iDim++) {
			total_index = iPoint*geometry->GetnDim() + iDim;
			rhs[total_index]  = 0.0;
			usol[total_index] = 0.0;
		}
	
	SetBoundaryDisplacements(geometry, config);
  
	if (config->GetHold_GridFixed())
		SetDomainDisplacements(geometry, config);
	
	NumError = config->GetGridDef_Error();
	if (NumError > MinLength) { 
		cout << "The numerical error is greater than the length of the smallest edge! MinLength: " << MinLength <<"."<<endl;	
		cin.get(); 
	}
    	
	StiffMatrix.CGSolution(rhs, usol, NumError, 999999, true, geometry, config);
	
  UpdateSpringGrid(geometry, config);
  
	if (UpdateGeo) {
		geometry->SetCG();
		geometry->SetControlVolume(config, UPDATE);
		geometry->SetBoundControlVolume(config, UPDATE);
	}
  
	Deallocate_StiffMatrix_Structure(geometry);

}

void CVolumetricMovement::TorsionalSpringMethod(CGeometry *geometry, CConfig *config, bool UpdateGeo) {
	
	Set2DMatrix_Structure(geometry);
	SetBoundary(geometry, config);
	SetSolution(geometry, config);
	UpdateGrid(geometry, config);
	
	if (UpdateGeo) {
		geometry->SetCG();
		geometry->SetControlVolume(config, UPDATE);
		geometry->SetBoundControlVolume(config, UPDATE);
	}
	
}

void CVolumetricMovement::SetRigidRotation(CGeometry *geometry, CConfig *config,
                                           unsigned short iZone, unsigned long iter) {
  
  int rank = MASTER_NODE;
#ifndef NO_MPI
	rank = MPI::COMM_WORLD.Get_rank();
#endif
  
  /*--- Local variables ---*/
	unsigned short iDim, nDim; 
	unsigned long iPoint;
	double r[3], rotCoord[3], *Coord, Center[3], Omega[3], Lref, dt;
  double rotMatrix[3][3] = {{0.0,0.0,0.0}, {0.0,0.0,0.0}, {0.0,0.0,0.0}};
	double dtheta, dphi, dpsi, cosTheta, sinTheta;
  double cosPhi, sinPhi, cosPsi, sinPsi;
  double DEG2RAD = PI_NUMBER/180.0;
  bool adjoint = config->IsAdjoint();

	/*--- Problem dimension and physical time step ---*/
	nDim = geometry->GetnDim();
	dt   = config->GetDelta_UnstTimeND();
  Lref = config->GetLength_Ref();

  /*--- For the unsteady adjoint, use reverse time ---*/
  if (adjoint) {
    /*--- Set the first adjoint mesh position to the final direct one ---*/
    if (iter == 0) dt = ((double)config->GetnExtIter()-1)*dt;
    /*--- Reverse the rotation direction for the adjoint ---*/
    else dt = -1.0*dt;
  } else {
    /*--- No rotation at all for the first direct solution ---*/
    if (iter == 0) dt = 0;
  }
  
	/*--- Center of rotation & angular velocity vector from config ---*/
  Center[0] = config->GetMotion_Origin_X(iZone);
  Center[1] = config->GetMotion_Origin_Y(iZone);
  Center[2] = config->GetMotion_Origin_Z(iZone);
  Omega[0]  = config->GetRotation_Rate_X(iZone);
  Omega[1]  = config->GetRotation_Rate_Y(iZone);
  Omega[2]  = config->GetRotation_Rate_Z(iZone);
  
	/*--- Compute delta change in the angle about the x, y, & z axes. ---*/
  
	dtheta = Omega[0]*dt;   
	dphi   = Omega[1]*dt; 
	dpsi   = Omega[2]*dt;
  
  if (rank == MASTER_NODE) {
    cout.precision(4);
		cout << "Delta rotation angles (about x, y, z axes): (";
    cout << dtheta/DEG2RAD << ", ";
    cout << dphi/DEG2RAD << ", ";
    cout << dpsi/DEG2RAD << ") degrees." << endl;
  }
  
	/*--- Store angles separately for clarity. Compute sines/cosines. ---*/
  
	cosTheta = cos(dtheta);  cosPhi = cos(dphi);  cosPsi = cos(dpsi);
	sinTheta = sin(dtheta);  sinPhi = sin(dphi);  sinPsi = sin(dpsi);
  
	/*--- Compute the rotation matrix. Note that the implicit
   ordering is rotation about the x-axis, y-axis, then z-axis. ---*/
  
	rotMatrix[0][0] = cosPhi*cosPsi;
	rotMatrix[1][0] = cosPhi*sinPsi;
	rotMatrix[2][0] = -sinPhi;
  
	rotMatrix[0][1] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;
	rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;
	rotMatrix[2][1] = sinTheta*cosPhi;
  
	rotMatrix[0][2] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
	rotMatrix[1][2] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
	rotMatrix[2][2] = cosTheta*cosPhi;
  
	/*--- Loop over and rotate each node in the volume mesh ---*/
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
    
    /*--- Coordinates of the current point ---*/
    Coord = geometry->node[iPoint]->GetCoord();
    
    /*--- Calculate non-dim. position from rotation center ---*/
    for (iDim = 0; iDim < nDim; iDim++)
      r[iDim] = (Coord[iDim]-Center[iDim])/Lref;
    if (nDim == 2) r[nDim] = 0.0;
    
    /*--- Compute transformed point coordinates ---*/
    rotCoord[0] = rotMatrix[0][0]*r[0] 
                + rotMatrix[0][1]*r[1] 
                + rotMatrix[0][2]*r[2] + Center[0];
    
    rotCoord[1] = rotMatrix[1][0]*r[0] 
                + rotMatrix[1][1]*r[1] 
                + rotMatrix[1][2]*r[2] + Center[1];
    
    rotCoord[2] = rotMatrix[2][0]*r[0] 
                + rotMatrix[2][1]*r[1] 
                + rotMatrix[2][2]*r[2] + Center[2];
    
    /*--- Store new node location ---*/
    for (iDim = 0; iDim < nDim; iDim++) {
      geometry->node[iPoint]->SetCoord(iDim,rotCoord[iDim]);
    }
  }
  
	/*--- After moving all nodes, update geometry class ---*/
	geometry->SetCG();
	geometry->SetControlVolume(config, UPDATE);
	geometry->SetBoundControlVolume(config, UPDATE);

}

void CVolumetricMovement::SetRigidPitching(CGeometry *geometry, CConfig *config, unsigned short iZone, unsigned long iter) {
  
  int rank = MASTER_NODE;
#ifndef NO_MPI
	rank = MPI::COMM_WORLD.Get_rank();
#endif
  
  /*--- Local variables ---*/
  double r[3], rotCoord[3],*Coord, Center[3], Omega[3], Ampl[3], Phase[3];
  double Lref, deltaT;
  double rotMatrix[3][3] = {{0.0,0.0,0.0}, {0.0,0.0,0.0}, {0.0,0.0,0.0}};
  double dtheta, dphi, dpsi, cosTheta, sinTheta;
  double cosPhi, sinPhi, cosPsi, sinPsi;
  double time_new, time_old;
  double DEG2RAD = PI_NUMBER/180.0;
  unsigned short iDim;
  unsigned short nDim = geometry->GetnDim();
  unsigned long iPoint;
  bool time_spectral = (config->GetUnsteady_Simulation() == TIME_SPECTRAL);
  bool adjoint = config->IsAdjoint();
	
  /*--- Retrieve values from the config file ---*/
  deltaT = config->GetDelta_UnstTimeND(); 
  Lref   = config->GetLength_Ref();

  /*--- For time-spectral, motion is the same in each zone (at each instance). ---*/
  if (time_spectral) {
	  iZone = ZONE_0;
  }

  /*--- Pitching origin, frequency, and amplitude from config. ---*/	
  Center[0] = config->GetMotion_Origin_X(iZone);
  Center[1] = config->GetMotion_Origin_Y(iZone);
  Center[2] = config->GetMotion_Origin_Z(iZone);
  Omega[0]  = config->GetPitching_Omega_X(iZone)/config->GetOmega_Ref();
  Omega[1]  = config->GetPitching_Omega_Y(iZone)/config->GetOmega_Ref();
  Omega[2]  = config->GetPitching_Omega_Z(iZone)/config->GetOmega_Ref();
  Ampl[0]   = config->GetPitching_Ampl_X(iZone)*DEG2RAD;
  Ampl[1]   = config->GetPitching_Ampl_Y(iZone)*DEG2RAD;
  Ampl[2]   = config->GetPitching_Ampl_Z(iZone)*DEG2RAD;
  Phase[0]   = config->GetPitching_Phase_X(iZone)*DEG2RAD;
  Phase[1]   = config->GetPitching_Phase_Y(iZone)*DEG2RAD;
  Phase[2]   = config->GetPitching_Phase_Z(iZone)*DEG2RAD;

  if (time_spectral) {    
	  /*--- period of oscillation & compute time interval using nTimeInstances ---*/
	  double Omega_mag = sqrt(pow(Omega[0],2)+pow(Omega[1],2)+pow(Omega[2],2));
	  double period = 2*PI_NUMBER/Omega_mag;
	  deltaT = period/(double)(config->GetnTimeInstances());
  }

  /*--- Compute delta time based on physical time step ---*/
  if (adjoint) {
    /*--- For the unsteady adjoint, we integrate backwards through
     physical time, so perform mesh motion in reverse. ---*/ 
    unsigned long nFlowIter  = config->GetnExtIter();
    unsigned long directIter = nFlowIter - iter - 1;
    time_new = static_cast<double>(directIter)*deltaT;
    time_old = time_new;
    if (iter != 0) time_old = (static_cast<double>(directIter)+1.0)*deltaT;
  } else {
    /*--- Forward time for the direct problem ---*/
    time_new = static_cast<double>(iter)*deltaT;
    if (time_spectral) {
    	/*--- For time-spectral, begin movement from the zero position ---*/
    	time_old = 0.0;
    } else {
    	time_old = time_new;
    	if (iter != 0) time_old = (static_cast<double>(iter)-1.0)*deltaT;
    }
  }
  
	/*--- Compute delta change in the angle about the x, y, & z axes. ---*/
  
	dtheta = -Ampl[0]*(sin(Omega[0]*time_new + Phase[0]) - sin(Omega[0]*time_old + Phase[0]));
	dphi   = -Ampl[1]*(sin(Omega[1]*time_new + Phase[1]) - sin(Omega[1]*time_old + Phase[1]));
	dpsi   = -Ampl[2]*(sin(Omega[2]*time_new + Phase[2]) - sin(Omega[2]*time_old + Phase[2]));
  
  if (rank == MASTER_NODE) {
    cout.precision(4);
		cout << "New pitching angles (about x, y, z axes): (";
    cout << Ampl[0]*sin(Omega[0]*time_new + Phase[0])/DEG2RAD << ", ";
    cout << Ampl[1]*sin(Omega[1]*time_new + Phase[1])/DEG2RAD << ", ";
    cout << Ampl[2]*sin(Omega[2]*time_new + Phase[2])/DEG2RAD << ") ";
    cout << "degrees." << endl;
  }
  
	/*--- Store angles separately for clarity. Compute sines/cosines. ---*/
  
	cosTheta = cos(dtheta);  cosPhi = cos(dphi);  cosPsi = cos(dpsi);
	sinTheta = sin(dtheta);  sinPhi = sin(dphi);  sinPsi = sin(dpsi);
  
	/*--- Compute the rotation matrix. Note that the implicit
   ordering is rotation about the x-axis, y-axis, then z-axis. ---*/
  
	rotMatrix[0][0] = cosPhi*cosPsi;
	rotMatrix[1][0] = cosPhi*sinPsi;
	rotMatrix[2][0] = -sinPhi;
  
	rotMatrix[0][1] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;
	rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;
	rotMatrix[2][1] = sinTheta*cosPhi;
  
	rotMatrix[0][2] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
	rotMatrix[1][2] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
	rotMatrix[2][2] = cosTheta*cosPhi;
  
	/*--- Loop over and rotate each node in the volume mesh ---*/
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
    
    /*--- Coordinates of the current point ---*/
    Coord = geometry->node[iPoint]->GetCoord();
    
    /*--- Calculate non-dim. position from rotation center ---*/
    for (iDim = 0; iDim < nDim; iDim++)
      r[iDim] = (Coord[iDim]-Center[iDim])/Lref;
    if (nDim == 2) r[nDim] = 0.0;
    
    /*--- Compute transformed point coordinates ---*/
    rotCoord[0] = rotMatrix[0][0]*r[0] 
                + rotMatrix[0][1]*r[1] 
                + rotMatrix[0][2]*r[2] + Center[0];
    
    rotCoord[1] = rotMatrix[1][0]*r[0] 
                + rotMatrix[1][1]*r[1] 
                + rotMatrix[1][2]*r[2] + Center[1];
    
    rotCoord[2] = rotMatrix[2][0]*r[0] 
                + rotMatrix[2][1]*r[1] 
                + rotMatrix[2][2]*r[2] + Center[2];
    
    /*--- Store new node location ---*/
    for (iDim = 0; iDim < nDim; iDim++) {
      geometry->node[iPoint]->SetCoord(iDim,rotCoord[iDim]);
    }
  }
  
	/*--- After moving all nodes, update geometry class ---*/
	geometry->SetCG();
	geometry->SetControlVolume(config, UPDATE);
	geometry->SetBoundControlVolume(config, UPDATE);
  
}

void CVolumetricMovement::SetRigidPlunging(CGeometry *geometry, CConfig *config, unsigned short iZone, unsigned long iter) {
  
  int rank = MASTER_NODE;
#ifndef NO_MPI
	rank = MPI::COMM_WORLD.Get_rank();
#endif
  
  /*--- Local variables ---*/
  double deltaX[3], newCoord[3], Center[3], *Coord, Omega[3], Ampl[3], Lref;
  double deltaT, time_new, time_old;
  unsigned short iDim, nDim = geometry->GetnDim();
  unsigned long iPoint;
  bool time_spectral = (config->GetUnsteady_Simulation() == TIME_SPECTRAL);
  bool adjoint = config->IsAdjoint();
	
  /*--- Retrieve values from the config file ---*/
  deltaT = config->GetDelta_UnstTimeND();
  Lref   = config->GetLength_Ref();
  
  /*--- For time-spectral, motion is the same in each zone (at each instance). ---*/
  if (time_spectral) {
	  iZone = ZONE_0;
  }
  
  /*--- Plunging frequency and amplitude from config. ---*/
  Center[0] = config->GetMotion_Origin_X(iZone);
  Center[1] = config->GetMotion_Origin_Y(iZone);
  Center[2] = config->GetMotion_Origin_Z(iZone);
  Omega[0]  = config->GetPlunging_Omega_X(iZone)/config->GetOmega_Ref();
  Omega[1]  = config->GetPlunging_Omega_Y(iZone)/config->GetOmega_Ref();
  Omega[2]  = config->GetPlunging_Omega_Z(iZone)/config->GetOmega_Ref();
  Ampl[0]   = config->GetPlunging_Ampl_X(iZone)/Lref;
  Ampl[1]   = config->GetPlunging_Ampl_Y(iZone)/Lref;
  Ampl[2]   = config->GetPlunging_Ampl_Z(iZone)/Lref;
  
  if (time_spectral) {
	  /*--- period of oscillation & compute time interval using nTimeInstances ---*/
	  double Omega_mag = sqrt(pow(Omega[0],2)+pow(Omega[1],2)+pow(Omega[2],2));
	  double period = 2*PI_NUMBER/Omega_mag;
	  deltaT = period/(double)(config->GetnTimeInstances());
  }
  
  /*--- Compute delta time based on physical time step ---*/
  if (adjoint) {
    /*--- For the unsteady adjoint, we integrate backwards through
     physical time, so perform mesh motion in reverse. ---*/
    unsigned long nFlowIter  = config->GetnExtIter();
    unsigned long directIter = nFlowIter - iter - 1;
    time_new = static_cast<double>(directIter)*deltaT;
    time_old = time_new;
    if (iter != 0) time_old = (static_cast<double>(directIter)+1.0)*deltaT;
  } else {
    /*--- Forward time for the direct problem ---*/
    time_new = static_cast<double>(iter)*deltaT;
    if (time_spectral) {
    	/*--- For time-spectral, begin movement from the zero position ---*/
    	time_old = 0.0;
    } else {
    	time_old = time_new;
    	if (iter != 0) time_old = (static_cast<double>(iter)-1.0)*deltaT;
    }
  }
  
	/*--- Compute delta change in the position in the x, y, & z directions. ---*/
	deltaX[0] = -Ampl[0]*(sin(Omega[0]*time_new) - sin(Omega[0]*time_old));
	deltaX[1] = -Ampl[1]*(sin(Omega[1]*time_new) - sin(Omega[1]*time_old));
	deltaX[2] = -Ampl[2]*(sin(Omega[2]*time_new) - sin(Omega[2]*time_old));
  
  if (rank == MASTER_NODE) {
    cout.precision(4);
		cout << "Delta plunging increments (dx, dy, dz): (";
    cout << deltaX[0] << ", ";
    cout << deltaX[1] << ", ";
    cout << deltaX[2] << ")." << endl;
  }
  
	/*--- Loop over and move each node in the volume mesh ---*/
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
    
    /*--- Coordinates of the current point ---*/
    Coord = geometry->node[iPoint]->GetCoord();
    
    /*--- Increment the node position using the delta values. ---*/
    for (iDim = 0; iDim < nDim; iDim++)
      newCoord[iDim] = Coord[iDim] + deltaX[iDim];
    
    /*--- Store new node location ---*/
    for (iDim = 0; iDim < nDim; iDim++) {
      geometry->node[iPoint]->SetCoord(iDim,newCoord[iDim]);
    }
  }
  
  /*--- Set the mesh motion center to the new location after
   incrementing the position with the rigid translation. This
   new loation will be used for subsequent pitching/rotation.---*/
  config->SetMotion_Origin_X(iZone,Center[0]+deltaX[0]);
  config->SetMotion_Origin_Y(iZone,Center[1]+deltaX[1]);
  config->SetMotion_Origin_Z(iZone,Center[2]+deltaX[2]);
  
	/*--- After moving all nodes, update geometry class ---*/
	geometry->SetCG();
	geometry->SetControlVolume(config, UPDATE);
	geometry->SetBoundControlVolume(config, UPDATE);
  
}

CSurfaceMovement::CSurfaceMovement(void) : CGridMovement() {
	nChunk = 0;
	ChunkDefinition = false;
}

CSurfaceMovement::~CSurfaceMovement(void) {}

void CSurfaceMovement::CopyBoundary(CGeometry *geometry, CConfig *config) {
	unsigned short iMarker;
	unsigned long iVertex, iPoint;
	double *Coord;
	
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
		for(iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {	
			iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
			Coord = geometry->node[iPoint]->GetCoord();
			geometry->vertex[iMarker][iVertex]->SetCoord(Coord);
		}
}

void CSurfaceMovement::SetParametricCoord(CGeometry *geometry, CConfig *config, CFreeFormChunk *chunk, unsigned short iChunk) {
	unsigned short iMarker, iDim;
	unsigned long iVertex, iPoint;
	double *car_coord, *car_coord_new, *par_coord, guess[3], max_diff, 
	my_max_diff = 0.0, diff;
	
#ifndef NO_MPI
	int rank = MPI::COMM_WORLD.Get_rank();
#else
	int rank = MASTER_NODE;
#endif
	
	guess[0] = 0.5; guess[1] = 0.5; guess[2] = 0.5;
		
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
		if (config->GetMarker_All_Moving(iMarker) == YES)
			for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
				car_coord = geometry->vertex[iMarker][iVertex]->GetCoord();
				iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
				
				/*--- If the point is inside the FFD, compute the value of the parametric coordinate ---*/
				if (chunk->GetPointFFD(geometry, config, iPoint)) {
					
					/*--- Find the parametric coordinate ---*/
					par_coord = chunk->GetParametricCoord_Iterative(car_coord, guess, 1E-10, 99999);
					
					/*--- If the parametric coordinates are in (0,1) the point belongs to the chunk ---*/
					if (((par_coord[0] >= - EPS) && (par_coord[0] <= 1.0 + EPS)) && 
							((par_coord[1] >= - EPS) && (par_coord[1] <= 1.0 + EPS)) && 
							((par_coord[2] >= - EPS) && (par_coord[2] <= 1.0 + EPS))) {
						
						/*--- Set the value of the parametric coordinate ---*/
						chunk->Set_MarkerIndex(iMarker);
						chunk->Set_VertexIndex(iVertex);
						chunk->Set_PointIndex(iPoint);
						chunk->Set_ParametricCoord(par_coord);
						chunk->Set_CartesianCoord(car_coord);						
						
						/*--- Compute the cartesian coordinates using the parametric coordinates 
						 to check that everithing is right ---*/
						car_coord_new = chunk->EvalCartesianCoord(par_coord);
						
						/*--- Compute max difference between original value and the recomputed value ---*/
						diff = 0.0; 
						for (iDim = 0; iDim < geometry->GetnDim(); iDim++)
							diff += (car_coord_new[iDim]-car_coord[iDim])*(car_coord_new[iDim]-car_coord[iDim]);
						diff = sqrt(diff);
						my_max_diff = max(my_max_diff, diff);
						
						guess[0] = par_coord[0]; guess[1] = par_coord[1]; guess[2] = par_coord[2];
					}
				}
			}
		
#ifndef NO_MPI
	MPI::COMM_WORLD.Allreduce(&my_max_diff, &max_diff, 1, MPI::DOUBLE, MPI::MAX); 	
#else
	max_diff = my_max_diff;
#endif
	
	if (rank == MASTER_NODE) 
		cout << "Compute parametric coord      | FFD box: " << chunk->GetTag() << ". Max diff: " << max_diff <<"."<< endl;
	
}

void CSurfaceMovement::SetParametricCoordCP(CGeometry *geometry, CConfig *config, CFreeFormChunk *ChunkParent, CFreeFormChunk *ChunkChild) {
	unsigned short iOrder, jOrder, kOrder;
	double *car_coord, *par_coord, guess[3];

#ifndef NO_MPI
	int rank = MPI::COMM_WORLD.Get_rank();
#else
	int rank = MASTER_NODE;
#endif
	
	for (iOrder = 0; iOrder < ChunkChild->GetlOrder(); iOrder++)
		for (jOrder = 0; jOrder < ChunkChild->GetmOrder(); jOrder++)
			for (kOrder = 0; kOrder < ChunkChild->GetnOrder(); kOrder++) {
				car_coord = ChunkChild->GetCoordControlPoints(iOrder, jOrder, kOrder);
				par_coord = ChunkParent->GetParametricCoord_Iterative(car_coord, guess, 1E-10, 99999);
				ChunkChild->SetParCoordControlPoints(par_coord, iOrder, jOrder, kOrder);
			}

	if (rank == MASTER_NODE)
		cout << "Compute parametric coord (CP) | FFD parent box: " << ChunkParent->GetTag() << ". FFD child box: " << ChunkChild->GetTag() <<"."<< endl;


}

void CSurfaceMovement::GetCartesianCoordCP(CGeometry *geometry, CConfig *config, CFreeFormChunk *ChunkParent, CFreeFormChunk *ChunkChild) {
	unsigned short iOrder, jOrder, kOrder, iDim;
	double *car_coord, *par_coord;
	
#ifndef NO_MPI
	int rank = MPI::COMM_WORLD.Get_rank();
#else
	int rank = MASTER_NODE;
#endif
		
	for (iOrder = 0; iOrder < ChunkChild->GetlOrder(); iOrder++)
		for (jOrder = 0; jOrder < ChunkChild->GetmOrder(); jOrder++)
			for (kOrder = 0; kOrder < ChunkChild->GetnOrder(); kOrder++) {
				par_coord = ChunkChild->GetParCoordControlPoints(iOrder, jOrder, kOrder);
				
				/*--- Clip the value of the parametric coordinates (just in case)  ---*/
				for (iDim = 0; iDim < 3; iDim++) {
					if (par_coord[iDim] >= 1.0) par_coord[iDim] = 1.0;
					if (par_coord[iDim] <= 0.0) par_coord[iDim] = 0.0;
				}

				car_coord = ChunkParent->EvalCartesianCoord(par_coord);
				ChunkChild->SetCoordControlPoints(car_coord, iOrder, jOrder, kOrder);
			}
	
	if (rank == MASTER_NODE)
		cout << "Update cartesian coord (CP)   | FFD parent box: " << ChunkParent->GetTag() << ". FFD child box: " << ChunkChild->GetTag() <<"."<< endl;

}


void CSurfaceMovement::UpdateParametricCoord(CGeometry *geometry, CConfig *config, CFreeFormChunk *chunk, unsigned short iChunk) {
	unsigned short iMarker, iDim;
	unsigned long iVertex, iPoint, iSurfacePoints;
	double car_coord[3], *car_coord_new, *car_coord_old, *par_coord, *var_coord, guess[3], max_diff, 
	my_max_diff = 0.0, diff;
	
#ifndef NO_MPI
	int rank = MPI::COMM_WORLD.Get_rank();
#else
	int rank = MASTER_NODE;
#endif
			
	/*--- Recompute the parametric coordinates ---*/
	for (iSurfacePoints = 0; iSurfacePoints < chunk->GetnSurfacePoint(); iSurfacePoints++) {
		
		/*--- Get the marker of the surface point ---*/
		iMarker = chunk->Get_MarkerIndex(iSurfacePoints);
		
		if (config->GetMarker_All_Moving(iMarker) == YES) {
			
			/*--- Get the vertex of the surface point ---*/
			iVertex = chunk->Get_VertexIndex(iSurfacePoints);
			iPoint = chunk->Get_PointIndex(iSurfacePoints);
	
			/*--- Get the parametric and cartesians coordinates of the 
			 surface point (they don't mach) ---*/
			par_coord = chunk->Get_ParametricCoord(iSurfacePoints);
			
			/*--- Compute and set the cartesian coord using the variation computed 
			 with the previous deformation ---*/
			var_coord = geometry->vertex[iMarker][iVertex]->GetVarCoord();
			car_coord_old = geometry->node[iPoint]->GetCoord();
			for (iDim = 0; iDim < 3; iDim++)
				car_coord[iDim] = car_coord_old[iDim] + var_coord[iDim];
			chunk->Set_CartesianCoord(car_coord, iSurfacePoints);

			/*--- Find the parametric coordinate using as guess the previous value ---*/	
			guess[0] = par_coord[0]; guess[1] = par_coord[1]; guess[2] = par_coord[2];
			par_coord = chunk->GetParametricCoord_Iterative(car_coord, guess, 1E-10, 99999);
					
			/*--- Set the new value of the parametric coordinates ---*/
			chunk->Set_ParametricCoord(par_coord, iSurfacePoints);
			
			/*--- Compute the cartesian coordinates using the parametric coordinates 
			 to check that everithing is right ---*/
			car_coord_new = chunk->EvalCartesianCoord(par_coord);
			
			/*--- Compute max difference between original value and the recomputed value ---*/
			diff = 0.0; 
			for (iDim = 0; iDim < geometry->GetnDim(); iDim++)
				diff += (car_coord_new[iDim]-car_coord[iDim])*(car_coord_new[iDim]-car_coord[iDim]);
			diff = sqrt(diff);
			my_max_diff = max(my_max_diff, diff);
				
		}
	}
		
#ifndef NO_MPI
	MPI::COMM_WORLD.Allreduce(&my_max_diff, &max_diff, 1, MPI::DOUBLE, MPI::MAX); 	
#else
	max_diff = my_max_diff;
#endif
	
	if (rank == MASTER_NODE) 
		cout << "Update parametric coord       | FFD box: " << chunk->GetTag() << ". Max diff: " << max_diff <<"."<< endl;
	
}

void CSurfaceMovement::SetCartesianCoord(CGeometry *geometry, CConfig *config, CFreeFormChunk *chunk, unsigned short iChunk) {
	double *car_coord_old, *car_coord_new, diff, my_max_diff = 0.0, max_diff,
	*par_coord, VarCoord[3];
	unsigned short iMarker, iDim;
	unsigned long iVertex, iPoint, iSurfacePoints;
	
#ifndef NO_MPI
	int rank = MPI::COMM_WORLD.Get_rank();
#else
	int rank = MASTER_NODE;
#endif
	
	/*--- Recompute the cartesians coordinates ---*/
	for (iSurfacePoints = 0; iSurfacePoints < chunk->GetnSurfacePoint(); iSurfacePoints++) {
		
		/*--- Get the marker of the surface point ---*/
		iMarker = chunk->Get_MarkerIndex(iSurfacePoints);
		
		if (config->GetMarker_All_Moving(iMarker) == YES) {
			
			/*--- Get the vertex of the surface point ---*/
			iVertex = chunk->Get_VertexIndex(iSurfacePoints);
			iPoint = chunk->Get_PointIndex(iSurfacePoints);

			/*--- Set to zero the variation of the coordinates ---*/
			for (iDim = 0; iDim < 3; iDim++) VarCoord[iDim] = 0.0;
			geometry->vertex[iMarker][iVertex]->SetVarCoord(VarCoord);

			/*--- Get the parametric coordinate of the surface point ---*/
			par_coord = chunk->Get_ParametricCoord(iSurfacePoints);
			
			/*--- Compute the new cartesian coordinate, and set the value in 
			 the chunk structure ---*/
			car_coord_new = chunk->EvalCartesianCoord(par_coord);
			chunk->Set_CartesianCoord(car_coord_new, iSurfacePoints);
			
			/*--- Get the original cartesian coordinates of the surface point ---*/
			car_coord_old = geometry->node[iPoint]->GetCoord();

			/*--- Set the value of the variation of the coordinates ---*/
			for (iDim = 0; iDim < 3; iDim++) {
				VarCoord[iDim] = car_coord_new[iDim] - car_coord_old[iDim];
				if (fabs(VarCoord[iDim]) < EPS) VarCoord[iDim] = 0.0;
			}
			
			diff = sqrt((car_coord_new[0]-car_coord_old[0])*(car_coord_new[0]-car_coord_old[0]) +
									(car_coord_new[1]-car_coord_old[1])*(car_coord_new[1]-car_coord_old[1]) +
									(car_coord_new[2]-car_coord_old[2])*(car_coord_new[2]-car_coord_old[2]));
			
			my_max_diff = max(my_max_diff, diff);
			
			/*--- Set the variation of the coordinates ---*/
			geometry->vertex[iMarker][iVertex]->SetVarCoord(VarCoord);
			
		}
	}
		
#ifndef NO_MPI
	MPI::COMM_WORLD.Allreduce(&my_max_diff, &max_diff, 1, MPI::DOUBLE, MPI::MAX); 	
#else
	max_diff = my_max_diff;
#endif
	
	if (rank == MASTER_NODE) 
		cout << "Update cartesian coord        | FFD box: " << chunk->GetTag() << ". Max diff: " << max_diff <<"."<< endl;
	
}

void CSurfaceMovement::SetFFDCPChange(CGeometry *geometry, CConfig *config, CFreeFormChunk *chunk, unsigned short iChunk, 
																			unsigned short iDV, bool ResetDef) {
	
	double movement[3], Ampl_old, Ampl_new, Ampl;	
	unsigned short design_chunk, index[3];
		
	design_chunk = int(config->GetParamDV(iDV,0));
	
	if (design_chunk == iChunk) {
		
		Ampl_old = config->GetDV_Value_Old(iDV);
		Ampl_new = config->GetDV_Value_New(iDV);
		Ampl = Ampl_new-Ampl_old;	
		
		index[0] = int(config->GetParamDV(iDV,1));
		index[1] = int(config->GetParamDV(iDV,2)); 
		index[2] = int(config->GetParamDV(iDV,3));
		
		movement[0] = config->GetParamDV(iDV,4)*Ampl; 
		movement[1] = config->GetParamDV(iDV,5)*Ampl; 
		movement[2] = config->GetParamDV(iDV,6)*Ampl;
		
		if (ResetDef == true) chunk->SetOriginalControlPoints();
		chunk->SetControlPoints(index, movement);
		
	}
		
}

void CSurfaceMovement::SetFFDCamber(CGeometry *geometry, CConfig *config, CFreeFormChunk *chunk, unsigned short iChunk, 
																		unsigned short iDV, bool ResetDef) {
	double Ampl_old, Ampl_new, Ampl, movement[3];	
	unsigned short design_chunk, index[3], kIndex;
	
	design_chunk = int(config->GetParamDV(iDV,0));
	
	if (design_chunk == iChunk) {
		
		/*--- Compute the variation of the design variable ---*/
		for (kIndex = 0; kIndex < 2; kIndex++) {
						
			Ampl_old = config->GetDV_Value_Old(iDV);
			Ampl_new = config->GetDV_Value_New(iDV);
			Ampl = Ampl_new-Ampl_old;	
			
			design_chunk = int(config->GetParamDV(iDV,0));
			if (design_chunk > nChunk) { cout <<"The chunk ID is bigger than the number of chunks!!"<< endl; exit(1); }
			
			index[0] = int(config->GetParamDV(iDV,1));
			index[1] = int(config->GetParamDV(iDV,2)); 
			index[2] = kIndex;
			
			movement[0] = 0.0; movement[1] = 0.0; 
			if (kIndex == 0) movement[2] = Ampl;
			else movement[2] = Ampl;
			
			if (ResetDef == true) chunk->SetOriginalControlPoints();
			chunk->SetControlPoints(index, movement);
		}
		
	}
	
}

void CSurfaceMovement::SetFFDThickness(CGeometry *geometry, CConfig *config, CFreeFormChunk *chunk, unsigned short iChunk, 
																			 unsigned short iDV, bool ResetDef) {
	double Ampl_old, Ampl_new, Ampl, movement[3];	
	unsigned short design_chunk, index[3], kIndex;
		
	design_chunk = int(config->GetParamDV(iDV,0));
	
	if (design_chunk == iChunk) {
		
		/*--- Compute the variation of the design variable ---*/
		for (kIndex = 0; kIndex < 2; kIndex++) {
			
			Ampl_old = config->GetDV_Value_Old(iDV);
			Ampl_new = config->GetDV_Value_New(iDV);
			Ampl = Ampl_new-Ampl_old;	
			
			design_chunk = int(config->GetParamDV(iDV,0));
			
			index[0] = int(config->GetParamDV(iDV,1));
			index[1] = int(config->GetParamDV(iDV,2)); 
			index[2] = kIndex;
			
			movement[0] = 0.0; movement[1] = 0.0; 
			if (kIndex == 0) movement[2] = -Ampl;
			else movement[2] = Ampl;
			
			if (ResetDef == true) chunk->SetOriginalControlPoints();
			chunk->SetControlPoints(index, movement);
		}
		
	}
	
}

void CSurfaceMovement::SetFFDVolume(CGeometry *geometry, CConfig *config, CFreeFormChunk *chunk, unsigned short iChunk, 
																			 unsigned short iDV, bool ResetDef) {
	double Ampl_old, Ampl_new, Ampl, movement[3]; 
	unsigned short design_chunk, index[3];
			
	design_chunk = int(config->GetParamDV(iDV,0));
	
	if (design_chunk == iChunk) {
		
		/*--- Compute the variation of the design variable ---*/
		Ampl_old = config->GetDV_Value_Old(iDV);
		Ampl_new = config->GetDV_Value_New(iDV);
		Ampl = Ampl_new-Ampl_old;	
				
		index[0] = int(config->GetParamDV(iDV,1));
		index[1] = int(config->GetParamDV(iDV,2)); 
		index[2] = 0;
		
		movement[0] = 0.0; movement[1] = 0.0; 
		movement[2] = Ampl;
		
		if (ResetDef == true) chunk->SetOriginalControlPoints();
		chunk->SetControlPoints(index, movement);
		
	}
	
}


void CSurfaceMovement::SetFFDDihedralAngle(CGeometry *geometry, CConfig *config, CFreeFormChunk *chunk, unsigned short iChunk, 
																					 unsigned short iDV, bool ResetDef) {
	unsigned short iOrder, jOrder, kOrder, design_chunk, index[3];
	double movement[3];
			
	design_chunk = int(config->GetParamDV(iDV,0));
	
	if (design_chunk == iChunk) {
		
		/*--- The angle of rotation. ---*/
		double theta_old = config->GetDV_Value_Old(iDV)*PI_NUMBER/180.0;
		double theta_new = config->GetDV_Value_New(iDV)*PI_NUMBER/180.0;
		double theta = theta_new-theta_old;
		
		/*--- Change the value of the control point if move is true ---*/
		for (iOrder = 0; iOrder < chunk->GetlOrder(); iOrder++)
			for (jOrder = 0; jOrder < chunk->GetmOrder(); jOrder++)
				for (kOrder = 0; kOrder < chunk->GetnOrder(); kOrder++) {
					index[0] = iOrder; index[1] = jOrder; index[2] = kOrder;
					double *coord = chunk->GetCoordControlPoints(iOrder, jOrder, kOrder);
					movement[0] = 0.0; movement[1] = 0.0; movement[2] = coord[1]*tan(theta);
					
					if (ResetDef == true) chunk->SetOriginalControlPoints();
					chunk->SetControlPoints(index, movement);
				}
		
	}

}

void CSurfaceMovement::SetFFDTwistAngle(CGeometry *geometry, CConfig *config, CFreeFormChunk *chunk, unsigned short iChunk, 
																				unsigned short iDV, bool ResetDef) {
	unsigned short iOrder, jOrder, kOrder;
	double  x, y, z, movement[3];
	unsigned short index[3], design_chunk;
	
	design_chunk = int(config->GetParamDV(iDV,0));
	
	if (design_chunk == iChunk) {
		
		/*--- xyz-coordinates of a point on the line of rotation. ---*/
		double a = config->GetParamDV(iDV,1);
		double b = config->GetParamDV(iDV,2);
		double c = config->GetParamDV(iDV,3);
		
    /*--- xyz-coordinate of the line's direction vector. ---*/
		double u = config->GetParamDV(iDV,4)-config->GetParamDV(iDV,1);
		double v = config->GetParamDV(iDV,5)-config->GetParamDV(iDV,2);
		double w = config->GetParamDV(iDV,6)-config->GetParamDV(iDV,3);
		
		/*--- The angle of rotation. ---*/
		double theta_old = config->GetDV_Value_Old(iDV)*PI_NUMBER/180.0;
		double theta_new = config->GetDV_Value_New(iDV)*PI_NUMBER/180.0;
		double theta = theta_new-theta_old;
		
		/*--- An intermediate value used in computations. ---*/
		double u2=u*u; double v2=v*v; double w2=w*w;     
		double l2 = u2 + v2 + w2; double l = sqrt(l2);
		double cosT; double sinT;  
		
		/*--- Change the value of the control point if move is true ---*/
		for (iOrder = 0; iOrder < chunk->GetlOrder(); iOrder++)
			for (jOrder = 0; jOrder < chunk->GetmOrder(); jOrder++)
				for (kOrder = 0; kOrder < chunk->GetnOrder(); kOrder++) {
					index[0] = iOrder; index[1] = jOrder; index[2] = kOrder;
					double *coord = chunk->GetCoordControlPoints(iOrder, jOrder, kOrder);
					x = coord[0]; y = coord[1]; z = coord[2];
					
					double factor = 0.0; 
					if ( z < config->GetParamDV(iDV,3) )
						factor = 0.0;
					if (( z >= config->GetParamDV(iDV,3)) && ( z <= config->GetParamDV(iDV,6)) )
						factor = (z-config->GetParamDV(iDV,3)) / (config->GetParamDV(iDV,6)-config->GetParamDV(iDV,3));
					if ( z > config->GetParamDV(iDV,6) )
						factor = 1.0;
					
					cosT = cos(theta*factor); 
					sinT = sin(theta*factor);  
					
					movement[0] = a*(v2 + w2) + u*(-b*v - c*w + u*x + v*y + w*z)
					+ (-a*(v2 + w2) + u*(b*v + c*w - v*y - w*z) + (v2 + w2)*x)*cosT
					+ l*(-c*v + b*w - w*y + v*z)*sinT;
					movement[0] = movement[0]/l2 - x;
					
					movement[1] = b*(u2 + w2) + v*(-a*u - c*w + u*x + v*y + w*z) 
					+ (-b*(u2 + w2) + v*(a*u + c*w - u*x - w*z) + (u2 + w2)*y)*cosT
					+ l*(c*u - a*w + w*x - u*z)*sinT;
					movement[1] = movement[1]/l2 - y;
					
					movement[2] = c*(u2 + v2) + w*(-a*u - b*v + u*x + v*y + w*z) 
					+ (-c*(u2 + v2) + w*(a*u + b*v - u*x - v*y) + (u2 + v2)*z)*cosT
					+ l*(-b*u + a*v - v*x + u*y)*sinT;
					movement[2] = movement[2]/l2 - z;
					
					if (ResetDef == true) chunk->SetOriginalControlPoints();
					chunk->SetControlPoints(index, movement);		
				}
		
	}
	
}


void CSurfaceMovement::SetFFDRotation(CGeometry *geometry, CConfig *config, CFreeFormChunk *chunk, unsigned short iChunk, 
																			unsigned short iDV, bool ResetDef) {
	unsigned short iOrder, jOrder, kOrder;
	double  movement[3], x, y, z;
	unsigned short index[3], design_chunk;
		
	design_chunk = int(config->GetParamDV(iDV,0));
	
	if (design_chunk == iChunk) {
		
		/*--- xyz-coordinates of a point on the line of rotation. ---*/
		double a = config->GetParamDV(0,1);
		double b = config->GetParamDV(0,2);
		double c = config->GetParamDV(0,3);
		
		/*--- xyz-coordinate of the line's direction vector. ---*/
		double u = config->GetParamDV(0,4)-config->GetParamDV(0,1);
		double v = config->GetParamDV(0,5)-config->GetParamDV(0,2);
		double w = config->GetParamDV(0,6)-config->GetParamDV(0,3);
		
		/*--- The angle of rotation. ---*/
		double theta_old = config->GetDV_Value_Old(0)*PI_NUMBER/180.0;
		double theta_new = config->GetDV_Value_New(0)*PI_NUMBER/180.0;
		double theta = theta_new-theta_old;
		
		/*--- An intermediate value used in computations. ---*/
		double u2=u*u; double v2=v*v; double w2=w*w;     
		double cosT = cos(theta); double sinT = sin(theta);  
		double l2 = u2 + v2 + w2; double l = sqrt(l2);
		
		/*--- Change the value of the control point if move is true ---*/
		for (iOrder = 0; iOrder < chunk->GetlOrder(); iOrder++)
			for (jOrder = 0; jOrder < chunk->GetmOrder(); jOrder++)
				for (kOrder = 0; kOrder < chunk->GetnOrder(); kOrder++) {
					index[0] = iOrder; index[1] = jOrder; index[2] = kOrder;
					double *coord = chunk->GetCoordControlPoints(iOrder, jOrder, kOrder);
					x = coord[0]; y = coord[1]; z = coord[2];
					movement[0] = a*(v2 + w2) + u*(-b*v - c*w + u*x + v*y + w*z)
					+ (-a*(v2 + w2) + u*(b*v + c*w - v*y - w*z) + (v2 + w2)*x)*cosT
					+ l*(-c*v + b*w - w*y + v*z)*sinT;
					movement[0] = movement[0]/l2 - x;
					
					movement[1] = b*(u2 + w2) + v*(-a*u - c*w + u*x + v*y + w*z) 
					+ (-b*(u2 + w2) + v*(a*u + c*w - u*x - w*z) + (u2 + w2)*y)*cosT
					+ l*(c*u - a*w + w*x - u*z)*sinT;
					movement[1] = movement[1]/l2 - y;
					
					movement[2] = c*(u2 + v2) + w*(-a*u - b*v + u*x + v*y + w*z) 
					+ (-c*(u2 + v2) + w*(a*u + b*v - u*x - v*y) + (u2 + v2)*z)*cosT
					+ l*(-b*u + a*v - v*x + u*y)*sinT;
					movement[2] = movement[2]/l2 - z;
					
					if (ResetDef == true) chunk->SetOriginalControlPoints();
					chunk->SetControlPoints(index, movement);		
				}
		
	}
	
}

void CSurfaceMovement::SetHicksHenne(CGeometry *boundary, CConfig *config, unsigned short iDV, bool ResetDef) {
	unsigned long iVertex, Point;
	unsigned short iMarker;
	double VarCoord[3], *Coord, *Normal, ek, fk, BumpSize = 1.0, BumpLoc = 1.5, xCoord;
	bool upper = true, double_surface = false;
		
	/*--- Reset airfoil deformation if first deformation ---*/
	if ((iDV == 0) || (ResetDef == true)) {
		for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
			for (iVertex = 0; iVertex < boundary->nVertex[iMarker]; iVertex++) {
				VarCoord[0] = 0.0; VarCoord[1] = 0.0; VarCoord[2] = 0.0;
				boundary->vertex[iMarker][iVertex]->SetVarCoord(VarCoord);
			}
	}
		
	/*--- Perform multiple airfoil deformation ---*/
	double Ampl_old = config->GetDV_Value_Old(iDV);
	double Ampl_new = config->GetDV_Value_New(iDV);
	double Ampl = Ampl_new-Ampl_old;
	double xk = config->GetParamDV(iDV,1);
	const double t2 = 3.0;	
	
	if (config->GetParamDV(iDV,0) == NO)  { upper = false; double_surface = true; }
	if (config->GetParamDV(iDV,0) == YES) { upper = true; double_surface = true; }
	
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
		for (iVertex = 0; iVertex < boundary->nVertex[iMarker]; iVertex++) {
			VarCoord[0] = 0.0; VarCoord[1] = 0.0; VarCoord[2] = 0.0;
				
			if (config->GetMarker_All_Moving(iMarker) == YES) {
	
				Point = boundary->vertex[iMarker][iVertex]->GetNode();
				Coord = boundary->vertex[iMarker][iVertex]->GetCoord();
				Normal = boundary->vertex[iMarker][iVertex]->GetNormal();

				/*--- Bump computation ---*/
				if (double_surface) {
					ek = log10(0.5)/log10(xk);
					fk = pow( sin( PI_NUMBER * pow(Coord[0],ek) ) , t2);
					/*--- Upper and lower surface ---*/
					if (( upper) && (Normal[1] > 0)) { VarCoord[1] =  Ampl*fk; }
					if ((!upper) && (Normal[1] < 0)) { VarCoord[1] = -Ampl*fk; }
				}
				else {
					xCoord = Coord[0] - BumpLoc;
					ek = log10(0.5)/log10(xk/BumpSize);
					fk = pow( sin( PI_NUMBER * pow(xCoord/BumpSize,ek)),t2);
					/*--- Only one surface ---*/
					if ((xCoord <= 0.0) || (xCoord >= BumpSize)) VarCoord[1] =  0.0;
					else VarCoord[1] =  Ampl*fk;
				}
			}
			boundary->vertex[iMarker][iVertex]->AddVarCoord(VarCoord);
		}
	
}

void CSurfaceMovement::SetDisplacement(CGeometry *boundary, CConfig *config, unsigned short iDV, bool ResetDef) {
	unsigned long iVertex;
	unsigned short iMarker;
	double VarCoord[3];
	double Ampl_old = config->GetDV_Value_Old(0);
	double Ampl_new = config->GetDV_Value_New(0);
	double Ampl = Ampl_new-Ampl_old;
	
	if (config->GetnDV() != 1) { cout << "This kind of design variable is not prepared for multiple desformations."; cin.get();	}
	
	double xDispl = config->GetParamDV(iDV,0);
	double yDispl = config->GetParamDV(iDV,1);
	double zDispl = 0;
	if (boundary->GetnDim() == 3) zDispl = config->GetParamDV(iDV,2);
	
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
		for (iVertex = 0; iVertex < boundary->nVertex[iMarker]; iVertex++) {
			VarCoord[0] = 0.0; VarCoord[1] = 0.0; VarCoord[2] = 0.0;
			if (config->GetMarker_All_Moving(iMarker) == YES) {
				VarCoord[0] = Ampl*xDispl;
				VarCoord[1] = Ampl*yDispl;
				if (boundary->GetnDim() == 3) VarCoord[2] = Ampl*zDispl;
			}
			boundary->vertex[iMarker][iVertex]->SetVarCoord(VarCoord);
		}
}

void CSurfaceMovement::SetRotation(CGeometry *boundary, CConfig *config, unsigned short iDV, bool ResetDef) {
	unsigned long iVertex;
	unsigned short iMarker;
	double VarCoord[3], *Coord;
	double  movement[3], x, y, z;
	
	if (config->GetnDV() != 1) { cout << "This kind of design variable is not prepared for multiple desformations."; cin.get();	}
	
	/*--- xyz-coordinates of a point on the line of rotation. */
	double a = config->GetParamDV(iDV,0);
	double b = config->GetParamDV(iDV,1);
	double c = 0.0;
	if (boundary->GetnDim() == 3) c = config->GetParamDV(0,2);
	
	/*--- xyz-coordinate of the line's direction vector. ---*/
	double u = config->GetParamDV(iDV,3)-config->GetParamDV(iDV,0);
	double v = config->GetParamDV(iDV,4)-config->GetParamDV(iDV,1);
	double w = 1.0;
	if (boundary->GetnDim() == 3) w = config->GetParamDV(iDV,5)-config->GetParamDV(iDV,2);
	
	/*--- The angle of rotation. ---*/
	double theta_old = config->GetDV_Value_Old(iDV)*PI_NUMBER/180.0;
	double theta_new = config->GetDV_Value_New(iDV)*PI_NUMBER/180.0;
	double theta = theta_new-theta_old;
	
	/*--- An intermediate value used in computations. ---*/
	double u2=u*u; double v2=v*v; double w2=w*w;     
	double cosT = cos(theta); double sinT = sin(theta);  
	double l2 = u2 + v2 + w2; double l = sqrt(l2);
	
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
		for (iVertex = 0; iVertex < boundary->nVertex[iMarker]; iVertex++) {
			VarCoord[0] = 0.0; VarCoord[1] = 0.0; VarCoord[2] = 0.0;
			if (config->GetMarker_All_Moving(iMarker) == YES) {
				Coord = boundary->vertex[iMarker][iVertex]->GetCoord();						
				x = Coord[0]; y = Coord[1]; z = Coord[2];
				
				movement[0] = a*(v2 + w2) + u*(-b*v - c*w + u*x + v*y + w*z)
				+ (-a*(v2 + w2) + u*(b*v + c*w - v*y - w*z) + (v2 + w2)*x)*cosT
				+ l*(-c*v + b*w - w*y + v*z)*sinT;
				movement[0] = movement[0]/l2 - x;
				
				movement[1] = b*(u2 + w2) + v*(-a*u - c*w + u*x + v*y + w*z) 
				+ (-b*(u2 + w2) + v*(a*u + c*w - u*x - w*z) + (u2 + w2)*y)*cosT
				+ l*(c*u - a*w + w*x - u*z)*sinT;
				movement[1] = movement[1]/l2 - y;
				
				movement[2] = c*(u2 + v2) + w*(-a*u - b*v + u*x + v*y + w*z) 
				+ (-c*(u2 + v2) + w*(a*u + b*v - u*x - v*y) + (u2 + v2)*z)*cosT
				+ l*(-b*u + a*v - v*x + u*y)*sinT;
				if (boundary->GetnDim() == 3) movement[2] = movement[2]/l2 - z;
				else movement[2] = 0.0;
				
				VarCoord[0] = movement[0];
				VarCoord[1] = movement[1];
				if (boundary->GetnDim() == 3) VarCoord[2] = movement[2];
				
			}
			boundary->vertex[iMarker][iVertex]->SetVarCoord(VarCoord);
		}	
}

void CSurfaceMovement::SetBoundary_Flutter2D(CGeometry *geometry, CConfig *config, 
                                             unsigned long iter) {
	
	double VarCoord[3], omega, w_red, deltaT, ampl, v_inf, *vel;
  double alpha, alpha_new, alpha_old, dx, dy;
  double time_new, time_old;
  double DEG2RAD = PI_NUMBER/180.0;
  unsigned short iDim, iMarker;
  unsigned short nDim = geometry->GetnDim();
  unsigned long iPoint, iVertex;
  bool adjoint = ((config->GetKind_Solver() == ADJ_EULER) || 
                  (config->GetKind_Solver() == ADJ_NAVIER_STOKES) ||
                  (config->GetKind_Solver() == ADJ_RANS));
#ifndef NO_MPI
	int rank = MPI::COMM_WORLD.Get_rank();
#else
	int rank = MASTER_NODE;
#endif
	
  /*--- Retrieve values from the config file ---*/
  deltaT    = config->GetDelta_UnstTimeND();
  vel       = config->GetVelocity_FreeStreamND();
  w_red     = config->GetReduced_Frequency();
  ampl      = config->GetPitching_Amplitude();
  
  /*--- Compute delta time based on physical time step ---*/
  if (adjoint) {
    /*--- For the unsteady adjoint, we integrate backwards through
     physical time, so perform mesh motion in reverse. ---*/ 
    unsigned long nFlowIter  = config->GetnExtIter();
    unsigned long directIter = nFlowIter - iter - 1;
    time_new = static_cast<double>(directIter)*deltaT;
    time_old = time_new;
    if (iter != 0) time_old = (static_cast<double>(directIter)+1.0)*deltaT;
  } else {
    /*--- Forward time for the direct problem ---*/
    time_new = static_cast<double>(iter)*deltaT;
    time_old = time_new;
    if (iter != 0) time_old = (static_cast<double>(iter)-1.0)*deltaT;
  }
	
  /*--- For now, hard code the origin and chord length. These can be
   inputs in the config file in the future. ---*/
  double x_origin = 0.248, y_origin = 0.0;
	double chord = 1.0;
  
  /*--- Compute the freestream velocity for use with the reduced frequency --*/
  v_inf = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    v_inf += vel[iDim]*vel[iDim];
  v_inf = sqrt(v_inf);
  
  /*--- Update the pitching angle at this time step. Flip sign for
   nose-up positive convention. ---*/
  omega     = 2.0*w_red*v_inf/chord;
  alpha_new = ampl*sin(omega*time_new);
  alpha_old = ampl*sin(omega*time_old);
  alpha     = -(1E-12 + (alpha_new - alpha_old))*DEG2RAD;
	
	if (rank == MASTER_NODE)
		cout << "New pitching angle (alpha): " << alpha_new << " degrees." << endl;
  
	/*--- Store movement and velocity of each node on the pitching surface ---*/
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_Moving(iMarker) == YES) {
      for(iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        dx = geometry->node[iPoint]->GetCoord(0) - x_origin;
        dy = geometry->node[iPoint]->GetCoord(1) - y_origin;
        VarCoord[0] = dx*cos(alpha) - dy*sin(alpha) - dx;
        VarCoord[1] = dx*sin(alpha) + dy*cos(alpha) - dy;
        VarCoord[2] = 0.0;
        /*--- Set position and velocity for this node ---*/
        geometry->vertex[iMarker][iVertex]->SetVarCoord(VarCoord);
      }
		}	
	}
  
}

void CSurfaceMovement::SetBoundary_Flutter3D(CGeometry *geometry, CConfig *config, 
                                             CFreeFormChunk **chunk, unsigned long iter) {
	
	double omega, w_red, deltaT, ampl, v_inf, *vel;
  double alpha, alpha_new, alpha_old;
  double time_new, time_old;
  unsigned short iDim;
  unsigned short nDim = geometry->GetnDim();
  bool adjoint = ((config->GetKind_Solver() == ADJ_EULER) || 
                  (config->GetKind_Solver() == ADJ_NAVIER_STOKES) ||
                  (config->GetKind_Solver() == ADJ_RANS));
#ifndef NO_MPI
	int rank = MPI::COMM_WORLD.Get_rank();
#else
	int rank = MASTER_NODE;
#endif
	
  /*--- Retrieve values from the config file ---*/
  deltaT = config->GetDelta_UnstTimeND();
  vel    = config->GetVelocity_FreeStreamND();
  w_red  = config->GetReduced_Frequency();
  ampl   = config->GetPitching_Amplitude();
  
  /*--- Compute delta time based on physical time step ---*/
  if (adjoint) {
    /*--- For the unsteady adjoint, we integrate backwards through
     physical time, so perform mesh motion in reverse. ---*/ 
    unsigned long nFlowIter  = config->GetnExtIter();
    unsigned long directIter = nFlowIter - iter - 1;
    time_new = static_cast<double>(directIter)*deltaT;
    time_old = time_new;
    if (iter != 0) time_old = (static_cast<double>(directIter)+1.0)*deltaT;
  } else {
    /*--- Forward time for the direct problem ---*/
    time_new = static_cast<double>(iter)*deltaT;
    time_old = time_new;
    if (iter != 0) time_old = (static_cast<double>(iter)-1.0)*deltaT;
  }
  
  /*--- Compute the freestream velocity for use with the reduced frequency --*/
  v_inf = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    v_inf += vel[iDim]*vel[iDim];
  v_inf = sqrt(v_inf);
  
  /*--- For now, hard code the origin and chord length. These can be
   inputs in the config file in the future. ---*/
	double chord = 1.0;
	
  /*--- Update the pitching angle at this time step. Flip sign for
   nose-up positive convention. ---*/
  omega     = 2.0*w_red*v_inf/chord;
  alpha_new = ampl*sin(omega*time_new);
  alpha_old = ampl*sin(omega*time_old);
  alpha     = (1E-10 + (alpha_new - alpha_old))*(-PI_NUMBER/180.0);
	
	if (rank == MASTER_NODE)
		cout << "New dihedral angle (alpha): " << alpha_new << " degrees." << endl;
	
	unsigned short iOrder, jOrder, kOrder;
	short iChunk;
	double movement[3];
	bool *move = new bool [nChunk];
	unsigned short *index = new unsigned short[3];
	
	move[0] = true; move[1] = true; move[2] = true;	
  
	/*--- Change the value of the control point if move is true ---*/
	for (iChunk = 0; iChunk < nChunk; iChunk++)
		if (move[iChunk])
			for (iOrder = 0; iOrder < chunk[iChunk]->GetlOrder(); iOrder++)
				for (jOrder = 0; jOrder < chunk[iChunk]->GetmOrder(); jOrder++)
					for (kOrder = 0; kOrder < chunk[iChunk]->GetnOrder(); kOrder++) {
						index[0] = iOrder; index[1] = jOrder; index[2] = kOrder;
						double *coord = chunk[iChunk]->GetCoordControlPoints(iOrder, jOrder, kOrder);
						movement[0] = 0.0; movement[1] = 0.0; movement[2] = coord[1]*tan(alpha);
						chunk[iChunk]->SetControlPoints(index, movement);
					}
	
	/*--- Recompute cartesian coordinates using the new control points position ---*/
	for (iChunk = 0; iChunk < nChunk; iChunk++)
		SetCartesianCoord(geometry, config, chunk[iChunk], iChunk);
	
}

void CSurfaceMovement::SetExternal_Deformation(CGeometry *geometry, CConfig *config, unsigned short iZone, unsigned long iter) {
  
  int rank = MASTER_NODE;
#ifndef NO_MPI
	rank = MPI::COMM_WORLD.Get_rank();
#endif
  
  /*--- Local variables ---*/
  
	unsigned short iDim, nDim; 
	unsigned long iPoint, flowIter = 0;
	double VarCoord[3], *Coord_Old = NULL, *Coord_New = NULL, Center[3], Lref;
  double NewCoord[3], rotMatrix[3][3] = {{0.0,0.0,0.0}, {0.0,0.0,0.0}, {0.0,0.0,0.0}};
  double r[3], rotCoord[3];
  unsigned long iVertex;
  unsigned short iMarker;
  char buffer[50];
  string motion_filename, UnstExt, text_line;
  ifstream motion_file;
  bool adjoint = ((config->GetKind_Solver() == ADJ_EULER) || 
                  (config->GetKind_Solver() == ADJ_NAVIER_STOKES) ||
                  (config->GetKind_Solver() == ADJ_RANS));
  
	/*--- Load stuff from config ---*/
  
	nDim = geometry->GetnDim();
  motion_filename = config->GetMotion_FileName();
  
  /*--- Center of rotation & angular velocity vector from config ---*/
	Lref   = config->GetLength_Ref();
  Center[0] = config->GetMotion_Origin_X(iZone);
  Center[1] = config->GetMotion_Origin_Y(iZone);
  Center[2] = config->GetMotion_Origin_Z(iZone);
  
  /*--- Set the extension for the correct unsteady mesh motion file ---*/  
  
  if (adjoint) {
    /*--- For the unsteady adjoint, we integrate backwards through
     physical time, so perform mesh motion in reverse. ---*/ 
    unsigned long nFlowIter = config->GetnExtIter() - 1;
    flowIter  = nFlowIter - iter;
    motion_filename.erase (motion_filename.end()-4, motion_filename.end());
    if ((int(flowIter) >= 0) && (int(flowIter) < 10)) sprintf (buffer, "_0000%d.dat", int(flowIter));
    if ((int(flowIter) >= 10) && (int(flowIter) < 100)) sprintf (buffer, "_000%d.dat", int(flowIter));
    if ((int(flowIter) >= 100) && (int(flowIter) < 1000)) sprintf (buffer, "_00%d.dat", int(flowIter));
    if ((int(flowIter) >= 1000) && (int(flowIter) < 10000)) sprintf (buffer, "_0%d.dat", int(flowIter));
    if  (int(flowIter) >= 10000) sprintf (buffer, "_%d.dat", int(flowIter));
    UnstExt = string(buffer);
    motion_filename.append(UnstExt);
  } else {
    /*--- Forward time for the direct problem ---*/
    flowIter = iter;
    motion_filename.erase (motion_filename.end()-4, motion_filename.end());
    if ((int(flowIter) >= 0) && (int(flowIter) < 10)) sprintf (buffer, "_0000%d.dat", int(flowIter));
    if ((int(flowIter) >= 10) && (int(flowIter) < 100)) sprintf (buffer, "_000%d.dat", int(flowIter));
    if ((int(flowIter) >= 100) && (int(flowIter) < 1000)) sprintf (buffer, "_00%d.dat", int(flowIter));
    if ((int(flowIter) >= 1000) && (int(flowIter) < 10000)) sprintf (buffer, "_0%d.dat", int(flowIter));
    if  (int(flowIter) >= 10000) sprintf (buffer, "_%d.dat", int(flowIter));
    UnstExt = string(buffer);
    motion_filename.append(UnstExt);
  }
  
  /*--- Open the motion file ---*/
  
  if (rank == MASTER_NODE)
    cout << "Reading in the arbitrary mesh motion from direct iteration " << flowIter << "." << endl;;
  motion_file.open(motion_filename.data(), ios::in);
  /*--- Throw error if there is no file ---*/
  if (motion_file.fail()) {
    cout << "There is no mesh motion file!" << endl;
    cout << "Press any key to exit..." << endl;
    cin.get(); exit(1);
  }
  
  /*--- Read in and store the new mesh node locations ---*/ 
  
  while (getline(motion_file,text_line)) {
    istringstream point_line(text_line);
    if (nDim == 2) point_line >> iPoint >> NewCoord[0] >> NewCoord[1];
    if (nDim == 3) point_line >> iPoint >> NewCoord[0] >> NewCoord[1] >> NewCoord[2];
    
#ifdef NO_MPI
    geometry->node[iPoint]->SetCoord_p1(NewCoord);
#else
    /*--- With MPI, need to match the local/global indices. ---*/
    unsigned long jPoint, GlobalIndex;
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if (config->GetMarker_All_Moving(iMarker) == YES) {
        for(iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          jPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          GlobalIndex = geometry->node[jPoint]->GetGlobalIndex();
          if (GlobalIndex == iPoint) {
            geometry->node[jPoint]->SetCoord_p1(NewCoord);
            break;
          }
        }
      }
    }
#endif
    
  }
  /*--- Close the restart file ---*/
  motion_file.close();
  
  /*--- If rotating as well, prepare the rotation matrix ---*/
  
  if (config->GetKind_GridMovement(iZone) == EXTERNAL_ROTATION) {
    
    /*--- Variables needed only for rotation ---*/
    
    double Omega[3], dt;
    double dtheta, dphi, dpsi, cosTheta, sinTheta;
    double cosPhi, sinPhi, cosPsi, sinPsi;
    
    /*--- Angular velocity vector from config ---*/
    
    dt = static_cast<double>(iter)*config->GetDelta_UnstTimeND();
    Omega[0]  = config->GetRotation_Rate_X(iZone);
    Omega[1]  = config->GetRotation_Rate_Y(iZone);
    Omega[2]  = config->GetRotation_Rate_Z(iZone);
    
    /*--- For the unsteady adjoint, use reverse time ---*/
    if (adjoint) {
      /*--- Set the first adjoint mesh position to the final direct one ---*/
      if (iter == 0) dt = ((double)config->GetnExtIter()-1) * dt;
      /*--- Reverse the rotation direction for the adjoint ---*/
      else dt = -1.0*dt;
    } else {
      /*--- No rotation at all for the first direct solution ---*/
      if (iter == 0) dt = 0;
    }
    
    /*--- Compute delta change in the angle about the x, y, & z axes. ---*/
    
    dtheta = Omega[0]*dt;   
    dphi   = Omega[1]*dt; 
    dpsi   = Omega[2]*dt;
    
    /*--- Store angles separately for clarity. Compute sines/cosines. ---*/
    
    cosTheta = cos(dtheta);  cosPhi = cos(dphi);  cosPsi = cos(dpsi);
    sinTheta = sin(dtheta);  sinPhi = sin(dphi);  sinPsi = sin(dpsi);
    
    /*--- Compute the rotation matrix. Note that the implicit
     ordering is rotation about the x-axis, y-axis, then z-axis. ---*/
    
    rotMatrix[0][0] = cosPhi*cosPsi;
    rotMatrix[1][0] = cosPhi*sinPsi;
    rotMatrix[2][0] = -sinPhi;
    
    rotMatrix[0][1] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;
    rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;
    rotMatrix[2][1] = sinTheta*cosPhi;
    
    rotMatrix[0][2] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
    rotMatrix[1][2] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
    rotMatrix[2][2] = cosTheta*cosPhi;
    
  }
  
  /*--- Loop through to find only moving surface markers ---*/
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_Moving(iMarker) == YES) {
      
      /*--- Loop over all surface points for this marker ---*/
      
      for(iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        
        /*--- Get current and new coordinates from file ---*/
        
        Coord_Old = geometry->node[iPoint]->GetCoord();
        Coord_New = geometry->node[iPoint]->GetCoord_p1();
        
        /*--- If we're also rotating, multiply each point by the
         rotation matrix. It is assumed that the coordinates in
         Coord_Old have already been rotated using SetRigid_Rotation(). ---*/
        
        if (config->GetKind_GridMovement(iZone) == EXTERNAL_ROTATION) {
          
          /*--- Calculate non-dim. position from rotation center ---*/
          
          for (iDim = 0; iDim < nDim; iDim++)
            r[iDim] = (Coord_New[iDim]-Center[iDim])/Lref;
          if (nDim == 2) r[nDim] = 0.0;
          
          /*--- Compute transformed point coordinates ---*/
          
          rotCoord[0] = rotMatrix[0][0]*r[0] 
                      + rotMatrix[0][1]*r[1] 
                      + rotMatrix[0][2]*r[2] + Center[0];
          
          rotCoord[1] = rotMatrix[1][0]*r[0] 
                      + rotMatrix[1][1]*r[1] 
                      + rotMatrix[1][2]*r[2] + Center[1];
          
          rotCoord[2] = rotMatrix[2][0]*r[0] 
                      + rotMatrix[2][1]*r[1] 
                      + rotMatrix[2][2]*r[2] + Center[2];
          
          /*--- Copy rotated coords back to original array for consistency ---*/
          for (iDim = 0; iDim < nDim; iDim++)
            Coord_New[iDim] = rotCoord[iDim];
        }
        
        /*--- Calculate delta change in the x, y, & z directions ---*/
        for (iDim = 0; iDim < nDim; iDim++)
          VarCoord[iDim] = (Coord_New[iDim]-Coord_Old[iDim])/Lref;
        if (nDim == 2) VarCoord[nDim] = 0.0;

        /*--- Set position changes to be applied by the spring analogy ---*/
        geometry->vertex[iMarker][iVertex]->SetVarCoord(VarCoord);
        
      }
    }	
  }
}

void CSurfaceMovement::SetNACA_4Digits(CGeometry *boundary, CConfig *config) {
	unsigned long iVertex, Point;
	unsigned short iMarker;
	double VarCoord[3], *Coord, *Normal, Ycurv, Yesp;

	if (config->GetnDV() != 1) { cout << "This kind of design variable is not prepared for multiple deformations."; cin.get();	}

	double Ya = config->GetParamDV(0,0) / 100.0; /*--- Maximum camber as a fraction of the chord 
					(100 m is the first of the four digits) ---*/
	double Xa = config->GetParamDV(0,1) / 10.0; /*--- Location of maximum camber as a fraction of 
					the chord (10 p is the second digit in the NACA xxxx description) ---*/
	double t = config->GetParamDV(0,2) / 100.0; /*--- Maximum thickness as a fraction of the
					  chord (so 100 t gives the last two digits in 
					  the NACA 4-digit denomination) ---*/
		
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
		for (iVertex = 0; iVertex < boundary->nVertex[iMarker]; iVertex++) {
			VarCoord[0] = 0.0; VarCoord[1] = 0.0; VarCoord[2] = 0.0;
			if (config->GetMarker_All_Moving(iMarker) == YES) {
				Point = boundary->vertex[iMarker][iVertex]->GetNode();
				Coord = boundary->vertex[iMarker][iVertex]->GetCoord();
				Normal = boundary->vertex[iMarker][iVertex]->GetNormal();
				
				if  (Coord[0] < Xa) Ycurv = (2.0*Xa*Coord[0]-pow(Coord[0],2.0))*(Ya/pow(Xa,2.0));
				else Ycurv = ((1.0-2.0*Xa)+2.0*Xa*Coord[0]-pow(Coord[0],2.0))*(Ya/pow((1.0-Xa), 2.0));
				
				Yesp = t*(1.4845*sqrt(Coord[0])-0.6300*Coord[0]-1.7580*pow(Coord[0],2.0)+
						  1.4215*pow(Coord[0],3.0)-0.518*pow(Coord[0],4.0));
				
				if (Normal[1] > 0) VarCoord[1] =  (Ycurv + Yesp) - Coord[1];
				if (Normal[1] < 0) VarCoord[1] =  (Ycurv - Yesp) - Coord[1];

			}
			boundary->vertex[iMarker][iVertex]->SetVarCoord(VarCoord);
		}
}

void CSurfaceMovement::SetParabolic(CGeometry *boundary, CConfig *config) {
	unsigned long iVertex, Point;
	unsigned short iMarker;
	double VarCoord[3], *Coord, *Normal;
	
	if (config->GetnDV() != 1) { cout << "This kind of design variable is not prepared for multiple desformations."; cin.get();	}
	
	double c = config->GetParamDV(0,0); /*--- Center of the parabola ---*/
	double t = config->GetParamDV(0,1) / 100.0; /*--- Thickness of the parabola ---*/
	
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
		for (iVertex = 0; iVertex < boundary->nVertex[iMarker]; iVertex++) {
			VarCoord[0] = 0.0; VarCoord[1] = 0.0; VarCoord[2] = 0.0;
			if (config->GetMarker_All_Moving(iMarker) == YES) {
				Point = boundary->vertex[iMarker][iVertex]->GetNode();
				Coord = boundary->vertex[iMarker][iVertex]->GetCoord();
				Normal = boundary->vertex[iMarker][iVertex]->GetNormal();
				
				if (Normal[1] > 0) {
					VarCoord[1] =  t*(Coord[0]*Coord[0]-Coord[0])/(2.0*(c*c-c)) - Coord[1];
				}
				if (Normal[1] < 0) {
					VarCoord[1] =  t*(Coord[0]-Coord[0]*Coord[0])/(2.0*(c*c-c)) - Coord[1];
				}
			}
			boundary->vertex[iMarker][iVertex]->SetVarCoord(VarCoord);
		}
}

void CSurfaceMovement::SetObstacle(CGeometry *boundary, CConfig *config) {
	unsigned long iVertex, Point;
	unsigned short iMarker;
	double VarCoord[3], *Coord, xCoord;
	
	if (config->GetnDV() != 1) { cout << "This kind of design variable is not prepared for multiple desformations."; cin.get();	}
	
	double H = config->GetParamDV(0,0); /*--- Non-dimensionalized height of the obstacle ---*/
	double L = config->GetParamDV(0,1); /*--- Non-dimensionalized length of the obstacle ---*/
	double xOffSet = 0; /*--- x offset ---*/

	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
		for (iVertex = 0; iVertex < boundary->nVertex[iMarker]; iVertex++) {
			VarCoord[0] = 0.0; VarCoord[1] = 0.0; VarCoord[2] = 0.0;
			if (config->GetMarker_All_Moving(iMarker) == YES) {
				Point = boundary->vertex[iMarker][iVertex]->GetNode();
				Coord = boundary->vertex[iMarker][iVertex]->GetCoord();
				xCoord = Coord[0]-xOffSet;
				if ((xCoord > 0) && (xCoord < L))
					VarCoord[1] = (27.0/4.0)*(H/(L*L*L))*xCoord*(xCoord-L)*(xCoord-L);
				else 
					VarCoord[1] = 0.0;
			}
			boundary->vertex[iMarker][iVertex]->SetVarCoord(VarCoord);
		}
}

void CSurfaceMovement::SetStretch(CGeometry *boundary, CConfig *config) {
	unsigned long iVertex, Point;
	unsigned short iMarker;
	double VarCoord[3], *Coord;
	
	if (config->GetnDV() != 1) { cout << "This kind of design variable is not prepared for multiple desformations."; cin.get();	}
	
	double End = config->GetParamDV(0,1);
	
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
		for (iVertex = 0; iVertex < boundary->nVertex[iMarker]; iVertex++) {
			VarCoord[0] = 0.0; VarCoord[1] = 0.0; VarCoord[2] = 0.0;
			if (config->GetMarker_All_Moving(iMarker) == YES) {
				Point = boundary->vertex[iMarker][iVertex]->GetNode();
				Coord = boundary->vertex[iMarker][iVertex]->GetCoord();
				VarCoord[0] = End - Coord[0];
			}
			boundary->vertex[iMarker][iVertex]->SetVarCoord(VarCoord);
		}
}

void CSurfaceMovement::ReadFFDInfo(CConfig *config, CGeometry *geometry, CFreeFormChunk **chunk, 
								   string val_mesh_filename) {
	string text_line, iTag;
	ifstream mesh_file;
	double coord[3];
	unsigned short degree[3], iChunk, iCornerPoints, iControlPoints, iMarker, iDegree, jDegree, kDegree, iChar;
	unsigned long iSurfacePoints, iPoint, jPoint, iVertex, nPoint, iElem = 0, nElem;
	unsigned short LevelChunk, nParentChunk, iParentChunk, nChildChunk, iChildChunk;

#ifndef NO_MPI
	int rank = MPI::COMM_WORLD.Get_rank();
#else
	int rank = MASTER_NODE;
#endif
	
	char *cstr = new char [val_mesh_filename.size()+1];
	strcpy (cstr, val_mesh_filename.c_str());
	
	mesh_file.open(cstr, ios::in);
	if (mesh_file.fail()) {
		cout << "There is no geometry file (ReadFFDInfo)!!" << endl;
		cout << "Press any key to exit..." << endl;
		cin.get();
		exit(1);
	}
			
	while (getline (mesh_file,text_line)) {
		
		/*--- Read the inner elements ---*/
		string::size_type position = text_line.find ("NELEM=",0);
		if (position != string::npos) {
			text_line.erase (0,6); nElem = atoi(text_line.c_str());
			while (iElem < nElem) {
				getline(mesh_file,text_line);
				iElem++; 
			}
		}
		
		/*--- Read the inner points ---*/
		position = text_line.find ("NPOINT=",0);
		if (position != string::npos) {
			text_line.erase (0,6); nPoint = atoi(text_line.c_str());
			while (iPoint < nPoint) {
				getline(mesh_file,text_line);
				iPoint++; 
			}
		}
		
		position = text_line.find ("NCHUNK=",0);
		if (position != string::npos) {
			text_line.erase (0,7);
			nChunk = atoi(text_line.c_str());
			if (rank == MASTER_NODE) cout << nChunk << " Free Form Deformation (FFD) chunks." << endl;
			unsigned short *nCornerPoints = new unsigned short[nChunk];
			unsigned short *nControlPoints = new unsigned short[nChunk];
			unsigned long *nSurfacePoints = new unsigned long[nChunk];
			
			getline (mesh_file,text_line);
			text_line.erase (0,7); 
			nLevel = atoi(text_line.c_str());
			if (rank == MASTER_NODE) cout << nLevel << " Free Form Deformation (FFD) nested levels." << endl;

			for (iChunk = 0 ; iChunk < nChunk; iChunk++) {
				
				/*--- Read the name of the FFD box ---*/
				getline (mesh_file,text_line);
				text_line.erase (0,10); 
				
				/*--- Remove extra data from the chunk name ---*/
				string::size_type position;
				for (iChar = 0; iChar < 20; iChar++) {
					position = text_line.find( " ", 0 );
					if(position != string::npos) text_line.erase (position,1);
					position = text_line.find( "\r", 0 );
					if(position != string::npos) text_line.erase (position,1);
					position = text_line.find( "\n", 0 );
					if(position != string::npos) text_line.erase (position,1);
				}
				
				string TagChunk = text_line.c_str();
				if (rank == MASTER_NODE) cout << "FFD box tag: " << TagChunk <<". ";

				/*--- Read the level of the FFD box ---*/
				getline (mesh_file,text_line);
				text_line.erase (0,12);
				LevelChunk = atoi(text_line.c_str());
				if (rank == MASTER_NODE) cout << "FFD box level: " << LevelChunk <<". ";
				
				/*--- Read the degree of the FFD box ---*/
				getline (mesh_file,text_line);
				text_line.erase (0,15); degree[0] = atoi(text_line.c_str());
				getline (mesh_file,text_line);
				text_line.erase (0,15); degree[1] = atoi(text_line.c_str());
				getline (mesh_file,text_line);
				text_line.erase (0,15); degree[2] = atoi(text_line.c_str());
				if (rank == MASTER_NODE) cout << "Degrees: " << degree[0] <<", " << degree[1] <<", "<< degree[2] <<". "<< endl;
				chunk[iChunk] = new CFreeFormChunk(int(degree[0]), int(degree[1]), int(degree[2]));				
				chunk[iChunk]->SetTag(TagChunk); chunk[iChunk]->SetLevel(LevelChunk);

				/*--- Read the number of parents boxes ---*/
				getline (mesh_file,text_line);
				text_line.erase (0,14);
				nParentChunk = atoi(text_line.c_str());
				if (rank == MASTER_NODE) cout << "Number of parent boxes: " << nParentChunk <<". ";
				for (iParentChunk = 0; iParentChunk < nParentChunk; iParentChunk++) {
					getline(mesh_file, text_line);
					
					/*--- Remove extra data from the chunk name ---*/
					string::size_type position;
					for (iChar = 0; iChar < 20; iChar++) {
						position = text_line.find( " ", 0 );
						if(position != string::npos) text_line.erase (position,1);
						position = text_line.find( "\r", 0 );
						if(position != string::npos) text_line.erase (position,1);
						position = text_line.find( "\n", 0 );
						if(position != string::npos) text_line.erase (position,1);
					}
					
					string ParentChunk = text_line.c_str();
					chunk[iChunk]->SetParentChunk(ParentChunk);
				}
				
				/*--- Read the number of children boxes ---*/
				getline (mesh_file,text_line);
				text_line.erase (0,15);
				nChildChunk = atoi(text_line.c_str());
				if (rank == MASTER_NODE) cout << "Number of child boxes: " << nChildChunk <<"." << endl;
				for (iChildChunk = 0; iChildChunk < nChildChunk; iChildChunk++) {
					getline(mesh_file, text_line);
					
					/*--- Remove extra data from the chunk name ---*/
					string::size_type position;
					for (iChar = 0; iChar < 20; iChar++) {
						position = text_line.find( " ", 0 );
						if(position != string::npos) text_line.erase (position,1);
						position = text_line.find( "\r", 0 );
						if(position != string::npos) text_line.erase (position,1);
						position = text_line.find( "\n", 0 );
						if(position != string::npos) text_line.erase (position,1);
					}
					
					string ChildChunk = text_line.c_str();
					chunk[iChunk]->SetChildChunk(ChildChunk);
				}
								
				/*--- Read the number of the corner points ---*/
				getline (mesh_file,text_line);
				text_line.erase (0,20); nCornerPoints[iChunk] = atoi(text_line.c_str());
				if (rank == MASTER_NODE) cout << "Corner points: " << nCornerPoints[iChunk] <<". ";
				
				/*--- Read the coordinates of the corner points ---*/
				for (iCornerPoints = 0; iCornerPoints < nCornerPoints[iChunk]; iCornerPoints++) {
					getline(mesh_file,text_line); istringstream chunk_line(text_line);
					chunk_line >> coord[0]; chunk_line >> coord[1]; chunk_line >> coord[2];
					chunk[iChunk]->SetCoordCornerPoints(coord, iCornerPoints);
				}
				
				/*--- Read the number of the control points ---*/
				getline (mesh_file,text_line);
				text_line.erase (0,21); nControlPoints[iChunk] = atoi(text_line.c_str());
				if (rank == MASTER_NODE) cout << "Control points: " << nControlPoints[iChunk] <<". ";
				
				/*--- Method to identify if there is a chunk definition ---*/
				if (nControlPoints[iChunk] != 0) ChunkDefinition = true;

				/*--- Read the coordinates of the control points ---*/
				for (iControlPoints = 0; iControlPoints < nControlPoints[iChunk]; iControlPoints++) {
					getline(mesh_file,text_line); istringstream chunk_line(text_line);
					chunk_line >> iDegree; chunk_line >> jDegree; chunk_line >> kDegree; 
					chunk_line >> coord[0]; chunk_line >> coord[1]; chunk_line >> coord[2]; 
					chunk[iChunk]->SetCoordControlPoints(coord, iDegree, jDegree, kDegree); 
				}
				
				getline (mesh_file,text_line);
				text_line.erase (0,21); nSurfacePoints[iChunk] = atoi(text_line.c_str());
				
				unsigned long my_nSurfPoints = nSurfacePoints[iChunk];
				unsigned long nSurfPoints = 0;
				
#ifndef NO_MPI
				MPI::COMM_WORLD.Allreduce(&my_nSurfPoints, &nSurfPoints, 1, MPI::UNSIGNED_LONG, MPI::SUM); 
#else
				nSurfPoints = my_nSurfPoints;
#endif
				
				if (rank == MASTER_NODE) cout << "Surface points: " << nSurfPoints <<"."<<endl;
												
				for (iSurfacePoints = 0; iSurfacePoints < nSurfacePoints[iChunk]; iSurfacePoints++) {
					getline(mesh_file,text_line); istringstream chunk_line(text_line);
					chunk_line >> iTag; chunk_line >> iPoint;
					iMarker = config->GetTag_Marker_All(iTag);
					chunk_line >> coord[0]; chunk_line >> coord[1]; chunk_line >> coord[2];
					for(iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
						jPoint =  geometry->vertex[iMarker][iVertex]->GetNode();
						if (iPoint == jPoint) {
							chunk[iChunk]->Set_MarkerIndex(iMarker);
							chunk[iChunk]->Set_VertexIndex(iVertex);
							chunk[iChunk]->Set_PointIndex(iPoint);
							chunk[iChunk]->Set_ParametricCoord(coord);
							chunk[iChunk]->Set_CartesianCoord(geometry->node[iPoint]->GetCoord());
						}
					}
				}
			}
			
			delete [] nCornerPoints;
			delete [] nControlPoints;
			delete [] nSurfacePoints;		
		}
	}
	mesh_file.close();
	if (nChunk == 0) {
		cout <<"There is no FFD box definition. Just in case, review the .su2 file" << endl;
	}

}

void CSurfaceMovement::WriteFFDInfo(CGeometry *geometry, CConfig *config, CFreeFormChunk **chunk, string val_mesh_filename) {
	ofstream mesh_file;
	unsigned short iOrder, jOrder, kOrder, iChunk, iCornerPoints, iMarker, iParentChunk, iChildChunk;
	unsigned long iVertex, iPoint, iSurfacePoints;
	char *cstr = new char [val_mesh_filename.size()+1];
	strcpy (cstr, val_mesh_filename.c_str());
	
	mesh_file.precision(15);
	mesh_file.open(cstr, ios::out | ios::app);
	
	mesh_file << "NCHUNK= " << nChunk << endl;
	mesh_file << "NLEVEL= " << nLevel << endl;
	
	for (iChunk = 0 ; iChunk < nChunk; iChunk++) {
		
		mesh_file << "CHUNK_TAG= " << chunk[iChunk]->GetTag() << endl;
		mesh_file << "CHUNK_LEVEL= " << chunk[iChunk]->GetLevel() << endl;

		mesh_file << "CHUNK_DEGREE_I= " << chunk[iChunk]->GetlOrder()-1 << endl;
		mesh_file << "CHUNK_DEGREE_J= " << chunk[iChunk]->GetmOrder()-1 << endl;
		mesh_file << "CHUNK_DEGREE_K= " << chunk[iChunk]->GetnOrder()-1 << endl;
		
		mesh_file << "CHUNK_PARENTS= " << chunk[iChunk]->GetnParentChunk() << endl;
		for (iParentChunk = 0; iParentChunk < chunk[iChunk]->GetnParentChunk(); iParentChunk++)
			mesh_file << chunk[iChunk]->GetParentChunkTag(iParentChunk) << endl;
		mesh_file << "CHUNK_CHILDREN= " << chunk[iChunk]->GetnChildChunk() << endl;
		for (iChildChunk = 0; iChildChunk < chunk[iChunk]->GetnChildChunk(); iChildChunk++)
			mesh_file << chunk[iChunk]->GetChildChunkTag(iChildChunk) << endl;
		
		mesh_file << "CHUNK_CORNER_POINTS= " << chunk[iChunk]->GetnCornerPoints() << endl;
		for (iCornerPoints = 0; iCornerPoints < chunk[iChunk]->GetnCornerPoints(); iCornerPoints++) {
			double *coord = chunk[iChunk]->GetCoordCornerPoints(iCornerPoints);
			mesh_file << coord[0] << "\t" << coord[1] << "\t" << coord[2] << endl;
		}

		/*--- No FFD definition ---*/
		if (chunk[iChunk]->GetnControlPoints() == 0) {
			mesh_file << "CHUNK_CONTROL_POINTS= 0" << endl;
			mesh_file << "CHUNK_SURFACE_POINTS= 0" << endl;				
		}
		else {
			mesh_file << "CHUNK_CONTROL_POINTS= " << chunk[iChunk]->GetnControlPoints() << endl;
			for (iOrder = 0; iOrder < chunk[iChunk]->GetlOrder(); iOrder++)
				for (jOrder = 0; jOrder < chunk[iChunk]->GetmOrder(); jOrder++)
					for (kOrder = 0; kOrder < chunk[iChunk]->GetnOrder(); kOrder++) {
						double *coord = chunk[iChunk]->GetCoordControlPoints(iOrder, jOrder, kOrder);
						mesh_file << iOrder << "\t" << jOrder << "\t" << kOrder << "\t" << coord[0] << "\t" << coord[1] << "\t" << coord[2] << endl;
					}
			
			mesh_file << "CHUNK_SURFACE_POINTS= " << chunk[iChunk]->GetnSurfacePoint() << endl;
			for (iSurfacePoints = 0; iSurfacePoints < chunk[iChunk]->GetnSurfacePoint(); iSurfacePoints++) {
				iMarker = chunk[iChunk]->Get_MarkerIndex(iSurfacePoints);
				iVertex = chunk[iChunk]->Get_VertexIndex(iSurfacePoints);
				iPoint = chunk[iChunk]->Get_PointIndex(iSurfacePoints);
				double *parcoord = chunk[iChunk]->Get_ParametricCoord(iSurfacePoints);
				mesh_file << scientific << config->GetMarker_All_Tag(iMarker) << "\t" << iPoint << "\t" << parcoord[0] << "\t" << parcoord[1] << "\t" << parcoord[2] << endl;
			}		
			
		}
		
	}
	mesh_file.close();
}

void CSurfaceMovement::WriteFFDInfo(CGeometry *geometry, CGeometry *domain, CConfig *config, CFreeFormChunk **chunk, string val_mesh_filename) {
	ofstream mesh_file;
	unsigned short iOrder, jOrder, kOrder, iChunk, iCornerPoints, iMarker, iParentChunk, iChildChunk;
	unsigned long iVertex, iPoint, iSurfacePoints;
	char *cstr = new char [val_mesh_filename.size()+1];
	strcpy (cstr, val_mesh_filename.c_str());
	
	mesh_file.precision(15);
	mesh_file.open(cstr, ios::out | ios::app);
	
	mesh_file << "NCHUNK= " << nChunk << endl;
	mesh_file << "NLEVEL= " << nLevel << endl;

	for (iChunk = 0 ; iChunk < nChunk; iChunk++) {
		
		mesh_file << "CHUNK_TAG= " << chunk[iChunk]->GetTag() << endl;
		mesh_file << "CHUNK_LEVEL= " << chunk[iChunk]->GetLevel() << endl;

		mesh_file << "CHUNK_DEGREE_I= " << chunk[iChunk]->GetlOrder()-1 << endl;
		mesh_file << "CHUNK_DEGREE_J= " << chunk[iChunk]->GetmOrder()-1 << endl;
		mesh_file << "CHUNK_DEGREE_K= " << chunk[iChunk]->GetnOrder()-1 << endl;
		
		mesh_file << "CHUNK_PARENTS= " << chunk[iChunk]->GetnParentChunk() << endl;
		for (iParentChunk = 0; iParentChunk < chunk[iChunk]->GetnParentChunk(); iParentChunk++)
			mesh_file << chunk[iChunk]->GetParentChunkTag(iParentChunk) << endl;
		
		mesh_file << "CHUNK_CHILDREN= " << chunk[iChunk]->GetnChildChunk() << endl;
		for (iChildChunk = 0; iChildChunk < chunk[iChunk]->GetnChildChunk(); iChildChunk++)
			mesh_file << chunk[iChunk]->GetChildChunkTag(iChildChunk) << endl;
		
		mesh_file << "CHUNK_CORNER_POINTS= " << chunk[iChunk]->GetnCornerPoints() << endl;
		for (iCornerPoints = 0; iCornerPoints < chunk[iChunk]->GetnCornerPoints(); iCornerPoints++) {
			double *coord = chunk[iChunk]->GetCoordCornerPoints(iCornerPoints);
			mesh_file << coord[0] << "\t" << coord[1] << "\t" << coord[2] << endl;
		}
		
		/*--- No FFD definition ---*/
		if (chunk[iChunk]->GetnControlPoints() == 0) {
			mesh_file << "CHUNK_CONTROL_POINTS= 0" << endl;
			mesh_file << "CHUNK_SURFACE_POINTS= 0" << endl;				
		}
		else {
			mesh_file << "CHUNK_CONTROL_POINTS= " << chunk[iChunk]->GetnControlPoints() << endl;
			for (iOrder = 0; iOrder < chunk[iChunk]->GetlOrder(); iOrder++)
				for (jOrder = 0; jOrder < chunk[iChunk]->GetmOrder(); jOrder++)
					for (kOrder = 0; kOrder < chunk[iChunk]->GetnOrder(); kOrder++) {
						double *coord = chunk[iChunk]->GetCoordControlPoints(iOrder, jOrder, kOrder);
						mesh_file << iOrder << "\t" << jOrder << "\t" << kOrder << "\t" << coord[0] << "\t" << coord[1] << "\t" << coord[2] << endl;
					}
			
			/*--- Compute the number of points on the new surfaces ---*/
			unsigned long nSurfacePoint = 0;
			for (iSurfacePoints = 0; iSurfacePoints < chunk[iChunk]->GetnSurfacePoint(); iSurfacePoints++) {
				iPoint = chunk[iChunk]->Get_PointIndex(iSurfacePoints);
				if (domain->GetGlobal_to_Local_Point(iPoint) != -1) nSurfacePoint++;
			}
			
			mesh_file << "CHUNK_SURFACE_POINTS= " << nSurfacePoint << endl;
			for (iSurfacePoints = 0; iSurfacePoints < chunk[iChunk]->GetnSurfacePoint(); iSurfacePoints++) {
				iMarker = chunk[iChunk]->Get_MarkerIndex(iSurfacePoints);
				iVertex = chunk[iChunk]->Get_VertexIndex(iSurfacePoints);
				iPoint = chunk[iChunk]->Get_PointIndex(iSurfacePoints);
				if (domain->GetGlobal_to_Local_Point(iPoint) != -1) {
					double *parcoord = chunk[iChunk]->Get_ParametricCoord(iSurfacePoints);
					mesh_file << scientific << config->GetMarker_All_Tag(iMarker) << "\t" << domain->GetGlobal_to_Local_Point(iPoint) << "\t" << parcoord[0] << "\t" << parcoord[1] << "\t" << parcoord[2] << endl;
				}
			}
			
		}
	}
	mesh_file.close();
}

CFreeFormChunk::CFreeFormChunk(void) : CGridMovement() { }

CFreeFormChunk::CFreeFormChunk(unsigned short val_lDegree, unsigned short val_mDegree, unsigned short val_nDegree) : CGridMovement() {
	unsigned short iCornerPoints, iOrder, jOrder, kOrder, iDim;
	
	/*--- Only for 3D problems and FFD with Hexahedron ---*/
	nDim = 3;
	nCornerPoints = 8;
	
	/*--- Allocate Corners points ---*/
	Coord_Corner_Points = new double* [nCornerPoints];
	for (iCornerPoints = 0; iCornerPoints < nCornerPoints; iCornerPoints++)
		Coord_Corner_Points[iCornerPoints] = new double [nDim];
	
	param_coord = new double[nDim]; param_coord_ = new double[nDim];
	cart_coord = new double[nDim]; cart_coord_ = new double[nDim];
	gradient = new double[nDim];

	lDegree = val_lDegree; lOrder = lDegree+1;
	mDegree = val_mDegree; mOrder = mDegree+1;
	nDegree = val_nDegree; nOrder = nDegree+1;
	nControlPoints = lOrder*mOrder*nOrder;
	
	Coord_Control_Points = new double*** [lOrder];
	ParCoord_Control_Points = new double*** [lOrder];
	Coord_Control_Points_Copy = new double*** [lOrder];
	for (iOrder = 0; iOrder < lOrder; iOrder++) {
		Coord_Control_Points[iOrder] = new double** [mOrder];
		ParCoord_Control_Points[iOrder] = new double** [mOrder];
		Coord_Control_Points_Copy[iOrder] = new double** [mOrder];
		for (jOrder = 0; jOrder < mOrder; jOrder++) {
			Coord_Control_Points[iOrder][jOrder] = new double* [nOrder];
			ParCoord_Control_Points[iOrder][jOrder] = new double* [nOrder];
			Coord_Control_Points_Copy[iOrder][jOrder] = new double* [nOrder];
			for (kOrder = 0; kOrder < nOrder; kOrder++) {
				Coord_Control_Points[iOrder][jOrder][kOrder] = new double [nDim];
				ParCoord_Control_Points[iOrder][jOrder][kOrder] = new double [nDim];
				Coord_Control_Points_Copy[iOrder][jOrder][kOrder] = new double [nDim];
			}
		}
	}
	
	/*--- Zero-initialization ---*/
	for (iOrder = 0; iOrder < lOrder; iOrder++) 
		for (jOrder = 0; jOrder < mOrder; jOrder++) 
			for (kOrder = 0; kOrder < nOrder; kOrder++)
				for (iDim = 0; iDim < nDim; iDim++)
					Coord_Control_Points[iOrder][jOrder][kOrder][iDim] = 0.0;
}

CFreeFormChunk::~CFreeFormChunk(void) {
	unsigned short iOrder, jOrder, kOrder, iCornerPoints;
	
	for (iOrder = 0; iOrder < lOrder; iOrder++) 
		for (jOrder = 0; jOrder < mOrder; jOrder++) 
			for (kOrder = 0; kOrder < nOrder; kOrder++) {
				delete [] Coord_Control_Points[iOrder][jOrder][kOrder];
				delete [] ParCoord_Control_Points[iOrder][jOrder][kOrder];
				delete [] Coord_Control_Points_Copy[iOrder][jOrder][kOrder];
			}
	delete [] Coord_Control_Points;
	delete [] ParCoord_Control_Points;
	delete [] Coord_Control_Points_Copy;

	delete [] param_coord;
	delete [] cart_coord;
	delete [] gradient;
	
	for (iCornerPoints = 0; iCornerPoints < nCornerPoints; iCornerPoints++)
		delete [] Coord_Corner_Points[iCornerPoints];
	delete [] Coord_Corner_Points;
}

void  CFreeFormChunk::SetUnitCornerPoints(void) {
	double coord [3];
	
	coord [0] = 0.0; coord [1] = 0.0; coord [2] = 0.0; this->SetCoordCornerPoints(coord,0);
	coord [0] = 1.0; coord [1] = 0.0; coord [2] = 0.0; this->SetCoordCornerPoints(coord,1);
	coord [0] = 1.0; coord [1] = 1.0; coord [2] = 0.0; this->SetCoordCornerPoints(coord,2);
	coord [0] = 0.0; coord [1] = 1.0; coord [2] = 0.0; this->SetCoordCornerPoints(coord,3);
	coord [0] = 0.0; coord [1] = 0.0; coord [2] = 1.0; this->SetCoordCornerPoints(coord,4);
	coord [0] = 1.0; coord [1] = 0.0; coord [2] = 1.0; this->SetCoordCornerPoints(coord,5);
	coord [0] = 1.0; coord [1] = 1.0; coord [2] = 1.0; this->SetCoordCornerPoints(coord,6);
	coord [0] = 0.0; coord [1] = 1.0; coord [2] = 1.0; this->SetCoordCornerPoints(coord,7);
}

void CFreeFormChunk::SetControlPoints_Parallelepiped (void) {
	unsigned short iDim, iDegree, jDegree, kDegree;
	
	/*--- Set base control points according to the notation of Vtk for hexahedrons ---*/
	for (iDim = 0; iDim < nDim; iDim++) {
		Coord_Control_Points	[0]			[0]			[0]			[iDim]	= Coord_Corner_Points[0][iDim];
		Coord_Control_Points	[lOrder-1]	[0]			[0]			[iDim]	= Coord_Corner_Points[1][iDim];
		Coord_Control_Points	[lOrder-1]	[mOrder-1]	[0]			[iDim]	= Coord_Corner_Points[2][iDim];
		Coord_Control_Points	[0]			[mOrder-1]	[0]			[iDim]	= Coord_Corner_Points[3][iDim];
		Coord_Control_Points	[0]			[0]			[nOrder-1]	[iDim]	= Coord_Corner_Points[4][iDim];
		Coord_Control_Points	[lOrder-1]	[0]			[nOrder-1]	[iDim]	= Coord_Corner_Points[5][iDim];
		Coord_Control_Points	[lOrder-1]	[mOrder-1]	[nOrder-1]	[iDim]	= Coord_Corner_Points[6][iDim];
		Coord_Control_Points	[0]			[mOrder-1]	[nOrder-1]	[iDim]	= Coord_Corner_Points[7][iDim];
	}
	
	/*--- Fill the rest of the cubic matrix of control points with uniform spacing (parallelepiped) ---*/
	for (iDegree = 0; iDegree <= lDegree; iDegree++)
		for (jDegree = 0; jDegree <= mDegree; jDegree++)
			for (kDegree = 0; kDegree <= nDegree; kDegree++) {
				Coord_Control_Points[iDegree][jDegree][kDegree][0] = Coord_Corner_Points[0][0] 
				+ double(iDegree)/double(lDegree)*(Coord_Corner_Points[1][0]-Coord_Corner_Points[0][0]);
				Coord_Control_Points[iDegree][jDegree][kDegree][1] = Coord_Corner_Points[0][1] 
				+ double(jDegree)/double(mDegree)*(Coord_Corner_Points[3][1]-Coord_Corner_Points[0][1]);				
				Coord_Control_Points[iDegree][jDegree][kDegree][2] = Coord_Corner_Points[0][2] 
				+ double(kDegree)/double(nDegree)*(Coord_Corner_Points[4][2]-Coord_Corner_Points[0][2]);
			}
}

void CFreeFormChunk::SetSupportCP(CFreeFormChunk *chunk) {
	unsigned short iDim, iOrder, jOrder, kOrder;
	unsigned short lOrder = chunk->GetlOrder();
	unsigned short mOrder = chunk->GetmOrder();
	unsigned short nOrder = chunk->GetnOrder();
	
	Coord_SupportCP = new double*** [lOrder];
	for (iOrder = 0; iOrder < lOrder; iOrder++) {
		Coord_SupportCP[iOrder] = new double** [mOrder];
		for (jOrder = 0; jOrder < mOrder; jOrder++) {
			Coord_SupportCP[iOrder][jOrder] = new double* [nOrder];
			for (kOrder = 0; kOrder < nOrder; kOrder++)
				Coord_SupportCP[iOrder][jOrder][kOrder] = new double [nDim];
		}
	}
	
	/*--- Set base support control points according to the notation of Vtk for hexahedrons ---*/
	for (iDim = 0; iDim < nDim; iDim++) {
		Coord_SupportCP	[0]			[0]			[0]			[iDim]	= Coord_Corner_Points[0][iDim];
		Coord_SupportCP	[lOrder-1]	[0]			[0]			[iDim]	= Coord_Corner_Points[1][iDim];
		Coord_SupportCP	[lOrder-1]	[mOrder-1]	[0]			[iDim]	= Coord_Corner_Points[2][iDim];
		Coord_SupportCP	[0]			[mOrder-1]	[0]			[iDim]	= Coord_Corner_Points[3][iDim];
		Coord_SupportCP	[0]			[0]			[nOrder-1]	[iDim]	= Coord_Corner_Points[4][iDim];
		Coord_SupportCP	[lOrder-1]	[0]			[nOrder-1]	[iDim]	= Coord_Corner_Points[5][iDim];
		Coord_SupportCP	[lOrder-1]	[mOrder-1]	[nOrder-1]	[iDim]	= Coord_Corner_Points[6][iDim];
		Coord_SupportCP	[0]			[mOrder-1]	[nOrder-1]	[iDim]	= Coord_Corner_Points[7][iDim];
	}
	
	/*--- Fill the rest of the cubic matrix of support control points with uniform spacing  ---*/
	for (iOrder = 0; iOrder < lOrder; iOrder++)
		for (jOrder = 0; jOrder < mOrder; jOrder++)
			for (kOrder = 0; kOrder < nOrder; kOrder++) {
				Coord_SupportCP[iOrder][jOrder][kOrder][0] = Coord_Corner_Points[0][0] 
				+ double(iOrder)/double(lOrder-1)*(Coord_Corner_Points[1][0]-Coord_Corner_Points[0][0]);
				Coord_SupportCP[iOrder][jOrder][kOrder][1] = Coord_Corner_Points[0][1] 
				+ double(jOrder)/double(mOrder-1)*(Coord_Corner_Points[3][1]-Coord_Corner_Points[0][1]);				
				Coord_SupportCP[iOrder][jOrder][kOrder][2] = Coord_Corner_Points[0][2] 
				+ double(kOrder)/double(nOrder-1)*(Coord_Corner_Points[4][2]-Coord_Corner_Points[0][2]);
			}
}

void CFreeFormChunk::SetSupportCPChange(CFreeFormChunk *chunk) {
	unsigned short iDim, iOrder, jOrder, kOrder;
	double movement[3], *car_coord_old, *car_coord_new, *par_coord;
	unsigned short lOrder = chunk->GetlOrder();
	unsigned short mOrder = chunk->GetmOrder();
	unsigned short nOrder = chunk->GetnOrder();
	unsigned short *index = new unsigned short[nDim];

	double ****param_Coord_SupportCP = new double*** [lOrder];
	for (iOrder = 0; iOrder < lOrder; iOrder++) {
		param_Coord_SupportCP[iOrder] = new double** [mOrder];
		for (jOrder = 0; jOrder < mOrder; jOrder++) {
			param_Coord_SupportCP[iOrder][jOrder] = new double* [nOrder];
			for (kOrder = 0; kOrder < nOrder; kOrder++)
				param_Coord_SupportCP[iOrder][jOrder][kOrder] = new double [nDim];
		}
	}
	
	for (iOrder = 0; iOrder < lOrder; iOrder++)
		for (jOrder = 0; jOrder < mOrder; jOrder++)
			for (kOrder = 0; kOrder < nOrder; kOrder++)
				for (iDim = 0; iDim < nDim; iDim++)
					param_Coord_SupportCP[iOrder][jOrder][kOrder][iDim] = 
					Coord_SupportCP[iOrder][jOrder][kOrder][iDim];
	
	for (iDim = 0; iDim < nDim; iDim++) {
		Coord_Control_Points[0][0][0][iDim]	= chunk->GetCoordCornerPoints(iDim, 0);
		Coord_Control_Points[1][0][0][iDim]	= chunk->GetCoordCornerPoints(iDim, 1);
		Coord_Control_Points[1][1][0][iDim]	= chunk->GetCoordCornerPoints(iDim, 2);
		Coord_Control_Points[0][1][0][iDim]	= chunk->GetCoordCornerPoints(iDim, 3);
		Coord_Control_Points[0][0][1][iDim]	= chunk->GetCoordCornerPoints(iDim, 4);
		Coord_Control_Points[1][0][1][iDim]	= chunk->GetCoordCornerPoints(iDim, 5);
		Coord_Control_Points[1][1][1][iDim]	= chunk->GetCoordCornerPoints(iDim, 6);
		Coord_Control_Points[0][1][1][iDim]	= chunk->GetCoordCornerPoints(iDim, 7);
	}
	
	for (iOrder = 0; iOrder < chunk->GetlOrder(); iOrder++)
		for (jOrder = 0; jOrder < chunk->GetmOrder(); jOrder++)
			for (kOrder = 0; kOrder < chunk->GetnOrder(); kOrder++) {
				par_coord = param_Coord_SupportCP[iOrder][jOrder][kOrder];
				car_coord_new = EvalCartesianCoord(par_coord);
				car_coord_old = chunk->GetCoordControlPoints(iOrder, jOrder, kOrder);
				index[0] = iOrder; index[1] = jOrder; index[2] = kOrder;
				movement[0] = car_coord_new[0] - car_coord_old[0]; 
				movement[1] = car_coord_new[1] - car_coord_old[1]; 
				movement[2] = car_coord_new[2] - car_coord_old[2]; 
				chunk->SetControlPoints(index, movement);
			}
}

void CFreeFormChunk::SetParaView (char chunk_filename[200], bool new_file) {
	ofstream chunk_file;
	unsigned short iCornerPoints, iDim, iControlPoints, iDegree, jDegree, kDegree;
		
	if (new_file) chunk_file.open(chunk_filename, ios::out);
	else chunk_file.open(chunk_filename, ios::out | ios::app);
	
	chunk_file << "# vtk DataFile Version 2.0" << endl;
	chunk_file << "Visualization of the FFD box" << endl;
	chunk_file << "ASCII" << endl;
	chunk_file << "DATASET UNSTRUCTURED_GRID" << endl;
	chunk_file << "POINTS " << nCornerPoints + nControlPoints << " float" << endl;
	
	chunk_file.precision(15);
	
	for(iCornerPoints = 0; iCornerPoints < nCornerPoints; iCornerPoints++) {
		for(iDim = 0; iDim < nDim; iDim++)
			chunk_file << scientific << Coord_Corner_Points[iCornerPoints][iDim] << "\t";
		chunk_file << "\n";
	}
	for (iDegree = 0; iDegree <= lDegree; iDegree++)
		for (jDegree = 0; jDegree <= mDegree; jDegree++)
			for (kDegree = 0; kDegree <= nDegree; kDegree++) {
				for(iDim = 0; iDim < nDim; iDim++)
					chunk_file << scientific << Coord_Control_Points[iDegree][jDegree][kDegree][iDim] << "\t";
				chunk_file << "\n";
			}
	
	chunk_file << "CELLS " << 1 + nControlPoints << "\t" << (8+1) + (1+1) * nControlPoints << endl;
	chunk_file << "8 0 1 2 3 4 5 6 7" << endl;
	for (iControlPoints = 0; iControlPoints < nControlPoints; iControlPoints++)
		chunk_file << "1 " << iControlPoints + 8 << endl;
	
	chunk_file << "CELL_TYPES " << 1 + nControlPoints<< endl;
	chunk_file << "12" << endl;
	for (iControlPoints = 0; iControlPoints < nControlPoints; iControlPoints++)
		chunk_file << "1" << endl;
	
	chunk_file.close();
}

void CFreeFormChunk::SetTecplot(char chunk_filename[200], bool new_file) {
	ofstream chunk_file;
	unsigned short iDim, iDegree, jDegree, kDegree;
	
	if (new_file) {
		chunk_file.open(chunk_filename, ios::out);
		chunk_file << "TITLE = \"Visualization of the FFD box\"" << endl;
		chunk_file << "VARIABLES = \"x\", \"y\", \"z\"" << endl;
	}
	else chunk_file.open(chunk_filename, ios::out | ios::app);

	chunk_file << "ZONE I="<<lDegree+1<<", J="<<mDegree+1<<", K="<<nDegree+1<<", DATAPACKING=POINT" << endl;
	
	chunk_file.precision(15);
	
	for (kDegree = 0; kDegree <= nDegree; kDegree++)
		for (jDegree = 0; jDegree <= mDegree; jDegree++)
			for (iDegree = 0; iDegree <= lDegree; iDegree++) {
				for(iDim = 0; iDim < nDim; iDim++)
					chunk_file << scientific << Coord_Control_Points[iDegree][jDegree][kDegree][iDim] << "\t";
				chunk_file << "\n";
			}
		
	chunk_file.close();
}


double *CFreeFormChunk::GetParametricCoord_Analytical(double *cart_coord) {
	unsigned short iDim;
	double *e1, *e2, *e3, *e12, *e23, *e13, *p;
	
	/*--- Auxiliary Basis Vectors of the deformed chunk ---*/
	e1 = new double[3]; e2 = new double[3]; e3 = new double[3];
	for (iDim = 0; iDim < nDim; iDim++) {
		e1[iDim] = Coord_Corner_Points[1][iDim]-Coord_Corner_Points[0][iDim];
		e2[iDim] = Coord_Corner_Points[3][iDim]-Coord_Corner_Points[0][iDim];
		e3[iDim] = Coord_Corner_Points[4][iDim]-Coord_Corner_Points[0][iDim];
	}
	
	/*--- Respective Cross-Products ---*/
	e12 = new double[3]; e23 = new double[3]; e13 = new double[3];
	CrossProduct(e1,e2,e12);
	CrossProduct(e1,e3,e13);
	CrossProduct(e2,e3,e23);
	
	/*--- p is Tranlated vector from the origin ---*/
	p = new double[3];
	for (iDim = 0; iDim < nDim; iDim++)
		p[iDim] = cart_coord[iDim] - Coord_Corner_Points[0][iDim];
	
	param_coord[0] = DotProduct(e23,p)/DotProduct(e23,e1);
	param_coord[1] = DotProduct(e13,p)/DotProduct(e13,e2);
	param_coord[2] = DotProduct(e12,p)/DotProduct(e12,e3);
	
	delete [] e1;
  delete [] e2;
  delete [] e3;
  delete [] e12;
  delete [] e23;
  delete [] e13;
  delete [] p;
	
	return param_coord;
}

double *CFreeFormChunk::EvalCartesianCoord(double *param_coord) {
	unsigned short iDim, iDegree, jDegree, kDegree;
	
	for (iDim = 0; iDim < nDim; iDim++)
		cart_coord[iDim] = 0;
	
	for (iDegree = 0; iDegree <= lDegree; iDegree++)
		for (jDegree = 0; jDegree <= mDegree; jDegree++)
			for (kDegree = 0; kDegree <= nDegree; kDegree++)
				for (iDim = 0; iDim < nDim; iDim++) {
					cart_coord[iDim] += Coord_Control_Points[iDegree][jDegree][kDegree][iDim]
					* GetBernstein(lDegree, iDegree, param_coord[0])
					* GetBernstein(mDegree, jDegree, param_coord[1])
					* GetBernstein(nDegree, kDegree, param_coord[2]);
				}
	
	return cart_coord;
}

double CFreeFormChunk::GetBernstein(short val_n, short val_i, double val_t) {
	double value;

	if (val_i > val_n) { value = 0; return value; }
	if (val_i == 0) {
		if (val_t == 0) value = 1;
		else if (val_t == 1) value = 0;
		else value = Binomial(val_n,val_i)*(pow(val_t, val_i)) * pow(1.0 - val_t, val_n - val_i);
	}
	else if (val_i == val_n) {
		if (val_t == 0) value = 0;
		else if (val_t == 1) value = 1;
		else value = pow(val_t,val_n);
	}
	else value = Binomial(val_n,val_i)*(pow(val_t,val_i)) * pow(1.0-val_t, val_n - val_i);
	
	return value;
}

double CFreeFormChunk::GetBernsteinDerivative(short val_n, short val_i, 
											   double val_t, short val_order) {
	double value = 0.0;
	
	/*--- Verify this subroutine, it provides negative val_n, 
	 which is a wrong value for GetBernstein ---*/
	
	if (val_order == 0) { 
		value = GetBernstein(val_n, val_i, val_t); return value; 
	}
	
	if (val_i == 0) { 
		value = val_n*(-GetBernsteinDerivative(val_n-1, val_i, val_t, val_order-1)); return value; 
	}
	else {
		if (val_n == 0) { 
			value = val_t; return value; 
		}
		else {
			value = val_n*(GetBernsteinDerivative(val_n-1, val_i-1, val_t, val_order-1) - GetBernsteinDerivative(val_n-1, val_i, val_t, val_order-1));
			return value;
		}
	}

	return value;
}

double *CFreeFormChunk::GetGradient_Analytical(double *val_coord, double *xyz) {
	unsigned short iDim, jDim, lmn[3];
	
	/*--- Set the Degree of the Berstein polynomials ---*/
	lmn[0] = lDegree; lmn[1] = mDegree; lmn[2] = nDegree;
	
	for (iDim = 0; iDim < nDim; iDim++) gradient[iDim] = 0;
	
	for (iDim = 0; iDim < nDim; iDim++)
		for (jDim = 0; jDim < nDim; jDim++)
			gradient[jDim] += GetDerivative2(val_coord, iDim, xyz,  lmn) *  
			GetDerivative3(val_coord, iDim, jDim, lmn);
	
	return gradient;
}

double *CFreeFormChunk::GetGradient_Numerical(double *uvw, double *xyz) {
	double delta = 1E-6, parametric[3], *coord_eval, functional_plus, functional_minus;
	
	parametric[0] = uvw[0] + delta;
	parametric[1] = uvw[1]; 
	parametric[2] = uvw[2]; 
	coord_eval = EvalCartesianCoord(parametric);
	functional_plus = ((coord_eval[0]-xyz[0])*(coord_eval[0]-xyz[0]) + 
					   (coord_eval[1]-xyz[1])*(coord_eval[1]-xyz[1]) +
					   (coord_eval[2]-xyz[2])*(coord_eval[2]-xyz[2]));
	parametric[0] = uvw[0] - delta;
	parametric[1] = uvw[1]; 
	parametric[2] = uvw[2]; 
	coord_eval = EvalCartesianCoord(parametric);
	functional_minus = ((coord_eval[0]-xyz[0])*(coord_eval[0]-xyz[0]) + 
						(coord_eval[1]-xyz[1])*(coord_eval[1]-xyz[1]) +
						(coord_eval[2]-xyz[2])*(coord_eval[2]-xyz[2]));
	gradient[0] = 0.5*(functional_plus-functional_minus)/delta;
	
	parametric[0] = uvw[0];
	parametric[1] = uvw[1] + delta;
	parametric[2] = uvw[2]; 
	coord_eval = EvalCartesianCoord(parametric);
	functional_plus = ((coord_eval[0]-xyz[0])*(coord_eval[0]-xyz[0]) + 
					   (coord_eval[1]-xyz[1])*(coord_eval[1]-xyz[1]) +
					   (coord_eval[2]-xyz[2])*(coord_eval[2]-xyz[2]));
	parametric[0] = uvw[0];
	parametric[1] = uvw[1] - delta;
	parametric[2] = uvw[2]; 
	coord_eval = EvalCartesianCoord(parametric);
	functional_minus = ((coord_eval[0]-xyz[0])*(coord_eval[0]-xyz[0]) + 
						(coord_eval[1]-xyz[1])*(coord_eval[1]-xyz[1]) +
						(coord_eval[2]-xyz[2])*(coord_eval[2]-xyz[2]));
	gradient[1] = 0.5*(functional_plus-functional_minus)/delta;
	
	parametric[0] = uvw[0];
	parametric[1] = uvw[1]; 
	parametric[2] = uvw[2] + delta;
	coord_eval = EvalCartesianCoord(parametric);
	functional_plus = ((coord_eval[0]-xyz[0])*(coord_eval[0]-xyz[0]) + 
					   (coord_eval[1]-xyz[1])*(coord_eval[1]-xyz[1]) +
					   (coord_eval[2]-xyz[2])*(coord_eval[2]-xyz[2]));
	parametric[0] = uvw[0];
	parametric[1] = uvw[1]; 
	parametric[2] = uvw[2] - delta;
	coord_eval = EvalCartesianCoord(parametric);
	functional_minus = ((coord_eval[0]-xyz[0])*(coord_eval[0]-xyz[0]) + 
						(coord_eval[1]-xyz[1])*(coord_eval[1]-xyz[1]) +
						(coord_eval[2]-xyz[2])*(coord_eval[2]-xyz[2]));
	gradient[2] = 0.5*(functional_plus-functional_minus)/delta;
	
	return gradient;
}

double *CFreeFormChunk::GetParametricCoord_Iterative(double *xyz, double *guess, double tol, 
																										 unsigned long it_max) {
	double **Hessian, Indep_Term[3], under_relax = 1.0, MinNormError, NormError;
	unsigned short iDim, RandonCounter;
	unsigned long iter;
	
	/*--- Allocate the Hessian ---*/
	Hessian = new double* [nDim];
	for (iDim = 0; iDim < nDim; iDim++) {
		Hessian[iDim] = new double[nDim];
		param_coord[iDim] = guess[iDim];
		Indep_Term [iDim] = 0.0;
	}
	
	RandonCounter = 0; MinNormError = 1E6;
	
	for (iter = 0; iter < it_max; iter++) {
		
		/*--- The independent term of the solution of our system is -Gradient(sol_old) ---*/
		gradient = GetGradient_Analytical(param_coord, xyz);
		
		for (iDim = 0; iDim < nDim; iDim++) 
			Indep_Term[iDim] = -gradient[iDim];
						
		/*--- Relaxation of the Newton Method ---*/
		for (iDim = 0; iDim < nDim; iDim++) 
			Indep_Term[iDim] = under_relax * Indep_Term[iDim];
		
		/*--- Hessian = The Matrix of our system, getHessian(sol_old,xyz,...) ---*/
		GetHessian_Analytical(param_coord, xyz, Hessian);
		
		/*--- Gauss elimination algorithm. Solution will be stored on Indep_Term ---*/
		Gauss_Elimination(Hessian, Indep_Term, nDim);				
		
		/*--- Solution is in fact par_new-par_old; Must Update doing par_new=par_old + solution ---*/
		for (iDim = 0; iDim < nDim; iDim++) 
			param_coord[iDim] += Indep_Term[iDim];
		
		/*--- If the gradient is small, we have converged ---*/
		if ((abs(Indep_Term[0]) < tol) && (abs(Indep_Term[1]) < tol) && (abs(Indep_Term[2]) < tol))	break;
		NormError = sqrt(Indep_Term[0]*Indep_Term[0] + Indep_Term[1]*Indep_Term[1] + Indep_Term[2]*Indep_Term[2]);
		MinNormError = min(NormError, MinNormError);
		
		/*--- If we have no convergence with 50 iterations probably we are out of the chunk, then 
		 we try with a ramdom choice ---*/
		if (((iter % 50) == 0) && (iter != 0)) {
			RandonCounter++;
			param_coord[0] = double(rand())/double(RAND_MAX);
			param_coord[1] = double(rand())/double(RAND_MAX);
			param_coord[2] = double(rand())/double(RAND_MAX);
		}
		
		if (RandonCounter == 100) {
			cout << "I can not localize this point: " << xyz[0] <<" "<< xyz[1] <<" "<< xyz[2] <<". Min Error: "<< MinNormError <<"."<< endl;
			param_coord[0] = 0.0; param_coord[1] = 0.0; param_coord[2] = 0.0;
			Indep_Term[0] = 0.0; Indep_Term[1] = 0.0; Indep_Term[2] = 0.0;
			break;
		}
	}
	
	/*--- There is no convergence of the point inversion algorithm ---*/
	if ((abs(Indep_Term[0]) > tol) || (abs(Indep_Term[1]) > tol) || (abs(Indep_Term[2]) > tol)) 
		cout << "No Convergence Detected After " << iter << " Iterations" << endl;
	
	for (iDim = 0; iDim < nDim; iDim++) 
		delete [] Hessian[iDim];
	delete [] Hessian;
	
	/*--- Real Solution is now par_coord; Return it ---*/
	return param_coord;
}

unsigned short CFreeFormChunk::Binomial (unsigned short n, unsigned short m) {
	unsigned short result;

	if ( (m == 0) || (m == n) ) result = 1;
	else result = Factorial(n) / (Factorial(n-m)*Factorial(m));
		
	return result;
}

unsigned long CFreeFormChunk::BinomialOpt (unsigned long n, unsigned long m) {
	unsigned long b[100], i , j;
	if (n+1 > 100) cout << "ERROR!!! Increase the size of b in the BinomialOpt subroutine!" <<endl;
	
	b[0] = 1;
	for (i = 1; i <= n; ++i) {
		b[i] = 1;
		for(j = i-1U; j > 0; --j)
			b[j] += b[j-1U];
	}
	
	return b[m];
}

unsigned short CFreeFormChunk::Factorial (unsigned short n) {
	
	if ( n > 1 ) n = n*Factorial(n-1);
	if ( n == 0 ) n = 1;
	
	return n;
}


bool CFreeFormChunk::GetPointFFD(CGeometry *geometry, CConfig *config, unsigned long iPoint) {
	double *Coord;
	unsigned short iVar, jVar;
	bool Inside;
	
	unsigned short Index[5][7] = {
		{0, 1, 2, 5, 0, 1, 2},
		{0, 2, 7, 5, 0, 2, 7},
		{0, 2, 3, 7, 0, 2, 3},
		{0, 5, 7, 4, 0, 5, 7},
		{2, 7, 5, 6, 2, 7, 5}};
	
	Coord = geometry->node[iPoint]->GetCoord();
	
	/*--- 1st tetrahedron {V0, V1, V2, V5}
	 2nd tetrahedron {V0, V2, V7, V5}
	 3th tetrahedron {V0, V2, V3, V7}
	 4th tetrahedron {V0, V5, V7, V4}
	 5th tetrahedron {V2, V7, V5, V6} ---*/
	
	for (iVar = 0; iVar < 5; iVar++) {
		Inside = true;
		for (jVar = 0; jVar < 4; jVar++) {
			double Distance_Point = geometry->Point2Plane_Distance(Coord, 
																														 Coord_Corner_Points[Index[iVar][jVar+1]], 
																														 Coord_Corner_Points[Index[iVar][jVar+2]], 
																														 Coord_Corner_Points[Index[iVar][jVar+3]]);
			
			double Distance_Vertex = geometry->Point2Plane_Distance(Coord_Corner_Points[Index[iVar][jVar]], 
																															Coord_Corner_Points[Index[iVar][jVar+1]], 
																															Coord_Corner_Points[Index[iVar][jVar+2]], 
																															Coord_Corner_Points[Index[iVar][jVar+3]]);
			if (Distance_Point*Distance_Vertex < 0.0) Inside = false;					
		}
		if (Inside) break;
	}
	
	return Inside;

}

void CFreeFormChunk::SetDeformationZone(CGeometry *geometry, CConfig *config, unsigned short iChunk) {
	double *Coord;
	unsigned short iMarker, iVar, jVar;
	unsigned long iVertex, iPoint;
	bool Inside;
	
	unsigned short Index[5][7] = {
		{0, 1, 2, 5, 0, 1, 2},
		{0, 2, 7, 5, 0, 2, 7},
		{0, 2, 3, 7, 0, 2, 3},
		{0, 5, 7, 4, 0, 5, 7},
		{2, 7, 5, 6, 2, 7, 5}};
	
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
		if (config->GetMarker_All_Moving(iMarker) == YES)
			for(iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {	
				iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
				geometry->node[iPoint]->SetMove(false);
				
				Coord = geometry->node[iPoint]->GetCoord();
				
				/*--- 1st tetrahedron {V0, V1, V2, V5}
				 2nd tetrahedron {V0, V2, V7, V5}
				 3th tetrahedron {V0, V2, V3, V7}
				 4th tetrahedron {V0, V5, V7, V4}
				 5th tetrahedron {V2, V7, V5, V6} ---*/
				
				for (iVar = 0; iVar < 5; iVar++) {
					Inside = true;
					for (jVar = 0; jVar < 4; jVar++) {
						double Distance_Point = geometry->Point2Plane_Distance(Coord, 
																																	 Coord_Corner_Points[Index[iVar][jVar+1]], 
																																	 Coord_Corner_Points[Index[iVar][jVar+2]], 
																																	 Coord_Corner_Points[Index[iVar][jVar+3]]);				
						double Distance_Vertex = geometry->Point2Plane_Distance(Coord_Corner_Points[Index[iVar][jVar]], 
																																		Coord_Corner_Points[Index[iVar][jVar+1]], 
																																		Coord_Corner_Points[Index[iVar][jVar+2]], 
																																		Coord_Corner_Points[Index[iVar][jVar+3]]);
						if (Distance_Point*Distance_Vertex < 0.0) Inside = false;					
					}
					if (Inside) break;
				}
				
				if (Inside) {
					geometry->node[iPoint]->SetMove(true);
				}
				
			}
}

double CFreeFormChunk::GetDerivative1 (double *uvw, unsigned short val_diff, unsigned short *ijk, unsigned short *lmn) {
	unsigned short iDim;
	double value = GetBernsteinDerivative(lmn[val_diff], ijk[val_diff], uvw[val_diff], 1);
	
	for (iDim = 0; iDim < nDim; iDim++)
		if (iDim != val_diff)
			value *= GetBernstein(lmn[iDim], ijk[iDim], uvw[iDim]);
	
	return value;	
}

double CFreeFormChunk::GetDerivative2 (double *uvw, unsigned short dim, double *xyz, unsigned short *lmn) {
	
	unsigned short iDegree, jDegree, kDegree;
	double value = 0.0;
	
	for (iDegree = 0; iDegree <= lmn[0]; iDegree++)
		for (jDegree = 0; jDegree <= lmn[1]; jDegree++)
			for (kDegree = 0; kDegree <= lmn[2]; kDegree++)
				value += Coord_Control_Points[iDegree][jDegree][kDegree][dim] 
				* GetBernstein(lmn[0], iDegree, uvw[0])
				* GetBernstein(lmn[1], jDegree, uvw[1])
				* GetBernstein(lmn[2], kDegree, uvw[2]);
	
	return 2.0*(value - xyz[dim]);	
}

double CFreeFormChunk::GetDerivative3(double *uvw, unsigned short dim, unsigned short diff_this, unsigned short *lmn) {
	unsigned short iDegree, jDegree, kDegree, ijk[3];
	double value = 0;
	
	for (iDegree = 0; iDegree <= lmn[0]; iDegree++)
		for (jDegree = 0; jDegree <= lmn[1]; jDegree++)
			for (kDegree = 0; kDegree <= lmn[2]; kDegree++) {
				ijk[0] = iDegree; ijk[1] = jDegree; ijk[2] = kDegree;
				value += Coord_Control_Points[iDegree][jDegree][kDegree][dim] * 
				GetDerivative1(uvw, diff_this, ijk, lmn);
			}
	
	return value;
}

double CFreeFormChunk::GetDerivative4 (double *uvw, unsigned short val_diff, unsigned short val_diff2,
																			 unsigned short *ijk, unsigned short *lmn) {
	unsigned short iDim;
	double value;
	
	if (val_diff == val_diff2) {
		value = GetBernsteinDerivative(lmn[val_diff], ijk[val_diff], uvw[val_diff], 2);
		for (iDim = 0; iDim < nDim; iDim++)
			if (iDim != val_diff)
				value *= GetBernstein(lmn[iDim], ijk[iDim], uvw[iDim]);
	}
	else {
		value = GetBernsteinDerivative(lmn[val_diff],  ijk[val_diff],  uvw[val_diff], 1) *
		GetBernsteinDerivative(lmn[val_diff2], ijk[val_diff2], uvw[val_diff2], 1);
		for (iDim = 0; iDim < nDim; iDim++)
			if ((iDim != val_diff) && (iDim != val_diff2))
				value *= GetBernstein(lmn[iDim], ijk[iDim], uvw[iDim]);
	}
	
	return value;
}

double CFreeFormChunk::GetDerivative5(double *uvw, unsigned short dim, unsigned short diff_this, unsigned short diff_this_also, 
																			unsigned short *lmn) {
	
	unsigned short iDegree, jDegree, kDegree, ijk[3];
	double value = 0.0;
	
	for (iDegree = 0; iDegree <= lmn[0]; iDegree++)
		for (jDegree = 0; jDegree <= lmn[1]; jDegree++)
			for (kDegree = 0; kDegree <= lmn[2]; kDegree++) {
				ijk[0] = iDegree; ijk[1] = jDegree; ijk[2] = kDegree;
				value += Coord_Control_Points[iDegree][jDegree][kDegree][dim] *
				GetDerivative4(uvw, diff_this, diff_this_also, ijk, lmn);
			}
	
	return value;
}

void CFreeFormChunk::GetHessian_Analytical(double *uvw, double *xyz, double **val_Hessian) {
	
	unsigned short iDim, jDim;
	unsigned short l, m, n, lmn[3];
	
	/*--- Set the Degree of the Berstein polynomials ---*/
	lmn[0] = lDegree; lmn[1] = mDegree; lmn[2] = nDegree;
	
	/*--- Berstein polynomials degrees ---*/
	l = lmn[0]; m = lmn[1]; n = lmn[2];
	
	for (iDim = 0; iDim < nDim; iDim++)
		for (jDim = 0; jDim < nDim; jDim++)
			val_Hessian[iDim][jDim] = 0.0;
	
	/*--- Note that being all the functions linear combinations of polynomials, they are C^\infty,
	 and the Hessian will be symmetric; no need to compute the under-diagonal part, for example ---*/
	for (iDim = 0; iDim < nDim; iDim++) {
		val_Hessian[0][0] += 2.0 * GetDerivative3(uvw,iDim,0,lmn) * GetDerivative3(uvw,iDim,0,lmn) + 
		GetDerivative2(uvw,iDim,xyz,lmn) * GetDerivative5(uvw,iDim,0,0,lmn);
		
		val_Hessian[1][1] += 2.0 * GetDerivative3(uvw,iDim,1,lmn) * GetDerivative3(uvw,iDim,1,lmn) + 
		GetDerivative2(uvw,iDim,xyz,lmn) * GetDerivative5(uvw,iDim,1,1,lmn);
		
		val_Hessian[2][2] += 2.0 * GetDerivative3(uvw,iDim,2,lmn) * GetDerivative3(uvw,iDim,2,lmn) + 
		GetDerivative2(uvw,iDim,xyz,lmn) * GetDerivative5(uvw,iDim,2,2,lmn);
		
		val_Hessian[0][1] += 2.0 * GetDerivative3(uvw,iDim,0,lmn) * GetDerivative3(uvw,iDim,1,lmn) +
		GetDerivative2(uvw,iDim,xyz,lmn) * GetDerivative5(uvw,iDim,0,1,lmn);
		
		val_Hessian[0][2] += 2.0 * GetDerivative3(uvw,iDim,0,lmn) * GetDerivative3(uvw,iDim,2,lmn) +
		GetDerivative2(uvw,iDim,xyz,lmn) * GetDerivative5(uvw,iDim,0,2,lmn);
		
		val_Hessian[1][2] += 2.0 * GetDerivative3(uvw,iDim,1,lmn) * GetDerivative3(uvw,iDim,2,lmn) +
		GetDerivative2(uvw,iDim,xyz,lmn) * GetDerivative5(uvw,iDim,1,2,lmn);
	}
	
	val_Hessian[1][0] = val_Hessian[0][1];
	val_Hessian[2][0] = val_Hessian[0][2];
	val_Hessian[2][1] = val_Hessian[1][2];
}

void CFreeFormChunk::Gauss_Elimination(double** A, double* rhs, unsigned short nVar) {
	unsigned short jVar, kVar, iVar;
    double weight, aux;
	
	if (nVar == 1)
		rhs[0] /= (A[0][0]+EPS*EPS);
	else {
		/*--- Transform system in Upper Matrix ---*/
		for (iVar = 1; iVar < nVar; iVar++) {
			for (jVar = 0; jVar < iVar; jVar++) {
				weight = A[iVar][jVar]/(A[jVar][jVar]+EPS*EPS);
				for (kVar = jVar; kVar < nVar; kVar++)
					A[iVar][kVar] -= weight*A[jVar][kVar];
				rhs[iVar] -= weight*rhs[jVar];
			}
		}
		/*--- Backwards substitution ---*/
		rhs[nVar-1] = rhs[nVar-1]/(A[nVar-1][nVar-1]+EPS*EPS);
		for (short iVar = nVar-2; iVar >= 0; iVar--) {
			aux = 0;
			for (jVar = iVar+1; jVar < nVar; jVar++)
				aux += A[iVar][jVar]*rhs[jVar];
			rhs[iVar] = (rhs[iVar]-aux)/(A[iVar][iVar]+EPS*EPS);
			if (iVar == 0) break;
		}
	}
}
