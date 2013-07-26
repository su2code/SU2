/*!
 * \file numerics_direct_electric.cpp
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

CGalerkin_Flow::CGalerkin_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;
}

CGalerkin_Flow::~CGalerkin_Flow(void) { }

void CGalerkin_Flow::ComputeResidual(double **val_stiffmatrix_elem, CConfig *config) {
  
	double a[4], b[4], c[4], d[4], Area, B_Matrix[4][4];
	unsigned short iVar, jVar;
  
	if (nDim == 2) {
		for (unsigned short iDim = 0; iDim < nDim; iDim++) {
			a[iDim] = Coord_0[iDim]-Coord_2[iDim];
			b[iDim] = Coord_1[iDim]-Coord_2[iDim];
		}
    
		Area = 0.5*fabs(a[0]*b[1]-a[1]*b[0]);	/* Norm of the normal component of area, area = 1/2*cross(a,b) */
    
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

CSourcePieceWise_Elec::CSourcePieceWise_Elec(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;
	Ni_times_Nj = new double*[4];
	for (unsigned short iVar = 0; iVar < 4; iVar++) {
		Ni_times_Nj[iVar] = new double [4];
	}
  
	Ni_times_Nj[0][0] = 2.0/120*6;	Ni_times_Nj[0][1] = 1.0/120*6;	Ni_times_Nj[0][2] = 1.0/120*6;	Ni_times_Nj[0][3] = 1.0/120*6;
	Ni_times_Nj[1][0] = 1.0/120*6;	Ni_times_Nj[1][1] = 2.0/120*6;	Ni_times_Nj[1][2] = 1.0/120*6;	Ni_times_Nj[1][3] = 1.0/120*6;
	Ni_times_Nj[2][0] = 1.0/120*6;	Ni_times_Nj[2][1] = 1.0/120*6;	Ni_times_Nj[2][2] = 2.0/120*6;	Ni_times_Nj[2][3] = 1.0/120*6;
	Ni_times_Nj[3][0] = 1.0/120*6;	Ni_times_Nj[3][1] = 1.0/120*6;	Ni_times_Nj[3][2] = 1.0/120*6;	Ni_times_Nj[3][3] = 2.0/120*6;
  
}

CSourcePieceWise_Elec::~CSourcePieceWise_Elec(void) { }

void CSourcePieceWise_Elec::ComputeResidual_MacCormack(double *val_residual, CConfig *config) {
  
	if (config->GetKind_GasModel() == ARGON) {
		double Kai_n, Kai_np1 = 0.0, ne, np;
		double diff_ru_Elec_x, diff_ru_Pos_x, diff_ru_Elec_y, diff_ru_Pos_y;
		double alpha;
		double dt = TimeStep;
		double rho_Pos = 0.0, rho_Elec = 0.0, mass_Elec, mass_Pos;
		double a[4], b[4], Area = 0.0;
		if (nDim == 2) {
			for (unsigned short iDim = 0; iDim < nDim; iDim++) {
				a[iDim] = Coord_0[iDim]-Coord_2[iDim];
				b[iDim] = Coord_1[iDim]-Coord_2[iDim];
			}
			Area = 0.5*fabs(a[0]*b[1]-a[1]*b[0]);			/*--- Norm of the normal component of area, area = 1/2*cross(a,b) ---*/
      
			rho_Pos  = 1.0/3.0*(U_0[0] + U_1[0] + U_2[0]);
			rho_Elec = 1.0/3.0*(U_0[1] + U_1[1] + U_2[1]);
      
			/*--- Source q ---*/
			mass_Elec = config->GetParticle_Mass(2);
			mass_Pos  = config->GetParticle_Mass(1);
      
			ne = rho_Elec / mass_Elec;
			np = rho_Pos  / mass_Pos;
      
			Kai_n = ELECTRON_CHARGE/FREE_PERMITTIVITY * (ne - np );
			alpha = pow(dt*ELECTRON_CHARGE,2)/ FREE_PERMITTIVITY * (rho_Elec / (mass_Elec*mass_Elec) + rho_Pos / (mass_Pos*mass_Pos));
      
			diff_ru_Pos_x  = 1.0/3.0*(ConsVar_Grad_0[1][0] + ConsVar_Grad_1[1][0] + ConsVar_Grad_2[1][0]);
			diff_ru_Elec_x = 1.0/3.0*(ConsVar_Grad_0[2][0] + ConsVar_Grad_1[2][0] + ConsVar_Grad_2[2][0]);
			diff_ru_Pos_y  = 1.0/3.0*(ConsVar_Grad_0[1][1] + ConsVar_Grad_1[1][1] + ConsVar_Grad_2[1][1]);
			diff_ru_Elec_y = 1.0/3.0*(ConsVar_Grad_0[2][1] + ConsVar_Grad_1[2][1] + ConsVar_Grad_2[2][1]);
      
      
			Kai_np1 = 1.0/(1.0+alpha) * (Kai_n - dt*ELECTRON_CHARGE/FREE_PERMITTIVITY*((diff_ru_Elec_x+diff_ru_Elec_y)/mass_Elec - (diff_ru_Pos_x+diff_ru_Pos_y)/mass_Pos) );
      
			/*--- Residual = transpose(N) * source (q) * Area  ---*/
			val_residual[0] =-1.0/3.0 * Kai_np1 * Area;
			val_residual[1] =-1.0/3.0 * Kai_np1 * Area;
			val_residual[2] =-1.0/3.0 * Kai_np1 * Area;
      
      //            if (fabs(diff_ru_Elec_x) > 1E-5) {
      //                cout << " alpha = " << alpha << " dt = " << dt << " Kai_np1 = " << Kai_np1 << " Kai_n = " << Kai_n;
      //                cout << " n_rho = " << rho_Elec/mass_Elec - rho_Pos/mass_Pos << endl;
      //                cout << " diff_ru_Elec_x = " << diff_ru_Elec_x;
      //                cout << " diff_ru_Elec_y = " << diff_ru_Elec_y << endl;
      //                cout << " source = " << val_residual[0] << ", " << val_residual[1] << ", " << val_residual[2] << endl;
      //            }
      
		}
	}
	if (config->GetKind_GasModel() == AIR7 || config->GetKind_GasModel() == ARGON_SID || config->GetKind_GasModel() == O2 || config->GetKind_GasModel() == N2 || config->GetKind_GasModel() == AIR5) {
		double Chi_n_0, Chi_n_1, Chi_n_2;
		double Chi_np1_0 = 0.0, Chi_np1_1 = 0.0, Chi_np1_2 = 0.0;
		double ne_0, ne_1, ne_2;
		double np_0, np_1, np_2;
		double diff_ru_Neg_x_0 = 0.0, diff_ru_Neg_x_1 = 0.0, diff_ru_Neg_x_2 = 0.0;
		double diff_ru_Pos_x_0 = 0.0, diff_ru_Pos_x_1 = 0.0, diff_ru_Pos_x_2 = 0.0;
		double diff_ru_Neg_y_0 = 0.0, diff_ru_Neg_y_1 = 0.0, diff_ru_Neg_y_2 = 0.0;
		double diff_ru_Pos_y_0 = 0.0, diff_ru_Pos_y_1 = 0.0, diff_ru_Pos_y_2 = 0.0;
		double alpha_0, alpha_1, alpha_2;
		double dt = TimeStep;
		double rho_Pos_0 = 0.0, rho_Pos_1 = 0.0, rho_Pos_2 = 0.0;
		double rho_Neg_0 = 0.0, rho_Neg_1 = 0.0, rho_Neg_2 = 0.0;
		double mass_Pos_0 = 0.0, mass_Pos_1 = 0.0, mass_Pos_2 = 0.0;
		double massTotal_Pos_0 = 0.0, massTotal_Pos_1 = 0.0, massTotal_Pos_2 = 0.0;
		double mass_Neg_0 = 0.0, mass_Neg_1 = 0.0, mass_Neg_2 = 0.0;
		double massTotal_Neg_0 = 0.0, massTotal_Neg_1 = 0.0, massTotal_Neg_2 = 0.0;
		double a[4], b[4], Area = 0.0;
		unsigned short counterPos = 0, counterNeg = 0;
		unsigned short loc;
		unsigned short iSpecies;
    
		if (nDim == 2) {
			for (unsigned short iDim = 0; iDim < nDim; iDim++) {
				a[iDim] = Coord_0[iDim]-Coord_2[iDim];
				b[iDim] = Coord_1[iDim]-Coord_2[iDim];
			}
			Area = 0.5*fabs(a[0]*b[1]-a[1]*b[0]);			/*--- Norm of the normal component of area, area = 1/2*cross(a,b) ---*/
      
      
			for (iSpecies = 0; iSpecies < config->GetnSpecies(); iSpecies++) {
				if ( iSpecies < nDiatomics ) loc = (nDim+3)*iSpecies;
				else loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
				if (config->GetParticle_ChargeNumber(iSpecies) > 0) {
					rho_Pos_0 += U_0[loc+0];
					rho_Pos_1 += U_1[loc+0];
					rho_Pos_2 += U_2[loc+0];
					mass_Pos_0 += U_0[loc+0] * config->GetMolar_Mass(iSpecies);
					mass_Pos_1 += U_1[loc+0] * config->GetMolar_Mass(iSpecies);
					mass_Pos_2 += U_2[loc+0] * config->GetMolar_Mass(iSpecies);
					massTotal_Pos_0 += U_0[loc+0];
					massTotal_Pos_1 += U_1[loc+0];
					massTotal_Pos_2 += U_2[loc+0];
					diff_ru_Pos_x_0 += ConsVar_Grad_0[loc+1][0];
					diff_ru_Pos_x_1 += ConsVar_Grad_1[loc+1][0];
					diff_ru_Pos_x_2 += ConsVar_Grad_2[loc+1][0];
					diff_ru_Pos_y_0 += ConsVar_Grad_0[loc+2][1];
					diff_ru_Pos_y_1 += ConsVar_Grad_1[loc+2][1];
					diff_ru_Pos_y_2 += ConsVar_Grad_2[loc+2][1];
					counterPos++;
				} else if (config->GetParticle_ChargeNumber(iSpecies) < 0) {
					rho_Neg_0 += U_0[loc+0];
					rho_Neg_1 += U_1[loc+0];
					rho_Neg_2 += U_2[loc+0];
					mass_Neg_0 += U_0[loc+0] * config->GetMolar_Mass(iSpecies);
					mass_Neg_1 += U_1[loc+0] * config->GetMolar_Mass(iSpecies);
					mass_Neg_2 += U_2[loc+0] * config->GetMolar_Mass(iSpecies);
					massTotal_Neg_0 += U_0[loc+0];
					massTotal_Neg_1 += U_1[loc+0];
					massTotal_Neg_2 += U_2[loc+0];
					diff_ru_Neg_x_0 += ConsVar_Grad_0[loc+1][0];
					diff_ru_Neg_x_1 += ConsVar_Grad_1[loc+1][0];
					diff_ru_Neg_x_2 += ConsVar_Grad_2[loc+1][0];
					diff_ru_Neg_y_0 += ConsVar_Grad_0[loc+2][1];
					diff_ru_Neg_y_1 += ConsVar_Grad_1[loc+2][1];
					diff_ru_Neg_y_2 += ConsVar_Grad_2[loc+2][1];
					counterNeg++;
				}
			}
			rho_Pos_0 = rho_Pos_0/counterPos;
			rho_Pos_1 = rho_Pos_1/counterPos;
			rho_Pos_2 = rho_Pos_2/counterPos;
			rho_Neg_0 = rho_Neg_0/counterNeg;
			rho_Neg_1 = rho_Neg_1/counterNeg;
			rho_Neg_2 = rho_Neg_2/counterNeg;
      
			mass_Pos_0 = mass_Pos_0 / massTotal_Pos_0;
			mass_Pos_1 = mass_Pos_1 / massTotal_Pos_1;
			mass_Pos_2 = mass_Pos_2 / massTotal_Pos_2;
			mass_Neg_0 = mass_Neg_0 / massTotal_Neg_0;
			mass_Neg_1 = mass_Neg_1 / massTotal_Neg_1;
			mass_Neg_2 = mass_Neg_2 / massTotal_Neg_2;
      
			diff_ru_Pos_x_0 = diff_ru_Pos_x_0 / counterPos;
			diff_ru_Pos_x_1 = diff_ru_Pos_x_1 / counterPos;
			diff_ru_Pos_x_2 = diff_ru_Pos_x_2 / counterPos;
			diff_ru_Pos_y_0 = diff_ru_Pos_y_0 / counterPos;
			diff_ru_Pos_y_1 = diff_ru_Pos_y_1 / counterPos;
			diff_ru_Pos_y_2 = diff_ru_Pos_y_2 / counterPos;
      
			diff_ru_Neg_x_0 = diff_ru_Neg_x_0 / counterNeg;
			diff_ru_Neg_x_1 = diff_ru_Neg_x_1 / counterNeg;
			diff_ru_Neg_x_2 = diff_ru_Neg_x_2 / counterNeg;
			diff_ru_Neg_y_0 = diff_ru_Neg_y_0 / counterNeg;
			diff_ru_Neg_y_1 = diff_ru_Neg_y_1 / counterNeg;
			diff_ru_Neg_y_2 = diff_ru_Neg_y_2 / counterNeg;
      
      
			/*--- Source q ---*/
			ne_0 = rho_Neg_0 / mass_Neg_0;
			ne_1 = rho_Neg_1 / mass_Neg_1;
			ne_2 = rho_Neg_2 / mass_Neg_2;
			np_0 = rho_Pos_0 / mass_Pos_0;
			np_1 = rho_Pos_1 / mass_Pos_1;
			np_2 = rho_Pos_2 / mass_Pos_2;
      
			Chi_n_0 = ELECTRON_CHARGE/FREE_PERMITTIVITY * (ne_0 - np_0);
			Chi_n_1 = ELECTRON_CHARGE/FREE_PERMITTIVITY * (ne_1 - np_1);
			Chi_n_2 = ELECTRON_CHARGE/FREE_PERMITTIVITY * (ne_2 - np_2);
      
			alpha_0 = pow(dt*ELECTRON_CHARGE,2)/ FREE_PERMITTIVITY * (rho_Neg_0 / (mass_Neg_0*mass_Neg_0) + rho_Pos_0 / (mass_Pos_0*mass_Pos_0));
			alpha_1 = pow(dt*ELECTRON_CHARGE,2)/ FREE_PERMITTIVITY * (rho_Neg_1 / (mass_Neg_1*mass_Neg_1) + rho_Pos_1 / (mass_Pos_1*mass_Pos_1));
			alpha_2 = pow(dt*ELECTRON_CHARGE,2)/ FREE_PERMITTIVITY * (rho_Neg_2 / (mass_Neg_2*mass_Neg_2) + rho_Pos_2 / (mass_Pos_2*mass_Pos_2));
      
			Chi_np1_0 = 1.0/(1.0+alpha_0) * (Chi_n_0 - dt*ELECTRON_CHARGE/FREE_PERMITTIVITY*((diff_ru_Neg_x_0+diff_ru_Neg_y_0)/mass_Neg_0 - (diff_ru_Pos_x_0+diff_ru_Pos_y_0)/mass_Pos_0));
			Chi_np1_1 = 1.0/(1.0+alpha_1) * (Chi_n_1 - dt*ELECTRON_CHARGE/FREE_PERMITTIVITY*((diff_ru_Neg_x_1+diff_ru_Neg_y_1)/mass_Neg_1 - (diff_ru_Pos_x_1+diff_ru_Pos_y_1)/mass_Pos_1));
			Chi_np1_2 = 1.0/(1.0+alpha_2) * (Chi_n_2 - dt*ELECTRON_CHARGE/FREE_PERMITTIVITY*((diff_ru_Neg_x_2+diff_ru_Neg_y_2)/mass_Neg_2 - (diff_ru_Pos_x_2+diff_ru_Pos_y_2)/mass_Pos_2));
      
			/*--- Residual = transpose(N) * source (q) * Area  ---*/
			val_residual[0] =-1.0/3.0 * Chi_np1_0 * Area;
			val_residual[1] =-1.0/3.0 * Chi_np1_1 * Area;
			val_residual[2] =-1.0/3.0 * Chi_np1_2 * Area;
      
      
		}
    
		if (nDim == 3) {
			val_residual[0] = 0.0;
			val_residual[1] = 0.0;
			val_residual[2] = 0.0;
		}
	}
}

void CSourcePieceWise_Elec::ComputeResidual(double *val_residual, CConfig *config) {
  
	if (config->GetKind_GasModel() == ARGON) {
		double Kai_n, ne, np;
		double rho_Pos = 0.0, rho_Elec = 0.0, mass_Pos;
		double a[4], b[4], Area = 0.0, f[4];
		mass_Pos  = config->GetParticle_Mass(1);
		if (nDim == 2) {
			for (unsigned short iDim = 0; iDim < nDim; iDim++) {
				a[iDim] = Coord_0[iDim]-Coord_2[iDim];
				b[iDim] = Coord_1[iDim]-Coord_2[iDim];
			}
			Area = 0.5*fabs(a[0]*b[1]-a[1]*b[0]);			/*--- Norm of the normal component of area, area = 1/2*cross(a,b) ---*/
      
			rho_Pos  = 1.0/3.0*(U_0[0] + U_1[0] + U_2[0]) ;
			rho_Elec = 1.0/3.0*(U_0[1] + U_1[1] + U_2[1]) ;
      
			/*--- Source q ---*/
      
			ne = rho_Elec / ELECTRON_MASS;
			np = rho_Pos  / mass_Pos;
      
			Kai_n = ELECTRON_CHARGE/FREE_PERMITTIVITY * (ne - np );
      
			/*--- Residual = transpose(N) * source (q) * Area  ---*/
			val_residual[0] =-1.0/3.0 * Kai_n * Area;
			val_residual[1] =-1.0/3.0 * Kai_n * Area;
			val_residual[2] =-1.0/3.0 * Kai_n * Area;
		}
		if (nDim == 3) {
      
			f[0] = ELECTRON_CHARGE/FREE_PERMITTIVITY * ( U_0[1]/ELECTRON_MASS - U_0[0]/mass_Pos );
			f[1] = ELECTRON_CHARGE/FREE_PERMITTIVITY * ( U_1[1]/ELECTRON_MASS - U_1[0]/mass_Pos );
			f[2] = ELECTRON_CHARGE/FREE_PERMITTIVITY * ( U_2[1]/ELECTRON_MASS - U_2[0]/mass_Pos );
			f[3] = ELECTRON_CHARGE/FREE_PERMITTIVITY * ( U_3[1]/ELECTRON_MASS - U_3[0]/mass_Pos );
      
			/*--- Residual = transpose(N) * source (q) * Area  ---*/
			for (unsigned short iVar = 0; iVar < 4; iVar ++ ) {
				val_residual[iVar] = 0.0;
				for (unsigned short jVar = 0; jVar < 4; jVar ++ )
					val_residual[iVar] -= Ni_times_Nj[iVar][jVar] * f[jVar] * Volume;
			}
		}
	}
}
