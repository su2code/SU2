/*!
 * \file numerics_adjoint_mean.cpp
 * \brief This file contains all the convective term discretization.
 * \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
 * \version 2.0.7
 *
 * Stanford University Unstructured (SU2).
 * Copyright (C) 2012-2013 Aerospace Design Laboratory (ADL).
 *
 * SU2 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * SU2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SU2. If not, see <http://www.gnu.org/licenses/>.
 */

#include "../include/numerics_structure.hpp"
#include <limits>

CUpwRoe_AdjTNE2::CUpwRoe_AdjTNE2(unsigned short val_nDim, unsigned short val_nVar,
                                 CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
	implicit = (config->GetKind_TimeIntScheme_AdjTNE2() == EULER_IMPLICIT);
	rotating_frame = config->GetRotating_Frame();
	grid_movement = config->GetGrid_Movement();
  
	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;
  
	Residual_Roe = new double [nVar];
	RoeVelocity = new double [nDim];
	Velocity_i = new double [nDim];
	Velocity_j = new double [nDim];
	Lambda = new double [nVar];
	P_Tensor = new double* [nVar];
	invP_Tensor = new double* [nVar];
	Proj_flux_tensor_i = new double*[nVar];
	Proj_flux_tensor_j = new double*[nVar];
	Proj_ModJac_Tensor = new double*[nVar];
	for (iVar = 0; iVar < nVar; iVar++) {
		P_Tensor[iVar] = new double [nVar];
		invP_Tensor[iVar] = new double [nVar];
		Proj_flux_tensor_i[iVar] = new double[nVar];
		Proj_flux_tensor_j[iVar] = new double[nVar];
		Proj_ModJac_Tensor[iVar] = new double[nVar];
	}
  
}

CUpwRoe_AdjTNE2::~CUpwRoe_AdjTNE2(void) {
  
	delete [] Residual_Roe;
	delete [] RoeVelocity;
	delete [] Velocity_i;
	delete [] Velocity_j;
	delete [] Lambda;
	for (iVar = 0; iVar < nVar; iVar++) {
		delete [] P_Tensor[iVar];
		delete [] invP_Tensor[iVar];
		delete [] Proj_flux_tensor_i[iVar];
		delete [] Proj_flux_tensor_j[iVar];
		delete [] Proj_ModJac_Tensor[iVar];
	}
	delete [] P_Tensor;
	delete [] invP_Tensor;
	delete [] Proj_flux_tensor_i;
	delete [] Proj_flux_tensor_j;
	delete [] Proj_ModJac_Tensor;
  
}

void CUpwRoe_AdjTNE2::ComputeResidual (double *val_residual_i, double *val_residual_j,
                                       double **val_Jacobian_ii, double **val_Jacobian_ij,
                                       double **val_Jacobian_ji, double **val_Jacobian_jj,
                                       CConfig *config) {
  
	/*--- Compute the area ---*/
	area = 0;
	for (iDim = 0; iDim < nDim; iDim++)
		area += Normal[iDim]*Normal[iDim];
	area = sqrt(area);
	rarea = 1.0 / area;
  
	/*--- Components of the normal & unit normal vector of the current face ---*/
	Sx = Normal[0];
	Sy = Normal[1];
	Sz = 0.0; if (nDim == 3) Sz = Normal[2];
	nx = Sx * rarea;
	ny = Sy * rarea;
	nz = Sz * rarea;
  
	/*--- Flow variable states at point i (left, _l) and j (right, _r)---*/
	rho_l  = U_i[0]; rho_r  = U_j[0];
	u_l = U_i[1]/U_i[0]; v_l = U_i[2]/U_i[0]; w_l = 0.0;
	u_r = U_j[1]/U_j[0]; v_r = U_j[2]/U_j[0]; w_r = 0.0;
	if (nDim == 3) w_l = U_i[3]/U_i[0];
	if (nDim == 3) w_r = U_j[3]/U_j[0];
	h_l = Enthalpy_i; h_r = Enthalpy_j;
  
	/*--- One-half speed squared ---*/
	q_l = ONE2 * ((u_l*u_l) + (v_l*v_l) + (w_l*w_l));
	q_r = ONE2 * ((u_r*u_r) + (v_r*v_r) + (w_r*w_r));
  
	/*--- Projected velocity ---*/
	Q_l = (u_l * Sx) + (v_l * Sy) + (w_l * Sz);
	Q_r = (u_r * Sx) + (v_r * Sy) + (w_r * Sz);
  
	/*--- Mean adjoint variables ---*/
	psi1 = ONE2 * (Psi_i[0] + Psi_j[0]);
	psi2 = ONE2 * (Psi_i[1] + Psi_j[1]);
	psi3 = ONE2 * (Psi_i[2] + Psi_j[2]);
	psi4 = 0.0; if (nDim == 3) psi4 = ONE2 * (Psi_i[3] + Psi_j[3]);
	psi5 = ONE2 * (Psi_i[nVar-1] + Psi_j[nVar-1]);
  
	/*--- Left state ---*/
	l1psi = (Sx * psi2) + (Sy * psi3) + (Sz * psi4) + (Q_l * psi5);
	l2psi = psi1 + (u_l * psi2) + (v_l * psi3) + (w_l * psi4) + (h_l * psi5);
  
	val_residual_i[0] = Q_l * psi1 - l2psi * Q_l + l1psi * Gamma_Minus_One * q_l;
	val_residual_i[1] = Q_l * psi2 + l2psi * Sx  - l1psi * Gamma_Minus_One * u_l;
	val_residual_i[2] = Q_l * psi3 + l2psi * Sy  - l1psi * Gamma_Minus_One * v_l;
	if (nDim == 3) val_residual_i[3] = Q_l * psi4 + l2psi * Sz  - l1psi * Gamma_Minus_One * w_l;
	val_residual_i[nVar-1] = Q_l * psi5 + l1psi * Gamma_Minus_One;
  
	/*--- Right state ---*/
	l1psi = (Sx * psi2) + (Sy * psi3) + (Sz * psi4) + (Q_r * psi5);
	l2psi = psi1 + (u_r * psi2) + (v_r * psi3) + (w_r * psi4) + (h_r * psi5);
  
	val_residual_j[0] = -(Q_r * psi1 - l2psi * Q_r + l1psi * Gamma_Minus_One * q_r);
	val_residual_j[1] = -(Q_r * psi2 + l2psi * Sx  - l1psi * Gamma_Minus_One * u_r);
	val_residual_j[2] = -(Q_r * psi3 + l2psi * Sy  - l1psi * Gamma_Minus_One * v_r);
	if (nDim == 3) val_residual_j[3] = -(Q_r * psi4 + l2psi * Sz  - l1psi * Gamma_Minus_One * w_r);
	val_residual_j[nVar-1] = -(Q_r * psi5 + l1psi * Gamma_Minus_One);
  
  
	/*--- f_{roe} = P^{-T} |lambda| P^T \delta \psi ---*/
	psi1_l = Psi_i[0];
	psi2_l = Psi_i[1];
	psi3_l = Psi_i[2];
	psi4_l = 0.0; if (nDim == 3) psi4_l = Psi_i[3];
	psi5_l = Psi_i[nVar-1];
  
	psi1_r = Psi_j[0];
	psi2_r = Psi_j[1];
	psi3_r = Psi_j[2];
	psi4_r = 0.0; if (nDim == 3) psi4_r = Psi_j[3];
	psi5_r = Psi_j[nVar-1];
  
	/*--- Roe averaging ---*/
	rrho_l   = 1.0 / rho_l;
	weight   = sqrt(rho_r * rrho_l);
	rweight1 = 1.0 / (1.0 + weight);
	weight  *= rweight1;
  
	h = h_l * rweight1 + weight * h_r;
	u = u_l * rweight1 + weight * u_r;
	v = v_l * rweight1 + weight * v_r;
	w = w_l * rweight1 + weight * w_r;
  
	psi1 = ONE2 * (psi1_r - psi1_l);
	psi2 = ONE2 * (psi2_r - psi2_l);
	psi3 = ONE2 * (psi3_r - psi3_l);
	psi4 = ONE2 * (psi4_r - psi4_l);
	psi5 = ONE2 * (psi5_r - psi5_l);
  
	q2 = (u*u) + (v*v) + (w*w);
	Q  = (u * Sx) + (v * Sy) + (w * Sz);
	vn = nx * u   + ny * v   + nz * w;
	cc = Gamma_Minus_One * h - 0.5 * Gamma_Minus_One * q2;
	c  = sqrt(cc);
  
  
	/*--- Contribution to velocity projection due to grid movement ---*/
	if (grid_movement) {
		double ProjGridVel = 0.0;
		for (iDim = 0; iDim < nDim; iDim++)
			ProjGridVel += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*Normal[iDim];
		Q -= ProjGridVel;
	}
  
	/*--- Eigenvalues from the primal solution ---*/
	absQ  = fabs(Q);
	absQp = fabs(Q + c * area);
	absQm = fabs(Q - c * area);
  
	alpha  = ONE2 * Gamma_Minus_One * q2 / cc;
	beta_u = psi2 + u * psi5;
	beta_v = psi3 + v * psi5;
	beta_w = psi4 + w * psi5;
	eta    = Gamma_Minus_One / cc;
	l1psi  = (nx * psi2) + (ny * psi3) + (nz * psi4) + (vn * psi5);
	l2psi  = psi1 + (u * psi2) + (v * psi3) + (w * psi4) + (h * psi5);
	l1l2p  = (l2psi + c * l1psi) * absQp;
	l1l2m  = (l2psi - c * l1psi) * absQm;
  
	/*--- adjoint flux computation in the x,y and z coordinate system ---*/
	Residual_Roe[0] = ((1.0-alpha)*l2psi - (1.0-alpha)*cc/Gamma_Minus_One*psi5
                     - u*beta_u*(1.0-(nx*nx)) - v*beta_v*(1.0-(ny*ny))
                     - w*beta_w*(1.0-(nz*nz)) + ny*nz*(w*beta_v + v*beta_w)
                     + nx*nz*(w*beta_u + u*beta_w) + ny*nx*(v*beta_u + u*beta_v) ) * absQ
  - ONE2 / c * vn * (l1l2p - l1l2m) + ONE2 * alpha *  (l1l2p + l1l2m);
  
	Residual_Roe[1] = (l2psi*u*eta - u*psi5 + beta_u*(1.0-(nx*nx))
                     - nx*(beta_v*ny + beta_w*nz) ) * absQ + ONE2*nx/c  * (l1l2p - l1l2m )
  - ONE2*eta*u * (l1l2p + l1l2m );
  
	Residual_Roe[2] = (l2psi*v*eta - v*psi5 + beta_v*(1.0-(ny*ny))
                     - ny*(beta_w*nz + beta_u*nx) ) * absQ + ONE2*ny/c  * (l1l2p - l1l2m )
  - ONE2*eta*v * (l1l2p + l1l2m );
  
	if (nDim == 3) Residual_Roe[3] = (l2psi*w*eta - w*psi5 + beta_w*(1.0-(nz*nz)) - nz*(beta_u*nx + beta_v*ny) ) * absQ
    + ONE2*nz/c  * (l1l2p - l1l2m ) - ONE2*eta*w * (l1l2p + l1l2m );
  
	Residual_Roe[nVar-1] = (psi5 - l2psi*eta) * absQ + ONE2*eta*(l1l2p + l1l2m);
  
	for (iVar = 0; iVar < nVar; iVar++) {
		val_residual_i[iVar]   += Residual_Roe[iVar];
		val_residual_j[iVar]   -= Residual_Roe[iVar];
	}
  
	/*--- Flux contribution due to grid movement ---*/
	if (grid_movement) {
		double ProjGridVel = 0.0;
		for (iDim = 0; iDim < nDim; iDim++)
			ProjGridVel += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*Normal[iDim];
		for (iVar = 0; iVar < nVar; iVar++) {
			val_residual_i[iVar] -= ProjGridVel * 0.5*(Psi_i[iVar]+Psi_j[iVar]);
			val_residual_j[iVar] += ProjGridVel * 0.5*(Psi_i[iVar]+Psi_j[iVar]);
		}
	}
  
	/*--- Implicit Contributions ---*/
	if (implicit) {
    
		/*--- Prepare variables for use in matrix routines ---*/
		RoeDensity = U_i[0]*sqrt(U_j[0]/U_i[0]);
		RoeSoundSpeed = c;
		UnitNormal[0] = nx;  UnitNormal[1] = ny;  if (nDim == 3 ) UnitNormal[2] = nz;
		RoeVelocity[0]   = u;   RoeVelocity[1]   = v;   if (nDim == 3 ) RoeVelocity[2]   = w;
		Velocity_i[0]    = u_l; Velocity_i[1]    = v_l; if (nDim == 3 ) Velocity_i[2]    = w_l;
		Velocity_j[0]    = u_r; Velocity_j[1]    = v_r; if (nDim == 3 ) Velocity_j[2]    = w_r;
		Energy_i = U_i[nDim+1] / U_i[0]; Energy_j = U_j[nDim+1] / U_j[0];
    
		/*--- Jacobians of the inviscid flux, scaled by
		 0.5 because val_resconv ~ 0.5*(fc_i+fc_j)*Normal ---*/
		GetInviscidProjJac(Velocity_i, &Energy_i, Normal, 0.5, Proj_flux_tensor_i);
		GetInviscidProjJac(Velocity_j, &Energy_j, Normal, 0.5, Proj_flux_tensor_j);
    
		/*--- Compute P, inverse P, and store eigenvalues ---*/
		GetPMatrix_inv(&RoeDensity, RoeVelocity, &RoeSoundSpeed, UnitNormal, invP_Tensor);
		GetPMatrix(&RoeDensity, RoeVelocity, &RoeSoundSpeed, UnitNormal, P_Tensor);
    
		/*--- Flow eigenvalues ---*/
		for (iDim = 0; iDim < nDim; iDim++)
			Lambda[iDim] = absQ;
		Lambda[nVar-2] = absQp;
		Lambda[nVar-1] = absQm;
    
		/*--- Roe's Flux approximation ---*/
		for (iVar = 0; iVar < nVar; iVar++) {
			for (jVar = 0; jVar < nVar; jVar++) {
				Proj_ModJac_Tensor_ij = 0.0;
				/*--- Compute |Proj_ModJac_Tensor| = P x |Lambda| x inverse P ---*/
				for (kVar = 0; kVar < nVar; kVar++)
					Proj_ModJac_Tensor_ij += P_Tensor[iVar][kVar]*Lambda[kVar]*invP_Tensor[kVar][jVar];
				Proj_ModJac_Tensor[iVar][jVar] = 0.5*Proj_ModJac_Tensor_ij*area;
			}
		}
    
		/*--- Transpose the matrices and store the Jacobians. Note the negative
		 sign for the ji and jj Jacobians bc the normal direction is flipped. ---*/
		for (iVar = 0; iVar < nVar; iVar++) {
			for (jVar = 0; jVar < nVar; jVar++) {
				val_Jacobian_ii[jVar][iVar] = Proj_flux_tensor_i[iVar][jVar] - Proj_ModJac_Tensor[iVar][jVar];
				val_Jacobian_ij[jVar][iVar] = Proj_flux_tensor_i[iVar][jVar] + Proj_ModJac_Tensor[iVar][jVar];
				val_Jacobian_ji[jVar][iVar] = -(Proj_flux_tensor_j[iVar][jVar] - Proj_ModJac_Tensor[iVar][jVar]);
				val_Jacobian_jj[jVar][iVar] = -(Proj_flux_tensor_j[iVar][jVar] + Proj_ModJac_Tensor[iVar][jVar]);
			}
		}
    
		/*--- Jacobian contribution due to grid movement ---*/
		if (grid_movement) {
			double ProjGridVel = 0.0;
			for (iDim = 0; iDim < nDim; iDim++)
				ProjGridVel += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*Normal[iDim];
			for (iVar = 0; iVar < nVar; iVar++) {
				/*--- Adjust Jacobian main diagonal ---*/
				val_Jacobian_ii[iVar][iVar] -= 0.5*ProjGridVel;
				val_Jacobian_ij[iVar][iVar] -= 0.5*ProjGridVel;
				val_Jacobian_ji[iVar][iVar] += 0.5*ProjGridVel;
				val_Jacobian_jj[iVar][iVar] += 0.5*ProjGridVel;
			}
		}
    
	}
}


CCentJST_AdjTNE2::CCentJST_AdjTNE2(unsigned short val_nDim, unsigned short val_nVar,
                                   CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;
  
	grid_movement = config->GetGrid_Movement();
	rotating_frame = config->GetRotating_Frame();
  
	Diff_Psi = new double [nVar]; Diff_Lapl = new double [nVar];
	Und_Lapl_i = new double [nVar]; Und_Lapl_j = new double [nVar];
	Velocity_i = new double [nDim]; Velocity_j = new double [nDim];
	MeanPhi = new double [nDim];
  
	Param_p = 0.3;
	Param_Kappa_2 = config->GetKappa_2nd_AdjTNE2();
	Param_Kappa_4 = config->GetKappa_4th_AdjTNE2();
	implicit = (config->GetKind_TimeIntScheme_AdjTNE2() == EULER_IMPLICIT);
}

CCentJST_AdjTNE2::~CCentJST_AdjTNE2(void) {
  
	delete [] Diff_Psi; delete [] Diff_Lapl;
	delete [] Und_Lapl_i; delete [] Und_Lapl_j;
	delete [] Velocity_i; delete [] Velocity_j;
	delete [] MeanPhi;
}

void CCentJST_AdjTNE2::ComputeResidual (double *val_resconv_i, double *val_resvisc_i,
                                        double *val_resconv_j, double *val_resvisc_j,
                                        double **val_Jacobian_ii, double **val_Jacobian_ij,
                                        double **val_Jacobian_ji, double **val_Jacobian_jj,
                                        CConfig *config) {
  
	/*--- Mean Values ---*/
	MeanPsiRho =  0.5*(Psi_i[0]+Psi_j[0]);
	for (iDim = 0; iDim < nDim; iDim++)
		MeanPhi[iDim] =  0.5*(Psi_i[iDim+1]+Psi_j[iDim+1]);
	MeanPsiE =  0.5*(Psi_i[nVar-1]+Psi_j[nVar-1]);
  
	/*--- Point i convective residual evaluation ---*/
	ProjVelocity_i = 0; ProjPhi = 0; ProjPhi_Vel = 0; sq_vel = 0; Area = 0;
	for (iDim = 0; iDim < nDim; iDim++) {
		Velocity_i[iDim] = U_i[iDim+1] / U_i[0];
		ProjVelocity_i += Velocity_i[iDim]*Normal[iDim];
		ProjPhi += MeanPhi[iDim]*Normal[iDim];
		ProjPhi_Vel += MeanPhi[iDim]*Velocity_i[iDim];
		sq_vel += 0.5*Velocity_i[iDim]*Velocity_i[iDim];
		Area += Normal[iDim]*Normal[iDim];
	}
	Area = sqrt(Area);
	phis1 = ProjPhi + ProjVelocity_i*MeanPsiE;
	phis2 = MeanPsiRho + ProjPhi_Vel + Enthalpy_i*MeanPsiE;
  
	val_resconv_i[0] = ProjVelocity_i*MeanPsiRho - phis2*ProjVelocity_i + Gamma_Minus_One*phis1*sq_vel;
	for (iDim = 0; iDim < nDim; iDim++)
		val_resconv_i[iDim+1] = ProjVelocity_i*MeanPhi[iDim] + phis2*Normal[iDim] - Gamma_Minus_One*phis1*Velocity_i[iDim];
	val_resconv_i[nVar-1] = ProjVelocity_i*MeanPsiE + Gamma_Minus_One*phis1;
  
	/*--- Flux contributions due to grid movement at point i (TDE) ---*/
	if (grid_movement) {
		double ProjGridVel = 0.0;
		for (iDim = 0; iDim < nDim; iDim++)
			ProjGridVel += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*Normal[iDim];
		val_resconv_i[0] -= ProjGridVel*MeanPsiRho;
		for (iDim = 0; iDim < nDim; iDim++)
			val_resconv_i[iDim+1] -= ProjGridVel*MeanPhi[iDim];
		val_resconv_i[nVar-1] -= ProjGridVel*MeanPsiE;
	}
  
	/*--- Jacobians of the inviscid flux ---*/
	if (implicit) {
		val_Jacobian_ii[0][0] = 0.0;
		for (jDim = 0; jDim < nDim; jDim++)
			val_Jacobian_ii[0][jDim+1] = -0.5*ProjVelocity_i*Velocity_i[jDim] + Gamma_Minus_One*sq_vel*0.5*Normal[jDim];
		val_Jacobian_ii[0][nVar-1] = 0.5*ProjVelocity_i*(Gamma_Minus_One*sq_vel - Enthalpy_i);
		for (iDim = 0; iDim < nDim; iDim++) {
			val_Jacobian_ii[iDim+1][0] = 0.5*Normal[iDim];
			for (jDim = 0; jDim < nDim; jDim++)
				val_Jacobian_ii[iDim+1][jDim+1] = 0.5*Normal[iDim]*Velocity_i[jDim] - 0.5*Gamma_Minus_One*Velocity_i[iDim]*Normal[jDim];
			val_Jacobian_ii[iDim+1][iDim+1] += 0.5*ProjVelocity_i;
			val_Jacobian_ii[iDim+1][nVar-1] = 0.5*Enthalpy_i*Normal[iDim] - 0.5*Gamma_Minus_One*Velocity_i[iDim]*ProjVelocity_i;
		}
		val_Jacobian_ii[nVar-1][0] = 0;
		for (jDim = 0; jDim < nDim; jDim++)
			val_Jacobian_ii[nVar-1][jDim+1] = 0.5*Gamma_Minus_One*Normal[jDim];
		val_Jacobian_ii[nVar-1][nVar-1] = 0.5*Gamma*ProjVelocity_i;
    
		for (iVar = 0; iVar < nVar; iVar++)
			for (jVar = 0; jVar < nVar; jVar++)
				val_Jacobian_ij[iVar][jVar] = val_Jacobian_ii[iVar][jVar];
    
		/*--- Jacobian contributions due to grid movement at point i (TDE) ---*/
		if (grid_movement) {
			double ProjGridVel = 0.0;
			for (iDim = 0; iDim < nDim; iDim++)
				ProjGridVel += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*Normal[iDim];
			for (iVar = 0; iVar < nVar; iVar++) {
				val_Jacobian_ii[iVar][iVar] -= 0.5*ProjGridVel;
				val_Jacobian_ij[iVar][iVar] -= 0.5*ProjGridVel;
			}
		}
	}
  
  
	/*--- Point j convective residual evaluation ---*/
	ProjVelocity_j = 0; ProjPhi_Vel = 0; sq_vel = 0;
	for (iDim = 0; iDim < nDim; iDim++) {
		Velocity_j[iDim] = U_j[iDim+1] / U_j[0];
		ProjVelocity_j += Velocity_j[iDim]*Normal[iDim];
		ProjPhi_Vel += MeanPhi[iDim]*Velocity_j[iDim];
		sq_vel += 0.5*Velocity_j[iDim]*Velocity_j[iDim];
	}
  
	phis1 = ProjPhi + ProjVelocity_j*MeanPsiE;
	phis2 = MeanPsiRho + ProjPhi_Vel + Enthalpy_j*MeanPsiE;
  
	val_resconv_j[0] = -(ProjVelocity_j*MeanPsiRho - phis2*ProjVelocity_j + Gamma_Minus_One*phis1*sq_vel);
	for (iDim = 0; iDim < nDim; iDim++)
		val_resconv_j[iDim+1] = -(ProjVelocity_j*MeanPhi[iDim] + phis2*Normal[iDim] - Gamma_Minus_One*phis1*Velocity_j[iDim]);
	val_resconv_j[nVar-1] = -(ProjVelocity_j*MeanPsiE + Gamma_Minus_One*phis1);
  
	/*--- Flux contributions due to grid movement at point j (TDE) ---*/
	if (grid_movement) {
		double ProjGridVel = 0.0;
		for (iDim = 0; iDim < nDim; iDim++)
			ProjGridVel += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*Normal[iDim];
		val_resconv_j[0] += ProjGridVel*MeanPsiRho;
		for (iDim = 0; iDim < nDim; iDim++)
			val_resconv_j[iDim+1] += ProjGridVel*MeanPhi[iDim];
		val_resconv_j[nVar-1] += ProjGridVel*MeanPsiE;
	}
  
	/*--- Jacobians of the inviscid flux ---*/
	if (implicit) {
		val_Jacobian_jj[0][0] = 0.0;
		for (jDim = 0; jDim < nDim; jDim++)
			val_Jacobian_jj[0][jDim+1] = 0.5*ProjVelocity_j*Velocity_j[jDim] - Gamma_Minus_One*sq_vel*0.5*Normal[jDim];
		val_Jacobian_jj[0][nVar-1] = -0.5*ProjVelocity_j*(Gamma_Minus_One*sq_vel - Enthalpy_j);
		for (iDim = 0; iDim < nDim; iDim++) {
			val_Jacobian_jj[iDim+1][0] = -0.5*Normal[iDim];
			for (jDim = 0; jDim < nDim; jDim++)
				val_Jacobian_jj[iDim+1][jDim+1] = -0.5*Normal[iDim]*Velocity_j[jDim] + 0.5*Gamma_Minus_One*Velocity_j[iDim]*Normal[jDim];
			val_Jacobian_jj[iDim+1][iDim+1] -= 0.5*ProjVelocity_j;
			val_Jacobian_jj[iDim+1][nVar-1] = -0.5*Enthalpy_j*Normal[iDim] + 0.5*Gamma_Minus_One*Velocity_j[iDim]*ProjVelocity_j;
		}
		val_Jacobian_jj[nVar-1][0] = 0;
		for (jDim = 0; jDim < nDim; jDim++)
			val_Jacobian_jj[nVar-1][jDim+1] = -0.5*Gamma_Minus_One*Normal[jDim];
		val_Jacobian_jj[nVar-1][nVar-1] = -0.5*Gamma*ProjVelocity_j;
    
		for (iVar = 0; iVar < nVar; iVar++)
			for (jVar = 0; jVar < nVar; jVar++)
				val_Jacobian_ji[iVar][jVar] = val_Jacobian_jj[iVar][jVar];
    
		/*--- Jacobian contributions due to grid movement at point j (TDE) ---*/
		if (grid_movement) {
			double ProjGridVel = 0.0;
			for (iDim = 0; iDim < nDim; iDim++)
				ProjGridVel += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*Normal[iDim];
			for (iVar = 0; iVar < nVar; iVar++) {
				val_Jacobian_jj[iVar][iVar] += 0.5*ProjGridVel;
				val_Jacobian_ji[iVar][iVar] += 0.5*ProjGridVel;
			}
		}
	}
  
	/*--- Computes differences btw. variables and Laplacians ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		Diff_Lapl[iVar] = Und_Lapl_i[iVar]-Und_Lapl_j[iVar];
		Diff_Psi[iVar]  = Psi_i[iVar]-Psi_j[iVar];
	}
  
	/*--- Adjustment to projected velocity due to mesh motion (TDE) ---*/
	if (grid_movement) {
		double ProjGridVel_i = 0.0; double ProjGridVel_j = 0.0; double ProjGridVel = 0.0;
		for (iDim = 0; iDim < nDim; iDim++) {
			ProjGridVel += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*Normal[iDim];
			ProjGridVel_i += GridVel_i[iDim]*Normal[iDim];
			ProjGridVel_j += GridVel_j[iDim]*Normal[iDim];
		}
		ProjVelocity_i -= ProjGridVel;
		ProjVelocity_j += ProjGridVel;
	}
  
	/*--- Compute the spectral radius and stretching factor ---*/
	Local_Lambda_i = (fabs(ProjVelocity_i)+SoundSpeed_i*Area);
	Local_Lambda_j = (fabs(ProjVelocity_j)+SoundSpeed_j*Area);
	MeanLambda = 0.5*(Local_Lambda_i+Local_Lambda_j);
  
	Phi_i = pow(Lambda_i/(4.0*MeanLambda+EPS),Param_p);
	Phi_j = pow(Lambda_j/(4.0*MeanLambda+EPS),Param_p);
	StretchingFactor = 4.0*Phi_i*Phi_j/(Phi_i+Phi_j+EPS);
  
	double sc2 = 3.0*(double(Neighbor_i)+double(Neighbor_j))/(double(Neighbor_i)*double(Neighbor_j));
	sc4 = sc2*sc2/4.0;
	Epsilon_2 = Param_Kappa_2*0.5*(Sensor_i+Sensor_j)*sc2;
	Epsilon_4 = max(0.0, Param_Kappa_4-Epsilon_2)*sc4;
  
	/*--- Compute viscous residual 1st- & 3rd-order dissipation ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		Residual = (Epsilon_2*Diff_Psi[iVar]-Epsilon_4*Diff_Lapl[iVar])*StretchingFactor*MeanLambda;
		val_resvisc_i[iVar] = -Residual;
		val_resvisc_j[iVar] =  Residual;
		if (implicit) {
			val_Jacobian_ii[iVar][iVar] -= Epsilon_2 + double(Neighbor_i+1)*Epsilon_4*StretchingFactor*MeanLambda;
			val_Jacobian_ij[iVar][iVar] += Epsilon_2 + double(Neighbor_j+1)*Epsilon_4*StretchingFactor*MeanLambda;
			val_Jacobian_ji[iVar][iVar] += Epsilon_2 + double(Neighbor_i+1)*Epsilon_4*StretchingFactor*MeanLambda;
			val_Jacobian_jj[iVar][iVar] -= Epsilon_2 + double(Neighbor_j+1)*Epsilon_4*StretchingFactor*MeanLambda;
		}
	}
}


CCentLax_AdjTNE2::CCentLax_AdjTNE2(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
	Gamma = config->GetGamma();
	Gamma_Minus_One = Gamma - 1.0;
  
	Diff_Psi = new double [nVar]; 	MeanPhi = new double [nDim];
	Velocity_i = new double [nDim]; Velocity_j = new double [nDim];
  
	implicit = (config->GetKind_TimeIntScheme_AdjTNE2() == EULER_IMPLICIT);
	rotating_frame = config->GetRotating_Frame();
	grid_movement = config->GetGrid_Movement();
  
	Param_p = 0.3;
	Param_Kappa_0 = config->GetKappa_1st_AdjTNE2();
  
}

CCentLax_AdjTNE2::~CCentLax_AdjTNE2(void) {
  
	delete [] Diff_Psi; delete [] MeanPhi;
	delete [] Velocity_i; delete [] Velocity_j;
  
}

void CCentLax_AdjTNE2::ComputeResidual (double *val_resconv_i, double *val_resvisc_i, double *val_resconv_j, double *val_resvisc_j,
                                        double **val_Jacobian_ii, double **val_Jacobian_ij, double **val_Jacobian_ji, double **val_Jacobian_jj,
                                        CConfig *config) {
  
	/*--- Mean value of the adjoint variables ---*/
	MeanPsiRho =  0.5*(Psi_i[0]+Psi_j[0]);
	for (iDim = 0; iDim < nDim; iDim++)
		MeanPhi[iDim] =  0.5*(Psi_i[iDim+1]+Psi_j[iDim+1]);
	MeanPsiE =  0.5*(Psi_i[nVar-1]+Psi_j[nVar-1]);
  
	/*--- Evaluation at point i ---*/
	ProjVelocity_i = 0; ProjPhi = 0; ProjPhi_Vel = 0; sq_vel = 0; Area = 0;
	for (iDim = 0; iDim < nDim; iDim++) {
		Velocity_i[iDim] = U_i[iDim+1] / U_i[0];
		ProjVelocity_i += Velocity_i[iDim]*Normal[iDim];
		ProjPhi += MeanPhi[iDim]*Normal[iDim];
		ProjPhi_Vel += MeanPhi[iDim]*Velocity_i[iDim];
		sq_vel += 0.5*Velocity_i[iDim]*Velocity_i[iDim];
		Area += Normal[iDim]*Normal[iDim];
	}
	Area = sqrt(Area);
	phis1 = ProjPhi + ProjVelocity_i*MeanPsiE;
	phis2 = MeanPsiRho + ProjPhi_Vel + Enthalpy_i*MeanPsiE;
  
	/*--- Compute inviscid residual at point i ---*/
	val_resconv_i[0] = ProjVelocity_i*MeanPsiRho - phis2*ProjVelocity_i + Gamma_Minus_One*phis1*sq_vel;
	for (iDim = 0; iDim < nDim; iDim++)
		val_resconv_i[iDim+1] = ProjVelocity_i*MeanPhi[iDim] + phis2*Normal[iDim] - Gamma_Minus_One*phis1*Velocity_i[iDim];
	val_resconv_i[nVar-1] = ProjVelocity_i*MeanPsiE + Gamma_Minus_One*phis1;
  
	/*--- Flux contributions due to grid movement at point i (TDE) ---*/
	if (grid_movement) {
		double ProjGridVel = 0.0;
		for (iDim = 0; iDim < nDim; iDim++)
			ProjGridVel += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*Normal[iDim];
		val_resconv_i[0] -= ProjGridVel*MeanPsiRho;
		for (iDim = 0; iDim < nDim; iDim++)
			val_resconv_i[iDim+1] -= ProjGridVel*MeanPhi[iDim];
		val_resconv_i[nVar-1] -= ProjGridVel*MeanPsiE;
	}
  
	/*--- Inviscid contribution to the implicit part ---*/
	if (implicit) {
		val_Jacobian_ii[0][0] = 0.0;
		for (jDim = 0; jDim < nDim; jDim++)
			val_Jacobian_ii[0][jDim+1] = -0.5*ProjVelocity_i*Velocity_i[jDim] + Gamma_Minus_One*sq_vel*0.5*Normal[jDim];
		val_Jacobian_ii[0][nVar-1] = 0.5*ProjVelocity_i*(Gamma_Minus_One*sq_vel - Enthalpy_i);
		for (iDim = 0; iDim < nDim; iDim++) {
			val_Jacobian_ii[iDim+1][0] = 0.5*Normal[iDim];
			for (jDim = 0; jDim < nDim; jDim++)
				val_Jacobian_ii[iDim+1][jDim+1] = 0.5*Normal[iDim]*Velocity_i[jDim] - 0.5*Gamma_Minus_One*Velocity_i[iDim]*Normal[jDim];
			val_Jacobian_ii[iDim+1][iDim+1] += 0.5*ProjVelocity_i;
			val_Jacobian_ii[iDim+1][nVar-1] = 0.5*Enthalpy_i*Normal[iDim] - 0.5*Gamma_Minus_One*Velocity_i[iDim]*ProjVelocity_i;
		}
		val_Jacobian_ii[nVar-1][0] = 0;
		for (jDim = 0; jDim < nDim; jDim++)
			val_Jacobian_ii[nVar-1][jDim+1] = 0.5*Gamma_Minus_One*Normal[jDim];
		val_Jacobian_ii[nVar-1][nVar-1] = 0.5*Gamma*ProjVelocity_i;
    
		for (iVar = 0; iVar < nVar; iVar++)
			for (jVar = 0; jVar < nVar; jVar++)
				val_Jacobian_ij[iVar][jVar] = val_Jacobian_ii[iVar][jVar];
    
		/*--- Jacobian contributions due to grid movement at point i (TDE) ---*/
		if (grid_movement) {
			double ProjGridVel = 0.0;
			for (iDim = 0; iDim < nDim; iDim++)
				ProjGridVel += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*Normal[iDim];
			for (iVar = 0; iVar < nVar; iVar++) {
				val_Jacobian_ii[iVar][iVar] -= 0.5*ProjGridVel;
				val_Jacobian_ij[iVar][iVar] -= 0.5*ProjGridVel;
			}
		}
	}
  
	/*--- Evaluation at point j ---*/
	ProjVelocity_j = 0; ProjPhi_Vel = 0; sq_vel = 0;
	for (iDim = 0; iDim < nDim; iDim++) {
		Velocity_j[iDim] = U_j[iDim+1] / U_j[0];
		ProjVelocity_j += Velocity_j[iDim]*Normal[iDim];
		ProjPhi_Vel += MeanPhi[iDim]*Velocity_j[iDim];
		sq_vel += 0.5*Velocity_j[iDim]*Velocity_j[iDim];
	}
  
	phis1 = ProjPhi + ProjVelocity_j*MeanPsiE;
	phis2 = MeanPsiRho + ProjPhi_Vel + Enthalpy_j*MeanPsiE;
  
	/*--- Compute inviscid residual at point j ---*/
	val_resconv_j[0] = -(ProjVelocity_j*MeanPsiRho - phis2*ProjVelocity_j + Gamma_Minus_One*phis1*sq_vel);
	for (iDim = 0; iDim < nDim; iDim++)
		val_resconv_j[iDim+1] = -(ProjVelocity_j*MeanPhi[iDim] + phis2*Normal[iDim] - Gamma_Minus_One*phis1*Velocity_j[iDim]);
	val_resconv_j[nVar-1] = -(ProjVelocity_j*MeanPsiE + Gamma_Minus_One*phis1);
  
	/*--- Flux contributions due to grid movement at point j (TDE) ---*/
	if (grid_movement) {
		double ProjGridVel = 0.0;
		for (iDim = 0; iDim < nDim; iDim++)
			ProjGridVel += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*Normal[iDim];
		val_resconv_j[0] += ProjGridVel*MeanPsiRho;
		for (iDim = 0; iDim < nDim; iDim++)
			val_resconv_j[iDim+1] += ProjGridVel*MeanPhi[iDim];
		val_resconv_j[nVar-1] += ProjGridVel*MeanPsiE;
	}
  
	/*--- Inviscid contribution to the implicit part ---*/
	if (implicit) {
		val_Jacobian_jj[0][0] = 0.0;
		for (jDim = 0; jDim < nDim; jDim++)
			val_Jacobian_jj[0][jDim+1] = 0.5*ProjVelocity_j*Velocity_j[jDim] - Gamma_Minus_One*sq_vel*0.5*Normal[jDim];
		val_Jacobian_jj[0][nVar-1] = -0.5*ProjVelocity_j*(Gamma_Minus_One*sq_vel - Enthalpy_j);
		for (iDim = 0; iDim < nDim; iDim++) {
			val_Jacobian_jj[iDim+1][0] = -0.5*Normal[iDim];
			for (jDim = 0; jDim < nDim; jDim++)
				val_Jacobian_jj[iDim+1][jDim+1] = -0.5*Normal[iDim]*Velocity_j[jDim] + 0.5*Gamma_Minus_One*Velocity_j[iDim]*Normal[jDim];
			val_Jacobian_jj[iDim+1][iDim+1] -= 0.5*ProjVelocity_j;
			val_Jacobian_jj[iDim+1][nVar-1] = -0.5*Enthalpy_j*Normal[iDim] + 0.5*Gamma_Minus_One*Velocity_j[iDim]*ProjVelocity_j;
		}
		val_Jacobian_jj[nVar-1][0] = 0;
		for (jDim = 0; jDim < nDim; jDim++)
			val_Jacobian_jj[nVar-1][jDim+1] = -0.5*Gamma_Minus_One*Normal[jDim];
		val_Jacobian_jj[nVar-1][nVar-1] = -0.5*Gamma*ProjVelocity_j;
    
		for (iVar = 0; iVar < nVar; iVar++)
			for (jVar = 0; jVar < nVar; jVar++)
				val_Jacobian_ji[iVar][jVar] = val_Jacobian_jj[iVar][jVar];
    
		/*--- Jacobian contributions due to grid movement at point j (TDE) ---*/
		if (grid_movement) {
			double ProjGridVel = 0.0;
			for (iDim = 0; iDim < nDim; iDim++)
				ProjGridVel += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*Normal[iDim];
			for (iVar = 0; iVar < nVar; iVar++) {
				val_Jacobian_jj[iVar][iVar] += 0.5*ProjGridVel;
				val_Jacobian_ji[iVar][iVar] += 0.5*ProjGridVel;
			}
		}
	}
  
	/*--- Computes differences btw. variables ---*/
	for (iVar = 0; iVar < nVar; iVar++)
		Diff_Psi[iVar] = Psi_i[iVar]-Psi_j[iVar];
  
	/*--- Adjustment to projected velocity due to mesh motion (TDE) ---*/
	if (grid_movement) {
		double ProjGridVel_i = 0.0; double ProjGridVel_j = 0.0; double ProjGridVel = 0.0;
		for (iDim = 0; iDim < nDim; iDim++) {
			ProjGridVel += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*Normal[iDim];
			ProjGridVel_i += GridVel_i[iDim]*Normal[iDim];
			ProjGridVel_j += GridVel_j[iDim]*Normal[iDim];
		}
		ProjVelocity_i -= ProjGridVel;
		ProjVelocity_j += ProjGridVel;
	}
  
	/*--- Compute spectral radius ---*/
	Local_Lambda_i = (fabs(ProjVelocity_i)+SoundSpeed_i*Area);
	Local_Lambda_j = (fabs(ProjVelocity_j)+SoundSpeed_j*Area);
	MeanLambda = 0.5*(Local_Lambda_i+Local_Lambda_j);
  
	/*--- Compute streching factor ---*/
	Phi_i = pow(Lambda_i/(4.0*MeanLambda+EPS),Param_p);
	Phi_j = pow(Lambda_j/(4.0*MeanLambda+EPS),Param_p);
	StretchingFactor = 4.0*Phi_i*Phi_j/(Phi_i+Phi_j+EPS);
  
	sc2 = 3.0*(double(Neighbor_i)+double(Neighbor_j))/(double(Neighbor_i)*double(Neighbor_j));
	Epsilon_0 = Param_Kappa_0*sc2*double(nDim)/3.0;
  
	/*--- Artifical dissipation evaluation ---*/
	for (iVar = 0; iVar < nVar; iVar++) {
		Residual = Epsilon_0*StretchingFactor*MeanLambda*Diff_Psi[iVar];
		val_resvisc_i[iVar] = -Residual;
		val_resvisc_j[iVar] =  Residual;
	}
  
	/*--- Contribution to implicit part ---*/
	if (implicit) {
		for (iVar = 0; iVar < nVar; iVar++) {
			val_Jacobian_ii[iVar][iVar] -= Epsilon_0*StretchingFactor*MeanLambda;
			val_Jacobian_ij[iVar][iVar] += Epsilon_0*StretchingFactor*MeanLambda;
			val_Jacobian_ji[iVar][iVar] += Epsilon_0*StretchingFactor*MeanLambda;
			val_Jacobian_jj[iVar][iVar] -= Epsilon_0*StretchingFactor*MeanLambda;
		}
	}
  
}
