/*!
 * \file numerics_adjoint_mean.cpp
 * \brief This file contains the numerical methods for adjoint compressible flow.
 * \author F. Palacios, T. Economon
 * \version 6.2.0 "Falcon"
 *
 * The current SU2 release has been coordinated by the
 * SU2 International Developers Society <www.su2devsociety.org>
 * with selected contributions from the open-source community.
 *
 * The main research teams contributing to the current release are:
 *  - Prof. Juan J. Alonso's group at Stanford University.
 *  - Prof. Piero Colonna's group at Delft University of Technology.
 *  - Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *  - Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *  - Prof. Rafael Palacios' group at Imperial College London.
 *  - Prof. Vincent Terrapon's group at the University of Liege.
 *  - Prof. Edwin van der Weide's group at the University of Twente.
 *  - Lab. of New Concepts in Aeronautics at Tech. Institute of Aeronautics.
 *
 * Copyright 2012-2019, Francisco D. Palacios, Thomas D. Economon,
 *                      Tim Albring, and the SU2 contributors.
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

CUpwRoe_AdjFlow::CUpwRoe_AdjFlow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  implicit = (config->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);
  grid_movement = config->GetGrid_Movement();
  
  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;
  
  Residual_Roe = new su2double [nVar];
  RoeVelocity = new su2double [nDim];
  Velocity_i = new su2double [nDim];
  Velocity_j = new su2double [nDim];
  Lambda = new su2double [nVar];
  P_Tensor = new su2double* [nVar];
  invP_Tensor = new su2double* [nVar];
  ProjFlux_i = new su2double*[nVar];
  ProjFlux_j = new su2double*[nVar];
  Proj_ModJac_Tensor = new su2double*[nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    P_Tensor[iVar] = new su2double [nVar];
    invP_Tensor[iVar] = new su2double [nVar];
    ProjFlux_i[iVar] = new su2double[nVar];
    ProjFlux_j[iVar] = new su2double[nVar];
    Proj_ModJac_Tensor[iVar] = new su2double[nVar];
  }
  
}

CUpwRoe_AdjFlow::~CUpwRoe_AdjFlow(void) {
  
  delete [] Residual_Roe;
  delete [] RoeVelocity;
  delete [] Velocity_i;
  delete [] Velocity_j;
  delete [] Lambda;
  for (iVar = 0; iVar < nVar; iVar++) {
    delete [] P_Tensor[iVar];
    delete [] invP_Tensor[iVar];
    delete [] ProjFlux_i[iVar];
    delete [] ProjFlux_j[iVar];
    delete [] Proj_ModJac_Tensor[iVar];
  }
  delete [] P_Tensor;
  delete [] invP_Tensor;
  delete [] ProjFlux_i;
  delete [] ProjFlux_j;
  delete [] Proj_ModJac_Tensor;
  
}

void CUpwRoe_AdjFlow::ComputeResidual (su2double *val_residual_i, su2double *val_residual_j, su2double **val_Jacobian_ii,
                                       su2double **val_Jacobian_ij, su2double **val_Jacobian_ji, su2double **val_Jacobian_jj, CConfig *config) {
  
  /*--- Compute the area ---*/
  
  area = 0.0;
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
  
  rho_l  = V_i[nDim+2]; rho_r  = V_j[nDim+2];
  u_l = V_i[1]; v_l = V_i[2]; w_l = 0.0; if (nDim == 3) w_l = V_i[3];
  u_r = V_j[1]; v_r = V_j[2]; w_r = 0.0; if (nDim == 3) w_r = V_j[3];
  h_l = V_i[nDim+3]; h_r = V_j[nDim+3];
  
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
    su2double ProjGridVel = 0.0;
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
  
  /*--- adjoint flux computation in the x, y and z coordinate system ---*/
  
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
    su2double ProjGridVel = 0.0;
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
    
    RoeDensity = V_i[nDim+2]*sqrt(V_j[nDim+2]/V_i[nDim+2]);
    RoeSoundSpeed = c;
    UnitNormal[0] = nx;  UnitNormal[1] = ny;  if (nDim == 3 ) UnitNormal[2] = nz;
    RoeVelocity[0]   = u;   RoeVelocity[1]   = v;   if (nDim == 3 ) RoeVelocity[2]   = w;
    Velocity_i[0]    = u_l; Velocity_i[1]    = v_l; if (nDim == 3 ) Velocity_i[2]    = w_l;
    Velocity_j[0]    = u_r; Velocity_j[1]    = v_r; if (nDim == 3 ) Velocity_j[2]    = w_r;
    
    Pressure_i = V_i[nDim+1];
    Density_i = V_i[nDim+2];
    Enthalpy_i = V_i[nDim+3];
    Energy_i = Enthalpy_i - Pressure_i/Density_i;

    Pressure_j = V_i[nDim+1];
    Density_j = V_i[nDim+2];
    Enthalpy_j = V_i[nDim+3];
    Energy_j = Enthalpy_j - Pressure_j/Density_j;
    
    /*--- Jacobians of the inviscid flux, scaled by
     0.5 because val_resconv ~ 0.5*(fc_i+fc_j)*Normal ---*/
    
    GetInviscidProjJac(Velocity_i, &Energy_i, Normal, 0.5, ProjFlux_i);
    GetInviscidProjJac(Velocity_j, &Energy_j, Normal, 0.5, ProjFlux_j);
    
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
        val_Jacobian_ii[jVar][iVar] = ProjFlux_i[iVar][jVar] - Proj_ModJac_Tensor[iVar][jVar];
        val_Jacobian_ij[jVar][iVar] = ProjFlux_i[iVar][jVar] + Proj_ModJac_Tensor[iVar][jVar];
        val_Jacobian_ji[jVar][iVar] = -(ProjFlux_j[iVar][jVar] - Proj_ModJac_Tensor[iVar][jVar]);
        val_Jacobian_jj[jVar][iVar] = -(ProjFlux_j[iVar][jVar] + Proj_ModJac_Tensor[iVar][jVar]);
      }
    }
    
    /*--- Jacobian contribution due to grid movement ---*/
    
    if (grid_movement) {
      su2double ProjGridVel = 0.0;
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

CCentJST_AdjFlow::CCentJST_AdjFlow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;
  
  grid_movement = config->GetGrid_Movement();
  
  Diff_Psi = new su2double [nVar]; Diff_Lapl = new su2double [nVar];
  Velocity_i = new su2double [nDim]; Velocity_j = new su2double [nDim];
  MeanPhi = new su2double [nDim];
  
  Param_p = 0.3;
  Param_Kappa_2 = config->GetKappa_2nd_AdjFlow();
  Param_Kappa_4 = config->GetKappa_4th_AdjFlow();
  implicit = (config->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);
  
}

CCentJST_AdjFlow::~CCentJST_AdjFlow(void) {
  
  delete [] Diff_Psi; delete [] Diff_Lapl;
  delete [] Velocity_i; delete [] Velocity_j;
  delete [] MeanPhi;
}

void CCentJST_AdjFlow::ComputeResidual (su2double *val_resconv_i, su2double *val_resvisc_i, su2double *val_resconv_j, su2double *val_resvisc_j,
                                    su2double **val_Jacobian_ii, su2double **val_Jacobian_ij, su2double **val_Jacobian_ji, su2double **val_Jacobian_jj,
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
  
  /*--- Flux contributions due to grid movement at point i ---*/
  
  if (grid_movement) {
    su2double ProjGridVel = 0.0;
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
    
    /*--- Jacobian contributions due to grid movement at point i ---*/
    
    if (grid_movement) {
      su2double ProjGridVel = 0.0;
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
  
  /*--- Flux contributions due to grid motion at point j ---*/
  
  if (grid_movement) {
    su2double ProjGridVel = 0.0;
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
    
    /*--- Jacobian contributions due to grid motion at point j ---*/
    
    if (grid_movement) {
      su2double ProjGridVel = 0.0;
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

  /*--- Adjustment to projected velocity due to grid motion ---*/
  
  if (grid_movement) {
    su2double ProjGridVel = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      ProjGridVel += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*Normal[iDim];
    }
    ProjVelocity_i -= ProjGridVel;
    ProjVelocity_j += ProjGridVel;
  }
  
  /*--- Compute the spectral radius and stretching factor ---*/
  
  Local_Lambda_i = (fabs(ProjVelocity_i)+SoundSpeed_i*Area);
  Local_Lambda_j = (fabs(ProjVelocity_j)+SoundSpeed_j*Area);
  MeanLambda = 0.5*(Local_Lambda_i+Local_Lambda_j);
  
  Phi_i = pow(Lambda_i/(4.0*MeanLambda), Param_p);
  Phi_j = pow(Lambda_j/(4.0*MeanLambda), Param_p);
  StretchingFactor = 4.0*Phi_i*Phi_j/(Phi_i+Phi_j);
  
  su2double sc2 = 3.0*(su2double(Neighbor_i)+su2double(Neighbor_j))/(su2double(Neighbor_i)*su2double(Neighbor_j));
  sc4 = sc2*sc2/4.0;
  Epsilon_2 = Param_Kappa_2*0.5*(Sensor_i+Sensor_j)*sc2;
  Epsilon_4 = max(0.0, Param_Kappa_4-Epsilon_2)*sc4;
    
  /*--- Compute viscous residual 1st- & 3rd-order dissipation ---*/
  
  for (iVar = 0; iVar < nVar; iVar++) {
    Residual = (Epsilon_2*Diff_Psi[iVar]-Epsilon_4*Diff_Lapl[iVar])*StretchingFactor*MeanLambda;
    val_resvisc_i[iVar] = -Residual;
    val_resvisc_j[iVar] =  Residual;
    if (implicit) {
      val_Jacobian_ii[iVar][iVar] -= Epsilon_2 + su2double(Neighbor_i+1)*Epsilon_4*StretchingFactor*MeanLambda;
      val_Jacobian_ij[iVar][iVar] += Epsilon_2 + su2double(Neighbor_j+1)*Epsilon_4*StretchingFactor*MeanLambda;
      val_Jacobian_ji[iVar][iVar] += Epsilon_2 + su2double(Neighbor_i+1)*Epsilon_4*StretchingFactor*MeanLambda;
      val_Jacobian_jj[iVar][iVar] -= Epsilon_2 + su2double(Neighbor_j+1)*Epsilon_4*StretchingFactor*MeanLambda;
    }
  }
  
}

CCentLax_AdjFlow::CCentLax_AdjFlow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;
  
  Diff_Psi = new su2double [nVar];   MeanPhi = new su2double [nDim];
  Velocity_i = new su2double [nDim]; Velocity_j = new su2double [nDim];
  
  implicit = (config->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);

  grid_movement = config->GetGrid_Movement();
  
  Param_p = 0.3;
  Param_Kappa_0 = config->GetKappa_1st_AdjFlow();
  
}

CCentLax_AdjFlow::~CCentLax_AdjFlow(void) {
  
  delete [] Diff_Psi; delete [] MeanPhi;
  delete [] Velocity_i; delete [] Velocity_j;
  
}

void CCentLax_AdjFlow::ComputeResidual (su2double *val_resconv_i, su2double *val_resvisc_i, su2double *val_resconv_j, su2double *val_resvisc_j,
                                    su2double **val_Jacobian_ii, su2double **val_Jacobian_ij, su2double **val_Jacobian_ji, su2double **val_Jacobian_jj,
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

  /*--- Flux contributions due to grid motion at point i ---*/
  if (grid_movement) {
    su2double ProjGridVel = 0.0;
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

    /*--- Jacobian contributions due to grid motion at point i ---*/
    if (grid_movement) {
      su2double ProjGridVel = 0.0;
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
  
  /*--- Flux contributions due to grid movement at point j ---*/
  if (grid_movement) {
    su2double ProjGridVel = 0.0;
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
    
    /*--- Jacobian contributions due to grid movement at point j ---*/
    if (grid_movement) {
      su2double ProjGridVel = 0.0;
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
  
  /*--- Adjustment to projected velocity due to grid motion ---*/
  if (grid_movement) {
    su2double ProjGridVel = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      ProjGridVel += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*Normal[iDim];
    ProjVelocity_i -= ProjGridVel;
    ProjVelocity_j += ProjGridVel;
  }
  
  /*--- Compute spectral radius ---*/
  Local_Lambda_i = (fabs(ProjVelocity_i)+SoundSpeed_i*Area);
  Local_Lambda_j = (fabs(ProjVelocity_j)+SoundSpeed_j*Area);
  MeanLambda = 0.5*(Local_Lambda_i+Local_Lambda_j);
  
  /*--- Compute streching factor ---*/
  Phi_i = pow(Lambda_i/(4.0*MeanLambda), Param_p);
  Phi_j = pow(Lambda_j/(4.0*MeanLambda), Param_p);
  StretchingFactor = 4.0*Phi_i*Phi_j/(Phi_i+Phi_j);
  
  sc2 = 3.0*(su2double(Neighbor_i)+su2double(Neighbor_j))/(su2double(Neighbor_i)*su2double(Neighbor_j));
  Epsilon_0 = Param_Kappa_0*sc2*su2double(nDim)/3.0;
  
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

CAvgGrad_AdjFlow::CAvgGrad_AdjFlow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  unsigned short iDim;
  
  implicit = (config->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);
  
  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;
  
  Velocity_i = new su2double [nDim];
  Velocity_j = new su2double [nDim];
  Mean_Velocity = new su2double [nDim];
  Mean_GradPhi = new su2double* [nDim];
  for (iDim = 0; iDim < nDim; iDim++)
    Mean_GradPhi[iDim] = new su2double [nDim];
  Mean_GradPsiE = new su2double [nDim];
  Edge_Vector = new su2double [nDim];
  
}

CAvgGrad_AdjFlow::~CAvgGrad_AdjFlow(void) {
  delete [] Velocity_i;
  delete [] Velocity_j;
  delete [] Mean_Velocity;
  delete [] Edge_Vector;
  delete [] Mean_GradPsiE;
  for (unsigned short iDim = 0; iDim < nDim; iDim++)
    delete [] Mean_GradPhi[iDim];
}

void CAvgGrad_AdjFlow::ComputeResidual(su2double *val_residual_i, su2double *val_residual_j,
                                   su2double **val_Jacobian_ii, su2double **val_Jacobian_ij,
                                   su2double **val_Jacobian_ji, su2double **val_Jacobian_jj, CConfig *config) {
  unsigned short iDim, jDim;
  su2double sq_vel_i, ViscDens_i, XiDens_i;
  su2double sq_vel_j, ViscDens_j, XiDens_j;
  su2double dist_ij_2, dPhiE_dn;
  
  su2double Prandtl_Lam      = config->GetPrandtl_Lam();
  su2double Prandtl_Turb     = config->GetPrandtl_Turb();

  /*--- States in point i ---*/
  
  sq_vel_i = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Velocity_i[iDim] = V_i[iDim+1];
    sq_vel_i += 0.5*Velocity_i[iDim]*Velocity_i[iDim];
  }
  Pressure_i = V_i[nDim+1];
  Density_i = V_i[nDim+2];
  Enthalpy_i = V_i[nDim+3];
  SoundSpeed_i = sqrt(fabs(Pressure_i*Gamma/Density_i));
  
  /*--- Laminar and Eddy viscosity ---*/
  
  Laminar_Viscosity_i = V_i[nDim+5];
  Eddy_Viscosity_i = V_i[nDim+6];
  
  ViscDens_i = (Laminar_Viscosity_i + Eddy_Viscosity_i) / Density_i;
  XiDens_i = Gamma*(Laminar_Viscosity_i/Prandtl_Lam + Eddy_Viscosity_i/Prandtl_Turb) / Density_i;
  
  /*--- States in point j ---*/
  
  sq_vel_j = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Velocity_j[iDim] = V_j[iDim+1];
    sq_vel_j += 0.5*Velocity_j[iDim]*Velocity_j[iDim];
  }
  Pressure_j = V_j[nDim+1];
  Density_j = V_j[nDim+2];
  Enthalpy_j = V_j[nDim+3];
  SoundSpeed_j = sqrt(fabs(Pressure_j*Gamma/Density_j));
  
  /*--- Laminar and Eddy viscosity ---*/
  
  Laminar_Viscosity_j = V_j[nDim+5];
  Eddy_Viscosity_j = V_j[nDim+6];
  
  ViscDens_j = (Laminar_Viscosity_j + Eddy_Viscosity_j) / Density_j;
  XiDens_j = Gamma*(Laminar_Viscosity_j/Prandtl_Lam + Eddy_Viscosity_j/Prandtl_Turb) / Density_j;
  
  /*--- Compute vector going from iPoint to jPoint ---*/
  
  dist_ij_2 = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Edge_Vector[iDim] = Coord_j[iDim]-Coord_i[iDim];
    dist_ij_2 += Edge_Vector[iDim]*Edge_Vector[iDim];
  }
  
  /*--- Average of the derivatives of the adjoint variables ---*/
  
  for (iDim = 0; iDim < nDim; iDim++) {
    Mean_GradPsiE[iDim] =  0.5*(PsiVar_Grad_i[nVar-1][iDim]+PsiVar_Grad_j[nVar-1][iDim]);
    for (jDim = 0; jDim < nDim; jDim++)
      Mean_GradPhi[iDim][jDim] =  0.5*(PsiVar_Grad_i[iDim+1][jDim]+PsiVar_Grad_j[iDim+1][jDim]);
  }
  
  dPhiE_dn = 0;
  for (iDim = 0; iDim < nDim; iDim++)
    dPhiE_dn += Mean_GradPsiE[iDim]*Normal[iDim];
  
  /*--- Compute the viscous residual and jacobian ---*/
  
  GetAdjViscousFlux_Jac(Pressure_i, Pressure_j, Density_i, Density_j,
                        ViscDens_i, ViscDens_j, Velocity_i, Velocity_j, sq_vel_i, sq_vel_j,
                        XiDens_i, XiDens_j, Mean_GradPhi, Mean_GradPsiE,
                        dPhiE_dn, Normal, Edge_Vector, dist_ij_2, val_residual_i, val_residual_j,
                        val_Jacobian_ii, val_Jacobian_ij, val_Jacobian_ji, val_Jacobian_jj,
                        implicit);
  
}

CAvgGradCorrected_AdjFlow::CAvgGradCorrected_AdjFlow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  implicit = (config->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);
  
  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;
  
  Velocity_i = new su2double [nDim];
  Velocity_j = new su2double [nDim];
  Mean_Velocity = new su2double [nDim];
  
  Mean_GradPsiVar = new su2double* [nVar];
  for (unsigned short iVar = 0; iVar < nVar; iVar++)
    Mean_GradPsiVar[iVar] = new su2double [nDim];
  
  Edge_Vector = new su2double [nDim];
  Proj_Mean_GradPsiVar_Edge = new su2double [nVar];
  
  Mean_GradPhi = new su2double* [nDim];
  for (unsigned short iDim = 0; iDim < nDim; iDim++)
    Mean_GradPhi[iDim] = new su2double [nDim];
  Mean_GradPsiE = new su2double [nDim];
  
}

CAvgGradCorrected_AdjFlow::~CAvgGradCorrected_AdjFlow(void) {
  
  delete [] Velocity_i;
  delete [] Velocity_j;
  delete [] Mean_Velocity;
  delete [] Edge_Vector;
  delete [] Proj_Mean_GradPsiVar_Edge;
  
  for (unsigned short iVar = 0; iVar < nVar; iVar++)
    delete [] Mean_GradPsiVar[iVar];
  delete [] Mean_GradPsiVar;
  
  for (unsigned short iDim = 0; iDim < nDim; iDim++)
    delete [] Mean_GradPhi[iDim];
  delete [] Mean_GradPhi;
  delete [] Mean_GradPsiE;
  
}

void CAvgGradCorrected_AdjFlow::ComputeResidual(su2double *val_residual_i,
                                                su2double *val_residual_j,
                                                su2double **val_Jacobian_ii,
                                                su2double **val_Jacobian_ij,
                                            su2double **val_Jacobian_ji,
                                                su2double **val_Jacobian_jj,
                                                CConfig *config) {
  
  unsigned short iVar, iDim, jDim;
  su2double Density_i, sq_vel_i, Pressure_i, ViscDens_i, XiDens_i;
  su2double Density_j, sq_vel_j, Pressure_j, ViscDens_j, XiDens_j;
  su2double dist_ij_2, dPhiE_dn;
  
  su2double Prandtl_Lam  = config->GetPrandtl_Lam();
  su2double Prandtl_Turb = config->GetPrandtl_Turb();
  
  /*--- States in point i ---*/
  
  sq_vel_i = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Velocity_i[iDim] = V_i[iDim+1];
    sq_vel_i += 0.5*Velocity_i[iDim]*Velocity_i[iDim];
  }
  Pressure_i = V_i[nDim+1];
  Density_i  = V_i[nDim+2];
  Enthalpy_i = V_i[nDim+3];

  /*--- Laminar and Eddy viscosity ---*/
  
  Laminar_Viscosity_i = V_i[nDim+5];
  Eddy_Viscosity_i    = V_i[nDim+6];
  
  ViscDens_i = (Laminar_Viscosity_i + Eddy_Viscosity_i) / Density_i;
  XiDens_i   = Gamma*(Laminar_Viscosity_i/Prandtl_Lam +
                      Eddy_Viscosity_i/Prandtl_Turb) / Density_i;

  /*--- States in point j ---*/
  
  sq_vel_j = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Velocity_j[iDim] = V_j[iDim+1];
    sq_vel_j += 0.5*Velocity_j[iDim]*Velocity_j[iDim];
  }
  Pressure_j = V_j[nDim+1];
  Density_j  = V_j[nDim+2];
  Enthalpy_j = V_j[nDim+3];
  
  /*--- Laminar and Eddy viscosity ---*/
  
  Laminar_Viscosity_j = V_j[nDim+5];
  Eddy_Viscosity_j    = V_j[nDim+6];
  
  ViscDens_j = (Laminar_Viscosity_j + Eddy_Viscosity_j) / Density_j;
  XiDens_j   = Gamma*(Laminar_Viscosity_j/Prandtl_Lam +
                      Eddy_Viscosity_j/Prandtl_Turb) / Density_j;
  
  /*--- Compute vector going from iPoint to jPoint ---*/
  
  dist_ij_2 = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Mean_Velocity[iDim] = 0.5*(Velocity_i[iDim]+Velocity_j[iDim]);
    Edge_Vector[iDim] = Coord_j[iDim]-Coord_i[iDim];
    dist_ij_2 += Edge_Vector[iDim]*Edge_Vector[iDim];
  }
  
  /*--- Mean gradient approximation. Projection of the mean gradient in the direction of the edge, weiss correction ---*/
  
  for (iVar = 0; iVar < nVar; iVar++) {
    Proj_Mean_GradPsiVar_Edge[iVar] = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      Mean_GradPsiVar[iVar][iDim] = 0.5*(PsiVar_Grad_i[iVar][iDim] + PsiVar_Grad_j[iVar][iDim]);
      Proj_Mean_GradPsiVar_Edge[iVar] += Mean_GradPsiVar[iVar][iDim]*Edge_Vector[iDim];
    }
    for (iDim = 0; iDim < nDim; iDim++)
      Mean_GradPsiVar[iVar][iDim] -= (Proj_Mean_GradPsiVar_Edge[iVar] -
                                      (Psi_j[iVar]-Psi_i[iVar]))*Edge_Vector[iDim]/dist_ij_2;
  }
  
  /*--- Average of the derivatives of the adjoint variables ---*/
  
  for (iDim = 0; iDim < nDim; iDim++) {
    Mean_GradPsiE[iDim] = Mean_GradPsiVar[nVar-1][iDim];
    for (jDim = 0; jDim < nDim; jDim++)
      Mean_GradPhi[iDim][jDim] = Mean_GradPsiVar[iDim+1][jDim];
  }
  
  dPhiE_dn = 0;
  for (iDim = 0; iDim < nDim; iDim++)
    dPhiE_dn += Mean_GradPsiE[iDim]*Normal[iDim];
  
  /*--- Compute the viscous residual and jacobian ---*/
  
  GetAdjViscousFlux_Jac(Pressure_i, Pressure_j, Density_i, Density_j,
                        ViscDens_i, ViscDens_j, Velocity_i, Velocity_j, sq_vel_i, sq_vel_j,
                        XiDens_i, XiDens_j, Mean_GradPhi, Mean_GradPsiE,
                        dPhiE_dn, Normal, Edge_Vector, dist_ij_2, val_residual_i, val_residual_j,
                        val_Jacobian_ii, val_Jacobian_ij, val_Jacobian_ji, val_Jacobian_jj,
                        implicit);
  
}

CSourceViscous_AdjFlow::CSourceViscous_AdjFlow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  unsigned short iDim;
  
  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;
  
  Velocity = new su2double [nVar];
  GradDensity = new su2double [nDim];
  GradInvDensity = new su2double [nDim];
  dPoDensity2 = new su2double [nDim];
  alpha = new su2double [nDim];
  beta = new su2double [nDim];
  Sigma_5_vec = new su2double [nDim];
  
  GradVel_o_Rho = new su2double* [nDim];
  sigma = new su2double* [nDim];
  Sigma_phi = new su2double* [nDim];
  Sigma_5_Tensor = new su2double* [nDim];
  Sigma = new su2double* [nDim];
  
  for (iDim = 0; iDim < nDim; iDim++) {
    GradVel_o_Rho[iDim] = new su2double [nDim];
    sigma[iDim] = new su2double [nDim];
    Sigma_phi[iDim] = new su2double [nDim];
    Sigma_5_Tensor[iDim] = new su2double [nDim];
    Sigma[iDim] = new su2double [nDim];
  }
  
}

CSourceViscous_AdjFlow::~CSourceViscous_AdjFlow(void) {
  unsigned short iDim;
  
  for (iDim = 0; iDim < nDim; iDim++) {
    delete [] GradVel_o_Rho[iDim];
    delete [] sigma[iDim];
    delete [] Sigma_phi[iDim];
    delete [] Sigma_5_Tensor[iDim];
    delete [] Sigma[iDim];
  }
  
  delete [] GradVel_o_Rho;
  delete [] sigma;
  delete [] Sigma_phi;
  delete [] Sigma_5_Tensor;
  delete [] Sigma;
  
  delete [] Velocity;
  delete [] GradDensity;
  delete [] GradInvDensity;
  delete [] dPoDensity2;
  delete [] alpha;
  delete [] beta;
  delete [] Sigma_5_vec;
  
}

void CSourceViscous_AdjFlow::ComputeResidual (su2double *val_residual, CConfig *config) {
  
  unsigned short iDim, jDim;
  
//  su2double Temperature = V_i[0];
  su2double Pressure = V_i[nDim+1];
  su2double Density = V_i[nDim+2];
//  su2double Enthalpy = V_i[nDim+3];
  su2double Laminar_Viscosity = V_i[nDim+5];
  su2double Eddy_Viscosity = V_i[nDim+6];
  
//  su2double Energy = Enthalpy - Pressure/Density;
  su2double invDensity     = 1.0/Density;
  su2double invDensitysq   = invDensity*invDensity;
  su2double invDensitycube = invDensitysq*invDensity;
  su2double Prandtl_Lam      = config->GetPrandtl_Lam();
  su2double Prandtl_Turb     = config->GetPrandtl_Turb();
  su2double mu_tot_1 = Laminar_Viscosity + Eddy_Viscosity;
  su2double mu_tot_2 = Laminar_Viscosity/Prandtl_Lam + Eddy_Viscosity/Prandtl_Turb;
//  su2double Gas_Constant = config->GetGas_ConstantND();
  
  /*--- Required gradients of the flow variables, point j ---*/
  
  for (iDim = 0; iDim < nDim; iDim++) {
    
    /*--- Gradient density ---*/
    
    GradDensity[iDim] = PrimVar_Grad_i[nDim+2][iDim];
    
    /*--- Gradient (1/rho) ---*/
    
    GradInvDensity[iDim] = -GradDensity[iDim]*invDensitysq;
    
    /*--- Computation of the derivatives of P/(Density^2) ---*/
    
    dPoDensity2[iDim] = (PrimVar_Grad_i[nVar-1][iDim]*Density - 2.0*GradDensity[iDim]*Pressure)*invDensitycube;
    
    /*--- Abbreviations: alpha, beta, sigma_5_vec ---*/
    
    alpha[iDim] = Gamma*mu_tot_2*GradInvDensity[iDim];
    beta[iDim] = Gamma*mu_tot_2*dPoDensity2[iDim]/Gamma_Minus_One;
    Sigma_5_vec[iDim] = Gamma*mu_tot_2*PsiVar_Grad_i[nVar-1][iDim];
    
  }
  
  /*--- Definition of tensors and derivatives of velocity over density ---*/
  
  su2double div_vel = 0.0, div_phi = 0.0, vel_gradpsi5 = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    div_vel += PrimVar_Grad_i[iDim+1][iDim];
    div_phi += PsiVar_Grad_i[iDim+1][iDim];
    vel_gradpsi5 += V_i[iDim+1]*PsiVar_Grad_i[nVar-1][iDim];
    for (jDim = 0; jDim < nDim; jDim++) {
      sigma[iDim][jDim] = mu_tot_1*(PrimVar_Grad_i[iDim+1][jDim]+PrimVar_Grad_i[jDim+1][iDim]);
      Sigma_phi[iDim][jDim] = mu_tot_1*(PsiVar_Grad_i[iDim+1][jDim]+PsiVar_Grad_i[jDim+1][iDim]);
      Sigma_5_Tensor[iDim][jDim] = mu_tot_1*(V_i[jDim+1]*PsiVar_Grad_i[nVar-1][iDim]+V_i[iDim+1]*PsiVar_Grad_i[nVar-1][jDim]);
      GradVel_o_Rho[iDim][jDim] = (PrimVar_Grad_i[iDim+1][jDim]*Density - V_i[iDim+1]*GradDensity[jDim])*invDensitysq;
    }
  }
  
  for (iDim = 0; iDim < nDim; iDim++) {
    sigma[iDim][iDim] -= TWO3*mu_tot_1*div_vel;
    Sigma_phi[iDim][iDim] -= TWO3*mu_tot_1*div_phi;
    Sigma_5_Tensor[iDim][iDim] -= TWO3*mu_tot_1*vel_gradpsi5;
  }
  
  for (iDim = 0; iDim < nDim; iDim++) {
    for (jDim = 0; jDim < nDim; jDim++) {
      Sigma[iDim][jDim] = Sigma_phi[iDim][jDim] + Sigma_5_Tensor[iDim][jDim];
    }
  }
  
  /*--- Vector-Tensors products ---*/
  
  su2double gradT_gradpsi5 = 0.0, sigma_gradpsi = 0.0, vel_sigma_gradpsi5 = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    gradT_gradpsi5 += PrimVar_Grad_i[0][iDim]*PsiVar_Grad_i[nVar-1][iDim];
    for (jDim = 0; jDim < nDim; jDim++) {
      sigma_gradpsi += sigma[iDim][jDim]*PsiVar_Grad_i[jDim+1][iDim];
      vel_sigma_gradpsi5 += V_i[iDim+1]*sigma[iDim][jDim]*PsiVar_Grad_i[nVar-1][jDim];
    }
  }
  
  /*--- Residuals ---*/
  
  su2double alpha_gradpsi5 = 0.0, beta_gradpsi5 = 0.0, Sigma_gradvel_o_rho = 0.0, Sigma5_vel_gradvel = 0.0, sq_vel = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    alpha_gradpsi5 += alpha[iDim]*PsiVar_Grad_i[nVar-1][iDim];
    beta_gradpsi5 += beta[iDim]*PsiVar_Grad_i[nVar-1][iDim];
    for (jDim = 0; jDim < nDim; jDim++) {
      Sigma_gradvel_o_rho += Sigma[iDim][jDim]*GradVel_o_Rho[iDim][jDim];
      Sigma5_vel_gradvel += Sigma_5_vec[iDim]*(V_i[jDim+1]*PrimVar_Grad_i[jDim+1][iDim]);
    }
    sq_vel += V_i[iDim+1]*V_i[iDim+1];
  }
  
  val_residual[0] = (-vel_sigma_gradpsi5*invDensity + 0.5*sq_vel*alpha_gradpsi5 -
                     beta_gradpsi5) * Volume;
  for (iDim = 0; iDim < nDim; iDim++) {
    val_residual[iDim+1] = 0.0;
    for (jDim = 0; jDim < nDim; jDim++) {
      val_residual[iDim+1] += (sigma[iDim][jDim]*PsiVar_Grad_i[nVar-1][jDim]*invDensity - V_i[iDim+1]*alpha[jDim]*PsiVar_Grad_i[nVar-1][jDim]) * Volume;
    }
  }
  val_residual[nVar-1] = alpha_gradpsi5 * Volume;
  
//  val_residual[0] += (Sigma5_vel_gradvel*invDensity - Sigma_gradvel_o_rho) * Volume;
//  for (iDim = 0; iDim < nDim; iDim++) {
//    for (jDim = 0; jDim < nDim; jDim++) {
//      val_residual[iDim+1] += (Sigma[iDim][jDim]*GradInvDensity[jDim] -
//                              Sigma_5_vec[jDim]*PrimVar_Grad_i[iDim+1][jDim]*invDensity) * Volume;
//    }
//  }
  
  /*--- Laminar viscosity sensitivity for NS ---*/
  
  if (config->GetKind_Solver() != ADJ_RANS) {
    
//    su2double Temperature_Ref = config->GetTemperature_Ref();
//    su2double Temperature_Dim = Temperature*Temperature_Ref;
//    
//    su2double S = 0.0;
//    if (config->GetSystemMeasurements() == SI) { S = 110.4; }
//    if (config->GetSystemMeasurements() == US) { S = 198.72; }
//    su2double dVisc_T = ((Laminar_Viscosity)/(2.0*Temperature_Dim*(Temperature_Dim + S)))*(Temperature_Dim + 3.0*S)*Temperature_Ref;
//    
//    su2double Cp = (Gamma/Gamma_Minus_One)*Gas_Constant;
//    su2double kappa_psi = (sigma_gradpsi + vel_sigma_gradpsi5)/mu_tot_1;
//    su2double theta = (kappa_psi + Cp/Prandtl_Lam*gradT_gradpsi5)*dVisc_T*Gamma_Minus_One/(Gas_Constant*Density);
//    
//    val_residual[0] += (theta*(sq_vel-Energy))*Volume;
//    for (iDim = 0; iDim < nDim; iDim++)
//      val_residual[iDim+1] -= theta*V_i[iDim+1]*Volume;
//    val_residual[nVar-1] += theta*Volume;
    
  }
  
//  /*--- Coupling terms coming from the continuous adjoint turbulent equations ---*/
//  
//  if ((config->GetKind_Solver() == ADJ_RANS) && (!config->GetFrozen_Visc_Cont())) {
//    
//    /*--- Closure constants ---*/
//    
//    su2double cv1_3 = 7.1*7.1*7.1;
//    su2double k2 = 0.41*0.41;
//    su2double cb1 = 0.1355;
//    su2double cw2 = 0.3;
//    su2double cw3_6 = pow(2.0,6.0);
//    su2double sigma = 2./3.;
//    su2double cb2 = 0.622;
//    su2double cw1 = cb1/k2+(1+cb2)/sigma;
//    
//    su2double nu, Ji, Ji_2, Ji_3, fv1;
//    nu = Laminar_Viscosity/Density;
//    Ji = TurbVar_i[0]/nu;
//    Ji_2 = Ji*Ji;
//    Ji_3 = Ji_2*Ji;
//    fv1 = Ji_3/(Ji_3+cv1_3);
//    
//    /*--- Contributions due to variation of viscosities ---*/
//    
//    su2double Temperature_Ref = config->GetTemperature_Ref();
//    su2double Temperature_Dim = Temperature*Temperature_Ref;
//    
//    su2double S = 0.0;
//    if (config->GetSystemMeasurements() == SI) { S = 110.4; }
//    if (config->GetSystemMeasurements() == US) { S = 198.72; }
//    su2double dVisc_T = ((Laminar_Viscosity)/(2.0*Temperature_Dim*(Temperature_Dim + S)))*(Temperature_Dim + 3.0*S)*Temperature_Ref;
//    
//    su2double Cp = (Gamma/Gamma_Minus_One)*Gas_Constant;
//    su2double kappa_psi = (sigma_gradpsi + vel_sigma_gradpsi5)/mu_tot_1 + Cp/Prandtl_Turb*gradT_gradpsi5;
//    su2double cv1_const = 3.0*cv1_3/(Ji_3+cv1_3);
//    su2double theta = (kappa_psi*(1.0-Eddy_Viscosity/Laminar_Viscosity*cv1_const) -
//                    Cp/Prandtl_Turb*gradT_gradpsi5*(1.0-Prandtl_Turb/Prandtl_Lam))*dVisc_T*Gamma_Minus_One/(Gas_Constant*Density);
//    su2double xi = kappa_psi*(1.0+cv1_const)*Eddy_Viscosity/Density;
//    
//    val_residual[0] += (theta*(sq_vel-Energy) + xi)*Volume;
//    for (iDim = 0; iDim < nDim; iDim++)
//      val_residual[iDim+1] -= theta*V_i[iDim+1]*Volume;
//    val_residual[nVar-1] += theta*Volume;
//    
//    /*--- Coupling residuals ---*/
//    
//    if (dist_i > 0.0) {
//      su2double fv2, Omega, Shat, dist_0_2, one_o_oneplusJifv1;
//      su2double r, g, g_6, glim, fw;
//      su2double dfw_g, dg_r, dr_nuhat, dr_Shat;
//      su2double dShat_fv2, dfv2_fv1, dfv1_Ji, dJi_nu, dJi_nuhat, dfv2_Ji;
//      
//      /*--- Vorticity ---*/
//      Omega = (PrimVar_Grad_i[1][1]-PrimVar_Grad_i[2][0])*(PrimVar_Grad_i[1][1]-PrimVar_Grad_i[2][0]);
//      if (nDim == 3) Omega += (PrimVar_Grad_i[1][2]-PrimVar_Grad_i[3][0])*(PrimVar_Grad_i[1][2]-PrimVar_Grad_i[3][0]) +
//        (PrimVar_Grad_i[2][2]-PrimVar_Grad_i[3][1])*(PrimVar_Grad_i[2][2]-PrimVar_Grad_i[3][1]);
//      Omega = sqrt(Omega);
//      
//      dist_0_2 = dist_i*dist_i;
//      one_o_oneplusJifv1 = 1.0/(1.0+Ji*fv1);
//      fv2 = 1.0 - Ji*one_o_oneplusJifv1;
//      Shat = max(Omega + TurbVar_i[0]*fv2/(k2*dist_0_2), TURB_EPS);
//      
//      r = min(TurbVar_i[0]/(Shat*k2*dist_0_2), 10.);
//      g = r + cw2*(pow(r,6.)-r);
//      g_6 = pow(g,6.);
//      glim = pow((1+cw3_6)/(g_6+cw3_6),1./6.);
//      fw = g*glim;
//      
//      dfw_g  = glim*cw3_6/(g_6+cw3_6);
//      dg_r = 1.0 + cw2*(6.0*pow(r,5.0)-1.0);
//      dr_nuhat = 1.0/(Shat*k2*dist_0_2);
//      dr_Shat = -dr_nuhat*TurbVar_i[0]/Shat;
//      
//      dShat_fv2 = TurbVar_i[0]/(k2*dist_0_2);
//      dfv2_fv1 = Ji_2*one_o_oneplusJifv1*one_o_oneplusJifv1;
//      dfv1_Ji = 3.0*cv1_3*Ji_2/((Ji_3+cv1_3)*(Ji_3+cv1_3));
//      dJi_nuhat = 1.0/nu;
//      dJi_nu = -Ji/nu;
//      dfv2_Ji = -one_o_oneplusJifv1*one_o_oneplusJifv1;
//      
//      /*--- Terms 1 & 2: -Fcv\B7nabla(TurbPsi_i) - Fs\B7TurbPsi_i ---*/
//      
//      su2double gradTurbVar_gradTurbPsi = 0, vel_gradTurbPsi = 0;
//      for (iDim = 0; iDim < nDim; iDim++) {
//        gradTurbVar_gradTurbPsi += TurbVar_Grad_i[0][iDim]*TurbPsi_Grad_i[0][iDim];
//        vel_gradTurbPsi += V_i[iDim+1]*TurbPsi_Grad_i[0][iDim];
//      }
//      
//      su2double alpha_coeff = Gamma_Minus_One/(Gas_Constant*Density)*dVisc_T;
//      su2double beta_coeff = alpha_coeff*(sq_vel-Energy)-Laminar_Viscosity_i/Density;
//      su2double Fs_coeff = TurbPsi_i[0]*(cb1*TurbVar_i[0]-cw1*TurbVar_i[0]*TurbVar_i[0]/dist_0_2*dfw_g*dg_r*dr_Shat)*
//      dShat_fv2*(dfv2_Ji+dfv2_fv1*dfv1_Ji)*dJi_nu;
//      su2double Gamma = Fs_coeff - gradTurbVar_gradTurbPsi/sigma;
//      
//      val_residual[0] -= (Gamma*beta_coeff - TurbVar_i[0]*vel_gradTurbPsi)/Density*Volume;
//      for (iDim = 0; iDim < nDim; iDim++)
//        val_residual[iDim+1] += (Gamma*alpha_coeff*V_i[iDim+1] - TurbVar_i[0]*TurbPsi_Grad_i[0][iDim])/Density*Volume;
//      val_residual[nVar-1] -= (Gamma*alpha_coeff)/Density*Volume;
//      
//      // this should improve stability (when commented):
//      /*--- Terms 3: -partial{T^s}_GradVel x GradN ---*/
//      //      su2double Ms_coeff = (cb1*TurbVar_i[0]-cw1*TurbVar_i[0]*TurbVar_i[0]/dist_0_2*dfw_g*dg_r*dr_Shat);
//      //      Ms_coeff *= TurbPsi_i[0]/(Omega + TURB_EPS);
//      //
//      //      for (iDim = 0; iDim < nDim; iDim++) {
//      //        for (jDim = 0; jDim < nDim; jDim++) {
//      //          val_residual[0] += Ms_coeff*(PrimVar_Grad_i[iDim+1][jDim]-PrimVar_Grad_i[jDim+1][iDim])*
//      //          GradVel_o_Rho[iDim][jDim]*dV;
//      //          val_residual[iDim+1] -= Ms_coeff*(PrimVar_Grad_i[iDim+1][jDim]-PrimVar_Grad_i[jDim+1][iDim])*
//      //          GradInvDensity[jDim]*dV;
//      //        }
//      //      }
//      
//    }
//  }

}

CSourceConservative_AdjFlow::CSourceConservative_AdjFlow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;
  
  Velocity = new su2double [nDim];
  Residual_i = new su2double [nVar];
  Residual_j = new su2double [nVar];
  Mean_Residual = new su2double [nVar];
  
  Mean_PrimVar_Grad = new su2double* [nVar];
  for (unsigned short iVar = 0; iVar < nVar; iVar++)
    Mean_PrimVar_Grad[iVar] = new su2double [nDim];
}

CSourceConservative_AdjFlow::~CSourceConservative_AdjFlow(void) {
  delete [] Mean_Residual;
  delete [] Residual_j;
  delete [] Residual_i;
  delete [] Velocity;
  
  for (unsigned short iVar = 0; iVar < nVar; iVar++)
    delete [] Mean_PrimVar_Grad[iVar];
  delete [] Mean_PrimVar_Grad;
}

void CSourceConservative_AdjFlow::ComputeResidual (su2double *val_residual, CConfig *config) {
  unsigned short iDim, jDim, iVar;
  su2double rho, nu, Ji, fv1, fv2, Omega, Shat, dist_sq, Ji_2, Ji_3, one_o_oneplusJifv1;
  su2double r, g, g_6, glim, dfw_g, dg_r, dr_nuhat, dr_Shat, Ms_coeff, invOmega;
  
  su2double cv1_3 = 7.1*7.1*7.1;
  su2double k2 = 0.41*0.41;
  su2double cb1 = 0.1355;
  su2double cw2 = 0.3;
  su2double cw3_6 = pow(2.0,6.0);
  su2double sigma = 2./3.;
  su2double cb2 = 0.622;
  su2double cw1 = cb1/k2+(1+cb2)/sigma;
  
  for (iVar = 0; iVar < nVar; iVar++) {
    Residual_i[iVar] = 0.0;
    Residual_j[iVar] = 0.0;
  }
  
  /*--- iPoint ---*/
  
  /*--- Density and velocities ---*/
  
  rho = U_i[0];
  for (iDim = 0; iDim < nDim; iDim++)
    Velocity[iDim] = U_i[iDim+1]/rho;
  
  /*--- Vorticity ---*/
  
  Omega = (PrimVar_Grad_i[1][1]-PrimVar_Grad_i[2][0])*(PrimVar_Grad_i[1][1]-PrimVar_Grad_i[2][0]);
  if (nDim == 3) Omega += (PrimVar_Grad_i[1][2]-PrimVar_Grad_i[3][0])*(PrimVar_Grad_i[1][2]-PrimVar_Grad_i[3][0]) +
    (PrimVar_Grad_i[2][2]-PrimVar_Grad_i[3][1])*(PrimVar_Grad_i[2][2]-PrimVar_Grad_i[3][1]);
  Omega = sqrt(Omega);
  invOmega = 1.0/(Omega + TURB_EPS);
  
  /*--- Compute Ms_coeff -> coming from partial derivatives ---*/
  
  Ms_coeff = 0.0;
  if (dist_i > 0) {
    dist_sq = dist_i*dist_i;
    nu = Laminar_Viscosity_i/rho;
    Ji = TurbVar_i[0]/nu;
    Ji_2 = Ji*Ji;
    Ji_3 = Ji_2*Ji;
    fv1 = Ji_3/(Ji_3+cv1_3);
    one_o_oneplusJifv1 = 1.0/(1.0+Ji*fv1);
    fv2 = 1.0 - Ji*one_o_oneplusJifv1;
    Shat = max(Omega + TurbVar_i[0]*fv2/(k2*dist_sq), TURB_EPS);
    
    r = min(TurbVar_i[0]/(Shat*k2*dist_sq),10.);
    g = r + cw2*(pow(r,6.)-r);
    g_6 = pow(g,6.);
    glim = pow((1+cw3_6)/(g_6+cw3_6),1./6.);
    
    dfw_g  = glim*cw3_6/(g_6+cw3_6);
    dg_r = 1.0 + cw2*(6.0*pow(r,5.0)-1.0);
    dr_nuhat = 1.0/(Shat*k2*dist_sq);
    dr_Shat = -dr_nuhat*TurbVar_i[0]/Shat;
    
    Ms_coeff = (cb1*TurbVar_i[0]-cw1*TurbVar_i[0]*TurbVar_i[0]/dist_sq*dfw_g*dg_r*dr_Shat);
  }
  Ms_coeff *= TurbPsi_i[0]*invOmega/rho;
  
  /*--- Compute residual of iPoint ---*/
  
  for (iDim = 0; iDim < nDim; iDim++) {
    for (jDim = 0; jDim < nDim; jDim++) {
      Residual_i[0] -= Ms_coeff*(Velocity[jDim]*PrimVar_Grad_i[jDim+1][iDim]*Normal[iDim] -
                                 Velocity[jDim]*PrimVar_Grad_i[iDim+1][jDim]*Normal[iDim]);
      Residual_i[iDim+1] += Ms_coeff*(PrimVar_Grad_i[iDim+1][jDim]*Normal[jDim] -
                                      PrimVar_Grad_i[jDim+1][iDim]*Normal[jDim]);
    }
  }
  
  /*--- jPoint ---*/
  
  /*--- Density and velocities ---*/
  
  rho = U_j[0];
  for (iDim = 0; iDim < nDim; iDim++)
    Velocity[iDim] = U_j[iDim+1]/rho;
  
  /*--- Vorticity ---*/
  
  Omega = (PrimVar_Grad_j[1][1]-PrimVar_Grad_j[2][0])*(PrimVar_Grad_j[1][1]-PrimVar_Grad_j[2][0]);
  if (nDim == 3) Omega += (PrimVar_Grad_j[1][2]-PrimVar_Grad_j[3][0])*(PrimVar_Grad_j[1][2]-PrimVar_Grad_j[3][0]) +
    (PrimVar_Grad_j[2][2]-PrimVar_Grad_j[3][1])*(PrimVar_Grad_j[2][2]-PrimVar_Grad_j[3][1]);
  Omega = sqrt(Omega);
  invOmega = 1.0/(Omega + TURB_EPS);
  
  /*--- Compute Ms_coeff -> coming from partial derivatives ---*/
  
  Ms_coeff = 0.0;
  if (dist_j > 0) {
    dist_sq = dist_j*dist_j;
    nu = Laminar_Viscosity_j/rho;
    Ji = TurbVar_j[0]/nu;
    Ji_2 = Ji*Ji;
    Ji_3 = Ji_2*Ji;
    fv1 = Ji_3/(Ji_3+cv1_3);
    one_o_oneplusJifv1 = 1.0/(1.0+Ji*fv1);
    fv2 = 1.0 - Ji*one_o_oneplusJifv1;
    Shat = max(Omega + TurbVar_j[0]*fv2/(k2*dist_sq), TURB_EPS);
    
    r = min(TurbVar_j[0]/(Shat*k2*dist_sq),10.);
    g = r + cw2*(pow(r,6.)-r);
    g_6 = pow(g,6.);
    glim = pow((1+cw3_6)/(g_6+cw3_6),1./6.);
    
    dfw_g  = glim*cw3_6/(g_6+cw3_6);
    dg_r = 1.0 + cw2*(6.0*pow(r,5.0)-1.0);
    dr_nuhat = 1.0/(Shat*k2*dist_sq);
    dr_Shat = -dr_nuhat*TurbVar_j[0]/Shat;
    
    Ms_coeff = (cb1*TurbVar_j[0]-cw1*TurbVar_j[0]*TurbVar_j[0]/dist_sq*dfw_g*dg_r*dr_Shat);
  }
  Ms_coeff *= TurbPsi_j[0]*invOmega/rho;
  
  /*--- Compute residual of jPoint ---*/
  
  for (iDim = 0; iDim < nDim; iDim++) {
    for (jDim = 0; jDim < nDim; jDim++) {
      Residual_j[0] -= Ms_coeff*(Velocity[jDim]*PrimVar_Grad_j[jDim+1][iDim]*Normal[iDim] -
                                 Velocity[jDim]*PrimVar_Grad_j[iDim+1][jDim]*Normal[iDim]);
      Residual_j[iDim+1] += Ms_coeff*(PrimVar_Grad_j[iDim+1][jDim]*Normal[jDim] -
                                      PrimVar_Grad_j[jDim+1][iDim]*Normal[jDim]);
    }
  }
  
  /*--- Compute the mean residual ---*/
  
  for (iVar = 0; iVar < nVar; iVar++)
    val_residual[iVar] = 0.5*(Residual_i[iVar] + Residual_j[iVar]);
  
}

CSourceRotatingFrame_AdjFlow::CSourceRotatingFrame_AdjFlow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) { }

CSourceRotatingFrame_AdjFlow::~CSourceRotatingFrame_AdjFlow(void) { }

void CSourceRotatingFrame_AdjFlow::ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, CConfig *config) {
  
  unsigned short iDim, iVar, jVar;
  su2double Omega[3] = {0,0,0}, Phi[3] = {0,0,0};
  bool implicit = (config->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);

  /*--- Retrieve the angular velocity vector from config. ---*/

  for (iDim = 0; iDim < 3; iDim++){
    Omega[iDim] = config->GetRotation_Rate(iDim)/config->GetOmega_Ref();
  }
  
  /*--- Get the adjoint velocity vector at the current node. ---*/
  
  for (iDim = 0; iDim < nDim; iDim++)
    Phi[iDim] = Psi_i[iDim+1];
  
  /*--- Compute the source term as the Jacobian of the rotating frame
   source term multiplied by the adjoint state and the dual cell volume. ---*/
  
  if (nDim == 2) {
    val_residual[0] = 0.0;
    val_residual[1] =  Omega[2]*Phi[1]*Volume;
    val_residual[2] = -Omega[2]*Phi[0]*Volume;
    val_residual[3] = 0.0;
  } else {
    val_residual[0] = 0.0;
    val_residual[1] = (Omega[2]*Phi[1] - Omega[1]*Phi[2])*Volume;
    val_residual[2] = (Omega[0]*Phi[2] - Omega[2]*Phi[0])*Volume;
    val_residual[3] = (Omega[1]*Phi[0] - Omega[0]*Phi[1])*Volume;
    val_residual[4] = 0.0;
  }
  
  /*--- Calculate the source term Jacobian ---*/
  
  if (implicit) {
    for (iVar = 0; iVar < nVar; iVar++)
      for (jVar = 0; jVar < nVar; jVar++)
        val_Jacobian_i[iVar][jVar] = 0.0;
    if (nDim == 2) {
      val_Jacobian_i[1][2] =  Omega[2]*Volume;
      val_Jacobian_i[2][1] = -Omega[2]*Volume;
    } else {
      val_Jacobian_i[1][2] =  Omega[2]*Volume;
      val_Jacobian_i[1][3] = -Omega[1]*Volume;
      val_Jacobian_i[2][1] = -Omega[2]*Volume;
      val_Jacobian_i[2][3] =  Omega[0]*Volume;
      val_Jacobian_i[3][1] =  Omega[1]*Volume;
      val_Jacobian_i[3][2] = -Omega[0]*Volume;
    }
  }
  
}

CSourceAxisymmetric_AdjFlow::CSourceAxisymmetric_AdjFlow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) { }

CSourceAxisymmetric_AdjFlow::~CSourceAxisymmetric_AdjFlow(void) { }

void CSourceAxisymmetric_AdjFlow::ComputeResidual(su2double *val_residual, su2double **Jacobian_ii, CConfig *config) {
  
  su2double yinv;
  su2double Jacobian_Axisymmetric[4][4];
  
  if (Coord_i[1] > 0.0) yinv = 1.0/Coord_i[1];
  else yinv = 0.0;
  
  Jacobian_Axisymmetric[0][0] = 0;
  Jacobian_Axisymmetric[0][1] = 0;
  Jacobian_Axisymmetric[0][2] = 1.;
  Jacobian_Axisymmetric[0][3] = 0;
  
  Jacobian_Axisymmetric[1][0] = -U_i[1]*U_i[2]/(U_i[0]*U_i[0]);
  Jacobian_Axisymmetric[1][1] = U_i[2]/U_i[0];
  Jacobian_Axisymmetric[1][2] = U_i[1]/U_i[0];
  Jacobian_Axisymmetric[1][3] = 0;
  
  Jacobian_Axisymmetric[2][0] = -U_i[2]*U_i[2]/(U_i[0]*U_i[0]);
  Jacobian_Axisymmetric[2][1] = 0;
  Jacobian_Axisymmetric[2][2] = 2*U_i[2]/U_i[0];
  Jacobian_Axisymmetric[2][3] = 0;
  
  Jacobian_Axisymmetric[3][0] = -Gamma*U_i[2]*U_i[3]/(U_i[0]*U_i[0]) + (Gamma-1)*U_i[2]*(U_i[1]*U_i[1]+U_i[2]*U_i[2])/(U_i[0]*U_i[0]*U_i[0]);
  Jacobian_Axisymmetric[3][1] = -(Gamma-1)*U_i[2]*U_i[1]/(U_i[0]*U_i[0]);
  Jacobian_Axisymmetric[3][2] = Gamma*U_i[3]/U_i[0] - 1/2*(Gamma-1)*( (U_i[1]*U_i[1]+U_i[2]*U_i[2])/(U_i[0]*U_i[0]) + 2*U_i[2]*U_i[2]/(U_i[0]*U_i[0]) );
  Jacobian_Axisymmetric[3][3] = Gamma*U_i[2]/U_i[0];
  
  for (int iVar=0; iVar<4; iVar++)
    for (int jVar=0; jVar<4; jVar++)
      Jacobian_Axisymmetric[iVar][jVar] *= yinv*Volume;
  
  /* -- Residual = transpose(Jacobian) * psi --*/
  for (int iVar = 0; iVar < nVar; iVar++) {
    val_residual[iVar] = 0.0;
    for (int jVar = 0; jVar < nVar; jVar++) {
      val_residual[iVar] += Jacobian_Axisymmetric[jVar][iVar]*Psi_i[jVar];
      Jacobian_ii[iVar][jVar] = Jacobian_Axisymmetric[jVar][iVar];
    }
  }
}
