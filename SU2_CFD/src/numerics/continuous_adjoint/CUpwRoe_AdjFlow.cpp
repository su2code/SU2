/*!
 * \file CUpwRoe_AdjFlow.cpp
 * \brief Implementation of numerics class CUpwRoe_AdjFlow.
 * \author F. Palacios, T. Economon
 * \version 7.0.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2019, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../../include/numerics/continuous_adjoint/CUpwRoe_AdjFlow.hpp"

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
