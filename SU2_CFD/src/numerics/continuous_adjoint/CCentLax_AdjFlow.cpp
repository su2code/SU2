/*!
 * \file CCentLax_AdjFlow.cpp
 * \brief Implementation of numerics class CCentLax_AdjFlow.
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

#include "../../../include/numerics/continuous_adjoint/CCentLax_AdjFlow.hpp"

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
