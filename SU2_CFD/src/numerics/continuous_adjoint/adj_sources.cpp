/*!
 * \file adj_sources.cpp
 * \brief Implementation of adjoint source numerics classes.
 * \author F. Palacios, T. Economon
 * \version 8.0.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2023, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../../include/numerics/continuous_adjoint/adj_sources.hpp"

CSourceAxisymmetric_AdjFlow::CSourceAxisymmetric_AdjFlow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) { }

CSourceAxisymmetric_AdjFlow::~CSourceAxisymmetric_AdjFlow() = default;

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

CSourceConservative_AdjFlow::~CSourceConservative_AdjFlow() {
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
    Ji = ScalarVar_i[0]/nu;
    Ji_2 = Ji*Ji;
    Ji_3 = Ji_2*Ji;
    fv1 = Ji_3/(Ji_3+cv1_3);
    one_o_oneplusJifv1 = 1.0/(1.0+Ji*fv1);
    fv2 = 1.0 - Ji*one_o_oneplusJifv1;
    Shat = max(Omega + ScalarVar_i[0]*fv2/(k2*dist_sq), TURB_EPS);

    r = min(ScalarVar_i[0]/(Shat*k2*dist_sq),10.);
    g = r + cw2*(pow(r,6.)-r);
    g_6 = pow(g,6.);
    glim = pow((1+cw3_6)/(g_6+cw3_6),1./6.);

    dfw_g  = glim*cw3_6/(g_6+cw3_6);
    dg_r = 1.0 + cw2*(6.0*pow(r,5.0)-1.0);
    dr_nuhat = 1.0/(Shat*k2*dist_sq);
    dr_Shat = -dr_nuhat*ScalarVar_i[0]/Shat;

    Ms_coeff = (cb1*ScalarVar_i[0]-cw1*ScalarVar_i[0]*ScalarVar_i[0]/dist_sq*dfw_g*dg_r*dr_Shat);
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
    Ji = ScalarVar_j[0]/nu;
    Ji_2 = Ji*Ji;
    Ji_3 = Ji_2*Ji;
    fv1 = Ji_3/(Ji_3+cv1_3);
    one_o_oneplusJifv1 = 1.0/(1.0+Ji*fv1);
    fv2 = 1.0 - Ji*one_o_oneplusJifv1;
    Shat = max(Omega + ScalarVar_j[0]*fv2/(k2*dist_sq), TURB_EPS);

    r = min(ScalarVar_j[0]/(Shat*k2*dist_sq),10.);
    g = r + cw2*(pow(r,6.)-r);
    g_6 = pow(g,6.);
    glim = pow((1+cw3_6)/(g_6+cw3_6),1./6.);

    dfw_g  = glim*cw3_6/(g_6+cw3_6);
    dg_r = 1.0 + cw2*(6.0*pow(r,5.0)-1.0);
    dr_nuhat = 1.0/(Shat*k2*dist_sq);
    dr_Shat = -dr_nuhat*ScalarVar_j[0]/Shat;

    Ms_coeff = (cb1*ScalarVar_j[0]-cw1*ScalarVar_j[0]*ScalarVar_j[0]/dist_sq*dfw_g*dg_r*dr_Shat);
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

CSourceRotatingFrame_AdjFlow::CSourceRotatingFrame_AdjFlow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) :
                              CNumerics(val_nDim, val_nVar, config) { }

CSourceRotatingFrame_AdjFlow::~CSourceRotatingFrame_AdjFlow() = default;

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

CSourceViscous_AdjFlow::CSourceViscous_AdjFlow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) :
                        CNumerics(val_nDim, val_nVar, config) {
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

CSourceViscous_AdjFlow::~CSourceViscous_AdjFlow() {
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

  if (config->GetKind_Solver() != MAIN_SOLVER::ADJ_RANS) {

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
//  if ((config->GetKind_Solver() == MAIN_SOLVER::ADJ_RANS) && (!config->GetFrozen_Visc_Cont())) {
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
//    Ji = ScalarVar_i[0]/nu;
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
//      Shat = max(Omega + ScalarVar_i[0]*fv2/(k2*dist_0_2), TURB_EPS);
//
//      r = min(ScalarVar_i[0]/(Shat*k2*dist_0_2), 10.);
//      g = r + cw2*(pow(r,6.)-r);
//      g_6 = pow(g,6.);
//      glim = pow((1+cw3_6)/(g_6+cw3_6),1./6.);
//      fw = g*glim;
//
//      dfw_g  = glim*cw3_6/(g_6+cw3_6);
//      dg_r = 1.0 + cw2*(6.0*pow(r,5.0)-1.0);
//      dr_nuhat = 1.0/(Shat*k2*dist_0_2);
//      dr_Shat = -dr_nuhat*ScalarVar_i[0]/Shat;
//
//      dShat_fv2 = ScalarVar_i[0]/(k2*dist_0_2);
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
//        gradTurbVar_gradTurbPsi += ScalarVar_Grad_i[0][iDim]*TurbPsi_Grad_i[0][iDim];
//        vel_gradTurbPsi += V_i[iDim+1]*TurbPsi_Grad_i[0][iDim];
//      }
//
//      su2double alpha_coeff = Gamma_Minus_One/(Gas_Constant*Density)*dVisc_T;
//      su2double beta_coeff = alpha_coeff*(sq_vel-Energy)-Laminar_Viscosity_i/Density;
//      su2double Fs_coeff = TurbPsi_i[0]*(cb1*ScalarVar_i[0]-cw1*ScalarVar_i[0]*ScalarVar_i[0]/dist_0_2*dfw_g*dg_r*dr_Shat)*
//      dShat_fv2*(dfv2_Ji+dfv2_fv1*dfv1_Ji)*dJi_nu;
//      su2double Gamma = Fs_coeff - gradTurbVar_gradTurbPsi/sigma;
//
//      val_residual[0] -= (Gamma*beta_coeff - ScalarVar_i[0]*vel_gradTurbPsi)/Density*Volume;
//      for (iDim = 0; iDim < nDim; iDim++)
//        val_residual[iDim+1] += (Gamma*alpha_coeff*V_i[iDim+1] - ScalarVar_i[0]*TurbPsi_Grad_i[0][iDim])/Density*Volume;
//      val_residual[nVar-1] -= (Gamma*alpha_coeff)/Density*Volume;
//
//      // this should improve stability (when commented):
//      /*--- Terms 3: -partial{T^s}_GradVel x GradN ---*/
//      //      su2double Ms_coeff = (cb1*ScalarVar_i[0]-cw1*ScalarVar_i[0]*ScalarVar_i[0]/dist_0_2*dfw_g*dg_r*dr_Shat);
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

CSourceConservative_AdjTurb::CSourceConservative_AdjTurb(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {

  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;

}

CSourceConservative_AdjTurb::~CSourceConservative_AdjTurb() = default;

void CSourceConservative_AdjTurb::ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config) {

  /*--- SOURCE term  -->  \nabla ( \psi_\mu \B7 E^{s} )
   E^{s} = 2 c_{b2}/\sigma \nabla \hat{nu} ---*/

  unsigned short iDim;
  bool implicit = (config->GetKind_TimeIntScheme_AdjTurb() == EULER_IMPLICIT);

  su2double cb2 = 0.622;
  su2double sigma = 2./3.;
  su2double coeff = 2.0*cb2/sigma;
  su2double E_ij, proj_TurbVar_Grad_i, proj_TurbVar_Grad_j;

  E_ij = 0.0;  proj_TurbVar_Grad_i = 0.0; proj_TurbVar_Grad_j = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    proj_TurbVar_Grad_i += coeff*ScalarVar_Grad_i[0][iDim]*Normal[iDim];
    proj_TurbVar_Grad_j += coeff*ScalarVar_Grad_j[0][iDim]*Normal[iDim];
    E_ij += 0.5*(TurbPsi_i[0]*proj_TurbVar_Grad_i + TurbPsi_j[0]*proj_TurbVar_Grad_j);
  }

  val_residual[0] = E_ij;

  if (implicit) {
    val_Jacobian_i[0][0] = 0.5*proj_TurbVar_Grad_i;
    val_Jacobian_j[0][0] = 0.5*proj_TurbVar_Grad_j;
  }
}


CSourcePieceWise_AdjTurb::CSourcePieceWise_AdjTurb(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {

  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;

  Velocity = new su2double [nDim];
  tau = new su2double* [nDim];
  for (unsigned short iDim = 0; iDim < nDim; iDim++)
    tau[iDim] = new su2double [nDim];
}

CSourcePieceWise_AdjTurb::~CSourcePieceWise_AdjTurb() {
  delete [] Velocity;

  for (unsigned short iDim = 0; iDim < nDim; iDim++)
    delete [] tau[iDim];
  delete [] tau;
}

void CSourcePieceWise_AdjTurb::ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config) {
  unsigned short iDim, jDim;

  bool implicit = (config->GetKind_TimeIntScheme_AdjTurb() == EULER_IMPLICIT);
  su2double Prandtl_Turb = config->GetPrandtl_Turb();

  val_residual[0] = 0.0;
  if (implicit)
    val_Jacobian_i[0][0] = 0.0;

  if (dist_i > 0.0) {

    /*--- Computation of Vorticity and Divergence of velocity ---*/
    su2double div_vel = 0;
    for (iDim = 0; iDim < nDim; iDim++) {
      Velocity[iDim] = U_i[iDim+1]/U_i[0];
      div_vel += PrimVar_Grad_i[iDim+1][iDim];
    }

    su2double Vorticity = (PrimVar_Grad_i[2][0]-PrimVar_Grad_i[1][1])*(PrimVar_Grad_i[2][0]-PrimVar_Grad_i[1][1]);
    if (nDim == 3)
      Vorticity += ( (PrimVar_Grad_i[3][1]-PrimVar_Grad_i[2][2])*(PrimVar_Grad_i[3][1]-PrimVar_Grad_i[2][2]) +
                    (PrimVar_Grad_i[1][2]-PrimVar_Grad_i[3][0])*(PrimVar_Grad_i[1][2]-PrimVar_Grad_i[3][0]) );
    Vorticity = sqrt(Vorticity);

    /*--- FIRST PART: -Bs*TurbPsi_i ---*/
    /*--- CLOUSURE CONSTANTS ---*/
    su2double cv1 = 7.1;
    su2double cv1_3 = cv1*cv1*cv1;
    su2double k = 0.41;
    su2double k2 = k*k;
    su2double cb1 = 0.1355;
    su2double cw2 = 0.3;
    su2double cw3_6 = pow(2.0,6.0);
    su2double sigma = 2./3.;
    su2double cb2 = 0.622;
    su2double cw1 = cb1/k2+(1+cb2)/sigma;

    su2double nu, Ji, fv1, fv2, Shat, dist_0_2, Ji_2, Ji_3, one_o_oneplusJifv1;
    su2double r, g, g_6, glim, fw;
    su2double dTs_nuhat, dTs_Shat, dShat_nuhat, dTs_fw, dfw_g, dg_r, dr_nuhat, dr_Shat;
    su2double dShat_fv2, dfv2_fv1, dfv1_Ji, dJi_nuhat, dfv2_Ji;
    su2double Bs;

    dist_0_2 = dist_i*dist_i;
    nu = Laminar_Viscosity_i/U_i[0];
    Ji = ScalarVar_i[0]/nu;
    Ji_2 = Ji*Ji;
    Ji_3 = Ji_2*Ji;
    fv1 = Ji_3/(Ji_3+cv1_3);
    one_o_oneplusJifv1 = 1.0/(1.0+Ji*fv1);
    fv2 = 1.0 - Ji*one_o_oneplusJifv1;
    Shat = max(Vorticity + ScalarVar_i[0]*fv2/(k2*dist_0_2), TURB_EPS);

    //    r = ScalarVar_i[0]/(Shat*k2*dist_0_2);
    r = min(ScalarVar_i[0]/(Shat*k2*dist_0_2),10.);
    g = r + cw2*(pow(r,6.)-r);
    g_6 = pow(g,6.);
    glim = pow((1+cw3_6)/(g_6+cw3_6),1./6.);
    fw = g*glim;

    dTs_nuhat = cb1*Shat-2.0*cw1*fw*ScalarVar_i[0]/dist_0_2;
    dTs_Shat = cb1*ScalarVar_i[0];
    dTs_fw = -cw1*ScalarVar_i[0]*ScalarVar_i[0]/dist_0_2;
    dfw_g  = glim*cw3_6/(g_6+cw3_6);
    dg_r = 1.0 + cw2*(6.0*pow(r,5.0)-1.0);
    dr_nuhat = 1.0/(Shat*k2*dist_0_2);
    dr_Shat = -dr_nuhat*ScalarVar_i[0]/Shat;

    dShat_nuhat = fv2/(k2*dist_0_2);
    dShat_fv2 = ScalarVar_i[0]/(k2*dist_0_2);
    dfv2_fv1 = Ji_2*one_o_oneplusJifv1*one_o_oneplusJifv1;
    dfv1_Ji = 3.0*cv1_3*Ji_2/((Ji_3+cv1_3)*(Ji_3+cv1_3));
    dJi_nuhat = 1.0/nu;
    dfv2_Ji = -one_o_oneplusJifv1*one_o_oneplusJifv1;
    dShat_nuhat += dShat_fv2*(dfv2_fv1*dfv1_Ji+dfv2_Ji)*dJi_nuhat;

    Bs = dTs_nuhat;                       // nu_hat term
    Bs += dTs_Shat*dShat_nuhat;                 // S_hat term
    Bs += dTs_fw*dfw_g*dg_r*(dr_nuhat+dr_Shat*dShat_nuhat);   // fw terms

    val_residual[0] = -Bs*TurbPsi_i[0]*Volume;

    if (implicit)
      val_Jacobian_i[0][0] = -Bs*Volume;

    /*---SECOND PART: \partial_nu_hat mu^k F^{vk} cdot \grad Psi ---*/
    su2double dEddyVisc_nuhat;
    if (!config->GetFrozen_Visc_Cont())
      dEddyVisc_nuhat = U_i[0]*fv1*(1.0 + 3.0*cv1_3/(Ji_3+cv1_3));
    else
      dEddyVisc_nuhat = 0;

    for (iDim = 0; iDim < nDim; iDim++) {
      for (jDim = 0; jDim < nDim; jDim++)
        tau[iDim][jDim] = PrimVar_Grad_i[iDim+1][jDim] + PrimVar_Grad_i[jDim+1][iDim];
      tau[iDim][iDim] -= TWO3*div_vel;
    }

    su2double Gas_Constant = config->GetGas_ConstantND();
    su2double Cp = (Gamma/Gamma_Minus_One)*Gas_Constant;
    su2double tau_gradphi = 0.0, vel_tau_gradpsi5 = 0.0, gradT_gradpsi5 = 0.0;

    for (iDim = 0; iDim < nDim; iDim++) {
      gradT_gradpsi5 += PrimVar_Grad_i[0][iDim]*PsiVar_Grad_i[nVar-1][iDim];
      for (jDim = 0; jDim < nDim; jDim++) {
        tau_gradphi += tau[iDim][jDim]*PsiVar_Grad_i[iDim+1][jDim];
        vel_tau_gradpsi5 += Velocity[iDim]*tau[iDim][jDim]*PsiVar_Grad_i[nVar-1][jDim];
      }
    }
    val_residual[0] += (tau_gradphi + vel_tau_gradpsi5 + Cp/Prandtl_Turb*gradT_gradpsi5)*dEddyVisc_nuhat*Volume;

  }
}
