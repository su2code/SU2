/*!
 * \file CSourceViscous_AdjFlow.cpp
 * \brief Implementation of numerics class CSourceViscous_AdjFlow.
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

#include "../../../include/numerics/continuous_adjoint/CSourceViscous_AdjFlow.hpp"

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
