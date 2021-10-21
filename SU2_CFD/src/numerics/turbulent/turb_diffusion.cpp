/*!
 * \file turb_diffusion.cpp
 * \brief Implementation of numerics classes to compute viscous
 *        fluxes in turbulence problems.
 * \author F. Palacios, T. Economon
 * \version 7.2.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2021, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../../include/numerics/turbulent/turb_diffusion.hpp"


CAvgGrad_TurbSA::CAvgGrad_TurbSA(unsigned short val_nDim, unsigned short val_nVar,
                                 bool correct_grad, const CConfig* config) :
                 CAvgGrad_Scalar(val_nDim, val_nVar, correct_grad, config) { }

void CAvgGrad_TurbSA::ExtraADPreaccIn() { }

void CAvgGrad_TurbSA::FinishResidualCalc(const CConfig* config) {

  /*--- Compute mean effective viscosity ---*/

  su2double nu_i = Laminar_Viscosity_i/Density_i;
  su2double nu_j = Laminar_Viscosity_j/Density_j;
  su2double nu_e = 0.5*(nu_i+nu_j+ScalarVar_i[0]+ScalarVar_j[0]);

  Flux[0] = nu_e*Proj_Mean_GradScalarVar[0]/sigma;

  /*--- For Jacobians -> Use of TSL approx. to compute derivatives of the gradients ---*/

  if (implicit) {
    Jacobian_i[0][0] = (0.5*Proj_Mean_GradScalarVar[0]-nu_e*proj_vector_ij)/sigma;
    Jacobian_j[0][0] = (0.5*Proj_Mean_GradScalarVar[0]+nu_e*proj_vector_ij)/sigma;
  }

}

CAvgGrad_TurbSA_Neg::CAvgGrad_TurbSA_Neg(unsigned short val_nDim,
                                         unsigned short val_nVar,
                                         bool correct_grad,
                                         const CConfig* config) :
                     CAvgGrad_Scalar(val_nDim, val_nVar, correct_grad, config) { }

void CAvgGrad_TurbSA_Neg::ExtraADPreaccIn() { }

void CAvgGrad_TurbSA_Neg::FinishResidualCalc(const CConfig* config) {

  /*--- Compute mean effective viscosity ---*/

  su2double nu_i = Laminar_Viscosity_i/Density_i;
  su2double nu_j = Laminar_Viscosity_j/Density_j;

  su2double nu_ij = 0.5*(nu_i+nu_j);
  su2double nu_tilde_ij = 0.5*(ScalarVar_i[0]+ScalarVar_j[0]);

  su2double nu_e;

  if (nu_tilde_ij > 0.0) {
    nu_e = nu_ij + nu_tilde_ij;
  }
  else {
    su2double Xi = nu_tilde_ij/nu_ij;
    su2double fn = (cn1 + Xi*Xi*Xi)/(cn1 - Xi*Xi*Xi);
    nu_e = nu_ij + fn*nu_tilde_ij;
  }

  Flux[0] = nu_e*Proj_Mean_GradScalarVar_Normal[0]/sigma;

  /*--- For Jacobians -> Use of TSL approx. to compute derivatives of the gradients ---*/

  if (implicit) {
    Jacobian_i[0][0] = (0.5*Proj_Mean_GradScalarVar[0]-nu_e*proj_vector_ij)/sigma;
    Jacobian_j[0][0] = (0.5*Proj_Mean_GradScalarVar[0]+nu_e*proj_vector_ij)/sigma;
  }

}

CAvgGrad_TurbSST::CAvgGrad_TurbSST(unsigned short val_nDim,
                                   unsigned short val_nVar,
                                   const su2double *constants,
                                   bool correct_grad,
                                   const CConfig* config) :
  CAvgGrad_Scalar(val_nDim, val_nVar, correct_grad, config),
  sigma_k1(constants[0]),
  sigma_k2(constants[1]),
  sigma_om1(constants[2]),
  sigma_om2(constants[3]) {

}

void CAvgGrad_TurbSST::ExtraADPreaccIn() {
  AD::SetPreaccIn(F1_i); AD::SetPreaccIn(F1_j);
}

void CAvgGrad_TurbSST::FinishResidualCalc(const CConfig* config) {

  su2double sigma_kine_i, sigma_kine_j, sigma_omega_i, sigma_omega_j;
  su2double diff_i_kine, diff_i_omega, diff_j_kine, diff_j_omega;

  /*--- Compute the blended constant for the viscous terms ---*/
  sigma_kine_i = F1_i*sigma_k1 + (1.0 - F1_i)*sigma_k2;
  sigma_kine_j = F1_j*sigma_k1 + (1.0 - F1_j)*sigma_k2;
  sigma_omega_i = F1_i*sigma_om1 + (1.0 - F1_i)*sigma_om2;
  sigma_omega_j = F1_j*sigma_om1 + (1.0 - F1_j)*sigma_om2;

  /*--- Compute mean effective dynamic viscosity ---*/
  diff_i_kine = Laminar_Viscosity_i + sigma_kine_i*Eddy_Viscosity_i;
  diff_j_kine = Laminar_Viscosity_j + sigma_kine_j*Eddy_Viscosity_j;
  diff_i_omega = Laminar_Viscosity_i + sigma_omega_i*Eddy_Viscosity_i;
  diff_j_omega = Laminar_Viscosity_j + sigma_omega_j*Eddy_Viscosity_j;

  su2double diff_kine = 0.5*(diff_i_kine + diff_j_kine);
  su2double diff_omega = 0.5*(diff_i_omega + diff_j_omega);

  Flux[0] = diff_kine*Proj_Mean_GradScalarVar[0];
  Flux[1] = diff_omega*Proj_Mean_GradScalarVar[1];

  /*--- For Jacobians -> Use of TSL (Thin Shear Layer) approx. to compute derivatives of the gradients ---*/
  if (implicit) {
    su2double proj_on_rho = proj_vector_ij/Density_i;

    Jacobian_i[0][0] = -diff_kine*proj_on_rho;  Jacobian_i[0][1] = 0.0;
    Jacobian_i[1][0] = 0.0;                     Jacobian_i[1][1] = -diff_omega*proj_on_rho;

    proj_on_rho = proj_vector_ij/Density_j;

    Jacobian_j[0][0] = diff_kine*proj_on_rho;   Jacobian_j[0][1] = 0.0;
    Jacobian_j[1][0] = 0.0;                     Jacobian_j[1][1] = diff_omega*proj_on_rho;
  }

}
