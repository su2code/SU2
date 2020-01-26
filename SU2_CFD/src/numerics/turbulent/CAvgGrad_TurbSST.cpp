/*!
 * \file CAvgGrad_TurbSST.cpp
 * \brief Implementation of numerics class CAvgGrad_TurbSST.
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

#include "../../../include/numerics/turbulent/CAvgGrad_TurbSST.hpp"

CAvgGrad_TurbSST::CAvgGrad_TurbSST(unsigned short val_nDim,
                                   unsigned short val_nVar,
                                   const su2double *constants,
                                   bool correct_grad,
                                   CConfig *config)
 : CAvgGrad_Scalar(val_nDim, val_nVar, correct_grad, config) {

  sigma_k1  = constants[0];
  sigma_om1 = constants[2];
  sigma_k2  = constants[1];
  sigma_om2 = constants[3];

  F1_i = 0.0; F1_j = 0.0;
  diff_kine = 0.0;
  diff_omega = 0.0;

}

CAvgGrad_TurbSST::~CAvgGrad_TurbSST(void) {
}

void CAvgGrad_TurbSST::ExtraADPreaccIn() {
  AD::SetPreaccIn(F1_i); AD::SetPreaccIn(F1_j);
}

void CAvgGrad_TurbSST::FinishResidualCalc(su2double *val_residual, su2double **Jacobian_i, su2double **Jacobian_j, CConfig *config) {

  su2double sigma_kine_i, sigma_kine_j, sigma_omega_i, sigma_omega_j;
  su2double diff_i_kine, diff_i_omega, diff_j_kine, diff_j_omega;

  /*--- Compute the blended constant for the viscous terms ---*/
  sigma_kine_i  = F1_i*sigma_k1 + (1.0 - F1_i)*sigma_k2;
  sigma_kine_j  = F1_j*sigma_k1 + (1.0 - F1_j)*sigma_k2;
  sigma_omega_i = F1_i*sigma_om1 + (1.0 - F1_i)*sigma_om2;
  sigma_omega_j = F1_j*sigma_om1 + (1.0 - F1_j)*sigma_om2;

  /*--- Compute mean effective viscosity ---*/
  diff_i_kine  = Laminar_Viscosity_i + sigma_kine_i*Eddy_Viscosity_i;
  diff_j_kine  = Laminar_Viscosity_j + sigma_kine_j*Eddy_Viscosity_j;
  diff_i_omega = Laminar_Viscosity_i + sigma_omega_i*Eddy_Viscosity_i;
  diff_j_omega = Laminar_Viscosity_j + sigma_omega_j*Eddy_Viscosity_j;

  diff_kine  = 0.5*(diff_i_kine + diff_j_kine);    // Could instead use weighted average!
  diff_omega = 0.5*(diff_i_omega + diff_j_omega);

  val_residual[0] = diff_kine*Proj_Mean_GradTurbVar[0];
  val_residual[1] = diff_omega*Proj_Mean_GradTurbVar[1];

  /*--- For Jacobians -> Use of TSL approx. to compute derivatives of the gradients ---*/
  if (implicit) {
    Jacobian_i[0][0] = -diff_kine*proj_vector_ij/Density_i;    Jacobian_i[0][1] = 0.0;
    Jacobian_i[1][0] = 0.0;                      Jacobian_i[1][1] = -diff_omega*proj_vector_ij/Density_i;

    Jacobian_j[0][0] = diff_kine*proj_vector_ij/Density_j;     Jacobian_j[0][1] = 0.0;
    Jacobian_j[1][0] = 0.0;                      Jacobian_j[1][1] = diff_omega*proj_vector_ij/Density_j;
  }

}
