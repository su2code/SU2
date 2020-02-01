/*!
 * \file turb_diffusion.cpp
 * \brief Implementation of numerics classes to compute viscous
 *        fluxes in turbulence problems.
 * \author F. Palacios, T. Economon
 * \version 7.0.1 "Blackbird"
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

#include "../../../include/numerics/turbulent/turb_diffusion.hpp"

CAvgGrad_Scalar::CAvgGrad_Scalar(unsigned short val_nDim,
                                     unsigned short val_nVar,
                                     bool correct_grad,
                                     CConfig *config)
    : CNumerics(val_nDim, val_nVar, config), correct_gradient(correct_grad) {

  implicit = (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT);
  incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);

  Edge_Vector = new su2double [nDim];
  Proj_Mean_GradTurbVar_Normal = new su2double [nVar];
  Proj_Mean_GradTurbVar_Edge = new su2double [nVar];
  Proj_Mean_GradTurbVar = new su2double [nVar];
  Mean_GradTurbVar = new su2double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++)
    Mean_GradTurbVar[iVar] = new su2double [nDim];

}

CAvgGrad_Scalar::~CAvgGrad_Scalar(void) {

  delete [] Edge_Vector;
  delete [] Proj_Mean_GradTurbVar_Normal;
  delete [] Proj_Mean_GradTurbVar_Edge;
  delete [] Proj_Mean_GradTurbVar;
  for (iVar = 0; iVar < nVar; iVar++)
    delete [] Mean_GradTurbVar[iVar];
  delete [] Mean_GradTurbVar;

}

void CAvgGrad_Scalar::ComputeResidual(su2double *val_residual,
                                        su2double **Jacobian_i,
                                        su2double **Jacobian_j,
                                        CConfig *config) {

  AD::StartPreacc();
  AD::SetPreaccIn(Coord_i, nDim); AD::SetPreaccIn(Coord_j, nDim);
  AD::SetPreaccIn(Normal, nDim);
  AD::SetPreaccIn(TurbVar_Grad_i, nVar, nDim);
  AD::SetPreaccIn(TurbVar_Grad_j, nVar, nDim);
  if (correct_gradient) {
    AD::SetPreaccIn(TurbVar_i, nVar); AD::SetPreaccIn(TurbVar_j ,nVar);
  }
  ExtraADPreaccIn();

  if (incompressible) {
    AD::SetPreaccIn(V_i, nDim+6); AD::SetPreaccIn(V_j, nDim+6);

    Density_i = V_i[nDim+2];            Density_j = V_j[nDim+2];
    Laminar_Viscosity_i = V_i[nDim+4];  Laminar_Viscosity_j = V_j[nDim+4];
    Eddy_Viscosity_i = V_i[nDim+5];     Eddy_Viscosity_j = V_j[nDim+5];
  }
  else {
    AD::SetPreaccIn(V_i, nDim+7); AD::SetPreaccIn(V_j, nDim+7);

    Density_i = V_i[nDim+2];            Density_j = V_j[nDim+2];
    Laminar_Viscosity_i = V_i[nDim+5];  Laminar_Viscosity_j = V_j[nDim+5];
    Eddy_Viscosity_i = V_i[nDim+6];     Eddy_Viscosity_j = V_j[nDim+6];
  }

  /*--- Compute vector going from iPoint to jPoint ---*/

  dist_ij_2 = 0; proj_vector_ij = 0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Edge_Vector[iDim] = Coord_j[iDim]-Coord_i[iDim];
    dist_ij_2 += Edge_Vector[iDim]*Edge_Vector[iDim];
    proj_vector_ij += Edge_Vector[iDim]*Normal[iDim];
  }
  if (dist_ij_2 == 0.0) proj_vector_ij = 0.0;
  else proj_vector_ij = proj_vector_ij/dist_ij_2;

  /*--- Mean gradient approximation ---*/
  for (iVar = 0; iVar < nVar; iVar++) {
    Proj_Mean_GradTurbVar_Normal[iVar] = 0.0;
    Proj_Mean_GradTurbVar_Edge[iVar] = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      Mean_GradTurbVar[iVar][iDim] = 0.5*(TurbVar_Grad_i[iVar][iDim] +
                                          TurbVar_Grad_j[iVar][iDim]);
      Proj_Mean_GradTurbVar_Normal[iVar] += Mean_GradTurbVar[iVar][iDim] *
                                            Normal[iDim];
      if (correct_gradient)
        Proj_Mean_GradTurbVar_Edge[iVar] += Mean_GradTurbVar[iVar][iDim]*Edge_Vector[iDim];
    }
    Proj_Mean_GradTurbVar[iVar] = Proj_Mean_GradTurbVar_Normal[iVar];
    if (correct_gradient) {
      Proj_Mean_GradTurbVar[iVar] -= Proj_Mean_GradTurbVar_Edge[iVar]*proj_vector_ij -
      (TurbVar_j[iVar]-TurbVar_i[iVar])*proj_vector_ij;
    }
  }

  FinishResidualCalc(val_residual, Jacobian_i, Jacobian_j, config);

  AD::SetPreaccOut(val_residual, nVar);
  AD::EndPreacc();

}

CAvgGrad_TurbSA::CAvgGrad_TurbSA(unsigned short val_nDim,
                                 unsigned short val_nVar, bool correct_grad,
                                 CConfig *config)
   : CAvgGrad_Scalar(val_nDim, val_nVar, correct_grad, config), sigma(2./3.) {
}

CAvgGrad_TurbSA::~CAvgGrad_TurbSA(void) {
}

void CAvgGrad_TurbSA::ExtraADPreaccIn() {
}

void CAvgGrad_TurbSA::FinishResidualCalc(su2double *val_residual, su2double **Jacobian_i, su2double **Jacobian_j, CConfig *config) {

  /*--- Compute mean effective viscosity ---*/

  nu_i = Laminar_Viscosity_i/Density_i;
  nu_j = Laminar_Viscosity_j/Density_j;
  nu_e = 0.5*(nu_i+nu_j+TurbVar_i[0]+TurbVar_j[0]);

  val_residual[0] = nu_e*Proj_Mean_GradTurbVar[0]/sigma;

  /*--- For Jacobians -> Use of TSL approx. to compute derivatives of the gradients ---*/

  if (implicit) {
    Jacobian_i[0][0] = (0.5*Proj_Mean_GradTurbVar[0]-nu_e*proj_vector_ij)/sigma;
    Jacobian_j[0][0] = (0.5*Proj_Mean_GradTurbVar[0]+nu_e*proj_vector_ij)/sigma;
  }

}


CAvgGrad_TurbSA_Neg::CAvgGrad_TurbSA_Neg(unsigned short val_nDim,
                                         unsigned short val_nVar,
                                         bool correct_grad,
                                         CConfig *config)
    : CAvgGrad_Scalar(val_nDim, val_nVar, correct_grad, config),
      sigma(2./3.), cn1(16.0), fn(0.0) {
}

CAvgGrad_TurbSA_Neg::~CAvgGrad_TurbSA_Neg(void) {
}

void CAvgGrad_TurbSA_Neg::ExtraADPreaccIn() {
}

void CAvgGrad_TurbSA_Neg::FinishResidualCalc(su2double *val_residual,
                                                   su2double **Jacobian_i,
                                                   su2double **Jacobian_j,
                                                   CConfig *config) {

  /*--- Compute mean effective viscosity ---*/

  nu_i = Laminar_Viscosity_i/Density_i;
  nu_j = Laminar_Viscosity_j/Density_j;

  nu_ij = 0.5*(nu_i+nu_j);
  nu_tilde_ij = 0.5*(TurbVar_i[0]+TurbVar_j[0]);

  Xi = nu_tilde_ij/nu_ij;

  if (nu_tilde_ij > 0.0) {
    nu_e = nu_ij + nu_tilde_ij;
  }
  else {
    fn = (cn1 + Xi*Xi*Xi)/(cn1 - Xi*Xi*Xi);
    nu_e = nu_ij + fn*nu_tilde_ij;
  }

  val_residual[0] = nu_e*Proj_Mean_GradTurbVar_Normal[0]/sigma;

  /*--- For Jacobians -> Use of TSL approx. to compute derivatives of the gradients ---*/

  if (implicit) {
    Jacobian_i[0][0] = (0.5*Proj_Mean_GradTurbVar[0]-nu_e*proj_vector_ij)/sigma;
    Jacobian_j[0][0] = (0.5*Proj_Mean_GradTurbVar[0]+nu_e*proj_vector_ij)/sigma;
  }

}

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

  diff_kine  = 0.5*(diff_i_kine + diff_j_kine);
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
