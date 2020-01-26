/*!
 * \file CAvgGrad_AdjFlow.cpp
 * \brief Implementation of numerics class CAvgGrad_AdjFlow.
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

#include "../../../include/numerics/continuous_adjoint/CAvgGrad_AdjFlow.hpp"

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
