/*!
 * \file CSourceConservative_AdjFlow.cpp
 * \brief Implementation of numerics class CSourceConservative_AdjFlow.
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

#include "../../../include/numerics/continuous_adjoint/CSourceConservative_AdjFlow.hpp"

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
