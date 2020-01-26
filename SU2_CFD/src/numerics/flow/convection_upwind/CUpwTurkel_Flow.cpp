/*!
 * \file CUpwTurkel_Flow.cpp
 * \brief Implementation of numerics class CUpwTurkel_Flow.
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

#include "../../../../include/numerics/flow/convection_upwind/CUpwTurkel_Flow.hpp"

CUpwTurkel_Flow::CUpwTurkel_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {

  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  /* A grid is defined as dynamic if there's rigid grid movement or grid deformation AND the problem is time domain */
  dynamic_grid = config->GetDynamic_Grid();

  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;

  Beta_min = config->GetminTurkelBeta();
  Beta_max = config->GetmaxTurkelBeta();

  Diff_U = new su2double [nVar];
  Velocity_i = new su2double [nDim];
  Velocity_j = new su2double [nDim];
  RoeVelocity = new su2double [nDim];
  ProjFlux_i = new su2double [nVar];
  ProjFlux_j = new su2double [nVar];
  Lambda = new su2double [nVar];
  Epsilon = new su2double [nVar];
  absPeJac = new su2double* [nVar];
  invRinvPe = new su2double* [nVar];
  R_Tensor  = new su2double* [nVar];
  Matrix    = new su2double* [nVar];
  Art_Visc  = new su2double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    absPeJac[iVar] = new su2double [nVar];
    invRinvPe[iVar] = new su2double [nVar];
    Matrix[iVar] = new su2double [nVar];
    Art_Visc[iVar] = new su2double [nVar];
    R_Tensor[iVar] = new su2double [nVar];
  }
}

CUpwTurkel_Flow::~CUpwTurkel_Flow(void) {

  delete [] Diff_U;
  delete [] Velocity_i;
  delete [] Velocity_j;
  delete [] RoeVelocity;
  delete [] ProjFlux_i;
  delete [] ProjFlux_j;
  delete [] Lambda;
  delete [] Epsilon;
  for (iVar = 0; iVar < nVar; iVar++) {
    delete [] absPeJac[iVar];
    delete [] invRinvPe[iVar];
    delete [] Matrix[iVar];
    delete [] Art_Visc[iVar];
    delete [] R_Tensor[iVar];
  }
  delete [] Matrix;
  delete [] Art_Visc;
  delete [] absPeJac;
  delete [] invRinvPe;
  delete [] R_Tensor;

}

void CUpwTurkel_Flow::ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config) {

  su2double U_i[5] = {0.0,0.0,0.0,0.0,0.0}, U_j[5] = {0.0,0.0,0.0,0.0,0.0};

  /*--- Face area (norm or the normal vector) ---*/

  Area = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    Area += Normal[iDim]*Normal[iDim];
  Area = sqrt(Area);

  /*-- Unit Normal ---*/

  for (iDim = 0; iDim < nDim; iDim++)
    UnitNormal[iDim] = Normal[iDim]/Area;

  /*--- Primitive variables at point i ---*/

  for (iDim = 0; iDim < nDim; iDim++)
    Velocity_i[iDim] = V_i[iDim+1];
  Pressure_i = V_i[nDim+1];
  Density_i = V_i[nDim+2];
  Enthalpy_i = V_i[nDim+3];
  Energy_i = Enthalpy_i - Pressure_i/Density_i;
  SoundSpeed_i = sqrt(fabs(Pressure_i*Gamma/Density_i));

  /*--- Primitive variables at point j ---*/

  for (iDim = 0; iDim < nDim; iDim++)
    Velocity_j[iDim] = V_j[iDim+1];
  Pressure_j = V_j[nDim+1];
  Density_j = V_j[nDim+2];
  Enthalpy_j = V_j[nDim+3];
  Energy_j = Enthalpy_j - Pressure_j/Density_j;
  SoundSpeed_j = sqrt(fabs(Pressure_j*Gamma/Density_j));

  /*--- Recompute conservative variables ---*/

  U_i[0] = Density_i; U_j[0] = Density_j;
  for (iDim = 0; iDim < nDim; iDim++) {
    U_i[iDim+1] = Density_i*Velocity_i[iDim]; U_j[iDim+1] = Density_j*Velocity_j[iDim];
  }
  U_i[nDim+1] = Density_i*Energy_i; U_j[nDim+1] = Density_j*Energy_j;

  /*--- Roe-averaged variables at interface between i & j ---*/

  R = sqrt(fabs(Density_j/Density_i));
  RoeDensity = R*Density_i;
  sq_vel = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    RoeVelocity[iDim] = (R*Velocity_j[iDim]+Velocity_i[iDim])/(R+1);
    sq_vel += RoeVelocity[iDim]*RoeVelocity[iDim];
  }
  RoeEnthalpy = (R*Enthalpy_j+Enthalpy_i)/(R+1);
  RoeSoundSpeed = sqrt(fabs((Gamma-1)*(RoeEnthalpy-0.5*sq_vel)));
  RoePressure = RoeDensity/Gamma*RoeSoundSpeed*RoeSoundSpeed;

  /*--- Compute ProjFlux_i ---*/
  GetInviscidProjFlux(&Density_i, Velocity_i, &Pressure_i, &Enthalpy_i, Normal, ProjFlux_i);

  /*--- Compute ProjFlux_j ---*/
  GetInviscidProjFlux(&Density_j, Velocity_j, &Pressure_j, &Enthalpy_j, Normal, ProjFlux_j);

  ProjVelocity = 0.0; ProjVelocity_i = 0.0; ProjVelocity_j = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    ProjVelocity   += RoeVelocity[iDim]*UnitNormal[iDim];
    ProjVelocity_i += Velocity_i[iDim]*UnitNormal[iDim];
    ProjVelocity_j += Velocity_j[iDim]*UnitNormal[iDim];
  }

  /*--- Projected velocity adjustment due to mesh motion ---*/
  if (dynamic_grid) {
    su2double ProjGridVel = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      ProjGridVel   += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*UnitNormal[iDim];
    }
    ProjVelocity   -= ProjGridVel;
    ProjVelocity_i -= ProjGridVel;
    ProjVelocity_j -= ProjGridVel;
  }

  /*--- First few flow eigenvalues of A.Normal with the normal---*/
  for (iDim = 0; iDim < nDim; iDim++)
    Lambda[iDim] = ProjVelocity;

  local_Mach = sqrt(sq_vel)/RoeSoundSpeed;
  Beta      = max(Beta_min, min(local_Mach, Beta_max));
  Beta2      = Beta*Beta;

  one_m_Betasqr        = 1.0 - Beta2;  // 1-Beta*Beta
  one_p_Betasqr        = 1.0 + Beta2;  // 1+Beta*Beta
  sqr_one_m_Betasqr_Lam1 = pow((one_m_Betasqr*Lambda[0]),2); // [(1-Beta^2)*Lambda[0]]^2
  sqr_two_Beta_c_Area    = pow(2.0*Beta*RoeSoundSpeed*Area,2); // [2*Beta*c*Area]^2

  /*--- The rest of the flow eigenvalues of preconditioned matrix---*/
  Lambda[nVar-2] = 0.5 * ( one_p_Betasqr*Lambda[0] + sqrt( sqr_one_m_Betasqr_Lam1 + sqr_two_Beta_c_Area));
  Lambda[nVar-1] = 0.5 * ( one_p_Betasqr*Lambda[0] - sqrt( sqr_one_m_Betasqr_Lam1 + sqr_two_Beta_c_Area));

  s_hat = 1.0/Area * (Lambda[nVar-1] - Lambda[0]*Beta2);
  r_hat = 1.0/Area * (Lambda[nVar-2] - Lambda[0]*Beta2);
  t_hat = 0.5/Area * (Lambda[nVar-1] - Lambda[nVar-2]);
  rhoB2a2 = RoeDensity*Beta2*RoeSoundSpeed*RoeSoundSpeed;

  /*--- Diference variables iPoint and jPoint and absolute value of the eigen values---*/
  for (iVar = 0; iVar < nVar; iVar++) {
    Diff_U[iVar] = U_j[iVar]-U_i[iVar];
    Lambda[iVar] = fabs(Lambda[iVar]);
  }

  /*--- Compute the absolute Preconditioned Jacobian in entropic Variables (do it with the Unitary Normal) ---*/
  GetPrecondJacobian(Beta2, r_hat, s_hat, t_hat, rhoB2a2, Lambda, UnitNormal, absPeJac);

  /*--- Compute the matrix from entropic variables to conserved variables ---*/
  GetinvRinvPe(Beta2, RoeEnthalpy, RoeSoundSpeed, RoeDensity, RoeVelocity, invRinvPe);

  /*--- Compute the matrix from entropic variables to conserved variables ---*/
  GetRMatrix(RoePressure, RoeSoundSpeed, RoeDensity, RoeVelocity, R_Tensor);

  if (implicit) {
    /*--- Jacobians of the inviscid flux, scaled by
     0.5 because val_residual ~ 0.5*(fc_i+fc_j)*Normal ---*/
    GetInviscidProjJac(Velocity_i, &Energy_i, Normal, 0.5, val_Jacobian_i);
    GetInviscidProjJac(Velocity_j, &Energy_j, Normal, 0.5, val_Jacobian_j);
  }

  for (iVar = 0; iVar < nVar; iVar ++) {
    for (jVar = 0; jVar < nVar; jVar ++) {
      Matrix[iVar][jVar] = 0.0;
      for (kVar = 0; kVar < nVar; kVar++)
        Matrix[iVar][jVar]  += absPeJac[iVar][kVar]*R_Tensor[kVar][jVar];
    }
  }

  for (iVar = 0; iVar < nVar; iVar ++) {
    for (jVar = 0; jVar < nVar; jVar ++) {
      Art_Visc[iVar][jVar] = 0.0;
      for (kVar = 0; kVar < nVar; kVar++)
        Art_Visc[iVar][jVar]  += invRinvPe[iVar][kVar]*Matrix[kVar][jVar];
    }
  }

  /*--- Roe's Flux approximation ---*/
  for (iVar = 0; iVar < nVar; iVar++) {
    val_residual[iVar] = 0.5*(ProjFlux_i[iVar]+ProjFlux_j[iVar]);
    for (jVar = 0; jVar < nVar; jVar++) {
      val_residual[iVar] -= 0.5*Art_Visc[iVar][jVar]*Diff_U[jVar];
      if (implicit) {
        val_Jacobian_i[iVar][jVar] += 0.5*Art_Visc[iVar][jVar];
        val_Jacobian_j[iVar][jVar] -= 0.5*Art_Visc[iVar][jVar];
      }
    }
  }

  /*--- Contributions due to mesh motion---*/
  if (dynamic_grid) {
    ProjVelocity = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      ProjVelocity += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*UnitNormal[iDim];
    for (iVar = 0; iVar < nVar; iVar++) {
      val_residual[iVar] -= ProjVelocity * 0.5*(U_i[iVar]+U_j[iVar]);
      /*--- Implicit terms ---*/
      if (implicit) {
        val_Jacobian_i[iVar][iVar] -= 0.5*ProjVelocity;
        val_Jacobian_j[iVar][iVar] -= 0.5*ProjVelocity;
      }
    }
  }

}
