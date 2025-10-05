/*!
 * \file roe.cpp
 * \brief Implementations of Roe-type schemes.
 * \author F. Palacios, T. Economon
 * \version 8.3.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2025, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../../../include/numerics/flow/convection/roe.hpp"
#include "../../../../../Common/include/toolboxes/geometry_toolbox.hpp"

CUpwRoeBase_Flow::CUpwRoeBase_Flow(unsigned short val_nDim, unsigned short val_nVar, unsigned short val_nPrimVar,
                                   unsigned short val_nPrimVarGrad, const CConfig* config, bool val_low_dissipation)
    : CNumerics(val_nDim, val_nVar, val_nPrimVar, val_nPrimVarGrad, config) {

  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  /* A grid is defined as dynamic if there's rigid grid movement or grid deformation AND the problem is time domain */
  dynamic_grid = config->GetDynamic_Grid();
  kappa = config->GetRoe_Kappa();

  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;

  roe_low_dissipation = val_low_dissipation;

  Flux = new su2double [nVar];
  Diff_U = new su2double [nVar];
  ProjFlux_i = new su2double [nVar];
  ProjFlux_j = new su2double [nVar];
  Conservatives_i = new su2double [nVar];
  Conservatives_j = new su2double [nVar];
  Lambda = new su2double [nVar];
  P_Tensor = new su2double* [nVar];
  invP_Tensor = new su2double* [nVar];
  Jacobian_i = new su2double* [nVar];
  Jacobian_j = new su2double* [nVar];
  for (unsigned short iVar = 0; iVar < nVar; iVar++) {
    P_Tensor[iVar] = new su2double [nVar];
    invP_Tensor[iVar] = new su2double [nVar];
    Jacobian_i[iVar] = new su2double [nVar];
    Jacobian_j[iVar] = new su2double [nVar];
  }

}

CUpwRoeBase_Flow::~CUpwRoeBase_Flow() {

  delete [] Flux;
  delete [] Diff_U;
  delete [] ProjFlux_i;
  delete [] ProjFlux_j;
  delete [] Conservatives_i;
  delete [] Conservatives_j;
  delete [] Lambda;
  for (unsigned short iVar = 0; iVar < nVar; iVar++) {
    delete [] P_Tensor[iVar];
    delete [] invP_Tensor[iVar];
    delete [] Jacobian_i[iVar];
    delete [] Jacobian_j[iVar];
  }
  delete [] P_Tensor;
  delete [] invP_Tensor;
  delete [] Jacobian_i;
  delete [] Jacobian_j;

}

void CUpwRoeBase_Flow::FinalizeResidual(su2double *val_residual, su2double **val_Jacobian_i,
                                        su2double **val_Jacobian_j, const CConfig* config) {
/*---
 CUpwRoeBase_Flow::ComputeResidual initializes the residual (flux) and its Jacobians with the standard Roe averaging
 fc_{1/2} = kappa*(fc_i+fc_j)*Normal. It then calls this method, which derived classes specialize, to account for
 the dissipation part.
---*/
}

CNumerics::ResidualType<> CUpwRoeBase_Flow::ComputeResidual(const CConfig* config) {

  implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);

  unsigned short iVar, jVar, iDim;
  su2double ProjGridVel = 0.0, Energy_i, Energy_j;

  AD::StartPreacc();
  AD::SetPreaccIn(V_i, nDim+4); AD::SetPreaccIn(V_j, nDim+4); AD::SetPreaccIn(Normal, nDim);
  if (dynamic_grid) {
    AD::SetPreaccIn(GridVel_i, nDim); AD::SetPreaccIn(GridVel_j, nDim);
  }
  if (roe_low_dissipation){
    AD::SetPreaccIn(Sensor_i); AD::SetPreaccIn(Sensor_j);
    AD::SetPreaccIn(Dissipation_i); AD::SetPreaccIn(Dissipation_j);
  }

  /*--- Face area (norm of the normal vector) and unit normal ---*/

  Area = GeometryToolbox::Norm(nDim, Normal);

  for (iDim = 0; iDim < nDim; iDim++)
    UnitNormal[iDim] = Normal[iDim]/Area;

  /*--- Primitive variables at point i ---*/

  for (iDim = 0; iDim < nDim; iDim++)
    Velocity_i[iDim] = V_i[iDim+1];
  Pressure_i = V_i[nDim+1];
  Density_i  = V_i[nDim+2];
  Enthalpy_i = V_i[nDim+3];
  Energy_i = Enthalpy_i - Pressure_i/Density_i;

  /*--- Primitive variables at point j ---*/

  for (iDim = 0; iDim < nDim; iDim++)
    Velocity_j[iDim] = V_j[iDim+1];
  Pressure_j = V_j[nDim+1];
  Density_j  = V_j[nDim+2];
  Enthalpy_j = V_j[nDim+3];
  Energy_j = Enthalpy_j - Pressure_j/Density_j;

  /*--- Compute variables that are common to the derived schemes ---*/

  /*--- Roe-averaged variables at interface between i & j ---*/

  su2double R = sqrt(fabs(Density_j/Density_i));
  RoeDensity = R*Density_i;
  su2double sq_vel = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    RoeVelocity[iDim] = (R*Velocity_j[iDim]+Velocity_i[iDim])/(R+1);
    sq_vel += RoeVelocity[iDim]*RoeVelocity[iDim];
  }
  RoeEnthalpy = (R*Enthalpy_j+Enthalpy_i)/(R+1);
  RoeSoundSpeed2 = (Gamma-1)*(RoeEnthalpy-0.5*sq_vel);

  /*--- Negative RoeSoundSpeed^2, the jump variables is too large, clear fluxes and exit. ---*/

  if (RoeSoundSpeed2 <= 0.0) {
    for (iVar = 0; iVar < nVar; iVar++) {
      Flux[iVar] = 0.0;
      if (implicit){
        for (jVar = 0; jVar < nVar; jVar++) {
          Jacobian_i[iVar][jVar] = 0.0;
          Jacobian_j[iVar][jVar] = 0.0;
        }
      }
    }
    AD::SetPreaccOut(Flux, nVar);
    AD::EndPreacc();

    return ResidualType<>(Flux, Jacobian_i, Jacobian_j);
  }

  RoeSoundSpeed = sqrt(RoeSoundSpeed2);

  /*--- P tensor ---*/

  GetPMatrix(&RoeDensity, RoeVelocity, &RoeSoundSpeed, UnitNormal, P_Tensor);

  /*--- Projected velocity adjusted for mesh motion ---*/

  ProjVelocity = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    ProjVelocity += RoeVelocity[iDim]*UnitNormal[iDim];

  if (dynamic_grid) {
    for (iDim = 0; iDim < nDim; iDim++)
      ProjGridVel += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*UnitNormal[iDim];
    ProjVelocity -= ProjGridVel;
  }

  /*--- Flow eigenvalues ---*/

  for (iDim = 0; iDim < nDim; iDim++)
    Lambda[iDim] = ProjVelocity;

  Lambda[nVar-2] = ProjVelocity + RoeSoundSpeed;
  Lambda[nVar-1] = ProjVelocity - RoeSoundSpeed;

  /*--- Apply Mavriplis' entropy correction to eigenvalues ---*/

  su2double MaxLambda = fabs(ProjVelocity) + RoeSoundSpeed;

  for (iVar = 0; iVar < nVar; iVar++)
    Lambda[iVar] = max(fabs(Lambda[iVar]), config->GetEntropyFix_Coeff()*MaxLambda);

  /*--- Reconstruct conservative variables ---*/

  Conservatives_i[0] = Density_i;
  Conservatives_j[0] = Density_j;

  for (iDim = 0; iDim < nDim; iDim++) {
    Conservatives_i[iDim+1] = Density_i*Velocity_i[iDim];
    Conservatives_j[iDim+1] = Density_j*Velocity_j[iDim];
  }
  Conservatives_i[nDim+1] = Density_i*Energy_i;
  Conservatives_j[nDim+1] = Density_j*Energy_j;

  /*--- Compute left and right fluxes ---*/

  GetInviscidProjFlux(&Density_i, Velocity_i, &Pressure_i, &Enthalpy_i, Normal, ProjFlux_i);
  GetInviscidProjFlux(&Density_j, Velocity_j, &Pressure_j, &Enthalpy_j, Normal, ProjFlux_j);

  /*--- Initialize residual (flux) and Jacobians ---*/

  for (iVar = 0; iVar < nVar; iVar++)
    Flux[iVar] = 0.5*(ProjFlux_i[iVar]+ProjFlux_j[iVar]);

  if (implicit) {
    GetInviscidProjJac(Velocity_i, &Energy_i, Normal, 0.5, Jacobian_i);
    GetInviscidProjJac(Velocity_j, &Energy_j, Normal, 0.5, Jacobian_j);
  }

  /*--- Finalize in children class ---*/

  FinalizeResidual(Flux, Jacobian_i, Jacobian_j, config);

  /*--- Correct for grid motion ---*/

  if (dynamic_grid) {
    for (iVar = 0; iVar < nVar; iVar++) {
      Flux[iVar] -= ProjGridVel*Area * 0.5*(Conservatives_i[iVar]+Conservatives_j[iVar]);

      if (implicit) {
        Jacobian_i[iVar][iVar] -= 0.5*ProjGridVel*Area;
        Jacobian_j[iVar][iVar] -= 0.5*ProjGridVel*Area;
      }
    }
  }

  AD::SetPreaccOut(Flux, nVar);
  AD::EndPreacc();

  return ResidualType<>(Flux, Jacobian_i, Jacobian_j);

}

CUpwRoe_Flow::CUpwRoe_Flow(unsigned short val_nDim, unsigned short val_nVar, unsigned short val_nPrimVar,
                           unsigned short val_nPrimVarGrad, const CConfig* config, bool val_low_dissipation)
    : CUpwRoeBase_Flow(val_nDim, val_nVar, val_nPrimVar, val_nPrimVarGrad, config, val_low_dissipation) {}

void CUpwRoe_Flow::FinalizeResidual(su2double *val_residual, su2double **val_Jacobian_i,
                                    su2double **val_Jacobian_j, const CConfig* config) {

  unsigned short iVar, jVar, kVar;

  /*--- Compute inverse P tensor ---*/
  GetPMatrix_inv(&RoeDensity, RoeVelocity, &RoeSoundSpeed, UnitNormal, invP_Tensor);

  /*--- Diference between conservative variables at jPoint and iPoint ---*/
  for (iVar = 0; iVar < nVar; iVar++)
    Diff_U[iVar] = Conservatives_j[iVar]-Conservatives_i[iVar];

  /*--- Low dissipation formulation ---*/
  if (roe_low_dissipation)
    Dissipation_ij = GetRoe_Dissipation(Dissipation_i, Dissipation_j, Sensor_i, Sensor_j, config);
  else
    Dissipation_ij = 1.0;

  /*--- Standard Roe "dissipation" ---*/

  for (iVar = 0; iVar < nVar; iVar++) {
    for (jVar = 0; jVar < nVar; jVar++) {
      /*--- Compute |Proj_ModJac_Tensor| = P x |Lambda| x inverse P ---*/
      su2double Proj_ModJac_Tensor_ij = 0.0;
      for (kVar = 0; kVar < nVar; kVar++)
        Proj_ModJac_Tensor_ij += P_Tensor[iVar][kVar]*Lambda[kVar]*invP_Tensor[kVar][jVar];

      /*--- Update residual and Jacobians ---*/
      val_residual[iVar] -= (1.0-kappa)*Proj_ModJac_Tensor_ij*Diff_U[jVar]*Area*Dissipation_ij;

      if(implicit){
        val_Jacobian_i[iVar][jVar] += (1.0-kappa)*Proj_ModJac_Tensor_ij*Area;
        val_Jacobian_j[iVar][jVar] -= (1.0-kappa)*Proj_ModJac_Tensor_ij*Area;
      }
    }
  }

  CorrectResidual(val_residual);

}

CUpwL2Roe_Flow::CUpwL2Roe_Flow(unsigned short val_nDim, unsigned short val_nVar, unsigned short val_nPrimVar,
                               unsigned short val_nPrimVarGrad, const CConfig* config)
    : CUpwRoeBase_Flow(val_nDim, val_nVar, val_nPrimVar, val_nPrimVarGrad, config, false) {}

void CUpwL2Roe_Flow::FinalizeResidual(su2double *val_residual, su2double **val_Jacobian_i,
                                      su2double **val_Jacobian_j, const CConfig* config) {

  /*--- L2Roe: a low dissipation version of Roe's approximate Riemann solver for low Mach numbers. IJNMF 2015 ---*/

  unsigned short iVar, jVar, kVar, iDim;

  /*--- Clamped Mach number ---*/

  su2double M_i = 0.0, M_j = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    M_i += Velocity_i[iDim]*Velocity_i[iDim];
    M_j += Velocity_j[iDim]*Velocity_j[iDim];
  }
  M_i = sqrt(M_i / fabs(Pressure_i*Gamma/Density_i));
  M_j = sqrt(M_j / fabs(Pressure_j*Gamma/Density_j));

  su2double zeta = max(0.05,min(max(M_i,M_j),1.0));

  /*--- Compute wave amplitudes (characteristics) ---*/

  su2double proj_delta_vel = 0.0, delta_vel[3] = {0.0};
  for (iDim = 0; iDim < nDim; iDim++) {
    delta_vel[iDim] = Velocity_j[iDim] - Velocity_i[iDim];
    proj_delta_vel += delta_vel[iDim]*UnitNormal[iDim];
  }
  proj_delta_vel *= zeta;
  su2double delta_p = Pressure_j - Pressure_i;
  su2double delta_rho = Density_j - Density_i;

  su2double delta_wave[5] = {0.0};
  if (nDim == 2) {
    delta_wave[0] = delta_rho - delta_p/RoeSoundSpeed2;
    delta_wave[1] = (UnitNormal[1]*delta_vel[0]-UnitNormal[0]*delta_vel[1])*zeta;
    delta_wave[2] = proj_delta_vel + delta_p/(RoeDensity*RoeSoundSpeed);
    delta_wave[3] = -proj_delta_vel + delta_p/(RoeDensity*RoeSoundSpeed);
  } else {
    delta_wave[0] = delta_rho - delta_p/RoeSoundSpeed2;
    delta_wave[1] = (UnitNormal[0]*delta_vel[2]-UnitNormal[2]*delta_vel[0])*zeta;
    delta_wave[2] = (UnitNormal[1]*delta_vel[0]-UnitNormal[0]*delta_vel[1])*zeta;
    delta_wave[3] = proj_delta_vel + delta_p/(RoeDensity*RoeSoundSpeed);
    delta_wave[4] = -proj_delta_vel + delta_p/(RoeDensity*RoeSoundSpeed);
  }

  /*--- Update residual ---*/

  for (iVar = 0; iVar < nVar; iVar++)
    for (kVar = 0; kVar < nVar; kVar++)
      val_residual[iVar] -= (1.0-kappa)*Lambda[kVar]*delta_wave[kVar]*P_Tensor[iVar][kVar]*Area;

  if (!implicit) return;

  /*--- If implicit use the Jacobians of the standard Roe scheme as an approximation ---*/

  GetPMatrix_inv(&RoeDensity, RoeVelocity, &RoeSoundSpeed, UnitNormal, invP_Tensor);

  for (iVar = 0; iVar < nVar; iVar++) {
    for (jVar = 0; jVar < nVar; jVar++) {
      /*--- Compute |Proj_ModJac_Tensor| = P x |Lambda| x inverse P ---*/
      su2double Proj_ModJac_Tensor_ij = 0.0;
      for (kVar = 0; kVar < nVar; kVar++)
        Proj_ModJac_Tensor_ij += P_Tensor[iVar][kVar]*Lambda[kVar]*invP_Tensor[kVar][jVar];

      val_Jacobian_i[iVar][jVar] += (1.0-kappa)*Proj_ModJac_Tensor_ij*Area;
      val_Jacobian_j[iVar][jVar] -= (1.0-kappa)*Proj_ModJac_Tensor_ij*Area;
    }
  }

}

CUpwLMRoe_Flow::CUpwLMRoe_Flow(unsigned short val_nDim, unsigned short val_nVar, unsigned short val_nPrimVar,
                               unsigned short val_nPrimVarGrad, const CConfig* config) : CUpwRoeBase_Flow(val_nDim, val_nVar, val_nPrimVar, val_nPrimVarGrad, config, false) {}

void CUpwLMRoe_Flow::FinalizeResidual(su2double *val_residual, su2double **val_Jacobian_i,
                                      su2double **val_Jacobian_j, const CConfig* config) {

  /*--- Rieper, A low-Mach number fix for Roe's approximate Riemman Solver, JCP 2011 ---*/

  unsigned short iVar, jVar, kVar, iDim;

  /*--- Clamped Mach number ---*/

  su2double M_i = 0.0, M_j = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    M_i += Velocity_i[iDim]*Velocity_i[iDim];
    M_j += Velocity_j[iDim]*Velocity_j[iDim];
  }
  M_i = sqrt(M_i / fabs(Pressure_i*Gamma/Density_i));
  M_j = sqrt(M_j / fabs(Pressure_j*Gamma/Density_j));

  su2double zeta = max(0.05,min(max(M_i,M_j),1.0));

  /*--- Compute wave amplitudes (characteristics) ---*/

  su2double proj_delta_vel = 0.0, delta_vel[3] = {0.0};
  for (iDim = 0; iDim < nDim; iDim++) {
    delta_vel[iDim] = Velocity_j[iDim] - Velocity_i[iDim];
    proj_delta_vel += delta_vel[iDim]*UnitNormal[iDim];
  }
  proj_delta_vel *= zeta;
  su2double delta_p = Pressure_j - Pressure_i;
  su2double delta_rho = Density_j - Density_i;

  su2double delta_wave[5] = {0.0};
  if (nDim == 2) {
    delta_wave[0] = delta_rho - delta_p/RoeSoundSpeed2;
    delta_wave[1] = (UnitNormal[1]*delta_vel[0]-UnitNormal[0]*delta_vel[1]);
    delta_wave[2] = proj_delta_vel + delta_p/(RoeDensity*RoeSoundSpeed);
    delta_wave[3] = -proj_delta_vel + delta_p/(RoeDensity*RoeSoundSpeed);
  } else {
    delta_wave[0] = delta_rho - delta_p/RoeSoundSpeed2;
    delta_wave[1] = (UnitNormal[0]*delta_vel[2]-UnitNormal[2]*delta_vel[0]);
    delta_wave[2] = (UnitNormal[1]*delta_vel[0]-UnitNormal[0]*delta_vel[1]);
    delta_wave[3] = proj_delta_vel + delta_p/(RoeDensity*RoeSoundSpeed);
    delta_wave[4] = -proj_delta_vel + delta_p/(RoeDensity*RoeSoundSpeed);
  }

  /*--- Update residual ---*/

  for (iVar = 0; iVar < nVar; iVar++)
    for (kVar = 0; kVar < nVar; kVar++)
      val_residual[iVar] -= (1.0-kappa)*Lambda[kVar]*delta_wave[kVar]*P_Tensor[iVar][kVar]*Area;

  if (!implicit) return;

  /*--- If implicit use the Jacobians of the standard Roe scheme as an approximation ---*/

  GetPMatrix_inv(&RoeDensity, RoeVelocity, &RoeSoundSpeed, UnitNormal, invP_Tensor);

  for (iVar = 0; iVar < nVar; iVar++) {
    for (jVar = 0; jVar < nVar; jVar++) {
      /*--- Compute |Proj_ModJac_Tensor| = P x |Lambda| x inverse P ---*/
      su2double Proj_ModJac_Tensor_ij = 0.0;
      for (kVar = 0; kVar < nVar; kVar++)
        Proj_ModJac_Tensor_ij += P_Tensor[iVar][kVar]*Lambda[kVar]*invP_Tensor[kVar][jVar];

      val_Jacobian_i[iVar][jVar] += (1.0-kappa)*Proj_ModJac_Tensor_ij*Area;
      val_Jacobian_j[iVar][jVar] -= (1.0-kappa)*Proj_ModJac_Tensor_ij*Area;
    }
  }

}

CUpwTurkel_Flow::CUpwTurkel_Flow(unsigned short val_nDim, unsigned short val_nVar, unsigned short val_nPrimVar,
                                 unsigned short val_nPrimVarGrad, const CConfig* config) : CNumerics(val_nDim, val_nVar, val_nPrimVar, val_nPrimVarGrad, config) {

  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  /* A grid is defined as dynamic if there's rigid grid movement or grid deformation AND the problem is time domain */
  dynamic_grid = config->GetDynamic_Grid();

  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;

  Beta_min = config->GetminTurkelBeta();
  Beta_max = config->GetmaxTurkelBeta();

  Flux = new su2double [nVar];
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
  Jacobian_i = new su2double* [nVar];
  Jacobian_j = new su2double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    absPeJac[iVar] = new su2double [nVar];
    invRinvPe[iVar] = new su2double [nVar];
    Matrix[iVar] = new su2double [nVar];
    Art_Visc[iVar] = new su2double [nVar];
    R_Tensor[iVar] = new su2double [nVar];
    Jacobian_i[iVar] = new su2double [nVar];
    Jacobian_j[iVar] = new su2double [nVar];
  }
}

CUpwTurkel_Flow::~CUpwTurkel_Flow() {

  delete [] Flux;
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
    delete [] Jacobian_i[iVar];
    delete [] Jacobian_j[iVar];
  }
  delete [] Matrix;
  delete [] Art_Visc;
  delete [] absPeJac;
  delete [] invRinvPe;
  delete [] R_Tensor;
  delete [] Jacobian_i;
  delete [] Jacobian_j;

}

CNumerics::ResidualType<> CUpwTurkel_Flow::ComputeResidual(const CConfig* config) {

  implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);

  su2double U_i[5] = {0.0}, U_j[5] = {0.0};

  /*--- Face area (norm or the normal vector) ---*/

  Area = GeometryToolbox::Norm(nDim, Normal);

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
  Beta       = max(Beta_min, min(local_Mach, Beta_max));
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
     0.5 because Flux ~ 0.5*(fc_i+fc_j)*Normal ---*/
    GetInviscidProjJac(Velocity_i, &Energy_i, Normal, 0.5, Jacobian_i);
    GetInviscidProjJac(Velocity_j, &Energy_j, Normal, 0.5, Jacobian_j);
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
    Flux[iVar] = 0.5*(ProjFlux_i[iVar]+ProjFlux_j[iVar]);
    for (jVar = 0; jVar < nVar; jVar++) {
      Flux[iVar] -= 0.5*Art_Visc[iVar][jVar]*Diff_U[jVar];
      if (implicit) {
        Jacobian_i[iVar][jVar] += 0.5*Art_Visc[iVar][jVar];
        Jacobian_j[iVar][jVar] -= 0.5*Art_Visc[iVar][jVar];
      }
    }
  }

  /*--- Contributions due to mesh motion---*/
  if (dynamic_grid) {
    ProjVelocity = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      ProjVelocity += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*UnitNormal[iDim];
    for (iVar = 0; iVar < nVar; iVar++) {
      Flux[iVar] -= ProjVelocity * 0.5*(U_i[iVar]+U_j[iVar]);
      /*--- Implicit terms ---*/
      if (implicit) {
        Jacobian_i[iVar][iVar] -= 0.5*ProjVelocity;
        Jacobian_j[iVar][iVar] -= 0.5*ProjVelocity;
      }
    }
  }

  return ResidualType<>(Flux, Jacobian_i, Jacobian_j);

}

CUpwGeneralRoe_Flow::CUpwGeneralRoe_Flow(unsigned short val_nDim, unsigned short val_nVar, unsigned short val_nPrimVar,
                                         unsigned short val_nPrimVarGrad, const CConfig* config) : CNumerics(val_nDim, val_nVar, val_nPrimVar, val_nPrimVarGrad, config) {

  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  /* A grid is defined as dynamic if there's rigid grid movement or grid deformation AND the problem is time domain */
  dynamic_grid = config->GetDynamic_Grid();
  kappa = config->GetRoe_Kappa();


  Flux = new su2double [nVar];
  Diff_U = new su2double [nVar];
  Velocity_i = new su2double [nDim];
  Velocity_j = new su2double [nDim];
  RoeVelocity = new su2double [nDim];
  delta_vel  = new su2double [nDim];
  delta_wave = new su2double [nVar];
  ProjFlux_i = new su2double [nVar];
  ProjFlux_j = new su2double [nVar];
  Lambda = new su2double [nVar];
  Epsilon = new su2double [nVar];
  P_Tensor = new su2double* [nVar];
  invP_Tensor = new su2double* [nVar];
  Jacobian_i = new su2double* [nVar];
  Jacobian_j = new su2double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    P_Tensor[iVar] = new su2double [nVar];
    invP_Tensor[iVar] = new su2double [nVar];
    Jacobian_i[iVar] = new su2double [nVar];
    Jacobian_j[iVar] = new su2double [nVar];
  }
}

CUpwGeneralRoe_Flow::~CUpwGeneralRoe_Flow() {

  delete [] Flux;
  delete [] Diff_U;
  delete [] Velocity_i;
  delete [] Velocity_j;
  delete [] RoeVelocity;
  delete [] delta_vel;
  delete [] delta_wave;
  delete [] ProjFlux_i;
  delete [] ProjFlux_j;
  delete [] Lambda;
  delete [] Epsilon;
  for (iVar = 0; iVar < nVar; iVar++) {
    delete [] P_Tensor[iVar];
    delete [] invP_Tensor[iVar];
    delete [] Jacobian_i[iVar];
    delete [] Jacobian_j[iVar];
  }
  delete [] P_Tensor;
  delete [] invP_Tensor;
  delete [] Jacobian_i;
  delete [] Jacobian_j;

}

CNumerics::ResidualType<> CUpwGeneralRoe_Flow::ComputeResidual(const CConfig* config) {

  AD::StartPreacc();
  AD::SetPreaccIn(V_i, nDim+4); AD::SetPreaccIn(V_j, nDim+4); AD::SetPreaccIn(Normal, nDim);
  AD::SetPreaccIn(S_i, 2); AD::SetPreaccIn(S_j, 2);
  if (dynamic_grid) {
    AD::SetPreaccIn(GridVel_i, nDim); AD::SetPreaccIn(GridVel_j, nDim);
  }
  su2double U_i[5] = {0.0,0.0,0.0,0.0,0.0}, U_j[5] = {0.0,0.0,0.0,0.0,0.0};

  /*--- Face area (norm or the normal vector) ---*/

  Area = GeometryToolbox::Norm(nDim, Normal);

  /*-- Unit Normal ---*/

  for (iDim = 0; iDim < nDim; iDim++)
    UnitNormal[iDim] = Normal[iDim]/Area;

  /*--- Primitive variables at point i ---*/

  Velocity2_i = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Velocity_i[iDim] = V_i[iDim+1];
    Velocity2_i += Velocity_i[iDim]*Velocity_i[iDim];
  }

  Pressure_i = V_i[nDim+1];
  Density_i = V_i[nDim+2];
  Enthalpy_i = V_i[nDim+3];
  Energy_i = Enthalpy_i - Pressure_i/Density_i;
  StaticEnthalpy_i = Enthalpy_i - 0.5*Velocity2_i;
  StaticEnergy_i = StaticEnthalpy_i - Pressure_i/Density_i;

  Kappa_i = S_i[1]/Density_i;
  Chi_i = S_i[0] - Kappa_i*StaticEnergy_i;
  SoundSpeed_i = sqrt(Chi_i + StaticEnthalpy_i*Kappa_i);

  /*--- Primitive variables at point j ---*/

  Velocity2_j = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Velocity_j[iDim] = V_j[iDim+1];
    Velocity2_j += Velocity_j[iDim]*Velocity_j[iDim];
  }

  Pressure_j = V_j[nDim+1];
  Density_j = V_j[nDim+2];
  Enthalpy_j = V_j[nDim+3];
  Energy_j = Enthalpy_j - Pressure_j/Density_j;

  StaticEnthalpy_j = Enthalpy_j - 0.5*Velocity2_j;
  StaticEnergy_j = StaticEnthalpy_j - Pressure_j/Density_j;

  Kappa_j = S_j[1]/Density_j;
  Chi_j = S_j[0] - Kappa_j*StaticEnergy_j;
  SoundSpeed_j = sqrt(Chi_j + StaticEnthalpy_j*Kappa_j);

  /*--- Recompute conservative variables ---*/

  U_i[0] = Density_i; U_j[0] = Density_j;
  for (iDim = 0; iDim < nDim; iDim++) {
    U_i[iDim+1] = Density_i*Velocity_i[iDim]; U_j[iDim+1] = Density_j*Velocity_j[iDim];
  }
  U_i[nDim+1] = Density_i*Energy_i; U_j[nDim+1] = Density_j*Energy_j;

  /*--- Roe-averaged variables at interface between i & j ---*/

  ComputeRoeAverage();

  if (RoeSoundSpeed2 <= 0.0) {
    for (iVar = 0; iVar < nVar; iVar++) {
      Flux[iVar] = 0.0;
      for (jVar = 0; jVar < nVar; jVar++) {
      Jacobian_i[iVar][iVar] = 0.0;
      Jacobian_j[iVar][iVar] = 0.0;
      }
    }
    AD::SetPreaccOut(Flux, nVar);
    AD::EndPreacc();

    return ResidualType<>(Flux, Jacobian_i, Jacobian_j);
  }

  RoeSoundSpeed = sqrt(RoeSoundSpeed2);

  /*--- Compute ProjFlux_i ---*/
  GetInviscidProjFlux(&Density_i, Velocity_i, &Pressure_i, &Enthalpy_i, Normal, ProjFlux_i);

  /*--- Compute ProjFlux_j ---*/
  GetInviscidProjFlux(&Density_j, Velocity_j, &Pressure_j, &Enthalpy_j, Normal, ProjFlux_j);

  /*--- Compute P and Lambda (do it with the Normal) ---*/

  GetPMatrix(&RoeDensity, RoeVelocity, &RoeSoundSpeed, &RoeEnthalpy, &RoeChi, &RoeKappa, UnitNormal, P_Tensor);

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

  /*--- Flow eigenvalues and entropy correctors ---*/
  for (iDim = 0; iDim < nDim; iDim++)
    Lambda[iDim] = ProjVelocity;

  Lambda[nVar-2] = ProjVelocity + RoeSoundSpeed;
  Lambda[nVar-1] = ProjVelocity - RoeSoundSpeed;

  /*--- Compute absolute value with Mavriplis' entropy correction ---*/

  MaxLambda = fabs(ProjVelocity) + RoeSoundSpeed;
  Delta = config->GetEntropyFix_Coeff();

  for (iVar = 0; iVar < nVar; iVar++) {
    Lambda[iVar] = max(fabs(Lambda[iVar]), Delta*MaxLambda);
   }

//  /*--- Harten and Hyman (1983) entropy correction ---*/
//  for (iDim = 0; iDim < nDim; iDim++)
//    Epsilon[iDim] = 4.0*max(0.0, max(Lambda[iDim]-ProjVelocity_i, ProjVelocity_j-Lambda[iDim]));
//
//  Epsilon[nVar-2] = 4.0*max(0.0, max(Lambda[nVar-2]-(ProjVelocity_i+SoundSpeed_i),(ProjVelocity_j+SoundSpeed_j)-Lambda[nVar-2]));
//  Epsilon[nVar-1] = 4.0*max(0.0, max(Lambda[nVar-1]-(ProjVelocity_i-SoundSpeed_i),(ProjVelocity_j-SoundSpeed_j)-Lambda[nVar-1]));
//
//  for (iVar = 0; iVar < nVar; iVar++)
//    if ( fabs(Lambda[iVar]) < Epsilon[iVar] )
//      Lambda[iVar] = (Lambda[iVar]*Lambda[iVar] + Epsilon[iVar]*Epsilon[iVar])/(2.0*Epsilon[iVar]);
//    else
//      Lambda[iVar] = fabs(Lambda[iVar]);

//  for (iVar = 0; iVar < nVar; iVar++)
//    Lambda[iVar] = fabs(Lambda[iVar]);

  if (!implicit) {

    /*--- Compute wave amplitudes (characteristics) ---*/
    proj_delta_vel = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      delta_vel[iDim] = Velocity_j[iDim] - Velocity_i[iDim];
      proj_delta_vel += delta_vel[iDim]*Normal[iDim];
    }
    delta_p = Pressure_j - Pressure_i;
    delta_rho = Density_j - Density_i;
    proj_delta_vel = proj_delta_vel/Area;

    if (nDim == 2) {
      delta_wave[0] = delta_rho - delta_p/(RoeSoundSpeed*RoeSoundSpeed);
      delta_wave[1] = UnitNormal[1]*delta_vel[0]-UnitNormal[0]*delta_vel[1];
      delta_wave[2] = proj_delta_vel + delta_p/(RoeDensity*RoeSoundSpeed);
      delta_wave[3] = -proj_delta_vel + delta_p/(RoeDensity*RoeSoundSpeed);
    } else {
      delta_wave[0] = delta_rho - delta_p/(RoeSoundSpeed*RoeSoundSpeed);
      delta_wave[1] = UnitNormal[0]*delta_vel[2]-UnitNormal[2]*delta_vel[0];
      delta_wave[2] = UnitNormal[1]*delta_vel[0]-UnitNormal[0]*delta_vel[1];
      delta_wave[3] = proj_delta_vel + delta_p/(RoeDensity*RoeSoundSpeed);
      delta_wave[4] = -proj_delta_vel + delta_p/(RoeDensity*RoeSoundSpeed);
    }

    /*--- Roe's Flux approximation ---*/
    for (iVar = 0; iVar < nVar; iVar++) {
      Flux[iVar] = 0.5*(ProjFlux_i[iVar]+ProjFlux_j[iVar]);
      for (jVar = 0; jVar < nVar; jVar++)
        Flux[iVar] -= 0.5*Lambda[jVar]*delta_wave[jVar]*P_Tensor[iVar][jVar]*Area;
    }

    /*--- Flux contribution due to grid motion ---*/
    if (dynamic_grid) {
      ProjVelocity = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        ProjVelocity += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*Normal[iDim];
      for (iVar = 0; iVar < nVar; iVar++) {
        Flux[iVar] -= ProjVelocity * 0.5*(U_i[iVar]+U_j[iVar]);
      }
    }
  }
  else {

    /*--- Compute inverse P ---*/

    GetPMatrix_inv(invP_Tensor, &RoeDensity, RoeVelocity, &RoeSoundSpeed, &RoeChi , &RoeKappa, UnitNormal);

     /*--- Jacobians of the inviscid flux, scaled by
      0.5 because val_resconv ~ 0.5*(fc_i+fc_j)*Normal ---*/

    GetInviscidProjJac(Velocity_i, &Enthalpy_i, &Chi_i, &Kappa_i, Normal, 0.5, Jacobian_i);

    GetInviscidProjJac(Velocity_j, &Enthalpy_j, &Chi_j, &Kappa_j, Normal, 0.5, Jacobian_j);


    /*--- Diference variables iPoint and jPoint ---*/
    for (iVar = 0; iVar < nVar; iVar++)
      Diff_U[iVar] = U_j[iVar]-U_i[iVar];

    /*--- Roe's Flux approximation ---*/
    for (iVar = 0; iVar < nVar; iVar++) {
      Flux[iVar] = 0.5*(ProjFlux_i[iVar]+ProjFlux_j[iVar]);
      for (jVar = 0; jVar < nVar; jVar++) {
        Proj_ModJac_Tensor_ij = 0.0;

        /*--- Compute |Proj_ModJac_Tensor| = P x |Lambda| x inverse P ---*/

        for (kVar = 0; kVar < nVar; kVar++)
          Proj_ModJac_Tensor_ij += P_Tensor[iVar][kVar]*Lambda[kVar]*invP_Tensor[kVar][jVar];

        Flux[iVar] -= (1.0-kappa)*Proj_ModJac_Tensor_ij*Diff_U[jVar]*Area;
        Jacobian_i[iVar][jVar] += (1.0-kappa)*Proj_ModJac_Tensor_ij*Area;
        Jacobian_j[iVar][jVar] -= (1.0-kappa)*Proj_ModJac_Tensor_ij*Area;
      }
    }

    /*--- Jacobian contributions due to grid motion ---*/
    if (dynamic_grid) {
      ProjVelocity = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        ProjVelocity += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*Normal[iDim];
      for (iVar = 0; iVar < nVar; iVar++) {
        Flux[iVar] -= ProjVelocity * 0.5*(U_i[iVar]+U_j[iVar]);
        /*--- Implicit terms ---*/
        Jacobian_i[iVar][iVar] -= 0.5*ProjVelocity;
        Jacobian_j[iVar][iVar] -= 0.5*ProjVelocity;
      }
    }

  }

  AD::SetPreaccOut(Flux, nVar);
  AD::EndPreacc();

  return ResidualType<>(Flux, Jacobian_i, Jacobian_j);

}

void CUpwGeneralRoe_Flow::ComputeRoeAverage() {

  //su2double delta_rhoStaticEnergy, err_P, s, D;
  // su2double tol = 10-6;

  R = sqrt(fabs(Density_j/Density_i));
  RoeDensity = R*Density_i;
  sq_vel = 0;  for (iDim = 0; iDim < nDim; iDim++) {
    RoeVelocity[iDim] = (R*Velocity_j[iDim]+Velocity_i[iDim])/(R+1);
    sq_vel += RoeVelocity[iDim]*RoeVelocity[iDim];
  }

  RoeEnthalpy = (R*Enthalpy_j+Enthalpy_i)/(R+1);
  delta_rho = Density_j - Density_i;
  delta_p = Pressure_j - Pressure_i;
  RoeKappa = 0.5*(Kappa_i + Kappa_j);
  RoeKappa = (Kappa_i + Kappa_j + 4*RoeKappa)/6;
  RoeChi = 0.5*(Chi_i + Chi_j);
  RoeChi = (Chi_i + Chi_j + 4*RoeChi)/6;

//  RoeKappaStaticEnthalpy = 0.5*(StaticEnthalpy_i*Kappa_i + StaticEnthalpy_j*Kappa_j);
//  RoeKappaStaticEnthalpy = (StaticEnthalpy_i*Kappa_i + StaticEnthalpy_j*Kappa_j + 4*RoeKappaStaticEnthalpy)/6;
//  s = RoeChi + RoeKappaStaticEnthalpy;
//  D = s*s*delta_rho*delta_rho + delta_p*delta_p;
//  delta_rhoStaticEnergy = Density_j*StaticEnergy_j - Density_i*StaticEnergy_i;
//  err_P = delta_p - RoeChi*delta_rho - RoeKappa*delta_rhoStaticEnergy;
//
//
//  if (abs((D - delta_p*err_P)/Density_i)>1e-3 && abs(delta_rho/Density_i)>1e-3 && s/Density_i > 1e-3) {
//
//    RoeKappa = (D*RoeKappa)/(D - delta_p*err_P);
//    RoeChi = (D*RoeChi+ s*s*delta_rho*err_P)/(D - delta_p*err_P);
//
//  }

  RoeSoundSpeed2 = RoeChi + RoeKappa*(RoeEnthalpy-0.5*sq_vel);

}


CUpwRoeBase_FlowNew::CUpwRoeBase_FlowNew(unsigned short val_nDim, unsigned short val_nVar, unsigned short val_nPrimVar,
                                         unsigned short val_nPrimVarGrad, const CConfig* config, bool val_low_dissipation) :
                                                                                                                             CNumerics(val_nDim, val_nVar, val_nPrimVar, val_nPrimVarGrad, config) {

  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  /* A grid is defined as dynamic if there's rigid grid movement or grid deformation AND the problem is time domain */
  dynamic_grid = config->GetDynamic_Grid();
  kappa = config->GetRoe_Kappa(); // 1 is unstable

  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;

  roe_low_dissipation = val_low_dissipation;

  Flux = new su2double [nVar];
  Diff_U = new su2double [nVar];
  ProjFlux_i = new su2double [nVar];
  ProjFlux_j = new su2double [nVar];
  Conservatives_i = new su2double [nVar];
  Conservatives_j = new su2double [nVar];
  Lambda = new su2double [nVar];
  P_Tensor = new su2double* [nVar];
  invP_Tensor = new su2double* [nVar];
  Jacobian_i = new su2double* [nVar];
  Jacobian_j = new su2double* [nVar];
  for (unsigned short iVar = 0; iVar < nVar; iVar++) {
    P_Tensor[iVar] = new su2double [nVar];
    invP_Tensor[iVar] = new su2double [nVar];
    Jacobian_i[iVar] = new su2double [nVar];
    Jacobian_j[iVar] = new su2double [nVar];
  }

}

CUpwRoeBase_FlowNew::~CUpwRoeBase_FlowNew() {

  delete [] Flux;
  delete [] Diff_U;
  delete [] ProjFlux_i;
  delete [] ProjFlux_j;
  delete [] Conservatives_i;
  delete [] Conservatives_j;
  delete [] Lambda;
  for (unsigned short iVar = 0; iVar < nVar; iVar++) {
    delete [] P_Tensor[iVar];
    delete [] invP_Tensor[iVar];
    delete [] Jacobian_i[iVar];
    delete [] Jacobian_j[iVar];
  }
  delete [] P_Tensor;
  delete [] invP_Tensor;
  delete [] Jacobian_i;
  delete [] Jacobian_j;

}

void CUpwRoeBase_FlowNew::FinalizeResidual(su2double *val_residual, su2double **val_Jacobian_i,
                                           su2double **val_Jacobian_j, const CConfig* config) {
  /*---
   CUpwRoeBase_Flow::ComputeResidual initializes the residual (flux) and its Jacobians with the standard Roe averaging
   fc_{1/2} = kappa*(fc_i+fc_j)*Normal. It then calls this method, which derived classes specialize, to account for
   the dissipation part.
  ---*/
}

CNumerics::ResidualType<> CUpwRoeBase_FlowNew::ComputeResidual(const CConfig* config) {

  implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);

  unsigned short iVar, jVar, iDim, jDim;
  su2double ProjGridVel = 0.0, Energy_i, Energy_j;

  AD::StartPreacc();
  AD::SetPreaccIn(V_i, nDim+4); AD::SetPreaccIn(V_j, nDim+4); AD::SetPreaccIn(Normal, nDim);
  if (dynamic_grid) {
    AD::SetPreaccIn(GridVel_i, nDim); AD::SetPreaccIn(GridVel_j, nDim);
  }
  if (roe_low_dissipation){
    AD::SetPreaccIn(Sensor_i); AD::SetPreaccIn(Sensor_j);
    AD::SetPreaccIn(Dissipation_i); AD::SetPreaccIn(Dissipation_j);
  }

  /*--- Face area (norm or the normal vector) and unit normal ---*/

  Area = GeometryToolbox::Norm(nDim, Normal);

  for (iDim = 0; iDim < nDim; iDim++)
    UnitNormal[iDim] = Normal[iDim]/Area;

  /*--- Primitive variables at point i ---*/

  for (iDim = 0; iDim < nDim; iDim++)
    Velocity_i[iDim] = V_i[iDim+1];
  Pressure_i = V_i[nDim+1];
  Density_i  = V_i[nDim+2];
  Enthalpy_i = V_i[nDim+3];
  Energy_i = Enthalpy_i - Pressure_i/Density_i;

  /*--- Primitive variables at point j ---*/

  for (iDim = 0; iDim < nDim; iDim++)
    Velocity_j[iDim] = V_j[iDim+1];
  Pressure_j = V_j[nDim+1];
  Density_j  = V_j[nDim+2];
  Enthalpy_j = V_j[nDim+3];
  Energy_j = Enthalpy_j - Pressure_j/Density_j;

  /*--- Compute variables that are common to the derived schemes ---*/

  /*--- Roe-averaged variables at interface between i & j ---*/

  su2double RR = sqrt(fabs(Density_j/Density_i));
  RoeDensity = RR*Density_i;
  su2double sq_vel = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    RoeVelocity[iDim] = (RR*Velocity_j[iDim]+Velocity_i[iDim])/(RR+1);
    sq_vel += RoeVelocity[iDim]*RoeVelocity[iDim];
  }
  RoeEnthalpy = (RR*Enthalpy_j+Enthalpy_i)/(RR+1);
  RoeSoundSpeed2 = (Gamma-1)*(RoeEnthalpy-0.5*sq_vel);

  /*--- Negative RoeSoundSpeed^2, the jump variables is too large, clear fluxes and exit. ---*/

  if (RoeSoundSpeed2 <= 0.0) {
    for (iVar = 0; iVar < nVar; iVar++) {
      Flux[iVar] = 0.0;
      if (implicit){
        for (jVar = 0; jVar < nVar; jVar++) {
          Jacobian_i[iVar][jVar] = 0.0;
          Jacobian_j[iVar][jVar] = 0.0;
        }
      }
    }
    AD::SetPreaccOut(Flux, nVar);
    AD::EndPreacc();

    return ResidualType<>(Flux, Jacobian_i, Jacobian_j);
  }

  RoeSoundSpeed = sqrt(RoeSoundSpeed2);

  /*--- Projected velocity adjusted for mesh motion ---*/

  ProjVelocity = GeometryToolbox::DotProduct(nDim,RoeVelocity,UnitNormal);

  if (dynamic_grid) {
    for (iDim = 0; iDim < nDim; iDim++)
      ProjGridVel += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*UnitNormal[iDim];
    ProjVelocity -= ProjGridVel;
  }

  /*--- Reconstruct conservative variables ---*/

  Conservatives_i[0] = Density_i;
  Conservatives_j[0] = Density_j;

  for (iDim = 0; iDim < nDim; iDim++) {
    Conservatives_i[iDim+1] = Density_i*Velocity_i[iDim];
    Conservatives_j[iDim+1] = Density_j*Velocity_j[iDim];
  }
  Conservatives_i[nDim+1] = Density_i*Energy_i;
  Conservatives_j[nDim+1] = Density_j*Energy_j;



  const auto* vel = RoeVelocity;
  const auto& rho = RoeDensity;
  const auto& H = RoeEnthalpy;
  const auto& qn = ProjVelocity;
  const auto& a = RoeSoundSpeed;
  const auto& a2 = RoeSoundSpeed2;
  const auto& q2 = GeometryToolbox::SquaredNorm(nDim,RoeVelocity);


  const auto* velL = Velocity_i;
  const auto& rhoL = Density_i;
  const auto& pL = Pressure_i;
  const auto& HL = Enthalpy_i;
  const auto& aL = SoundSpeed_i;
  const auto& q2L = GeometryToolbox::SquaredNorm(nDim,Velocity_i);
  const auto& qnL  = GeometryToolbox::DotProduct(nDim,Velocity_i,UnitNormal);


  const auto* velR = Velocity_j;
  const auto& rhoR = Density_j;
  const auto& pR = Pressure_j;
  const auto& HR = Enthalpy_j;
  const auto& aR = SoundSpeed_j;
  const auto& q2R = GeometryToolbox::SquaredNorm(nDim,Velocity_j);
  const auto& qnR  = GeometryToolbox::DotProduct(nDim,Velocity_j,UnitNormal);

  const auto* n = UnitNormal;

  //!Wave Strengths
  const auto drho = Density_j - Density_i; //!Density difference
  const auto dp   = Pressure_j - Pressure_i;   //!Pressure difference
  const auto dqn  = qnR - qnL;  //!Normal velocity difference

  LdU[0] = (dp - rho*a*dqn )/(2*a2); //!Left-moving acoustic wave strength
  LdU[1] =  drho - dp/(a2);            //!Entropy wave strength
  LdU[2] = (dp + rho*a*dqn )/(2*a2); //!Right-moving acoustic wave strength
  LdU[3] = rho;                         //!Shear wave strength (not really, just a factor)

  //!Absolute values of the wave Speeds

  ws[0] = abs(qn-a); //!Left-moving acoustic wave
  ws[1] = abs(qn);   //!Entropy wave
  ws[2] = abs(qn+a); //!Right-moving acoustic wave
  ws[3] = abs(qn);   //!Shear waves

  for(int iWave = 0; iWave < 4; ++iWave) ws_orig[iWave] = ws[iWave];


  //!Harten's entropy fix JCP(1983), 49, pp357-393 is applied to avoid vanishing
  //!wave speeds by making a parabolic fit near ws = 0. It is an entropy fix for
  //!nonlinear waves (avoids expnsion shock), but applied here as an eigenvalue
  //!limiting, for the pusrpose of avoiding 0 eigenvalues (wave speeds).

  //!Nonlinear fields
  const su2double elimc_nonlinear = 0.25;
  const su2double elimc_linear = 0.05;

  dws[0] = elimc_nonlinear *a;
  dws[2] = elimc_nonlinear *a;
  dws[1] = elimc_linear *a;
  dws[3] = elimc_linear *a;

  for(int k = 0; k < 4; ++k) {
    const su2double fix = ws[k] < dws[k];
    ws[k] = fix * (0.5 * (ws[k] * ws[k] / dws[k] + dws[k])) + (1-fix) * ws[k];
  }

  //!Right Eigenvectors
  //!Note: Two shear wave components are combined into 1, so that tangent vectors
  //!      are not required. And that's why there are only 4 vectors here.
  //!      See "I do like CFD, VOL.1" about how tangent vectors are eliminated.

  //! Left-moving acoustic wave
  R[0][0] = 1;
  for (iDim = 0; iDim < nDim; ++iDim) {
    R[iDim+1][0] = vel[iDim] - a*n[iDim];
  }
  R[nDim+1][0] = H - a*qn;

  //! Entropy wave
  R[0][1] = 1;
  for (iDim = 0; iDim < nDim; ++iDim)
    R[iDim+1][1] = vel[iDim];
  R[nDim+1][1] = 0.5*q2;

  //! Right-moving acoustic wave
  R[0][2] = 1;
  for (iDim = 0; iDim < nDim; ++iDim) {
    R[iDim+1][2] = vel[iDim] + a*n[iDim];
  }
  R[nDim+1][2] = H + a*qn;

  //! Two shear wave components combined into 1 (wave strength incorporated).

  for (iDim = 0; iDim < nDim; ++iDim) {
    dVel[iDim] = velR[iDim] - velL[iDim];
  }

  R[0][3] = 0;
  R[nDim+1][3] = - qn*dqn;
  for (iDim = 0; iDim < nDim; ++iDim) {
    R[iDim+1][3] = dVel[iDim] - dqn*n[iDim];
    R[nDim+1][3] += vel[iDim] * dVel[iDim];
  }

  //! We are now ready to compute the Jacobians for the Roe flux:
  //!
  //! Roe flux function -> Fn_Roe = 0.5*[Fn(ucL)+Fn(ucR)] - 0.5*sum_{k=1,4}|lambda_k|*(LdU)_k*r_k

  //!--------------------------------------------------------------------------------
  //! Part 1. Compute dFn_Roe/ducL
  //!
  //!  dFn_Roe/ducL =   d(0.5*Fn(ucL))/duL
  //!                 - 0.5*sum_{k=1,4} [d(|lambda_k|)/ducL]*(LdU)_k*r_k
  //!                 - 0.5*sum_{k=1,4} |lambda_k|*[d(LdU)_k/ducL]*r_k
  //!                 - 0.5*sum_{k=1,4} |lambda_k|*(LdU)_k*[dr_k/ducL]
  //!
  //!  So, we proceed as follows:
  //!
  //!  1.1 Compute                d(0.5*Fn(ucL))/duL
  //!  1.2 Compute various deriavives that will be used in the following steps.
  //!  1.3 Add the second term, - 0.5*sum_{k=1,4} [d(|lambda_k|)/ducL]*(LdU)_k*r_k
  //!  1.4 Add the  third term, - 0.5*sum_{k=1,4} |lambda_k|*[d(LdU)_k/ducL]*r_k
  //!  1.5 Add the fourth term, - 0.5*sum_{k=1,4} |lambda_k|*(LdU)_k*[dr_k/ducL]
  //!

  //!--------------------------------------
  //! 1.1 Compute "d(0.5*Fn(ucL))/ducL"
  //!
  //!     (See "I Do Like CFD, VOL.1", page 55, for the analytical Jacobian, dFn(u)/du)


  //!  1st column
  dFnducL[0][0] = 0;
  for (iDim = 0; iDim < nDim; ++iDim)
    dFnducL[iDim+1][0] = 0.5*(Gamma_Minus_One)*q2L*n[iDim]  - velL[iDim]*qnL;
  dFnducL[nDim+1][0] = 0.5*(Gamma_Minus_One)*q2L*qnL - HL*qnL;

  //!  2nd column
  for (iDim = 0; iDim < nDim; ++iDim) {
    dFnducL[0][iDim+1] = n[iDim];
    for (jDim = 0; jDim < nDim; ++jDim)
      dFnducL[jDim+1][iDim+1] = velL[jDim] * n[iDim] - (Gamma_Minus_One) * velL[iDim] * n[jDim];
    dFnducL[iDim+1][iDim+1] += qnL;
    dFnducL[nDim+1][iDim+1] = HL * n[iDim] - (Gamma_Minus_One) * velL[iDim] * qnL;
  }

  //!  5th column
  dFnducL[0][nDim+1] = 0;
  for (iDim = 0; iDim < nDim; ++iDim)
    dFnducL[iDim+1][nDim+1] = (Gamma_Minus_One)*n[iDim];
  dFnducL[nDim+1][nDim+1] =  Gamma*qnL;

  //!  Factor 1/2
  for (iVar = 0; iVar < nVar; ++iVar)
    for (jVar = 0; jVar < nVar; ++jVar)
      dFnducL[iVar][jVar] *= 0.5;

  //!--------------------------------------
  //! 1.2 Compute various deriavives that will be used in the following steps.
  //!     (See "I Do Like CFD, VOL.1" for details)

  //! dqn/ducL

  dqn_ducL[0] = -0.5*(qnL+qn) / (rhoL+rho);
  for (iDim = 0; iDim < nDim; ++iDim)
    dqn_ducL[iDim+1] =            n[iDim]  / (rhoL+rho);
  dqn_ducL[nDim+1] =          0;

  //! d(|qn|)/ducL

  su2double mask = qn < 0;
  su2double sigQN = mask * -1 + (1-mask) * 1;
  dabs_qn_ducL[0] = -0.5*sigQN*(qnL+qn) / (rhoL+rho);
  for (iDim = 0; iDim < nDim; ++iDim)
  dabs_qn_ducL[iDim+1] =       sigQN*     n[iDim]  / (rhoL+rho);
  dabs_qn_ducL[nDim+1] =  0;

  //! da/ducL

  const su2double velLdotVel = GeometryToolbox::DotProduct(nDim,velL,vel);
  da_ducL[0] =  0.5*(Gamma_Minus_One)/a*( 0.5*( velLdotVel + q2 )
                                        +  0.5*(HL-H) - aL*aL/(Gamma_Minus_One) + 0.5*(Gamma-2)*q2L )/(rhoL+rho);
  for (iDim = 0; iDim < nDim; ++iDim)
    da_ducL[iDim+1] = -0.5*(Gamma-1)*(vel[iDim]+(Gamma-1)*velL[iDim])/a  / (rhoL+rho);
  da_ducL[nDim+1] =  0.5*Gamma*(Gamma-1)/a               / (rhoL+rho);

  //! drho/ducL

  drho_ducL[0] = 0.5*rho/rhoL;
  for (iDim = 0; iDim < nDim; ++iDim) drho_ducL[iDim+1] = 0;
  drho_ducL[nDim+1] = 0;

  for (iDim = 0; iDim < nDim; ++iDim) {
    dvel_ducL[iDim][0] = -0.5*(velL[iDim]+vel[iDim]) / (rhoL+rho);
    for (jDim = 0; jDim < nDim; ++jDim) dvel_ducL[iDim][jDim+1] = 0;
    dvel_ducL[iDim][iDim+1] =          1 / (rhoL+rho);
    dvel_ducL[iDim][nDim+1] =  0;
  }

  //! dH/ducL

  dH_ducL[0] = ( 0.5*(HL-H) - aL*aL/(Gamma_Minus_One) + 0.5*(Gamma-2)*q2L ) / (rhoL+rho);
  for (iDim = 0; iDim < nDim; ++iDim)
  dH_ducL[iDim+1] = -Gamma_Minus_One*velL[iDim] / (rhoL+rho);
  dH_ducL[nDim+1] =              Gamma / (rhoL+rho);

  //! d(rhoR-rhoL)/ducL = - drhoL/ducL = - (drhoL/dWL)*(dWL/ducL) = -(1,0,0,0,0)*dW/dU

  ddrho_ducL[0] = - 1;
  for (iDim = 0; iDim < nDim; ++iDim) ddrho_ducL[iDim+1] = 0;
  ddrho_ducL[nDim+1] = 0;

  //! d(pR-pL)/ducL = - dpL/ducL = - (dpL/dWL)*(dWL/ducL) = -(0,0,0,0,1)*dW/dU

  ddp_ducL[0]   = - ( 0.5*(Gamma_Minus_One)*q2L );
  for (iDim = 0; iDim < nDim; ++iDim)
    ddp_ducL[iDim+1]   = Gamma_Minus_One*velL[iDim];
  ddp_ducL[nDim+1]   = - Gamma_Minus_One;

  //! d(qnR-qnL)/ducL = - dqnL/ducL = - (dqnL/dWL)*(dWL/ducL) = -(0,nx,ny,nz,0)*dW/dU
  ddqn_ducL[0]  = - (-qnL/rhoL);
  for (iDim = 0; iDim < nDim; ++iDim)
  ddqn_ducL[iDim+1]  = - (  n[iDim]/rhoL);
  ddqn_ducL[nDim+1]  = 0;

  //! d(uR-uL)/ducL = - duL/ducL = - (duL/dWL)*(dWL/ducL) = -(0,1,0,0,0)*dW/dU

  for (iDim = 0; iDim < nDim; ++iDim) {
    ddvel_ducL[iDim][0]   = velL[iDim]/rhoL;
    for (jDim = 0; jDim < nDim; ++jDim) ddvel_ducL[iDim][jDim+1] = 0;
    ddvel_ducL[iDim][iDim+1] = -1.0/rhoL;
    ddvel_ducL[iDim][nDim+1] = 0;
  }


  for(int k = 0; k < nVar; ++k)
    ddws1_ducL[k] = elimc_nonlinear*da_ducL[k];

  for(int k = 0; k < nVar; ++k)
    ddws3_ducL[k] = elimc_nonlinear*da_ducL[k];

  for(int k = 0; k < nVar; ++k)
    ddws2_ducL[k] = elimc_linear*da_ducL[k];

  for(int k = 0; k < nVar; ++k)
    ddws4_ducL[k] = elimc_linear*da_ducL[k];

  //!--------------------------------------
  //! 1.3 Differentiate the absolute values of the wave speeds, and
  //!     add the second term, - 0.5*sum_{k=1,4} [d(|lambda_k|)/ducL]*(LdU)_k*r_k

  //!  dws1_ducL = d(|qn-a|)/dUcL
  //!
  //!  Note on entropy fix:
  //!　 Let   dws1' = 0.5 * ( ws[1]*ws[1]/dws[1]+dws[1] ) <- Entropy fix
  ////! 　Then, dws1'/dUcL = (dws1'/dws1) * dws1_ducL
  //!                    = ws[1]/dws[1] * dws1_ducL

  //!  Absolute value

  mask = qn-a > 0;
  auto s = (mask) * 1 + (1-mask) * -1;
  for (iVar = 0; iVar < nVar; ++iVar) {
    dws1_ducL[iVar] = s * (dqn_ducL[iVar] - da_ducL[iVar]);
  }


  //!  Entropy fix/Eigenvalue limiting
  su2double fix = ws_orig[0] < dws[0];
  for (iVar = 0; iVar < nVar; ++iVar)
    dws1_ducL[iVar] = fix * (ws_orig[0]/dws[0] * dws1_ducL[iVar] + 0.5 * (-ws_orig[0]*ws_orig[0]/(dws[0]*dws[0]) + 1)*ddws1_ducL[iVar]) +
                      (1-fix) * dws1_ducL[iVar];


  for(iVar = 0; iVar < nVar; ++iVar)
    for(jVar = 0; jVar < nVar; ++jVar)
      dFnducL[iVar][jVar] -= 0.5 * dws1_ducL[jVar]*LdU[0]*R[iVar][0];


  //!  dws2_ducL = d(|qn|)/dUcL

  for(iVar = 0; iVar < nVar; ++iVar)
    dws2_ducL[iVar] = dabs_qn_ducL[iVar];

  //!  Entropy fix/Eigenvalue limiting
  fix = ws_orig[0] < dws[0];
  for (iVar = 0; iVar < nVar; ++iVar)
    dws2_ducL[iVar] = fix * (ws_orig[1]/dws[1] * dws2_ducL[iVar] + 0.5 * (-ws_orig[1]*ws_orig[1]/(dws[1]*dws[1]) + 1)*ddws2_ducL[iVar]) +
                      (1-fix) * dws2_ducL[iVar];


  for(iVar = 0; iVar < nVar; ++iVar)
    for(jVar = 0; jVar < nVar; ++jVar)
      dFnducL[iVar][jVar] -= 0.5 * dws2_ducL[jVar]*LdU[1]*R[iVar][1];


  //!  dws3_ducL = d(|qn+a|)/dUcL
  //!
  //!  Note on entropy fix:
  //!　 Let   dws3' = 0.5 * ( ws[3]*ws[3]/dws[3]+dws[3] ) <- Entropy fix
  //! 　Then, dws3'/dUcL = (dws3'/dws3) * dws3_ducL
  //!                    = ws[3]/dws[3] * dws3_ducL

  //!  Absolute value

  mask = qn+a > 0;
  s = mask * 1 + (1-mask) * -1;
  for (iVar = 0; iVar < nVar; ++iVar) {
    dws3_ducL[iVar] = s * (dqn_ducL[iVar] + da_ducL[iVar]);
  }


  //!  Entropy fix/Eigenvalue limiting
  fix = ws_orig[2] < dws[2];
  for (iVar = 0; iVar < nVar; ++iVar)
    dws3_ducL[iVar] = fix * (ws_orig[2]/dws[2] * dws3_ducL[iVar] + 0.5 * (-ws_orig[2]*ws_orig[2]/(dws[2]*dws[2]) + 1)*ddws3_ducL[iVar])
                      + (1-fix) * dws3_ducL[iVar];

  for(iVar = 0; iVar < nVar; ++iVar)
    for(jVar = 0; jVar < nVar; ++jVar)
      dFnducL[iVar][jVar] -= 0.5 * dws3_ducL[jVar]*LdU[2]*R[iVar][2];

  //!  dws4_ducL = d(|qn|)/dUcL = dws1_ducL


  for(iVar = 0; iVar < nVar; ++iVar)
    dws4_ducL[iVar] = dabs_qn_ducL[iVar];

  //!  Entropy fix/Eigenvalue limiting
  fix = ws_orig[3] < dws[3];
  for (iVar = 0; iVar < nVar; ++iVar)
    dws4_ducL[iVar] = fix * (ws_orig[3]/dws[3] * dws4_ducL[iVar] + 0.5 * (-ws_orig[3]*ws_orig[3]/(dws[3]*dws[3]) + 1)*ddws4_ducL[iVar]) +
                      (1-fix) * dws4_ducL[iVar];


  for(iVar = 0; iVar < nVar; ++iVar)
    for(jVar = 0; jVar < nVar; ++jVar)
      dFnducL[iVar][jVar] -= 0.5 * dws4_ducL[jVar]*LdU[3]*R[iVar][3];

  //!--------------------------------------
  //! 1.4 Differentiate the wave strength, and
  //!     add the third term, - 0.5*sum_{k=1,4} |lambda_k|*[d(LdU)_k/ducL]*r_k.

  //!  dLdU1_ducL = d( dp/(2a^2) - rho/(2a)*dqn )/dUcL
  //!
  //!             = dp*[d(1/(2a^2))/dUcL] - dqn*[d(rho/(2a))/dUcL]
  //!             + [d(dp)/dUcL]/(2a^2)   - [d(dqn)/dUcL]*rho/(2a)
  //!
  //!             = -a^(-3)*dp*[da/dUcL]  - dqn*( [drho/dUcL]/(2a) - rho/(2*a^2)*[da/dUcL] )
  //!             + [d(dp)/dUcL]/(2a^2)   - [d(dqn)/dUcL]*rho/(2a)
  //!
  //!             = ( -2*dp + rho*a*dqn )/(2a^3)*[da/dUcL] - dqn*[drho/dUcL]/(2a)
  //!             + [d(dp)/dUcL]/(2a^2)   - [d(dqn)/dUcL]*rho/(2a)


  for(iVar = 0; iVar < nVar; ++iVar) {
    dLdU1_ducL[iVar] = 0.5*(-2*dp+rho*a*dqn )/(a2*a) * (da_ducL[iVar])
                       - 0.5*dqn/a * (drho_ducL[iVar])
                       + 0.5*(ddp_ducL[iVar])/a2
                       - 0.5*rho*(ddqn_ducL[iVar])/a;
  }

  for(iVar = 0; iVar < nVar; ++iVar)
    for(jVar = 0; jVar < nVar; ++jVar)
      dFnducL[iVar][jVar] -= 0.5 * ws[0]*dLdU1_ducL[jVar]*R[iVar][0];

  //!  dLdU2_ducL = d( drho - dp/a^2 )/dUcL
  //!              = [d(drho)/dUcL] - [d(dp)/dUcL]/a^2 + 2*dp/a^3*[da/dUcL]


  for(iVar = 0; iVar < nVar; ++iVar)
    dLdU2_ducL[iVar] = (ddrho_ducL[iVar]) - (ddp_ducL[iVar])/a2 + 2*dp/(a2*a) * (da_ducL[iVar]);

  for(iVar = 0; iVar < nVar; ++iVar)
    for(jVar = 0; jVar < nVar; ++jVar)
      dFnducL[iVar][jVar] -= 0.5 * ws[1]*dLdU2_ducL[jVar]*R[iVar][1];

  //!  dLdU3_ducL = d( dp/(2a^2) + rho/(2a)*dqn )/dUcL
  //!
  //!             = dp*[d(1/(2a^2))/dUcL] + dqn*[d(rho/(2a))/dUcL]
  //!             + [d(dp)/dUcL]/(2a^2)   + [d(dqn)/dUcL]*rho/(2a)
  //!
  //!             = -a^(-3)*dp*[da/dUcL]  + dqn*( [drho/dUcL]/(2a) - rho/(2*a^2)*[da/dUcL] )
  //!             + [d(dp)/dUcL]/(2a^2)   + [d(dqn)/dUcL]*rho/(2a)
  //!
  //!             = ( -2*dp - rho*a*dqn )/(2a^3)*[da/dUcL]  + dqn*[drho/dUcL]/(2a)
  //!             + [d(dp)/dUcL]/(2a^2)   + [d(dqn)/dUcL]*rho/(2a)


  for(iVar = 0; iVar < nVar; ++iVar) {
    dLdU3_ducL[iVar] = 0.5 * (-2 * dp - rho * a * dqn) / (a2*a) * (da_ducL[iVar]) + 0.5 * dqn / a * (drho_ducL[iVar]) +
                       0.5 * (ddp_ducL[iVar]) / a2 + 0.5 * rho * (ddqn_ducL[iVar]) / a;
  }

  for(iVar = 0; iVar < nVar; ++iVar)
    for(jVar = 0; jVar < nVar; ++jVar)
      dFnducL[iVar][jVar] -= 0.5 * ws[2]*dLdU3_ducL[jVar]*R[iVar][2];

  //!  dLdU4_ducL = d(rho)/dUcL

  for(iVar = 0; iVar < 5; ++iVar) dLdU4_ducL[iVar] = drho_ducL[iVar];

  for(iVar = 0; iVar < nVar; ++iVar)
    for(jVar = 0; jVar < nVar; ++jVar)
      dFnducL[iVar][jVar] -= 0.5 * ws[3]*dLdU4_ducL[jVar]*R[iVar][3];

  //!--------------------------------------
  //! 1.5 Differentiate the right-eigenvectors, and
  //!     add the fourth term, - 0.5*sum_{k=1,4} |lambda_k|*(LdU)_k*[dr_k/ducL]

  //! dR1_ducL = dR(:,1)/dUcL
  //!
  //! Left-moving acoustic wave
  //!
  //!        Eigenvector -> Differentiated
  //!
  //!  R(1,1) = 1      -> 0
  //!  R(2,1) = u - a*nx -> du/dUcL - da/dUcL*nx
  //!  R(3,1) = v - a*ny -> dv/dUcL - da/dUcL*ny
  //!  R(4,1) = w - a*nz -> dw/dUcL - da/dUcL*nz
  //!  R(5,1) = H - a*qn -> dH/dUcL - da/dUcL*qn - dqn/dUcL*a

  for (iVar = 0; iVar < nVar; ++iVar) dR1_ducL[0][iVar] = 0;

  for (iDim = 0; iDim < nDim; ++iDim)
    for (iVar = 0; iVar < nVar; ++iVar)
      dR1_ducL[iDim+1][iVar] = (dvel_ducL[iDim][iVar]) - (da_ducL[iVar]) * n[iDim];

  for (iVar = 0; iVar < nVar; ++iVar) dR1_ducL[nDim+1][iVar] = (dH_ducL[iVar]) - (da_ducL[iVar]) * qn - (dqn_ducL[iVar]) * a;

  for(iVar = 0; iVar < nVar; ++iVar)
    for(jVar = 0; jVar < nVar; ++jVar)
      dFnducL[iVar][jVar] -= 0.5 * ws[0]*LdU[0]*dR1_ducL[iVar][jVar];

  //! dR2_ducL = dR(:,2)/dUcL
  //!
  //! Entropy wave
  //!
  //!                  Eigenvector -> Differentiated
  //!
  //!  R(1,2) = 1                -> 0
  //!  R(2,2) = u                  -> du/dUcL
  //!  R(3,2) = v                  -> dv/dUcL
  //!  R(4,2) = w                  -> dw/dUcL
  //!  R(5,2) = 0.5*(u*u+v*v+w*w) -> u*du/dUcL + v*dv/dUcL + w*dw/dUcL

  for (iVar = 0; iVar < nVar; ++iVar) dR2_ducL[0][iVar] = 0;

  for (iDim = 0; iDim < nDim; ++iDim) {
    for (iVar = 0; iVar < nVar; ++iVar) {
      dR2_ducL[iDim+1][iVar] = (dvel_ducL[iDim][iVar]);
    }
  }

  for (iVar = 0; iVar < nVar; ++iVar) {
    dR2_ducL[nDim+1][iVar] = 0;
    for (iDim = 0; iDim < nDim; ++iDim)
      dR2_ducL[nDim+1][iVar] += vel[iDim]*(dvel_ducL[iDim][iVar]);
  }

  for(iVar = 0; iVar < nVar; ++iVar)
    for(jVar = 0; jVar < nVar; ++jVar)
      dFnducL[iVar][jVar] -= 0.5 * ws[1]*LdU[1]*dR2_ducL[iVar][jVar];

  //! dR3_ducL = dR(:,3)/dUcL
  //!
  //! Right-moving acoustic wave
  //!
  //!        Eigenvector -> Differentiated
  //!
  //!  R(1,3) = 1      -> 0
  //!  R(2,3) = u + a*nx -> du/dUcL + da/dUcL*nx
  //!  R(3,3) = v + a*ny -> dv/dUcL + da/dUcL*ny
  //!  R(4,3) = w + a*nz -> dw/dUcL + da/dUcL*nz
  //!  R(5,3) = H + a*qn -> dH/dUcL + da/dUcL*qn + dqn/dUcL*a

  for (iVar = 0; iVar < nVar; ++iVar) dR3_ducL[0][iVar] = 0;

  for (iDim = 0; iDim < nDim; ++iDim)
    for (iVar = 0; iVar < nVar; ++iVar)
      dR3_ducL[iDim+1][iVar] = (dvel_ducL[iDim][iVar]) + (da_ducL[iVar]) * n[iDim];

  for (iVar = 0; iVar < nVar; ++iVar) dR3_ducL[nDim+1][iVar] = (dH_ducL[iVar]) + (da_ducL[iVar]) * qn + (dqn_ducL[iVar]) * a;

  for(iVar = 0; iVar < nVar; ++iVar)
    for(jVar = 0; jVar < nVar; ++jVar)
      dFnducL[iVar][jVar] -= 0.5 * ws[2]*LdU[2]*dR3_ducL[iVar][jVar];

  //! dR4_ducL = dR(:,4)/dUcL
  //! Two shear wave components combined into 1 (wave strength incorporated).
  //! So, it is not really an eigenvector.
  //!
  //!                  Combined vector -> Differentiated
  //!
  //!  R(1,4) = 0                   -> 0
  //!  R(2,4) = du - dqn*nx            -> d(du)/dUcL - d(dqn)/dUcL*nx
  //!  R(3,4) = dv - dqn*ny            -> d(dv)/dUcL - d(dqn)/dUcL*ny
  //!  R(4,4) = dw - dqn*nz            -> d(dw)/dUcL - d(dqn)/dUcL*nz
  //!  R(5,4) = u*du+v*dv+w*dw-qn*dqn  -> du/dUcL*du     + d(du)/dUcL*u
  //!                                   + dv/dUcL*dv     + d(dv)/dUcL*v
  //!                                   + dw/dUcL*dw     + d(dw)/dUcL*w
  //!                                   - d(qn)/dUcL*dqn - d(dqn)/dUcL*qn

  for (iVar = 0; iVar < nVar; ++iVar) dR4_ducL[0][iVar] = 0;

  for (iDim = 0; iDim < nDim; ++iDim) {
    for (iVar = 0; iVar < nVar; ++iVar) {
      dR4_ducL[iDim+1][iVar] = (ddvel_ducL[iDim][iVar]) - (ddqn_ducL[iVar])*n[iDim];
    }
  }

  for (iVar = 0; iVar < nVar; ++iVar) {
    dR4_ducL[nDim+1][iVar] = - dqn*( dqn_ducL[iVar]) - qn*(ddqn_ducL[iVar]);
    for (iDim = 0; iDim < nDim; ++iDim)
      dR4_ducL[nDim+1][iVar] += dVel[iDim]*( dvel_ducL[iDim][iVar]) + vel[iDim]*(ddvel_ducL[iDim][iVar]);
  }

  for(iVar = 0; iVar < nVar; ++iVar)
    for(jVar = 0; jVar < nVar; ++jVar)
      dFnducL[iVar][jVar] -= 0.5 * ws[3]*LdU[3]*dR4_ducL[iVar][jVar];


  //!--------------------------------------------------------------------------------
  //! Part 2. Compute dFn_Roe/ducR
  //!
  //!  dFn_Roe/ducR =   d(0.5*Fn(ucR))/duR
  //!                 - 0.5*sum_{k=1,4} [d(|lambda_k|)/ducR]*(LdU)_k*r_k
  //!                 - 0.5*sum_{k=1,4} |lambda_k|*[d(LdU)_k/ducR]*r_k
  //!                 - 0.5*sum_{k=1,4} |lambda_k|*(LdU)_k*[dr_k/ducR]
  //!
  //!  So, we proceed as follows:
  //!
  //!  1.1 Compute                d(0.5*Fn(ucR))/duR
  //!  1.2 Compute various deriavives that will be used in the following steps.
  //!  1.3 Add the second term, - 0.5*sum_{k=1,4} [d(|lambda_k|)/ducR]*(LdU)_k*r_k
  //!  1.4 Add the  third term, - 0.5*sum_{k=1,4} |lambda_k|*[d(LdU)_k/ducR]*r_k
  //!  1.5 Add the fourth term, - 0.5*sum_{k=1,4} |lambda_k|*(LdU)_k*[dr_k/ducR]
  //!

  //!--------------------------------------
  //! 2.1 Compute "d(0.5*Fn(ucR))/ducR"
  //!
  //!     (See "I Do Like CFD, VOL.1", page 55, for the analytical Jacobian, dFn(u)/du)

  //!  1st column

  dFnducR[0][0] = 0;
  for (iDim = 0; iDim < nDim; ++iDim)
    dFnducR[iDim+1][0] = 0.5*(Gamma_Minus_One)*q2R*n[iDim]  - velR[iDim]*qnR;
  dFnducR[nDim+1][0] = 0.5*(Gamma_Minus_One)*q2R*qnR - HR*qnR;

  //!  2nd column
  for (iDim = 0; iDim < nDim; ++iDim) {
    dFnducR[0][iDim+1] =     n[iDim];
    for (jDim = 0; jDim < nDim; ++jDim) {
      dFnducR[jDim+1][iDim+1] =  velR[jDim]*n[iDim] - (Gamma_Minus_One)*velR[iDim]*n[jDim];
    }
    dFnducR[iDim+1][iDim+1] += qnR;
    dFnducR[nDim+1][iDim+1] =  HR*n[iDim] - (Gamma_Minus_One)*velR[iDim]*qnR;
  }

  //!  5th column
  dFnducR[0][nDim+1] =  0;
  for (iDim = 0; iDim < nDim; ++iDim)
    dFnducR[iDim+1][nDim+1] = (Gamma_Minus_One)*n[iDim];
  dFnducR[nDim+1][nDim+1] =  Gamma*qnR;

  //!  Factor 1/2
  for(iVar = 0; iVar < nVar; ++iVar)
    for(jVar = 0; jVar < nVar; ++jVar)
      dFnducR[iVar][jVar] *= 0.5;

  //!--------------------------------------
  //! 2.2 Compute various deriavives that will be used in the following steps.
  //!     (See "I Do Like CFD, VOL.1" for details)

  //! dqn/ducR


  dqn_ducR[0] = -0.5*(qnR+qn) / (rhoR+rho);
  for (iDim = 0; iDim < nDim; ++iDim)
    dqn_ducR[iDim+1] =            n[iDim]  / (rhoR+rho);
  dqn_ducR[nDim+1] =          0;

  //! d(|qn|)/ducR

  mask = qn < 0;
  sigQN = mask * -1 + (1-mask) * 1;
  dabs_qn_ducR[0] = -0.5*sigQN*(qnR+qn) / (rhoR+rho);
  for (iDim = 0; iDim < nDim; ++iDim)
    dabs_qn_ducR[iDim+1] =       sigQN*     n[iDim]  / (rhoR+rho);
  dabs_qn_ducR[nDim+1] =  0;

  //! da/ducR

  const su2double velRdotVel = GeometryToolbox::DotProduct(nDim,velR,vel);

  da_ducR[0] =  0.5*(Gamma-1)/a*( 0.5*( velRdotVel + q2 )
                                        +  0.5*(HR-H) - aR*aR/(Gamma_Minus_One) + 0.5*(Gamma-2)*q2R )/(rhoR+rho);
  for (iDim = 0; iDim < nDim; ++iDim)
    da_ducR[iDim+1] = -0.5*(Gamma_Minus_One)*(vel[iDim]+(Gamma_Minus_One)*velR[iDim])/a  / (rhoR+rho);
  da_ducR[nDim+1] =  0.5*Gamma*(Gamma_Minus_One)/a               / (rhoR+rho);

  //! drho/ducR

  drho_ducR[0] = 0.5*rho/rhoR;
  for (iDim = 0; iDim < nDim; ++iDim) drho_ducR[iDim+1] = 0;
  drho_ducR[nDim+1] = 0;

  //! du/ducR

  for (iDim = 0; iDim < nDim; ++iDim) {
    dvel_ducR[iDim][0] = -0.5*(velR[iDim]+vel[iDim]) / (rhoR+rho);
    for (jDim = 0; jDim < nDim; ++jDim) dvel_ducR[iDim][jDim+1] = 0;
    dvel_ducR[iDim][iDim+1] = 1 / (rhoR+rho);
    dvel_ducR[iDim][nDim+1] = 0;
  }

  //! dH/ducR

  dH_ducR[0] = ( 0.5*(HR-H) - aR*aR/(Gamma_Minus_One) + 0.5*(Gamma-2)*q2R ) / (rhoR+rho);
  for (iDim = 0; iDim < nDim; ++iDim)
  dH_ducR[iDim+1] = -Gamma_Minus_One*velR[iDim] / (rhoR+rho);
  dH_ducR[nDim+1] =              Gamma / (rhoR+rho);

  //! d(rhoR-rhoR)/ducR = drhoR/ducR = (drhoR/dWR)*(dWR/ducR) = (1,0,0,0,0)*dW/dU

  ddrho_ducR[0] =  (  1 );
  for (iDim = 0; iDim < nDim; ++iDim) ddrho_ducR[iDim+1] = 0;
  ddrho_ducR[nDim+1] = 0;

  //! d(pR-pR)/ducR = dpR/ducR = (dpR/dWR)*(dWR/ducR) = (0,0,0,0,1)*dW/dU

  ddp_ducR[0]   =  ( 0.5*(Gamma_Minus_One)*q2R );
  for (iDim = 0; iDim < nDim; ++iDim)
    ddp_ducR[iDim+1]   =  (    - (Gamma_Minus_One)*velR[iDim]  );
  ddp_ducR[nDim+1]   =  (       Gamma_Minus_One      );

  //! d(qnR-qnR)/ducR = dqnR/ducR = (dqnR/dWR)*(dWR/ducR) = (0,nx,ny,nz,0)*dW/dU

  ddqn_ducR[0]  =  (-qnR/rhoR);
  for (iDim = 0; iDim < nDim; ++iDim)
    ddqn_ducR[iDim+1]  =  (  n[iDim]/rhoR);
  ddqn_ducR[nDim+1]  =  ( 0    );


  for (iDim = 0; iDim < nDim; ++iDim) {
    ddvel_ducR[iDim][0]   =  ( -velR[iDim]/rhoR);
    for (jDim = 0; jDim < nDim; ++jDim)
      ddvel_ducR[iDim][jDim+1] = 0;
    ddvel_ducR[iDim][iDim+1] = 1.0/rhoR;
    dvel_ducR[iDim][nDim+1] = 0;
  }



  for (iVar = 0; iVar < nVar; ++iVar) ddws1_ducR[iVar] = elimc_nonlinear*da_ducR[iVar];


  for (iVar = 0; iVar < nVar; ++iVar) ddws3_ducR[iVar] = elimc_nonlinear*da_ducR[iVar];


  for (iVar = 0; iVar < nVar; ++iVar) ddws2_ducR[iVar] = elimc_linear*da_ducR[iVar];


  for (iVar = 0; iVar < nVar; ++iVar) ddws4_ducR[iVar] = elimc_linear*da_ducR[iVar];

  //!--------------------------------------
  //! 2.3 Differentiate the absolute values of the wave speeds, and
  //!     add the second term, - 0.5*sum_{k=1,4} [d(|lambda_k|)/ducR]*(LdU)_k*r_k
  //
  //!  dws1_ducR = d(|qn-a|)/dUcR
  //!
  //!  Note on entropy fix:
  //!　 Let   dws1' = 0.5 * ( ws[1]*ws[1]/dws[1]+dws[1] ) <- Entropy fix
  //! 　Then, dws1'/dUcR = (dws1'/dws1) * dws1_ducR    + ( -ws[1]*ws[1]/dws[1]**2*ddws[1] + ddws[1] )
  //!                    = ws[1]/dws[1] * dws1_ducR

  //!  Absolute value


  mask = qn-a > 0;
  s = (mask) * 1 + (1-mask) * -1;
  for (iVar = 0; iVar < nVar; ++iVar)
    dws1_ducR[iVar] =  s * (dqn_ducR[iVar] - da_ducR[iVar]);


  //!  Entropy fix/Eigenvalue limiting

  fix = ws_orig[0] < dws[0];
  for (iVar = 0; iVar < nVar; ++iVar)
    dws1_ducR[iVar] = fix * (ws_orig[0]/dws[0] * dws1_ducR[iVar] + 0.5 * (-ws_orig[0]*ws_orig[0]/(dws[0]*dws[0]) + 1)*ddws1_ducR[iVar]) +
                      (1-fix) * dws1_ducR[iVar];

  for (iVar = 0; iVar < nVar; ++iVar)
    for (jVar = 0; jVar < nVar; ++jVar)
      dFnducR[iVar][jVar] -= 0.5 * dws1_ducR[jVar]*LdU[0]*R[iVar][0];

  //!  dws2_ducR = d(|qn|)/dUcR


  for(iVar = 0; iVar < 5; ++iVar) dws2_ducR[iVar] = dabs_qn_ducR[iVar];

  //!  Entropy fix/Eigenvalue limiting
  fix = ws_orig[1] < dws[1];
  for (iVar = 0; iVar < nVar; ++iVar)
    dws2_ducR[iVar] = fix * (ws_orig[1]/dws[1] * dws2_ducR[iVar] + 0.5 * (-ws_orig[1]*ws_orig[1]/(dws[1]*dws[1]) + 1)*ddws2_ducR[iVar]) +
                      (1-fix) * dws2_ducR[iVar];

  for (iVar = 0; iVar < nVar; ++iVar)
    for (jVar = 0; jVar < nVar; ++jVar)
      dFnducR[iVar][jVar] -= 0.5 * dws2_ducR[jVar]*LdU[1]*R[iVar][1];

  //!  dws3_ducR = d(|qn+a|)/dUcR
  //!
  //!  Note on entropy fix:
  //!　 Let   dws3' = 0.5 * ( ws[3]*ws[3]/dws[3]+dws[3] ) <- Entropy fix
  //! 　Then, dws3'/dUcR = (dws3'/dws3) * dws3_ducR
  //!                    = ws[3]/dws[3] * dws3_ducR

  //!  Absolute value

  mask = qn+a > 0;
  s = mask * 1 + (1-mask) * -1;
  for (iVar = 0; iVar < nVar; ++iVar)
    dws3_ducR[iVar] =  s * (dqn_ducR[iVar] + da_ducR[iVar]);


  //!  Entropy fix/Eigenvalue limiting

  fix = ws_orig[2] < dws[2];
  for (iVar = 0; iVar < nVar; ++iVar)
    dws3_ducR[iVar] = fix * (ws_orig[2]/dws[2] * dws3_ducR[iVar] + 0.5 * (-ws_orig[2]*ws_orig[2]/(dws[2]*dws[2]) + 1)*ddws3_ducR[iVar]) +
                      (1-fix) * dws3_ducR[iVar];

  for (iVar = 0; iVar < nVar; ++iVar)
    for (jVar = 0; jVar < nVar; ++jVar)
      dFnducR[iVar][jVar] -= 0.5 * dws3_ducR[jVar]*LdU[2]*R[iVar][2];

  //!  dws4_ducR = d(|qn|)/dUcR = dws1_ducR


  for(iVar = 0; iVar < 5; ++iVar) dws4_ducR[iVar] = dabs_qn_ducR[iVar];

  //!  Entropy fix/Eigenvalue limiting
  fix = ws_orig[3] < dws[3];
  for (iVar = 0; iVar < nVar; ++iVar)
    dws4_ducR[iVar] = fix * (ws_orig[3]/dws[3] * dws4_ducR[iVar] + 0.5 * (-ws_orig[3]*ws_orig[3]/(dws[3]*dws[3]) + 1)*ddws4_ducR[iVar]) +
                      (1-fix) * dws4_ducR[iVar];

  for (iVar = 0; iVar < nVar; ++iVar)
    for (jVar = 0; jVar < nVar; ++jVar)
      dFnducR[iVar][jVar] -= 0.5 * dws4_ducR[jVar]*LdU[3]*R[iVar][3];

  //!--------------------------------------
  //! 2.4 Differentiate the wave strength, and
  //!     add the third term, - 0.5*sum_{k=1,4} |lambda_k|*[d(LdU)_k/ducR]*r_k.
  //
  //!  dLdU1_ducR = d( dp/(2a^2) - rho/(2a)*dqn )/dUcR
  //!
  //!             = dp*[d(1/(2a^2))/dUcR] - dqn*[d(rho/(2a))/dUcR]
  //!             + [d(dp)/dUcR]/(2a^2)   - [d(dqn)/dUcR]*rho/(2a)
  //!
  //!             = -a^(-3)*dp*[da/dUcR]  - dqn*( [drho/dUcR]/(2a) - rho/(2*a^2)*[da/dUcR] )
  //!             + [d(dp)/dUcR]/(2a^2)   - [d(dqn)/dUcR]*rho/(2a)
  //!
  //!             = ( -2*dp + rho*a*dqn )/(2a^3)*[da/dUcR] - dqn*[drho/dUcR]/(2a)
  //!             + [d(dp)/dUcR]/(2a^2)   - [d(dqn)/dUcR]*rho/(2a)


  for (iVar = 0; iVar < nVar; ++iVar) {
    dLdU1_ducR[iVar] = 0.5 * (-2 * dp + rho * a * dqn) / (a2 * a) * (da_ducR[iVar]) -
                       0.5 * dqn / a * (drho_ducR[iVar]) + 0.5 * (ddp_ducR[iVar]) / a2 -
                       0.5 * rho * (ddqn_ducR[iVar]) / a;
  }

  for (iVar = 0; iVar < nVar; ++iVar)
    for (jVar = 0; jVar < nVar; ++jVar)
      dFnducR[iVar][jVar] -= 0.5 * ws[0]*dLdU1_ducR[jVar]*R[iVar][0];

  //!  dLdU2_ducR = d( drho - dp/a^2 )/dUcR
  //!              = [d(drho)/dUcR] - [d(dp)/dUcR]/a^2 + 2*dp/a^3*[da/dUcR]


  for (iVar = 0; iVar < nVar; ++iVar)
    dLdU2_ducR[iVar] = (ddrho_ducR[iVar]) - (ddp_ducR[iVar])/a2 + 2*dp/(a2*a) * (da_ducR[iVar]);

  for (iVar = 0; iVar < nVar; ++iVar)
    for (jVar = 0; jVar < nVar; ++jVar)
      dFnducR[iVar][jVar] = dFnducR[iVar][jVar] - 0.5 * ws[1]*dLdU2_ducR[jVar]*R[iVar][1];

  //!  dLdU3_ducR = d( dp/(2a^2) + rho/(2a)*dqn )/dUcR
  //!
  //!             = dp*[d(1/(2a^2))/dUcR] + dqn*[d(rho/(2a))/dUcR]
  //!             + [d(dp)/dUcR]/(2a^2)   + [d(dqn)/dUcR]*rho/(2a)
  //!
  //!             = -a^(-3)*dp*[da/dUcR]  + dqn*( [drho/dUcR]/(2a) - rho/(2*a^2)*[da/dUcR] )
  //!             + [d(dp)/dUcR]/(2a^2)   + [d(dqn)/dUcR]*rho/(2a)
  //!
  //!             = ( -2*dp - rho*a*dqn )/(2a^3)*[da/dUcR]  + dqn*[drho/dUcR]/(2a)
  //!             + [d(dp)/dUcR]/(2a^2)   + [d(dqn)/dUcR]*rho/(2a)


  for (iVar = 0; iVar < nVar; ++iVar) {
    dLdU3_ducR[iVar] = 0.5 * (-2 * dp - rho * a * dqn) / (a2 * a) * (da_ducR[iVar]) +
                       0.5 * dqn / a * (drho_ducR[iVar]) + 0.5 * (ddp_ducR[iVar]) / a2 +
                       0.5 * rho * (ddqn_ducR[iVar]) / a;
  }

  for (iVar = 0; iVar < nVar; ++iVar)
    for (jVar = 0; jVar < nVar; ++jVar)
      dFnducR[iVar][jVar] = dFnducR[iVar][jVar] - 0.5 * ws[2]*dLdU3_ducR[jVar]*R[iVar][2];

  //!  dLdU4_ducR = d(rho)/dUcR
  for(iVar = 0; iVar < 5; ++iVar) dLdU4_ducR[iVar] = drho_ducR[iVar];

  for (iVar = 0; iVar < nVar; ++iVar)
    for (jVar = 0; jVar < nVar; ++jVar)
      dFnducR[iVar][jVar] = dFnducR[iVar][jVar] - 0.5 * ws[3]*dLdU4_ducR[jVar]*R[iVar][3];

  //!--------------------------------------
  //! 2.5 Differentiate the right-eigenvectors, and
  //!     add the fourth term, - 0.5*sum_{k=1,4} |lambda_k|*(LdU)_k*[dr_k/ducL]
  //
  //! dR1_ducR = dR(:,1)/dUcR
  //!
  //! Reft-moving acoustic wave
  //!
  //!        Eigenvector -> Differentiated
  //!
  //!  R(1,1) = 1      -> 0
  //!  R(2,1) = u - a*nx -> du/dUcR - da/dUcR*nx
  //!  R(3,1) = v - a*ny -> dv/dUcR - da/dUcR*ny
  //!  R(4,1) = w - a*nz -> dw/dUcR - da/dUcR*nz
  //!  R(5,1) = H - a*qn -> dH/dUcR - da/dUcR*qn - dqn/dUcR*a

  for (iVar = 0; iVar < nVar; ++iVar) dR1_ducR[0][iVar] = 0;
  for (iDim = 0; iDim < nDim; ++iDim) {
    for (iVar = 0; iVar < nVar; ++iVar)
      dR1_ducR[iDim+1][iVar] = (dvel_ducR[iDim][iVar]) - (da_ducR[iVar]) * n[iDim];
  }
  for (iVar = 0; iVar < nVar; ++iVar)
    dR1_ducR[nDim+1][iVar] = (dH_ducR[iVar]) - (da_ducR[iVar]) * qn - (dqn_ducR[iVar]) * a;

  for (iVar = 0; iVar < nVar; ++iVar)
    for (jVar = 0; jVar < nVar; ++jVar)
      dFnducR[iVar][jVar] -= 0.5 * ws[0]*LdU[0]*dR1_ducR[iVar][jVar];

  //! dR2_ducR = dR(:,2)/dUcR
  //!
  //! Entropy wave
  //!
  //!                  Eigenvector -> Differentiated
  //!
  //!  R(1,2) = 1                -> 0
  //!  R(2,2) = u                  -> du/dUcR
  //!  R(3,2) = v                  -> dv/dUcR
  //!  R(4,2) = w                  -> dw/dUcR
  //!  R(5,2) = 0.5*(u*u+v*v+w*w) -> u*du/dUcR + v*dv/dUcR + w*dw/dUcR

  for (iVar = 0; iVar < nVar; ++iVar) dR2_ducR[0][iVar] = 0;

  for (iDim = 0; iDim < nDim; ++iDim) {
    for (iVar = 0; iVar < nVar; ++iVar)
      dR2_ducR[iDim+1][iVar] = (dvel_ducR[iDim][iVar]);
  }

  for (iVar = 0; iVar < nVar; ++iVar) {
    dR2_ducR[nDim+1][iVar] = 0;
    for (iDim = 0; iDim < nDim; ++iDim) {
      dR2_ducR[nDim+1][iVar] += vel[iDim]*(dvel_ducR[iDim][iVar]);
    }
  }


  for (iVar = 0; iVar < nVar; ++iVar)
    for (jVar = 0; jVar < nVar; ++jVar)
      dFnducR[iVar][jVar] -= 0.5 * ws[1]*LdU[1]*dR2_ducR[iVar][jVar];

  //! dR3_ducR = dR(:,3)/dUcR
  //!
  //! Right-moving acoustic wave
  //!
  //!        Eigenvector -> Differentiated
  //!
  //!  R(1,3) = 1      -> 0
  //!  R(2,3) = u + a*nx -> du/dUcR + da/dUcR*nx
  //!  R(3,3) = v + a*ny -> dv/dUcR + da/dUcR*ny
  //!  R(4,3) = w + a*nz -> dw/dUcR + da/dUcR*nz
  //!  R(5,3) = H + a*qn -> dH/dUcR + da/dUcR*qn + dqn/dUcR*a

  for (iVar = 0; iVar < nVar; ++iVar) dR3_ducR[0][iVar] = 0;
  for (iDim = 0; iDim < nDim; ++iDim) {
    for (iVar = 0; iVar < nVar; ++iVar)
      dR3_ducR[iDim+1][iVar] = (dvel_ducR[iDim][iVar]) + (da_ducR[iVar]) * n[iDim];
  }
  for (iVar = 0; iVar < nVar; ++iVar) dR3_ducR[nDim+1][iVar] = (dH_ducR[iVar]) + (da_ducR[iVar]) * qn + (dqn_ducR[iVar]) * a;

  for (iVar = 0; iVar < nVar; ++iVar)
    for (jVar = 0; jVar < nVar; ++jVar)
      dFnducR[iVar][jVar] -= 0.5 * ws[2]*LdU[2]*dR3_ducR[iVar][jVar];

  //! dR4_ducR = dR(:,4)/dUcR
  //! Two shear wave components combined into 1 (wave strength incorporated).
  //! So, it is not really an eigenvector.
  //!
  //!                  Combined vector -> Differentiated
  //!
  //!  R(1,4) = 0                   -> 0
  //!  R(2,4) = du - dqn*nx            -> d(du)/dUcR - d(dqn)/dUcR*nx
  //!  R(3,4) = dv - dqn*ny            -> d(dv)/dUcR - d(dqn)/dUcR*ny
  //!  R(4,4) = dw - dqn*nz            -> d(dw)/dUcR - d(dqn)/dUcR*nz
  //!  R(5,4) = u*du+v*dv+w*dw-qn*dqn  -> du/dUcR*du     + d(du)/dUcR*u
  //!                                   + dv/dUcR*dv     + d(dv)/dUcR*v
  //!                                   + dw/dUcR*dw     + d(dw)/dUcR*w
  //!                                   - d(qn)/dUcR*dqn - d(dqn)/dUcR*qn


  for (iVar = 0; iVar < nVar; ++iVar) dR4_ducR[0][iVar] = 0;

  for (iDim = 0; iDim < nDim; ++iDim) {
    for (iVar = 0; iVar < nVar; ++iVar)
      dR4_ducR[iDim+1][iVar] = (ddvel_ducR[iDim][iVar]) - (ddqn_ducR[iVar])*n[iDim];
  }

  for (iVar = 0; iVar < nVar; ++iVar) {
    dR4_ducR[nDim+1][iVar] = - dqn*( dqn_ducR[iVar]) - qn*(ddqn_ducR[iVar]);
    for (iDim = 0; iDim < nDim; ++iDim) {
      dR4_ducR[nDim+1][iVar] += dVel[iDim]*( dvel_ducR[iDim][iVar]) + vel[iDim]*(ddvel_ducR[iDim][iVar]);
    }
  }

  for (iVar = 0; iVar < nVar; ++iVar)
    for (jVar = 0; jVar < nVar; ++jVar)
      dFnducR[iVar][jVar] -= 0.5 * ws[3]*LdU[3]*dR4_ducR[iVar][jVar];

  for (iVar = 0; iVar < nVar; ++iVar)
    for (jVar = 0; jVar < nVar; ++jVar) {
      Jacobian_i[iVar][jVar] = dFnducL[iVar][jVar] * Area;
      Jacobian_j[iVar][jVar] = dFnducR[iVar][jVar] * Area;
    }


  //!Dissipation Term: |An|(UR-UL) = R|Lambda|L*dU = sum_k of [ ws[k] * R(:,k) * L*dU[k] ]

  for (iVar = 0; iVar < nVar; ++iVar) {
    diss[iVar] = 0;
    for (int iWave = 0; iWave < 4; ++iWave)
      diss[iVar] += ws[iWave]*LdU[iWave]*R[iVar][iWave];
  }

  //!Compute the physical flux: fL = Fn(UL) and fR = Fn(UR)

  su2double flux_i[5], flux_j[5];

  ProjFlux_i[0] = rhoL*qnL;
  for (iDim = 0; iDim < nDim; ++iDim) {
    ProjFlux_i[iDim+1] = rhoL*qnL * velL[iDim] + pL*n[iDim];
  }
  ProjFlux_i[nDim+1] = rhoL*qnL * HL;

  ProjFlux_j[0] = rhoR*qnR;
  for (iDim = 0; iDim < nDim; ++iDim) {
    ProjFlux_j[iDim+1] = rhoR*qnR * velR[iDim] + pR*n[iDim];
  }
  ProjFlux_j[nDim+1] = rhoR*qnR * HR;

  //! This is the numerical flux: Roe flux = 1/2 *[  Fn(UL)+Fn(UR) - |An|(UR-UL) ]
  for (iVar = 0; iVar < nVar; ++iVar) {
    Flux[iVar] = 0.5 * (ProjFlux_i[iVar] + ProjFlux_j[iVar] - diss[iVar]) * Area;
  }


  /*--- Finalize in children class ---*/

  FinalizeResidual(Flux, Jacobian_i, Jacobian_j, config);

  /*--- Correct for grid motion ---*/

  if (dynamic_grid) {
    for (iVar = 0; iVar < nVar; iVar++) {
      Flux[iVar] -= ProjGridVel*Area * 0.5*(Conservatives_i[iVar]+Conservatives_j[iVar]);

      if (implicit) {
        Jacobian_i[iVar][iVar] -= 0.5*ProjGridVel*Area;
        Jacobian_j[iVar][iVar] -= 0.5*ProjGridVel*Area;
      }
    }
  }

  AD::SetPreaccOut(Flux, nVar);
  AD::EndPreacc();

  return ResidualType<>(Flux, Jacobian_i, Jacobian_j);

}

CUpwRoe_FlowNew::CUpwRoe_FlowNew(unsigned short val_nDim, unsigned short val_nVar, unsigned short val_nPrimVar,
                                 unsigned short val_nPrimVarGrad, const CConfig* config, bool val_low_dissipation) :
                                           CUpwRoeBase_FlowNew(val_nDim, val_nVar, val_nPrimVar, val_nPrimVarGrad, config, val_low_dissipation) {}

void CUpwRoe_FlowNew::FinalizeResidual(su2double *val_residual, su2double **val_Jacobian_i,
                                       su2double **val_Jacobian_j, const CConfig* config) {
  CorrectResidual(val_residual);
}
