/*!
 * \file roe.cpp
 * \brief Implementations of Roe-type schemes.
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

#include "../../../../include/numerics/flow/convection/roe.hpp"
#include "../../../../../Common/include/toolboxes/geometry_toolbox.hpp"

CUpwRoeBase_Flow::CUpwRoeBase_Flow(unsigned short val_nDim, unsigned short val_nVar, const CConfig* config,
                                   bool val_low_dissipation) : CNumerics(val_nDim, val_nVar, config) {

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

CUpwRoe_Flow::CUpwRoe_Flow(unsigned short val_nDim, unsigned short val_nVar, const CConfig* config,
              bool val_low_dissipation) : CUpwRoeBase_Flow(val_nDim, val_nVar, config, val_low_dissipation) {}

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

}

CUpwL2Roe_Flow::CUpwL2Roe_Flow(unsigned short val_nDim, unsigned short val_nVar, const CConfig* config) :
                CUpwRoeBase_Flow(val_nDim, val_nVar, config, false) {}

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

CUpwLMRoe_Flow::CUpwLMRoe_Flow(unsigned short val_nDim, unsigned short val_nVar, const CConfig* config) :
                CUpwRoeBase_Flow(val_nDim, val_nVar, config, false) {}

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

CUpwTurkel_Flow::CUpwTurkel_Flow(unsigned short val_nDim, unsigned short val_nVar, const CConfig* config) : CNumerics(val_nDim, val_nVar, config) {

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

CUpwGeneralRoe_Flow::CUpwGeneralRoe_Flow(unsigned short val_nDim, unsigned short val_nVar, const CConfig* config) : CNumerics(val_nDim, val_nVar, config) {

  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  /* A grid is defined as dynamic if there's rigid grid movement or grid deformation AND the problem is time domain */
  dynamic_grid = config->GetDynamic_Grid();
  kappa = config->GetRoe_Kappa(); // 1 is unstable


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
