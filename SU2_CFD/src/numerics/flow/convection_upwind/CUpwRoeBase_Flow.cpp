/*!
 * \file CUpwRoeBase_Flow.cpp
 * \brief Implementation of numerics class CUpwRoeBase_Flow.
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

#include "../../../../include/numerics/flow/convection_upwind/CUpwRoeBase_Flow.hpp"

CUpwRoeBase_Flow::CUpwRoeBase_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config,
                                   bool val_low_dissipation) : CNumerics(val_nDim, val_nVar, config) {

  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  /* A grid is defined as dynamic if there's rigid grid movement or grid deformation AND the problem is time domain */
  dynamic_grid = config->GetDynamic_Grid();
  kappa = config->GetRoe_Kappa(); // 1 is unstable

  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;

  roe_low_dissipation = val_low_dissipation;

  Diff_U = new su2double [nVar];
  Velocity_i = new su2double [nDim];
  Velocity_j = new su2double [nDim];
  RoeVelocity = new su2double [nDim];
  ProjFlux_i = new su2double [nVar];
  ProjFlux_j = new su2double [nVar];
  Conservatives_i = new su2double [nVar];
  Conservatives_j = new su2double [nVar];
  Lambda = new su2double [nVar];
  P_Tensor = new su2double* [nVar];
  invP_Tensor = new su2double* [nVar];
  for (unsigned short iVar = 0; iVar < nVar; iVar++) {
    P_Tensor[iVar] = new su2double [nVar];
    invP_Tensor[iVar] = new su2double [nVar];
  }
}

CUpwRoeBase_Flow::~CUpwRoeBase_Flow(void) {

  delete [] Diff_U;
  delete [] Velocity_i;
  delete [] Velocity_j;
  delete [] RoeVelocity;
  delete [] ProjFlux_i;
  delete [] ProjFlux_j;
  delete [] Conservatives_i;
  delete [] Conservatives_j;
  delete [] Lambda;
  for (unsigned short iVar = 0; iVar < nVar; iVar++) {
    delete [] P_Tensor[iVar];
    delete [] invP_Tensor[iVar];
  }
  delete [] P_Tensor;
  delete [] invP_Tensor;

}

void CUpwRoeBase_Flow::FinalizeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config) {
/*---
 CUpwRoeBase_Flow::ComputeResidual initializes the residual (flux) and its Jacobians with the standard Roe averaging
 fc_{1/2} = kappa*(fc_i+fc_j)*Normal. It then calls this method, which derived classes specialize, to account for
 the dissipation part.
---*/
}

void CUpwRoeBase_Flow::ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config) {

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

  Area = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    Area += Normal[iDim]*Normal[iDim];
  Area = sqrt(Area);

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
      val_residual[iVar] = 0.0;
      if (implicit){
        for (jVar = 0; jVar < nVar; jVar++) {
          val_Jacobian_i[iVar][jVar] = 0.0;
          val_Jacobian_j[iVar][jVar] = 0.0;
        }
      }
    }
    AD::SetPreaccOut(val_residual, nVar);
    AD::EndPreacc();
    return;
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
    val_residual[iVar] = kappa*(ProjFlux_i[iVar]+ProjFlux_j[iVar]);

  if (implicit) {
    GetInviscidProjJac(Velocity_i, &Energy_i, Normal, kappa, val_Jacobian_i);
    GetInviscidProjJac(Velocity_j, &Energy_j, Normal, kappa, val_Jacobian_j);
  }

  /*--- Finalize in children class ---*/

  FinalizeResidual(val_residual, val_Jacobian_i, val_Jacobian_j, config);

  /*--- Correct for grid motion ---*/

  if (dynamic_grid) {
    for (iVar = 0; iVar < nVar; iVar++) {
      val_residual[iVar] -= ProjGridVel*Area * 0.5*(Conservatives_i[iVar]+Conservatives_j[iVar]);

      if (implicit) {
        val_Jacobian_i[iVar][iVar] -= 0.5*ProjGridVel*Area;
        val_Jacobian_j[iVar][iVar] -= 0.5*ProjGridVel*Area;
      }
    }
  }

  AD::SetPreaccOut(val_residual, nVar);
  AD::EndPreacc();

}
