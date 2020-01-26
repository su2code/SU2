/*!
 * \file CCentBase_Flow.cpp
 * \brief Implementation of numerics class CCentBase_Flow.
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

#include "../../../../include/numerics/flow/convection_centered/CCentBase_Flow.hpp"

CCentBase_Flow::CCentBase_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) :
                CNumerics(val_nDim, val_nVar, config) {

  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  /* A grid is defined as dynamic if there's rigid grid movement or grid deformation AND the problem is time domain */
  dynamic_grid = config->GetDynamic_Grid();
  fix_factor = config->GetCent_Jac_Fix_Factor();

  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;

  /*--- Allocate required structures ---*/
  Diff_U = new su2double [nVar];
  Diff_Lapl = new su2double [nVar];
  Velocity_i = new su2double [nDim];
  Velocity_j = new su2double [nDim];
  MeanVelocity = new su2double [nDim];
  ProjFlux = new su2double [nVar];
}

CCentBase_Flow::~CCentBase_Flow(void) {
  delete [] Diff_U;
  delete [] Diff_Lapl;
  delete [] Velocity_i;
  delete [] Velocity_j;
  delete [] MeanVelocity;
  delete [] ProjFlux;
}

void CCentBase_Flow::ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j,
                                    CConfig *config) {

  su2double U_i[5] = {0.0,0.0,0.0,0.0,0.0}, U_j[5] = {0.0,0.0,0.0,0.0,0.0};

  bool preacc = SetPreaccInVars();

  if (preacc) {
    AD::SetPreaccIn(Normal, nDim);
    AD::SetPreaccIn(V_i, nDim+5); AD::SetPreaccIn(V_j, nDim+5);
    AD::SetPreaccIn(Lambda_i);    AD::SetPreaccIn(Lambda_j);
    if (dynamic_grid) {
      AD::SetPreaccIn(GridVel_i, nDim); AD::SetPreaccIn(GridVel_j, nDim);
    }
  }

  /*--- Pressure, density, enthalpy, energy, and velocity at points i and j ---*/

  Pressure_i = V_i[nDim+1];                       Pressure_j = V_j[nDim+1];
  Density_i  = V_i[nDim+2];                       Density_j  = V_j[nDim+2];
  Enthalpy_i = V_i[nDim+3];                       Enthalpy_j = V_j[nDim+3];
  SoundSpeed_i = V_i[nDim+4];                     SoundSpeed_j = V_j[nDim+4];
  Energy_i = Enthalpy_i - Pressure_i/Density_i;   Energy_j = Enthalpy_j - Pressure_j/Density_j;

  sq_vel_i = 0.0; sq_vel_j = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Velocity_i[iDim] = V_i[iDim+1];
    Velocity_j[iDim] = V_j[iDim+1];
    sq_vel_i += 0.5*Velocity_i[iDim]*Velocity_i[iDim];
    sq_vel_j += 0.5*Velocity_j[iDim]*Velocity_j[iDim];
  }

  /*--- Recompute conservative variables ---*/

  U_i[0] = Density_i; U_j[0] = Density_j;
  for (iDim = 0; iDim < nDim; iDim++) {
    U_i[iDim+1] = Density_i*Velocity_i[iDim]; U_j[iDim+1] = Density_j*Velocity_j[iDim];
  }
  U_i[nDim+1] = Density_i*Energy_i; U_j[nDim+1] = Density_j*Energy_j;

  /*--- Compute mean values of the variables ---*/

  MeanDensity = 0.5*(Density_i+Density_j);
  MeanPressure = 0.5*(Pressure_i+Pressure_j);
  MeanEnthalpy = 0.5*(Enthalpy_i+Enthalpy_j);
  for (iDim = 0; iDim < nDim; iDim++)
    MeanVelocity[iDim] =  0.5*(Velocity_i[iDim]+Velocity_j[iDim]);
  MeanEnergy = 0.5*(Energy_i+Energy_j);

  /*--- Get projected flux tensor ---*/

  GetInviscidProjFlux(&MeanDensity, MeanVelocity, &MeanPressure, &MeanEnthalpy, Normal, ProjFlux);

  /*--- Residual of the inviscid flux ---*/

  for (iVar = 0; iVar < nVar; iVar++)
    val_residual[iVar] = ProjFlux[iVar];

  /*--- Jacobians of the inviscid flux, scale = 0.5 because val_residual ~ 0.5*(fc_i+fc_j)*Normal ---*/

  if (implicit) {
    GetInviscidProjJac(MeanVelocity, &MeanEnergy, Normal, 0.5, val_Jacobian_i);
    for (iVar = 0; iVar < nVar; iVar++)
      for (jVar = 0; jVar < nVar; jVar++)
        val_Jacobian_j[iVar][jVar] = val_Jacobian_i[iVar][jVar];
  }

  /*--- Adjustment due to grid motion ---*/

  if (dynamic_grid) {
    ProjGridVel = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      ProjGridVel += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*Normal[iDim];

    for (iVar = 0; iVar < nVar; iVar++) {
      val_residual[iVar] -= ProjGridVel * 0.5*(U_i[iVar] + U_j[iVar]);
      if (implicit) {
        val_Jacobian_i[iVar][iVar] -= 0.5*ProjGridVel;
        val_Jacobian_j[iVar][iVar] -= 0.5*ProjGridVel;
      }
    }
  }

  /*--- Compute the local spectral radius and the stretching factor ---*/

  ProjVelocity_i = 0.0; ProjVelocity_j = 0.0; Area = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    ProjVelocity_i += Velocity_i[iDim]*Normal[iDim];
    ProjVelocity_j += Velocity_j[iDim]*Normal[iDim];
    Area += Normal[iDim]*Normal[iDim];
  }
  Area = sqrt(Area);

  /*--- Adjustment due to mesh motion ---*/

  if (dynamic_grid) {
    ProjVelocity_i -= ProjGridVel;
    ProjVelocity_j -= ProjGridVel;
  }

  /*--- Dissipation term ---*/

  Local_Lambda_i = (fabs(ProjVelocity_i)+SoundSpeed_i*Area);
  Local_Lambda_j = (fabs(ProjVelocity_j)+SoundSpeed_j*Area);
  MeanLambda = 0.5*(Local_Lambda_i+Local_Lambda_j);

  Phi_i = pow(Lambda_i/(4.0*MeanLambda), Param_p);
  Phi_j = pow(Lambda_j/(4.0*MeanLambda), Param_p);
  StretchingFactor = 4.0*Phi_i*Phi_j/(Phi_i+Phi_j);

  /*--- Compute differences btw. conservative variables, with a correction for enthalpy ---*/

  for (iVar = 0; iVar < nVar-1; iVar++) {
    Diff_U[iVar] = U_i[iVar]-U_j[iVar];
  }
  Diff_U[nVar-1] = Density_i*Enthalpy_i-Density_j*Enthalpy_j;

  DissipationTerm(val_residual, val_Jacobian_i, val_Jacobian_j);

  if (preacc) {
    AD::SetPreaccOut(val_residual, nVar);
    AD::EndPreacc();
  }
}

void CCentBase_Flow::ScalarDissipationJacobian(su2double **val_Jacobian_i, su2double **val_Jacobian_j) {

  /*--- n-1 diagonal entries ---*/

  for (iVar = 0; iVar < (nVar-1); iVar++) {
    val_Jacobian_i[iVar][iVar] += fix_factor*cte_0;
    val_Jacobian_j[iVar][iVar] -= fix_factor*cte_1;
  }

  /*--- Last row of Jacobian_i ---*/

  val_Jacobian_i[nVar-1][0] += fix_factor*cte_0*Gamma_Minus_One*sq_vel_i;
  for (iDim = 0; iDim < nDim; iDim++)
    val_Jacobian_i[nVar-1][iDim+1] -= fix_factor*cte_0*Gamma_Minus_One*Velocity_i[iDim];
  val_Jacobian_i[nVar-1][nVar-1] += fix_factor*cte_0*Gamma;

  /*--- Last row of Jacobian_j ---*/

  val_Jacobian_j[nVar-1][0] -= fix_factor*cte_1*Gamma_Minus_One*sq_vel_j;
  for (iDim = 0; iDim < nDim; iDim++)
    val_Jacobian_j[nVar-1][iDim+1] += fix_factor*cte_1*Gamma_Minus_One*Velocity_j[iDim];
  val_Jacobian_j[nVar-1][nVar-1] -= fix_factor*cte_1*Gamma;

}
