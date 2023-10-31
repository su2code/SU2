/*!
 * \file fds.cpp
 * \brief Implementation of Flux-Difference-Splitting schemes.
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

#include "../../../../include/numerics/flow/convection/fds.hpp"
#include "../../../../../Common/include/toolboxes/geometry_toolbox.hpp"

CUpwFDSInc_Flow::CUpwFDSInc_Flow(unsigned short val_nDim, unsigned short val_nVar, const CConfig *config) : CNumerics(val_nDim, val_nVar, config) {

  implicit         = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  variable_density = (config->GetVariable_Density_Model());
  energy           = config->GetEnergy_Equation();
  /* A grid is defined as dynamic if there's rigid grid movement or grid deformation AND the problem is time domain */
  dynamic_grid = config->GetDynamic_Grid();

  Flux         = new su2double[nVar];
  Diff_V       = new su2double[nVar];
  Velocity_i   = new su2double[nDim];
  Velocity_j   = new su2double[nDim];
  MeanVelocity = new su2double[nDim];
  ProjFlux_i   = new su2double[nVar];
  ProjFlux_j   = new su2double[nVar];
  Lambda       = new su2double[nVar];
  Epsilon      = new su2double[nVar];
  Precon       = new su2double*[nVar];
  invPrecon_A  = new su2double*[nVar];
  Jacobian_i   = new su2double*[nVar];
  Jacobian_j   = new su2double*[nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    Precon[iVar]      = new su2double[nVar];
    invPrecon_A[iVar] = new su2double[nVar];
    Jacobian_i[iVar]  = new su2double[nVar];
    Jacobian_j[iVar]  = new su2double[nVar];
  }

}

CUpwFDSInc_Flow::~CUpwFDSInc_Flow() {

  delete [] Flux;
  delete [] Diff_V;
  delete [] Velocity_i;
  delete [] Velocity_j;
  delete [] MeanVelocity;
  delete [] ProjFlux_i;
  delete [] ProjFlux_j;
  delete [] Lambda;
  delete [] Epsilon;

  for (iVar = 0; iVar < nVar; iVar++) {
    delete [] Precon[iVar];
    delete [] invPrecon_A[iVar];
  }
  delete [] Precon;
  delete [] invPrecon_A;

  if (Jacobian_i != nullptr) {
    for (iVar = 0; iVar < nVar; iVar++) {
      delete [] Jacobian_i[iVar];
      delete [] Jacobian_j[iVar];
    }
    delete [] Jacobian_i;
    delete [] Jacobian_j;
  }

}

CNumerics::ResidualType<> CUpwFDSInc_Flow::ComputeResidual(const CConfig *config) {

  su2double U_i[5] = {0.0,0.0,0.0,0.0,0.0}, U_j[5] = {0.0,0.0,0.0,0.0,0.0};
  su2double ProjGridVel = 0.0;

  AD::StartPreacc();
  AD::SetPreaccIn(V_i, nDim+9); AD::SetPreaccIn(V_j, nDim+9); AD::SetPreaccIn(Normal, nDim);
  if (dynamic_grid) {
    AD::SetPreaccIn(GridVel_i, nDim);
    AD::SetPreaccIn(GridVel_j, nDim);
  }

  /*--- Face area (norm or the normal vector) ---*/

  Area = GeometryToolbox::Norm(nDim, Normal);

  /*--- Compute and unitary normal vector ---*/

  for (iDim = 0; iDim < nDim; iDim++) {
    UnitNormal[iDim] = Normal[iDim]/Area;
    if (fabs(UnitNormal[iDim]) < EPS) UnitNormal[iDim] = EPS;
  }

  /*--- Set primitive variables at points iPoint and jPoint ---*/

  Pressure_i    = V_i[0];             Pressure_j    = V_j[0];
  Temperature_i = V_i[nDim+1];        Temperature_j = V_j[nDim+1];
  DensityInc_i  = V_i[nDim+2];        DensityInc_j  = V_j[nDim+2];
  BetaInc2_i    = V_i[nDim+3];        BetaInc2_j    = V_j[nDim+3];
  Cp_i          = V_i[nDim+7];        Cp_j          = V_j[nDim+7];
  Enthalpy_i    = Cp_i*Temperature_i; Enthalpy_j    = Cp_j*Temperature_j;

  ProjVelocity = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Velocity_i[iDim]    = V_i[iDim+1];
    Velocity_j[iDim]    = V_j[iDim+1];
    MeanVelocity[iDim]  = 0.5*(Velocity_i[iDim] + Velocity_j[iDim]);
    ProjVelocity       += MeanVelocity[iDim]*Normal[iDim];
  }

  /*--- Projected velocity adjustment due to mesh motion ---*/

  if (dynamic_grid) {
    ProjGridVel = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      ProjGridVel   += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*Normal[iDim];
    }
    ProjVelocity   -= ProjGridVel;
  }

  /*--- Mean variables at points iPoint and jPoint ---*/

  MeanDensity     = 0.5*(DensityInc_i  + DensityInc_j);
  MeanPressure    = 0.5*(Pressure_i    + Pressure_j);
  MeanBetaInc2    = 0.5*(BetaInc2_i    + BetaInc2_j);
  MeanEnthalpy    = 0.5*(Enthalpy_i    + Enthalpy_j);
  MeanCp          = 0.5*(Cp_i          + Cp_j);
  MeanTemperature = 0.5*(Temperature_i + Temperature_j);

  /*--- Artificial sound speed based on eigs of preconditioned system ---*/

  MeanSoundSpeed = sqrt(MeanBetaInc2*Area*Area);

  /*--- We need the derivative of the equation of state to build the
   preconditioning matrix. For now, the only option is the ideal gas
   law, but in the future, dRhodT should be in the fluid model. ---*/

  MeandRhodT = 0.0; dRhodT_i = 0.0; dRhodT_j = 0.0;
  if (variable_density) {
    MeandRhodT = -MeanDensity/MeanTemperature;
    dRhodT_i   = -DensityInc_i/Temperature_i;
    dRhodT_j   = -DensityInc_j/Temperature_j;
  }

  /*--- Compute ProjFlux_i ---*/

  GetInviscidIncProjFlux(&DensityInc_i, Velocity_i, &Pressure_i, &BetaInc2_i, &Enthalpy_i, Normal, ProjFlux_i);

  /*--- Compute ProjFlux_j ---*/

  GetInviscidIncProjFlux(&DensityInc_j, Velocity_j, &Pressure_j, &BetaInc2_j, &Enthalpy_j, Normal, ProjFlux_j);

  /*--- Eigenvalues of the preconditioned system ---*/

  if (nDim == 2) {
    Lambda[0] = ProjVelocity;
    Lambda[1] = ProjVelocity;
    Lambda[2] = ProjVelocity - MeanSoundSpeed;
    Lambda[3] = ProjVelocity + MeanSoundSpeed;
  }
  if (nDim == 3) {
    Lambda[0] = ProjVelocity;
    Lambda[1] = ProjVelocity;
    Lambda[2] = ProjVelocity;
    Lambda[3] = ProjVelocity - MeanSoundSpeed;
    Lambda[4] = ProjVelocity + MeanSoundSpeed;
  }

  /*--- Absolute value of the eigenvalues ---*/

  for (iVar = 0; iVar < nVar; iVar++)
    Lambda[iVar] = fabs(Lambda[iVar]);

  /*--- Build the preconditioning matrix using mean values ---*/

  GetPreconditioner(&MeanDensity, MeanVelocity, &MeanBetaInc2, &MeanCp, &MeanTemperature, &MeandRhodT, Precon);

  /*--- Build the absolute value of the preconditioned Jacobian, i.e.,
   |A_precon| = P x |Lambda| x inv(P), where P diagonalizes the matrix
   inv(Precon) x dF/dV and Lambda is the diag. matrix of its eigenvalues. ---*/

  GetPreconditionedProjJac(&MeanDensity, Lambda, &MeanBetaInc2, UnitNormal, invPrecon_A);

  /*--- Difference of primitive variables at iPoint and jPoint ---*/

  Diff_V[0] = Pressure_j - Pressure_i;
  for (iDim = 0; iDim < nDim; iDim++)
    Diff_V[iDim+1] = Velocity_j[iDim] - Velocity_i[iDim];
  Diff_V[nDim+1] = Temperature_j - Temperature_i;

  /*--- Build the inviscid Jacobian w.r.t. the primitive variables ---*/

  if (implicit) {
    GetInviscidIncProjJac(&DensityInc_i, Velocity_i, &BetaInc2_i, &Cp_i, &Temperature_i, &dRhodT_i, Normal, 0.5, Jacobian_i);
    GetInviscidIncProjJac(&DensityInc_j, Velocity_j, &BetaInc2_j, &Cp_j, &Temperature_j, &dRhodT_j, Normal, 0.5, Jacobian_j);
  }

  /*--- Compute dissipation as Precon x |A_precon| x dV. If implicit,
   store Precon x |A_precon| from dissipation term. ---*/

  for (iVar = 0; iVar < nVar; iVar++) {
    Flux[iVar] = 0.5*(ProjFlux_i[iVar]+ProjFlux_j[iVar]);
    for (jVar = 0; jVar < nVar; jVar++) {
      Proj_ModJac_Tensor_ij = 0.0;
      for (kVar = 0; kVar < nVar; kVar++)
        Proj_ModJac_Tensor_ij += Precon[iVar][kVar]*invPrecon_A[kVar][jVar];
      Flux[iVar] -= 0.5*Proj_ModJac_Tensor_ij*Diff_V[jVar];
      if (implicit) {
        Jacobian_i[iVar][jVar] += 0.5*Proj_ModJac_Tensor_ij;
        Jacobian_j[iVar][jVar] -= 0.5*Proj_ModJac_Tensor_ij;
      }
    }
  }

  /*--- Corrections due to grid motion ---*/
  if (dynamic_grid) {

    /*--- Recompute conservative variables ---*/

    U_i[0] = DensityInc_i; U_j[0] = DensityInc_j;
    for (iDim = 0; iDim < nDim; iDim++) {
      U_i[iDim+1] = DensityInc_i*Velocity_i[iDim]; U_j[iDim+1] = DensityInc_j*Velocity_j[iDim];
    }
    U_i[nDim+1] = DensityInc_i*Enthalpy_i; U_j[nDim+1] = DensityInc_j*Enthalpy_j;

    ProjVelocity = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      ProjVelocity += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*Normal[iDim];

    /*--- Residual contributions ---*/
    for (iVar = 0; iVar < nVar; iVar++) {
      Flux[iVar] -= ProjVelocity * 0.5*(U_i[iVar]+U_j[iVar]);

      /*--- Jacobian contributions ---*/
      /*--- Implicit terms ---*/
      if (implicit) {
        for (iDim = 0; iDim < nDim; iDim++){
          Jacobian_i[iDim+1][iDim+1] -= 0.5*ProjVelocity*DensityInc_i;
          Jacobian_j[iDim+1][iDim+1] -= 0.5*ProjVelocity*DensityInc_j;
        }
        Jacobian_i[nDim+1][nDim+1] -= 0.5*ProjVelocity*DensityInc_i*Cp_i;
        Jacobian_j[nDim+1][nDim+1] -= 0.5*ProjVelocity*DensityInc_j*Cp_j;
      }
    }
  }

  if (!energy) {
    Flux[nDim+1] = 0.0;
    if (implicit) {
      for (iVar = 0; iVar < nVar; iVar++) {
        Jacobian_i[iVar][nDim+1] = 0.0;
        Jacobian_j[iVar][nDim+1] = 0.0;

        Jacobian_i[nDim+1][iVar] = 0.0;
        Jacobian_j[nDim+1][iVar] = 0.0;
      }
    }
  }

  AD::SetPreaccOut(Flux, nVar);
  AD::EndPreacc();

  return ResidualType<>(Flux, Jacobian_i, Jacobian_j);

}
