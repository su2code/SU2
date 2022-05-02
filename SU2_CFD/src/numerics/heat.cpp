/*!
 * \file heat.cpp
 * \brief Implementation of numerics classes for heat transfer.
 * \author F. Palacios, T. Economon
 * \version 7.3.1 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2022, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../include/numerics/heat.hpp"
#include "../../../Common/include/toolboxes/geometry_toolbox.hpp"

CCentSca_Heat::CCentSca_Heat(unsigned short val_nDim, unsigned short val_nVar, const CConfig *config) :
               CNumerics(val_nDim, val_nVar, config) {

  implicit = (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT);
  dynamic_grid = config->GetDynamic_Grid();
  Param_Kappa_4 = config->GetKappa_4th_Heat();

}

void CCentSca_Heat::ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i,
                                    su2double **val_Jacobian_j, CConfig *config) {
  AD::StartPreacc();
  AD::SetPreaccIn(V_i, nDim+3); AD::SetPreaccIn(V_j, nDim+3);
  AD::SetPreaccIn(Temp_i); AD::SetPreaccIn(Temp_j);
  AD::SetPreaccIn(Und_Lapl_i, nVar); AD::SetPreaccIn(Und_Lapl_j, nVar);
  AD::SetPreaccIn(Normal, nDim);
  if (dynamic_grid) {
    AD::SetPreaccIn(GridVel_i, nDim); AD::SetPreaccIn(GridVel_j, nDim);
  }

  /*--- Primitive variables at point i and j ---*/

  Pressure_i =    V_i[0];       Pressure_j = V_j[0];
  DensityInc_i =  V_i[nDim+2];  DensityInc_j = V_j[nDim+2];
  BetaInc2_i =    V_i[nDim+3];  BetaInc2_j = V_j[nDim+3];

  /*--- Projected velocities at the current edge ---*/

  su2double ProjVelocity_i = 0.0;
  su2double ProjVelocity_j = 0.0;

  if (dynamic_grid) {
    for (unsigned short iDim = 0; iDim < nDim; iDim++) {
      su2double Velocity_i = V_i[iDim+1] - GridVel_i[iDim];
      su2double Velocity_j = V_j[iDim+1] - GridVel_j[iDim];
      ProjVelocity_i += Velocity_i*Normal[iDim];
      ProjVelocity_j += Velocity_j*Normal[iDim];
    }
  }
  else {
    for (unsigned short iDim = 0; iDim < nDim; iDim++) {
      ProjVelocity_i += V_i[iDim+1]*Normal[iDim];
      ProjVelocity_j += V_j[iDim+1]*Normal[iDim];
    }
  }

  /*--- Computing the second order centered scheme part ---*/

  const su2double ProjVelocity = 0.5*(ProjVelocity_i+ProjVelocity_j);

  val_residual[0] = 0.5*(Temp_i + Temp_j)*ProjVelocity;

  if (implicit) {
    val_Jacobian_i[0][0] = 0.5*ProjVelocity;
    val_Jacobian_j[0][0] = 0.5*ProjVelocity;
  }

  /*--- Adding artificial dissipation to stabilize the centered scheme ---*/

  const su2double Diff_Lapl = Und_Lapl_i[0]-Und_Lapl_j[0];
  const su2double Area2 = GeometryToolbox::SquaredNorm(nDim, Normal);

  const su2double SoundSpeed_i = sqrt(pow(ProjVelocity_i,2) + (BetaInc2_i/DensityInc_i)*Area2);
  const su2double SoundSpeed_j = sqrt(pow(ProjVelocity_j,2) + (BetaInc2_j/DensityInc_j)*Area2);

  const su2double Local_Lambda_i = fabs(ProjVelocity_i)+SoundSpeed_i;
  const su2double Local_Lambda_j = fabs(ProjVelocity_j)+SoundSpeed_j;
  const su2double MeanLambda = 0.5*(Local_Lambda_i+Local_Lambda_j);

  val_residual[0] += -Param_Kappa_4*Diff_Lapl*MeanLambda;

  if (implicit) {
    const su2double cte_0 = Param_Kappa_4*su2double(Neighbor_i+1)*MeanLambda;
    const su2double cte_1 = Param_Kappa_4*su2double(Neighbor_j+1)*MeanLambda;

    val_Jacobian_i[0][0] += cte_0;
    val_Jacobian_j[0][0] -= cte_1;
  }

  AD::SetPreaccOut(val_residual[0]);
  AD::EndPreacc();

}

CUpwSca_Heat::CUpwSca_Heat(unsigned short val_nDim, unsigned short val_nVar, const CConfig *config) :
              CNumerics(val_nDim, val_nVar, config) {

  implicit = (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT);
  dynamic_grid = config->GetDynamic_Grid();

}

void CUpwSca_Heat::ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i,
                                   su2double **val_Jacobian_j, CConfig *config) {

  AD::StartPreacc();
  AD::SetPreaccIn(V_i, nDim+1); AD::SetPreaccIn(V_j, nDim+1);
  AD::SetPreaccIn(Temp_i); AD::SetPreaccIn(Temp_j);
  AD::SetPreaccIn(Normal, nDim);
  if (dynamic_grid) {
    AD::SetPreaccIn(GridVel_i, nDim); AD::SetPreaccIn(GridVel_j, nDim);
  }

  su2double q_ij = 0.0;

  if (dynamic_grid) {
    for (unsigned short iDim = 0; iDim < nDim; iDim++) {
      const su2double Velocity_i = V_i[iDim+1] - GridVel_i[iDim];
      const su2double Velocity_j = V_j[iDim+1] - GridVel_j[iDim];
      q_ij += 0.5*(Velocity_i+Velocity_j)*Normal[iDim];
    }
  }
  else {
    for (unsigned short iDim = 0; iDim < nDim; iDim++) {
      q_ij += 0.5*(V_i[iDim+1]+V_j[iDim+1])*Normal[iDim];
    }
  }

  const su2double a0 = 0.5*(q_ij+fabs(q_ij));
  const su2double a1 = 0.5*(q_ij-fabs(q_ij));

  val_residual[0] = a0*Temp_i+a1*Temp_j;

  if (implicit) {
    val_Jacobian_i[0][0] = a0;
    val_Jacobian_j[0][0] = a1;
  }

  AD::SetPreaccOut(val_residual[0]);
  AD::EndPreacc();

}

CAvgGrad_Heat::CAvgGrad_Heat(unsigned short val_nDim, unsigned short val_nVar,
                             const CConfig *config, bool correct_) :
               CNumerics(val_nDim, val_nVar, config) {

  implicit = (config->GetKind_TimeIntScheme_Heat() == EULER_IMPLICIT);
  correct = correct_;
}

void CAvgGrad_Heat::ComputeResidual(su2double *val_residual, su2double **Jacobian_i,
                                    su2double **Jacobian_j, CConfig *config) {
  constexpr int nVar = 1;

  AD::StartPreacc();
  AD::SetPreaccIn(Coord_i, nDim); AD::SetPreaccIn(Coord_j, nDim);
  AD::SetPreaccIn(Normal, nDim);
  AD::SetPreaccIn(Temp_i); AD::SetPreaccIn(Temp_j);
  AD::SetPreaccIn(ConsVar_Grad_i[0],nDim); AD::SetPreaccIn(ConsVar_Grad_j[0],nDim);
  AD::SetPreaccIn(Thermal_Diffusivity_i); AD::SetPreaccIn(Thermal_Conductivity_j);

  su2double NormalGrad[nVar], CorrectedGrad[nVar];

  auto proj_vector_ij = ComputeProjectedGradient(nDim, nVar, Normal, Coord_i, Coord_j, ConsVar_Grad_i,
                                                 ConsVar_Grad_j, correct, &Temp_i, &Temp_j,
                                                 NormalGrad, CorrectedGrad);

  const su2double Thermal_Diffusivity_Mean = 0.5*(Thermal_Diffusivity_i + Thermal_Diffusivity_j);

  val_residual[0] = Thermal_Diffusivity_Mean*CorrectedGrad[0];

  /*--- For Jacobians -> Use of TSL approx. to compute derivatives of the gradients ---*/
  if (implicit) {
    Jacobian_i[0][0] = -Thermal_Diffusivity_Mean*proj_vector_ij;
    Jacobian_j[0][0] = Thermal_Diffusivity_Mean*proj_vector_ij;
  }

  AD::SetPreaccOut(val_residual, nVar);
  AD::EndPreacc();

}
