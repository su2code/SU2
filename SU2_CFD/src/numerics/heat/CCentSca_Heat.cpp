/*!
 * \file CCentSca_Heat.cpp
 * \brief Implementation of numerics class CCentSca_Heat.
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

#include "../../../include/numerics/heat/CCentSca_Heat.hpp"

CCentSca_Heat::CCentSca_Heat(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {

  implicit        = (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT);
  /* A grid is defined as dynamic if there's rigid grid movement or grid deformation AND the problem is time domain */
  dynamic_grid = config->GetDynamic_Grid();

  MeanVelocity = new su2double [nDim];

  Laminar_Viscosity_i = config->GetViscosity_FreeStreamND();
  Laminar_Viscosity_j = config->GetViscosity_FreeStreamND();

  Param_Kappa_4 = config->GetKappa_4th_Heat();
  Diff_Lapl = new su2double [nVar];
}

CCentSca_Heat::~CCentSca_Heat(void) {

  delete [] MeanVelocity;
  delete [] Diff_Lapl;

}

void CCentSca_Heat::ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config) {

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

  ProjVelocity = 0.0; ProjVelocity_i = 0.0; ProjVelocity_j = 0.0; Area = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    ProjVelocity_i += V_i[iDim+1]*Normal[iDim];
    ProjVelocity_j += V_j[iDim+1]*Normal[iDim];
    Area += Normal[iDim]*Normal[iDim];
  }
  Area = sqrt(Area);

  /*--- Computing the second order centered scheme part ---*/

  ProjVelocity = 0.5*(ProjVelocity_i+ProjVelocity_j);

  val_residual[0] = 0.5*(Temp_i + Temp_j)*ProjVelocity;

  if (implicit) {
    val_Jacobian_i[0][0] = 0.5*ProjVelocity;
    val_Jacobian_j[0][0] = 0.5*ProjVelocity;
  }

  /*--- Adding artificial dissipation to stabilize the centered scheme ---*/

  Diff_Lapl[0] = Und_Lapl_i[0]-Und_Lapl_j[0];

  SoundSpeed_i = sqrt(ProjVelocity_i*ProjVelocity_i + (BetaInc2_i/DensityInc_i)*Area*Area);
  SoundSpeed_j = sqrt(ProjVelocity_j*ProjVelocity_j + (BetaInc2_j/DensityInc_j)*Area*Area);

  Local_Lambda_i = fabs(ProjVelocity_i)+SoundSpeed_i;
  Local_Lambda_j = fabs(ProjVelocity_j)+SoundSpeed_j;
  MeanLambda = 0.5*(Local_Lambda_i+Local_Lambda_j);

  val_residual[0] += -Param_Kappa_4*Diff_Lapl[0]*MeanLambda;

  if (implicit) {
    cte_0 = Param_Kappa_4*su2double(Neighbor_i+1)*MeanLambda;
    cte_1 = Param_Kappa_4*su2double(Neighbor_j+1)*MeanLambda;

    val_Jacobian_i[0][0] += cte_0;
    val_Jacobian_j[0][0] -= cte_1;
  }

  AD::SetPreaccOut(val_residual[0]);
  AD::EndPreacc();

}
