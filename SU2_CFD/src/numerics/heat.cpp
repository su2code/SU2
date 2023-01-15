/*!
 * \file heat.cpp
 * \brief Implementation of numerics classes for heat transfer.
 * \author F. Palacios, T. Economon
 * \version 7.5.0 "Blackbird"
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
