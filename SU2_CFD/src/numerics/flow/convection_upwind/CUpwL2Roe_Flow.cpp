/*!
 * \file CUpwL2Roe_Flow.cpp
 * \brief Implementation of numerics class CUpwL2Roe_Flow.
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

#include "../../../../include/numerics/flow/convection_upwind/CUpwL2Roe_Flow.hpp"

CUpwL2Roe_Flow::CUpwL2Roe_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) :
                CUpwRoeBase_Flow(val_nDim, val_nVar, config, false) {}

CUpwL2Roe_Flow::~CUpwL2Roe_Flow() {}

void CUpwL2Roe_Flow::FinalizeResidual(su2double *val_residual, su2double **val_Jacobian_i,
                                      su2double **val_Jacobian_j, CConfig *config) {

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

  su2double proj_delta_vel = 0.0, delta_vel[3];
  for (iDim = 0; iDim < nDim; iDim++) {
    delta_vel[iDim] = Velocity_j[iDim] - Velocity_i[iDim];
    proj_delta_vel += delta_vel[iDim]*UnitNormal[iDim];
  }
  proj_delta_vel *= zeta;
  su2double delta_p = Pressure_j - Pressure_i;
  su2double delta_rho = Density_j - Density_i;

  su2double delta_wave[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
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
