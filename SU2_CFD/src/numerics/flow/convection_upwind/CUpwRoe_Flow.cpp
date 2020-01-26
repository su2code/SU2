/*!
 * \file CUpwRoe_Flow.cpp
 * \brief Implementation of numerics class CUpwRoe_Flow.
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

#include "../../../../include/numerics/flow/convection_upwind/CUpwRoe_Flow.hpp"

CUpwRoe_Flow::CUpwRoe_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config,
              bool val_low_dissipation) : CUpwRoeBase_Flow(val_nDim, val_nVar, config, val_low_dissipation) {}

CUpwRoe_Flow::~CUpwRoe_Flow() {}

void CUpwRoe_Flow::FinalizeResidual(su2double *val_residual, su2double **val_Jacobian_i,
                                    su2double **val_Jacobian_j, CConfig *config) {

  unsigned short iVar, jVar, kVar;

  /*--- Compute inverse P tensor ---*/
  GetPMatrix_inv(&RoeDensity, RoeVelocity, &RoeSoundSpeed, UnitNormal, invP_Tensor);

  /*--- Diference between conservative variables at jPoint and iPoint ---*/
  for (iVar = 0; iVar < nVar; iVar++)
    Diff_U[iVar] = Conservatives_j[iVar]-Conservatives_i[iVar];

  /*--- Low dissipation formulation ---*/
  if (roe_low_dissipation)
    SetRoe_Dissipation(Dissipation_i, Dissipation_j, Sensor_i, Sensor_j, Dissipation_ij, config);
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
