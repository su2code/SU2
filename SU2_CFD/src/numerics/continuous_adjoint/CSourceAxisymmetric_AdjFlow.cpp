/*!
 * \file CSourceAxisymmetric_AdjFlow.cpp
 * \brief Implementation of numerics class CSourceAxisymmetric_AdjFlow.
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

#include "../../../include/numerics/continuous_adjoint/CSourceAxisymmetric_AdjFlow.hpp"

CSourceAxisymmetric_AdjFlow::CSourceAxisymmetric_AdjFlow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) { }

CSourceAxisymmetric_AdjFlow::~CSourceAxisymmetric_AdjFlow(void) { }

void CSourceAxisymmetric_AdjFlow::ComputeResidual(su2double *val_residual, su2double **Jacobian_ii, CConfig *config) {

  su2double yinv;
  su2double Jacobian_Axisymmetric[4][4];

  if (Coord_i[1] > 0.0) yinv = 1.0/Coord_i[1];
  else yinv = 0.0;

  Jacobian_Axisymmetric[0][0] = 0;
  Jacobian_Axisymmetric[0][1] = 0;
  Jacobian_Axisymmetric[0][2] = 1.;
  Jacobian_Axisymmetric[0][3] = 0;

  Jacobian_Axisymmetric[1][0] = -U_i[1]*U_i[2]/(U_i[0]*U_i[0]);
  Jacobian_Axisymmetric[1][1] = U_i[2]/U_i[0];
  Jacobian_Axisymmetric[1][2] = U_i[1]/U_i[0];
  Jacobian_Axisymmetric[1][3] = 0;

  Jacobian_Axisymmetric[2][0] = -U_i[2]*U_i[2]/(U_i[0]*U_i[0]);
  Jacobian_Axisymmetric[2][1] = 0;
  Jacobian_Axisymmetric[2][2] = 2*U_i[2]/U_i[0];
  Jacobian_Axisymmetric[2][3] = 0;

  Jacobian_Axisymmetric[3][0] = -Gamma*U_i[2]*U_i[3]/(U_i[0]*U_i[0]) + (Gamma-1)*U_i[2]*(U_i[1]*U_i[1]+U_i[2]*U_i[2])/(U_i[0]*U_i[0]*U_i[0]);
  Jacobian_Axisymmetric[3][1] = -(Gamma-1)*U_i[2]*U_i[1]/(U_i[0]*U_i[0]);
  Jacobian_Axisymmetric[3][2] = Gamma*U_i[3]/U_i[0] - 1/2*(Gamma-1)*( (U_i[1]*U_i[1]+U_i[2]*U_i[2])/(U_i[0]*U_i[0]) + 2*U_i[2]*U_i[2]/(U_i[0]*U_i[0]) );
  Jacobian_Axisymmetric[3][3] = Gamma*U_i[2]/U_i[0];

  for (int iVar=0; iVar<4; iVar++)
    for (int jVar=0; jVar<4; jVar++)
      Jacobian_Axisymmetric[iVar][jVar] *= yinv*Volume;

  /* -- Residual = transpose(Jacobian) * psi --*/
  for (int iVar = 0; iVar < nVar; iVar++) {
    val_residual[iVar] = 0.0;
    for (int jVar = 0; jVar < nVar; jVar++) {
      val_residual[iVar] += Jacobian_Axisymmetric[jVar][iVar]*Psi_i[jVar];
      Jacobian_ii[iVar][jVar] = Jacobian_Axisymmetric[jVar][iVar];
    }
  }
}
