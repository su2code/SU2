/*!
 * \file CSourceIncRotatingFrame_Flow.cpp
 * \brief Implementation of numerics class CSourceIncRotatingFrame_Flow.
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

#include "../../../../include/numerics/flow/sources/CSourceIncRotatingFrame_Flow.hpp"

CSourceIncRotatingFrame_Flow::CSourceIncRotatingFrame_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {

  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);

  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;

  /*--- Retrieve the angular velocity vector from config. ---*/
  for (unsigned short iDim = 0; iDim < 3; iDim++)
    Omega[iDim] = config->GetRotation_Rate(iDim)/config->GetOmega_Ref();

}

CSourceIncRotatingFrame_Flow::~CSourceIncRotatingFrame_Flow(void) { }

void CSourceIncRotatingFrame_Flow::ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, CConfig *config) {

  unsigned short iDim, iVar, jVar;
  su2double Momentum[3] = {0,0,0},
            Velocity_i[3] = {0,0,0};

  /*--- Primitive variables plus momentum at the node (point i) ---*/

  DensityInc_i  = V_i[nDim+2];

  for (iDim = 0; iDim < nDim; iDim++) {
    Velocity_i[iDim] = V_i[iDim+1];
    Momentum[iDim]   = DensityInc_i*Velocity_i[iDim];
  }

  /*--- Calculate rotating frame source term residual as ( Omega X Rho-U ) ---*/

  if (nDim == 2) {
    val_residual[0] = 0.0;
    val_residual[1] = (Omega[1]*Momentum[2] - Omega[2]*Momentum[1])*Volume;
    val_residual[2] = (Omega[2]*Momentum[0] - Omega[0]*Momentum[2])*Volume;
    val_residual[3] = 0.0;
  } else {
    val_residual[0] = 0.0;
    val_residual[1] = (Omega[1]*Momentum[2] - Omega[2]*Momentum[1])*Volume;
    val_residual[2] = (Omega[2]*Momentum[0] - Omega[0]*Momentum[2])*Volume;
    val_residual[3] = (Omega[0]*Momentum[1] - Omega[1]*Momentum[0])*Volume;
    val_residual[4] = 0.0;
  }

  /*--- Calculate the source term Jacobian ---*/

  if (implicit) {
    for (iVar = 0; iVar < nVar; iVar++)
      for (jVar = 0; jVar < nVar; jVar++)
        val_Jacobian_i[iVar][jVar] = 0.0;
    if (nDim == 2) {
      val_Jacobian_i[1][2] = -DensityInc_i*Omega[2]*Volume;
      val_Jacobian_i[2][1] =  DensityInc_i*Omega[2]*Volume;
    } else {
      val_Jacobian_i[1][2] = -DensityInc_i*Omega[2]*Volume;
      val_Jacobian_i[1][3] =  DensityInc_i*Omega[1]*Volume;
      val_Jacobian_i[2][1] =  DensityInc_i*Omega[2]*Volume;
      val_Jacobian_i[2][3] = -DensityInc_i*Omega[0]*Volume;
      val_Jacobian_i[3][1] = -DensityInc_i*Omega[1]*Volume;
      val_Jacobian_i[3][2] =  DensityInc_i*Omega[0]*Volume;
    }
  }

}
