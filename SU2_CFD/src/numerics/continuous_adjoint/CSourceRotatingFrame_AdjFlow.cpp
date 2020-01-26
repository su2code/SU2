/*!
 * \file CSourceRotatingFrame_AdjFlow.cpp
 * \brief Implementation of numerics class CSourceRotatingFrame_AdjFlow.
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

#include "../../../include/numerics/continuous_adjoint/CSourceRotatingFrame_AdjFlow.hpp"

CSourceRotatingFrame_AdjFlow::CSourceRotatingFrame_AdjFlow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) { }

CSourceRotatingFrame_AdjFlow::~CSourceRotatingFrame_AdjFlow(void) { }

void CSourceRotatingFrame_AdjFlow::ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, CConfig *config) {
  
  unsigned short iDim, iVar, jVar;
  su2double Omega[3] = {0,0,0}, Phi[3] = {0,0,0};
  bool implicit = (config->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);

  /*--- Retrieve the angular velocity vector from config. ---*/

  for (iDim = 0; iDim < 3; iDim++){
    Omega[iDim] = config->GetRotation_Rate(iDim)/config->GetOmega_Ref();
  }
  
  /*--- Get the adjoint velocity vector at the current node. ---*/
  
  for (iDim = 0; iDim < nDim; iDim++)
    Phi[iDim] = Psi_i[iDim+1];
  
  /*--- Compute the source term as the Jacobian of the rotating frame
   source term multiplied by the adjoint state and the dual cell volume. ---*/
  
  if (nDim == 2) {
    val_residual[0] = 0.0;
    val_residual[1] =  Omega[2]*Phi[1]*Volume;
    val_residual[2] = -Omega[2]*Phi[0]*Volume;
    val_residual[3] = 0.0;
  } else {
    val_residual[0] = 0.0;
    val_residual[1] = (Omega[2]*Phi[1] - Omega[1]*Phi[2])*Volume;
    val_residual[2] = (Omega[0]*Phi[2] - Omega[2]*Phi[0])*Volume;
    val_residual[3] = (Omega[1]*Phi[0] - Omega[0]*Phi[1])*Volume;
    val_residual[4] = 0.0;
  }
  
  /*--- Calculate the source term Jacobian ---*/
  
  if (implicit) {
    for (iVar = 0; iVar < nVar; iVar++)
      for (jVar = 0; jVar < nVar; jVar++)
        val_Jacobian_i[iVar][jVar] = 0.0;
    if (nDim == 2) {
      val_Jacobian_i[1][2] =  Omega[2]*Volume;
      val_Jacobian_i[2][1] = -Omega[2]*Volume;
    } else {
      val_Jacobian_i[1][2] =  Omega[2]*Volume;
      val_Jacobian_i[1][3] = -Omega[1]*Volume;
      val_Jacobian_i[2][1] = -Omega[2]*Volume;
      val_Jacobian_i[2][3] =  Omega[0]*Volume;
      val_Jacobian_i[3][1] =  Omega[1]*Volume;
      val_Jacobian_i[3][2] = -Omega[0]*Volume;
    }
  }
  
}
