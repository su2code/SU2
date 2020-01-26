/*!
 * \file CSourceIncAxisymmetric_Flow.cpp
 * \brief Implementation of numerics class CSourceIncAxisymmetric_Flow.
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

#include "../../../../include/numerics/flow/sources/CSourceIncAxisymmetric_Flow.hpp"

CSourceIncAxisymmetric_Flow::CSourceIncAxisymmetric_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {

  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  energy   = config->GetEnergy_Equation();
  viscous  = config->GetViscous();

}

CSourceIncAxisymmetric_Flow::~CSourceIncAxisymmetric_Flow(void) { }

void CSourceIncAxisymmetric_Flow::ComputeResidual(su2double *val_residual, su2double **Jacobian_i, CConfig *config) {

  su2double yinv, Velocity_i[3];
  unsigned short iDim, jDim, iVar, jVar;

  if (Coord_i[1] > EPS) {

    yinv = 1.0/Coord_i[1];

    /*--- Set primitive variables at points iPoint. ---*/

    Pressure_i    = V_i[0];
    Temp_i        = V_i[nDim+1];
    DensityInc_i  = V_i[nDim+2];
    BetaInc2_i    = V_i[nDim+3];
    Cp_i          = V_i[nDim+7];
    Enthalpy_i    = Cp_i*Temp_i;

    for (iDim = 0; iDim < nDim; iDim++)
      Velocity_i[iDim] = V_i[iDim+1];

    /*--- Inviscid component of the source term. ---*/

    val_residual[0] = yinv*Volume*DensityInc_i*Velocity_i[1];
    val_residual[1] = yinv*Volume*DensityInc_i*Velocity_i[0]*Velocity_i[1];
    val_residual[2] = yinv*Volume*DensityInc_i*Velocity_i[1]*Velocity_i[1];
    val_residual[3] = yinv*Volume*DensityInc_i*Enthalpy_i*Velocity_i[1];

    if (implicit) {

      Jacobian_i[0][0] = 0.0;
      Jacobian_i[0][1] = 0.0;
      Jacobian_i[0][2] = 1.0;
      Jacobian_i[0][3] = 0.0;

      Jacobian_i[1][0] = 0.0;
      Jacobian_i[1][1] = Velocity_i[1];
      Jacobian_i[1][2] = Velocity_i[0];
      Jacobian_i[1][3] = 0.0;

      Jacobian_i[2][0] = 0.0;
      Jacobian_i[2][1] = 0.0;
      Jacobian_i[2][2] = 2.0*Velocity_i[1];
      Jacobian_i[2][3] = 0.0;

      Jacobian_i[3][0] = 0.0;
      Jacobian_i[3][1] = 0.0;
      Jacobian_i[3][2] = Enthalpy_i;
      Jacobian_i[3][3] = Cp_i*Velocity_i[1];

      for (iVar=0; iVar < nVar; iVar++)
        for (jVar=0; jVar < nVar; jVar++)
          Jacobian_i[iVar][jVar] *= yinv*Volume*DensityInc_i;
      
    }

    /*--- Add the viscous terms if necessary. ---*/

    if (viscous) {

      Laminar_Viscosity_i    = V_i[nDim+4];
      Eddy_Viscosity_i       = V_i[nDim+5];
      Thermal_Conductivity_i = V_i[nDim+6];

      su2double total_viscosity, div_vel;

      total_viscosity = (Laminar_Viscosity_i + Eddy_Viscosity_i);

      /*--- The full stress tensor is needed for variable density ---*/

      div_vel = 0.0;
      for (iDim = 0 ; iDim < nDim; iDim++)
        div_vel += PrimVar_Grad_i[iDim+1][iDim];

      for (iDim = 0 ; iDim < nDim; iDim++)
        for (jDim = 0 ; jDim < nDim; jDim++)
          tau[iDim][jDim] = (total_viscosity*(PrimVar_Grad_i[jDim+1][iDim] +
                                              PrimVar_Grad_i[iDim+1][jDim] )
                             -TWO3*total_viscosity*div_vel*delta[iDim][jDim]);
      
      /*--- Viscous terms. ---*/

      val_residual[0] -= 0.0;
      val_residual[1] -= Volume*(yinv*tau[0][1] - TWO3*AuxVar_Grad_i[0]);
      val_residual[2] -= Volume*(yinv*2.0*total_viscosity*PrimVar_Grad_i[2][1] -
                                 yinv*yinv*2.0*total_viscosity*Velocity_i[1] -
                                 TWO3*AuxVar_Grad_i[1]);
      val_residual[3] -= Volume*yinv*Thermal_Conductivity_i*PrimVar_Grad_i[nDim+1][1];

    }
    
  } else {

    for (iVar=0; iVar < nVar; iVar++)
      val_residual[iVar] = 0.0;

    if (implicit) {
      for (iVar=0; iVar < nVar; iVar++) {
        for (jVar=0; jVar < nVar; jVar++)
          Jacobian_i[iVar][jVar] = 0.0;
      }
    }
    
  }

  if (!energy) {
    val_residual[nDim+1] = 0.0;
    if (implicit) {
      for (iVar = 0; iVar < nVar; iVar++) {
        Jacobian_i[iVar][nDim+1] = 0.0;
        Jacobian_i[nDim+1][iVar] = 0.0;
      }
    }
  }
}
