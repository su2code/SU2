/*!
 * \file CSourceAxisymmetric_Flow.cpp
 * \brief Implementation of numerics class CSourceAxisymmetric_Flow.
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

#include "../../../../include/numerics/flow/sources/CSourceAxisymmetric_Flow.hpp"

CSourceAxisymmetric_Flow::CSourceAxisymmetric_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;
  
}

CSourceAxisymmetric_Flow::~CSourceAxisymmetric_Flow(void) { }

void CSourceAxisymmetric_Flow::ComputeResidual(su2double *val_residual, su2double **Jacobian_i, CConfig *config) {
  
  su2double yinv, Pressure_i, Enthalpy_i, Velocity_i, sq_vel;
  unsigned short iDim, iVar, jVar;
  
  bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  
  if (Coord_i[1] > EPS) {
    
    yinv = 1.0/Coord_i[1];
    
    sq_vel = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      Velocity_i = U_i[iDim+1] / U_i[0];
      sq_vel += Velocity_i *Velocity_i;
    }
    
    Pressure_i = (Gamma-1.0)*U_i[0]*(U_i[nDim+1]/U_i[0]-0.5*sq_vel);
    Enthalpy_i = (U_i[nDim+1] + Pressure_i) / U_i[0];
    
    val_residual[0] = yinv*Volume*U_i[2];
    val_residual[1] = yinv*Volume*U_i[1]*U_i[2]/U_i[0];
    val_residual[2] = yinv*Volume*(U_i[2]*U_i[2]/U_i[0]);
    val_residual[3] = yinv*Volume*Enthalpy_i*U_i[2];
    
    if (implicit) {
      Jacobian_i[0][0] = 0.0;
      Jacobian_i[0][1] = 0.0;
      Jacobian_i[0][2] = 1.0;
      Jacobian_i[0][3] = 0.0;
      
      Jacobian_i[1][0] = -U_i[1]*U_i[2]/(U_i[0]*U_i[0]);
      Jacobian_i[1][1] = U_i[2]/U_i[0];
      Jacobian_i[1][2] = U_i[1]/U_i[0];
      Jacobian_i[1][3] = 0.0;
      
      Jacobian_i[2][0] = -U_i[2]*U_i[2]/(U_i[0]*U_i[0]);
      Jacobian_i[2][1] = 0.0;
      Jacobian_i[2][2] = 2*U_i[2]/U_i[0];
      Jacobian_i[2][3] = 0.0;
      
      Jacobian_i[3][0] = -Gamma*U_i[2]*U_i[3]/(U_i[0]*U_i[0]) + (Gamma-1)*U_i[2]*(U_i[1]*U_i[1]+U_i[2]*U_i[2])/(U_i[0]*U_i[0]*U_i[0]);
      Jacobian_i[3][1] = -(Gamma-1)*U_i[2]*U_i[1]/(U_i[0]*U_i[0]);
      Jacobian_i[3][2] = Gamma*U_i[3]/U_i[0] - 1/2*(Gamma-1)*( (U_i[1]*U_i[1]+U_i[2]*U_i[2])/(U_i[0]*U_i[0]) + 2*U_i[2]*U_i[2]/(U_i[0]*U_i[0]) );
      Jacobian_i[3][3] = Gamma*U_i[2]/U_i[0];
      
      for (iVar=0; iVar < nVar; iVar++)
        for (jVar=0; jVar < nVar; jVar++)
          Jacobian_i[iVar][jVar] *= yinv*Volume;
      
    }
    
  }
  
  else {
    
    for (iVar=0; iVar < nVar; iVar++)
      val_residual[iVar] = 0.0;
    
    if (implicit) {
      for (iVar=0; iVar < nVar; iVar++) {
        for (jVar=0; jVar < nVar; jVar++)
          Jacobian_i[iVar][jVar] = 0.0;
      }
    }
    
  }
  
}
