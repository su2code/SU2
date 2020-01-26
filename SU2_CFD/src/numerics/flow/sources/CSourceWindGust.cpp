/*!
 * \file CSourceWindGust.cpp
 * \brief Implementation of numerics class CSourceWindGust.
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

#include "../../../../include/numerics/flow/sources/CSourceWindGust.hpp"

CSourceWindGust::CSourceWindGust(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
}

CSourceWindGust::~CSourceWindGust(void) { }

void CSourceWindGust::ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, CConfig *config) {
  
  su2double u_gust, v_gust, du_gust_dx, du_gust_dy, du_gust_dt, dv_gust_dx, dv_gust_dy, dv_gust_dt, smx, smy, se, rho, u, v, p;
  unsigned short GustDir = config->GetGust_Dir(); //Gust direction
  
  u_gust = WindGust_i[0];
  v_gust = WindGust_i[1];
  
  if (GustDir == X_DIR) {
    du_gust_dx = WindGustDer_i[0];
    du_gust_dy = WindGustDer_i[1];
    du_gust_dt = WindGustDer_i[2];
    dv_gust_dx = 0.0;
    dv_gust_dy = 0.0;
    dv_gust_dt = 0.0;
  } else {
    du_gust_dx = 0.0;
    du_gust_dy = 0.0;
    du_gust_dt = 0.0;
    dv_gust_dx = WindGustDer_i[0];
    dv_gust_dy = WindGustDer_i[1];
    dv_gust_dt = WindGustDer_i[2];
    
  }
  
  /*--- Primitive variables at point i ---*/
  u = V_i[1];
  v = V_i[2];
  p = V_i[nDim+1];
  rho = V_i[nDim+2];
  
  /*--- Source terms ---*/
  smx = rho*(du_gust_dt + (u+u_gust)*du_gust_dx + (v+v_gust)*du_gust_dy);
  smy = rho*(dv_gust_dt + (u+u_gust)*dv_gust_dx + (v+v_gust)*dv_gust_dy);
  se = u*smx + v*smy + p*(du_gust_dx + dv_gust_dy);
  
  if (nDim == 2) {
    val_residual[0] = 0.0;
    val_residual[1] = smx*Volume;
    val_residual[2] = smy*Volume;
    val_residual[3] = se*Volume;
  } else {
    SU2_MPI::Error("You should only be in the gust source term in two dimensions", CURRENT_FUNCTION);
  }
  
  /*--- For now the source term Jacobian is just set to zero ---*/
  
  unsigned short iVar, jVar;
  bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  
  /*--- Calculate the source term Jacobian ---*/
  
  if (implicit) {
    for (iVar = 0; iVar < nVar; iVar++)
      for (jVar = 0; jVar < nVar; jVar++)
        val_Jacobian_i[iVar][jVar] = 0.0;
  }
}
