/*!
 * \file CUpwLin_TransLM.cpp
 * \brief Implementation of numerics class CUpwLin_TransLM.
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

#include "../../../include/numerics/transition/CUpwLin_TransLM.hpp"

CUpwLin_TransLM::CUpwLin_TransLM(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  implicit = (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT);
  incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  
  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;
  
  Velocity_i = new su2double [nDim];
  Velocity_j = new su2double [nDim];
  
}

CUpwLin_TransLM::~CUpwLin_TransLM(void) {
  delete [] Velocity_i;
  delete [] Velocity_j;
}

void CUpwLin_TransLM::ComputeResidual (su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config) {
  
  
  Density_i = U_i[0];
  Density_j = U_j[0];
  
  q_ij = 0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Velocity_i[iDim] = U_i[iDim+1]/Density_i;
    Velocity_j[iDim] = U_j[iDim+1]/Density_j;
    q_ij += 0.5*(Velocity_i[iDim]+Velocity_j[iDim])*Normal[iDim];
  }
  
  a0 = 0.5*(q_ij+fabs(q_ij));
  a1 = 0.5*(q_ij-fabs(q_ij));
  val_residual[0] = a0*TransVar_i[0]+a1*TransVar_j[0];
  val_residual[1] = a0*TransVar_i[1]+a1*TransVar_j[1];
  //  cout << "Velicity x: " << Velocity_i[0] << ", " << Velocity_j[0] << endl;
  //  cout << "Velicity y: " << Velocity_i[1] << ", " << Velocity_j[1] << endl;
  //  cout << "val_resid: " << val_residual[0] << ", " << val_residual[1] << endl;
  
  
  if (implicit) {
    val_Jacobian_i[0][0] = a0;
    val_Jacobian_i[1][1] = a0;
  }
}
