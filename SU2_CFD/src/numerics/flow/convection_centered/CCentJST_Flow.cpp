/*!
 * \file CCentJST_Flow.cpp
 * \brief Implementation of numerics class CCentJST_Flow.
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

#include "../../../../include/numerics/flow/convection_centered/CCentJST_Flow.hpp"

CCentJST_Flow::CCentJST_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) :
               CCentBase_Flow(val_nDim, val_nVar, config) {

  /*--- Artifical dissipation parameters ---*/
  Param_p = 0.3;
  Param_Kappa_2 = config->GetKappa_2nd_Flow();
  Param_Kappa_4 = config->GetKappa_4th_Flow();

}

CCentJST_Flow::~CCentJST_Flow(void) {

}

void CCentJST_Flow::DissipationTerm(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j) {

  /*--- Compute differences btw. Laplacians ---*/

  for (iVar = 0; iVar < nVar; iVar++) {
    Diff_Lapl[iVar] = Und_Lapl_i[iVar]-Und_Lapl_j[iVar];
  }

  /*--- Compute dissipation coefficients ---*/

  sc2 = 3.0*(su2double(Neighbor_i)+su2double(Neighbor_j))/(su2double(Neighbor_i)*su2double(Neighbor_j));
  sc4 = sc2*sc2/4.0;

  Epsilon_2 = Param_Kappa_2*0.5*(Sensor_i+Sensor_j)*sc2;
  Epsilon_4 = max(0.0, Param_Kappa_4-Epsilon_2)*sc4;

  /*--- Compute viscous part of the residual ---*/

  for (iVar = 0; iVar < nVar; iVar++)
    val_residual[iVar] += (Epsilon_2*Diff_U[iVar] - Epsilon_4*Diff_Lapl[iVar])*StretchingFactor*MeanLambda;

  /*--- Jacobian computation ---*/

  if (implicit) {

    cte_0 = (Epsilon_2 + Epsilon_4*su2double(Neighbor_i+1))*StretchingFactor*MeanLambda;
    cte_1 = (Epsilon_2 + Epsilon_4*su2double(Neighbor_j+1))*StretchingFactor*MeanLambda;

    ScalarDissipationJacobian(val_Jacobian_i, val_Jacobian_j);
  }
}

bool CCentJST_Flow::SetPreaccInVars(void) {
  AD::StartPreacc();
  AD::SetPreaccIn(Sensor_i);  AD::SetPreaccIn(Und_Lapl_i, nVar);
  AD::SetPreaccIn(Sensor_j);  AD::SetPreaccIn(Und_Lapl_j, nVar);
  return true;
}
