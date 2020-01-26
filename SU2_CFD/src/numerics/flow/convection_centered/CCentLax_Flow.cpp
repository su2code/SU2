/*!
 * \file CCentLax_Flow.cpp
 * \brief Implementation of numerics class CCentLax_Flow.
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

#include "../../../../include/numerics/flow/convection_centered/CCentLax_Flow.hpp"

CCentLax_Flow::CCentLax_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) :
               CCentBase_Flow(val_nDim, val_nVar, config) {

  /*--- Artifical dissipation parameters ---*/
  Param_p = 0.3;
  Param_Kappa_0 = config->GetKappa_1st_Flow();

}

CCentLax_Flow::~CCentLax_Flow(void) {

}

void CCentLax_Flow::DissipationTerm(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j) {

  /*--- Compute dissipation coefficient ---*/

  sc0 = 3.0*(su2double(Neighbor_i)+su2double(Neighbor_j))/(su2double(Neighbor_i)*su2double(Neighbor_j));
  Epsilon_0 = Param_Kappa_0*sc0*su2double(nDim)/3.0;
  
  /*--- Compute viscous part of the residual ---*/
  
  for (iVar = 0; iVar < nVar; iVar++)
    val_residual[iVar] += Epsilon_0*Diff_U[iVar]*StretchingFactor*MeanLambda;
  
  /*--- Jacobian computation ---*/

  if (implicit) {
    
    cte_0 = Epsilon_0*StretchingFactor*MeanLambda;
    cte_1 = cte_0;
    
    ScalarDissipationJacobian(val_Jacobian_i, val_Jacobian_j);
  }
}

bool CCentLax_Flow::SetPreaccInVars(void) {
  AD::StartPreacc();
  return true;
}
