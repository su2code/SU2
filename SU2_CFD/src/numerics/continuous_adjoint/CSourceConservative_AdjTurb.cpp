/*!
 * \file CSourceConservative_AdjTurb.cpp
 * \brief Implementation of numerics class CSourceConservative_AdjTurb.
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

#include "../../../include/numerics/continuous_adjoint/CSourceConservative_AdjTurb.hpp"

CSourceConservative_AdjTurb::CSourceConservative_AdjTurb(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {

  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;

}

CSourceConservative_AdjTurb::~CSourceConservative_AdjTurb(void) {
}

void CSourceConservative_AdjTurb::ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config) {

  /*--- SOURCE term  -->  \nabla ( \psi_\mu \B7 E^{s} )
   E^{s} = 2 c_{b2}/\sigma \nabla \hat{nu} ---*/

  unsigned short iDim;
  bool implicit = (config->GetKind_TimeIntScheme_AdjTurb() == EULER_IMPLICIT);

  su2double cb2 = 0.622;
  su2double sigma = 2./3.;
  su2double coeff = 2.0*cb2/sigma;
  su2double E_ij, proj_TurbVar_Grad_i, proj_TurbVar_Grad_j;

  E_ij = 0.0;  proj_TurbVar_Grad_i = 0.0; proj_TurbVar_Grad_j = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    proj_TurbVar_Grad_i += coeff*TurbVar_Grad_i[0][iDim]*Normal[iDim];
    proj_TurbVar_Grad_j += coeff*TurbVar_Grad_j[0][iDim]*Normal[iDim];
    E_ij += 0.5*(TurbPsi_i[0]*proj_TurbVar_Grad_i + TurbPsi_j[0]*proj_TurbVar_Grad_j);
  }

  val_residual[0] = E_ij;

  if (implicit) {
    val_Jacobian_i[0][0] = 0.5*proj_TurbVar_Grad_i;
    val_Jacobian_j[0][0] = 0.5*proj_TurbVar_Grad_j;
  }
}
