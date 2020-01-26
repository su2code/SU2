/*!
 * \file CUpwLin_AdjTurb.cpp
 * \brief Implementation of numerics class CUpwLin_AdjTurb.
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

#include "../../../include/numerics/continuous_adjoint/CUpwLin_AdjTurb.hpp"

CUpwLin_AdjTurb::CUpwLin_AdjTurb(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {

  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;

  Velocity_i = new su2double [nDim];
}

CUpwLin_AdjTurb::~CUpwLin_AdjTurb(void) {
  delete [] Velocity_i;
}

void CUpwLin_AdjTurb::ComputeResidual (su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config) {
  bool implicit = (config->GetKind_TimeIntScheme_AdjTurb() == EULER_IMPLICIT);

  /*--- Non-conservative term  -->  -\nabla \psi_\mu  B^{cv}
   B^{cv} = -v ---*/

  unsigned short iDim;
  su2double proj_conv_flux = 0;

  for (iDim = 0; iDim < nDim; iDim++) {
    Velocity_i[iDim] = U_i[iDim+1]/U_i[0];
    proj_conv_flux += Velocity_i[iDim]*Normal[iDim]; // projection of convective flux at iPoint
  }
  su2double psinu0 = TurbPsi_i[0];
  su2double psinu1;
  if (proj_conv_flux > 0)
    psinu1 = psinu0 + proj_conv_flux;
  else
    psinu1 = psinu0;

  val_residual[0] = 0.5*( proj_conv_flux*(psinu0+psinu1)-fabs(proj_conv_flux)*(psinu1-psinu0));
  if (implicit) {
    val_Jacobian_i[0][0] = 0.5*( proj_conv_flux + fabs(proj_conv_flux));
  }
}
