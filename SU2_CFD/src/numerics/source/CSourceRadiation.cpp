/*!
 * \file CSourceRadiation.cpp
 * \brief This file contains the source term integration from Radiative Heat Transfer.
 * \author Ruben Sanchez
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

#include "../../../include/numerics/source/CSourceRadiation.hpp"

CSourceRadiation::CSourceRadiation(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {

  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);

}

CSourceRadiation::~CSourceRadiation(void) {

}

void CSourceRadiation::ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, CConfig *config) {

  unsigned short iDim;

  /*--- Zero the continuity contribution ---*/

  val_residual[0] = 0.0;

  /*--- Zero the momentum contribution. ---*/

  for (iDim = 0; iDim < nDim; iDim++)
    val_residual[iDim+1] = 0.0;

  /*--- Set the energy contribution ---*/

  val_residual[nDim+1] = - RadVar_Source[0]*Volume;

  /*--- Set the energy contribution to the Jacobian ---*/

  if (implicit) {

    val_Jacobian_i[0][0] = 0.0;
    val_Jacobian_i[0][1] = 0.0;
    val_Jacobian_i[0][2] = 0.0;
    val_Jacobian_i[0][3] = 0.0;

    val_Jacobian_i[1][0] = 0.0;
    val_Jacobian_i[1][1] = 0.0;
    val_Jacobian_i[1][2] = 0.0;
    val_Jacobian_i[1][3] = 0.0;

    val_Jacobian_i[2][0] = 0.0;
    val_Jacobian_i[2][1] = 0.0;
    val_Jacobian_i[2][2] = 0.0;
    val_Jacobian_i[2][3] = 0.0;

    val_Jacobian_i[3][0] = 0.0;
    val_Jacobian_i[3][1] = 0.0;
    val_Jacobian_i[3][2] = 0.0;
    val_Jacobian_i[3][3] = - RadVar_Source[1]*Volume;

  }

}
