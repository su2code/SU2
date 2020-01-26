/*!
 * \file CSourceGravity.cpp
 * \brief Implementation of numerics class CSourceGravity.
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

#include "../../../../include/numerics/flow/sources/CSourceGravity.hpp"

CSourceGravity::CSourceGravity(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {

}

CSourceGravity::~CSourceGravity(void) { }

void CSourceGravity::ComputeResidual(su2double *val_residual, CConfig *config) {
  unsigned short iVar;

  for (iVar = 0; iVar < nVar; iVar++)
    val_residual[iVar] = 0.0;

  /*--- Evaluate the source term  ---*/
  val_residual[nDim] = Volume * U_i[0] * STANDARD_GRAVITY;

}
