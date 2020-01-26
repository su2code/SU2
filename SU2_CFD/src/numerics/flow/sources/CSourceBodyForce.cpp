/*!
 * \file CSourceBodyForce.cpp
 * \brief Implementation of numerics class CSourceBodyForce.
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

#include "../../../../include/numerics/flow/sources/CSourceBodyForce.hpp"

CSourceBodyForce::CSourceBodyForce(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {

  /*--- Store the pointer to the constant body force vector. ---*/

  Body_Force_Vector = new su2double[nDim];
  for (unsigned short iDim = 0; iDim < nDim; iDim++)
    Body_Force_Vector[iDim] = config->GetBody_Force_Vector()[iDim];

}

CSourceBodyForce::~CSourceBodyForce(void) {

  if (Body_Force_Vector != NULL) delete [] Body_Force_Vector;

}

void CSourceBodyForce::ComputeResidual(su2double *val_residual, CConfig *config) {
  
  unsigned short iDim;
  su2double Force_Ref = config->GetForce_Ref();
  
  /*--- Zero the continuity contribution ---*/
  
  val_residual[0] = 0.0;
  
  /*--- Momentum contribution ---*/
  
  for (iDim = 0; iDim < nDim; iDim++)
    val_residual[iDim+1] = -Volume * U_i[0] * Body_Force_Vector[iDim] / Force_Ref;
  
  /*--- Energy contribution ---*/
  
  val_residual[nDim+1] = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    val_residual[nDim+1] += -Volume * U_i[iDim+1] * Body_Force_Vector[iDim] / Force_Ref;
  
}
