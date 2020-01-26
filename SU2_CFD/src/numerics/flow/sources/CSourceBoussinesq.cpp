/*!
 * \file CSourceBoussinesq.cpp
 * \brief Implementation of numerics class CSourceBoussinesq.
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

#include "../../../../include/numerics/flow/sources/CSourceBoussinesq.hpp"

CSourceBoussinesq::CSourceBoussinesq(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {

  /*--- Store the pointer to the constant body force vector. ---*/

  Gravity_Vector = new su2double[nDim];
  for (unsigned short iDim = 0; iDim < nDim; iDim++)
    Gravity_Vector[iDim] = 0.0;

  /*--- Gravity is downward in y-dir for 2D and downward z-dir for 3D. ---*/

  Gravity_Vector[nDim-1] = -STANDARD_GRAVITY;

}

CSourceBoussinesq::~CSourceBoussinesq(void) {

  if (Gravity_Vector != NULL) delete [] Gravity_Vector;

}

void CSourceBoussinesq::ComputeResidual(su2double *val_residual, CConfig *config) {

  unsigned short iDim;
  su2double Force_Ref = config->GetForce_Ref();
  su2double T0        = config->GetTemperature_FreeStreamND();
  su2double Beta      = config->GetThermal_Expansion_CoeffND();

  /*--- Zero the continuity contribution ---*/

  val_residual[0] = 0.0;

  /*--- Momentum contribution. Note that this form assumes we have
   subtracted the operating density * gravity, i.e., removed the
   hydrostatic pressure component (important for pressure BCs). ---*/

  for (iDim = 0; iDim < nDim; iDim++)
    val_residual[iDim+1] = Volume * DensityInc_i * ( Beta * (U_i[nDim+1] - T0)) * Gravity_Vector[iDim] / Force_Ref;

  /*--- Zero the energy contribution ---*/

  val_residual[nDim+1] = 0.0;

}
