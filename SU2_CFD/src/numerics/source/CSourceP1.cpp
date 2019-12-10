/*!
 * \file CSourceP1.cpp
 * \brief Numerical methods for source term integration of the P1 radiation model.
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

#include "../../../include/numerics/CNumericsRadiation.hpp"
#include "../../../include/numerics/source/CSourceP1.hpp"

CSourceP1::CSourceP1(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumericsRadiation(val_nDim, val_nVar, config) {

}

CSourceP1::~CSourceP1(void) {

}

void CSourceP1::ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, CConfig *config) {

  unsigned short iDim;

  /*--- Retrieve the energy at the node i ---*/
  Energy_i = RadVar_i[0];

  /*--- Retrieve the temperature at the node i ---*/
  Temperature_i = V_i[nDim+1];

  /*--- Compute the blackbody intensity for gray media ---*/
  BlackBody_Intensity = 4.0*STEFAN_BOLTZMANN*pow(Temperature_i,4.0);

  /*--- Source term from black-body and energy contributions ---*/
  val_residual[0] = Absorption_Coeff * Volume * (BlackBody_Intensity - Energy_i);

  /*--- Contribution to the Jacobian ---*/
  if (implicit) {
    val_Jacobian_i[0][0] = - Absorption_Coeff * Volume;
  }

}
