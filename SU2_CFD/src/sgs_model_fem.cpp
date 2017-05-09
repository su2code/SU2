/*!
 * \file sgs_model_fem.cpp
 * \brief Main subroutines for calculating sub-grid scale models for LES with DG FEM framework.
 * \author J. Alonso, E. van der Weide, T. Economon, P. Urbanczyk
 * \version 5.0.0 "Raven"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
 *                      Dr. Thomas D. Economon (economon@stanford.edu).
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *
 * Copyright (C) 2012-2015 SU2, the open-source CFD code.
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

#include "../include/sgs_model.hpp"

CSGSModel::CSGSModel(void) {

  /*--- Initialization of member variables ---*/

}

CSGSModel::~CSGSModel(void) {

  /*--- Array deallocation ---*/

}

CSmagorinskyModel::CSmagorinskyModel(void) :
        CSGSModel() {

  const_smag = 0.1;
  filter_mult = 2.0;

}

CSmagorinskyModel::~CSmagorinskyModel(void) {

}

su2double CSmagorinskyModel::ComputeEddyViscosity(const unsigned short nDim,
    const su2double rho, const su2double velGrad[3][3],
    const su2double lenScale, const su2double distToWall) {

  /* Constant coefficient Smagorinsky SGS is calculated:
   * ( C_s * L_c )^2 * |S(x,t)|
   * C_s = Smagorinsky constant
   * L_c = Filter width
   * S(x,t) = Rate of Strain Tensor ( 1/2 [ du_i/dx_j + du_j/dx_i] )
   */

  const su2double filter_width = filter_mult * lenScale;

  /* Calculate the magnitude of the rate of strain tensor */

  su2double strain_rate_squared = 0.0;

  for (unsigned short int i = 0; i < nDim; ++i) {
    strain_rate_squared += 2 * velGrad[i][i] * velGrad[i][i];
    for (unsigned short int j = 0; j < nDim; ++j) {
      if (i != j) {
        strain_rate_squared += (velGrad[i][j] + velGrad[j][i])
                * (velGrad[i][j] + velGrad[j][i]);
      }
    }
  }

  su2double strain_rate = sqrt(strain_rate_squared);

  su2double sgs_viscosity = const_smag * const_smag * filter_width
      * filter_width * strain_rate;

  return (sgs_viscosity);

}

void CSmagorinskyModel::ComputeGradEddyViscosity(const unsigned short nDim,
    const su2double rho, const su2double rhoGrad[3],
    const su2double velGrad[3][3], const su2double velHess[3][3][3],
    const su2double lenScale, const su2double distToWall,
    su2double ViscosityTurbGrad[3]) {
  /* Not implemented yet. */
}

CWALEModel::CWALEModel(void) :
        CSGSModel() {

}

CWALEModel::~CWALEModel(void) {

}

su2double CWALEModel::ComputeEddyViscosity(const unsigned short nDim,
    const su2double rho, const su2double velGrad[3][3],
    const su2double lenScale, const su2double distToWall) {
  /* Not implemented yet. */
}

void CWALEModel::ComputeGradEddyViscosity(const unsigned short nDim,
    const su2double rho, const su2double rhoGrad[3],
    const su2double velGrad[3][3], const su2double velHess[3][3][3],
    const su2double lenScale, const su2double distToWall,
    su2double ViscosityTurbGrad[3]) {
  /* Not implemented yet. */
}
