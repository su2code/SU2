/*!
 * \file CSourceVolumetricHeat.cpp
 * \brief Numerical methods for volumetric heat source term integration.
 * \author Ruben Sanchez
 * \version 6.2.0 "Falcon"
 *
 * The current SU2 release has been coordinated by the
 * SU2 International Developers Society <www.su2devsociety.org>
 * with selected contributions from the open-source community.
 *
 * The main research teams contributing to the current release are:
 *  - Prof. Juan J. Alonso's group at Stanford University.
 *  - Prof. Piero Colonna's group at Delft University of Technology.
 *  - Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *  - Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *  - Prof. Rafael Palacios' group at Imperial College London.
 *  - Prof. Vincent Terrapon's group at the University of Liege.
 *  - Prof. Edwin van der Weide's group at the University of Twente.
 *  - Lab. of New Concepts in Aeronautics at Tech. Institute of Aeronautics.
 *
 * Copyright 2012-2018, Francisco D. Palacios, Thomas D. Economon,
 *                      Tim Albring, and the SU2 contributors.
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

#include "../../../include/numerics/source/CSourceVolumetricHeat.hpp"

CSourceVolumetricHeat::CSourceVolumetricHeat(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {

  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);

}

CSourceVolumetricHeat::~CSourceVolumetricHeat(void) {

}

void CSourceVolumetricHeat::ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, CConfig *config) {

  unsigned short iDim;

  /*--- Zero the continuity contribution ---*/

  val_residual[0] = 0.0;

  /*--- Zero the momentum contribution. ---*/

  for (iDim = 0; iDim < nDim; iDim++)
    val_residual[iDim+1] = 0.0;

  /*--- Set the energy contribution ---*/

  val_residual[nDim+1] = -1.0*config->GetHeatSource_Val()*Volume;

  /*--- Jacobian contribution is 0 as the heat source is constant ---*/

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
    val_Jacobian_i[3][3] = 0.0;

  }

}
