/*!
 * \file CHeatVariable.cpp
 * \brief Definition of the variables for heat equation problems.
 * \author F. Palacios, T. Economon
 * \version 8.0.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2023, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../include/variables/CHeatVariable.hpp"

CHeatVariable::CHeatVariable(su2double heat, unsigned long npoint, unsigned long ndim, unsigned long nvar, CConfig *config)
  : CScalarVariable(npoint, ndim, nvar, config) {

  const bool dual_time = ((config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_1ST) ||
                          (config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_2ND));

  /*--- Initialization of heat variable ---*/

  Solution = heat;
  Solution_Old = heat;

  /*--- Allocate residual structures ---*/

  Res_TruncError.resize(nPoint, nVar) = su2double(0.0);

  /*--- Only for residual smoothing (multigrid) ---*/

  for (unsigned long iMesh = 0; iMesh <= config->GetnMGLevels(); iMesh++) {
    if (config->GetMG_CorrecSmooth(iMesh) > 0) {
      Residual_Sum.resize(nPoint, nVar);
      Residual_Old.resize(nPoint, nVar);
      break;
    }
  }

  /*--- Allocate and initialize solution for dual time strategy ---*/
  if (dual_time) {
    Solution_time_n  = heat;
    Solution_time_n1 = heat;
  }

  Max_Lambda_Visc.resize(nPoint);

  if (config->GetMultizone_Problem())
    Set_BGSSolution_k();
}
