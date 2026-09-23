/*!
 * \file CTurbSAVariable.cpp
 * \brief Definition of the solution fields.
 * \author F. Palacios, A. Bueno
 * \version 8.5.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md)
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


#include "../../include/variables/CTurbSAVariable.hpp"

CTurbSAVariable::CTurbSAVariable(su2double val_nu_tilde, su2double val_muT, unsigned long npoint,
                                 unsigned long ndim, unsigned long nvar, CConfig *config) :
                 CTurbVariable(npoint, ndim, nvar, config) {

  /*--- Initialize solution (check if the Stochastic Backscatter Model is active) ---*/
  bool backscatter = config->GetSBSParam().StochasticBackscatter;
  if (!(backscatter && config->GetSBSParam().SBS_Ctau > 0.0)) {
    Solution_Old = Solution = val_nu_tilde;
  } else {
    for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++) {
      Solution_Old(iPoint, 0) = Solution(iPoint, 0) = val_nu_tilde;
      for (unsigned long iVar = 1; iVar < nVar; iVar++) {
        Solution_Old(iPoint, iVar) = Solution(iPoint, iVar) = 0.0;
      }
    }
  }

  muT.resize(nPoint) = val_muT;

  /*--- Allocate and initialize solution for the dual time strategy ---*/
  bool dual_time = ((config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_1ST) ||
                    (config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_2ND));

  if (dual_time) {
    Solution_time_n  = Solution;
    Solution_time_n1 = Solution;
  }


  if (config->GetKind_HybridRANSLES() != NO_HYBRIDRANSLES) {
    lesMode.resize(nPoint) = su2double(0.0);
    if (backscatter) {
      sbsInBox.resize(nPoint) = 1;
      stochSource.resize(nPoint, nDim) = su2double(0.0);
      stochSourceOld.resize(nPoint, nDim) = su2double(0.0);
      besselIntegral.resize(nPoint) = su2double(0.0);
    }
  }

}
