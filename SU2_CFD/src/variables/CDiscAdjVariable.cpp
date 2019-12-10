/*!
 * \file CDiscAdjVariable.cpp
 * \brief Main subroutines for the discrete adjoint variable structure.
 * \author T. Albring
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


#include "../../include/variables/CDiscAdjVariable.hpp"


CDiscAdjVariable::CDiscAdjVariable(const su2double* sol, unsigned long npoint, unsigned long ndim, unsigned long nvar, CConfig *config)
  : CVariable(npoint, ndim, nvar, config) {

  bool dual_time = (config->GetTime_Marching() == DT_STEPPING_1ST) ||
                   (config->GetTime_Marching() == DT_STEPPING_2ND);

  bool fsi = config->GetFSI_Simulation();

  if (dual_time) {
    DualTime_Derivative.resize(nPoint,nVar) = su2double(0.0);
    DualTime_Derivative_n.resize(nPoint,nVar) = su2double(0.0);

    Solution_time_n.resize(nPoint,nVar) = su2double(0.0);
    Solution_time_n1.resize(nPoint,nVar) = su2double(0.0);
  }

  Solution_Direct.resize(nPoint,nVar);
  Sensitivity.resize(nPoint,nDim) = su2double(0.0);

  for (unsigned long iPoint = 0; iPoint < nPoint; ++iPoint)
    for (unsigned long iVar = 0; iVar < nVar; ++iVar)
      Solution(iPoint,iVar) = sol[iVar];

  if (fsi) {
    Geometry_Direct.resize(nPoint,nDim) = su2double(0.0);
    Solution_Geometry.resize(nPoint,nDim) = su2double(1e-16);
    Solution_Geometry_Old.resize(nPoint,nDim) = su2double(0.0);
    Cross_Term_Derivative.resize(nPoint,nVar) = su2double(0.0);
    Geometry_CrossTerm_Derivative.resize(nPoint,nDim) = su2double(0.0);
    Geometry_CrossTerm_Derivative_Flow.resize(nPoint,nDim) = su2double(0.0);

    Solution_BGS.resize(nPoint,nVar) = su2double(0.0);
    Solution_Geometry_BGS_k.resize(nPoint,nDim) = su2double(0.0);
  }

  if (config->GetMultizone_Problem()) {
    External.resize(nPoint,nVar) = su2double(0.0);
    Solution_BGS.resize(nPoint,nVar) = su2double(0.0);
  }
}

void CDiscAdjVariable::Set_OldSolution_Geometry() { Solution_Geometry_Old = Solution_Geometry; }
