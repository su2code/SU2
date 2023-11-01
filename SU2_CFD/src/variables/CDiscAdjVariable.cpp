/*!
 * \file CDiscAdjVariable.cpp
 * \brief Main subroutines for the discrete adjoint variable structure.
 * \author T. Albring
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

#include "../../include/variables/CDiscAdjVariable.hpp"

CDiscAdjVariable::CDiscAdjVariable(const su2double* sol, unsigned long npoint, unsigned long ndim,
                                   unsigned long nvar, CConfig *config) :
  CVariable(npoint, ndim, nvar, config, true) {

  if (config->GetTime_Domain())
    DualTime_Derivative.resize(nPoint,nVar) = su2double(0.0);

  if (config->GetTime_Marching() != TIME_MARCHING::STEADY)
    DualTime_Derivative_n.resize(nPoint,nVar) = su2double(0.0);

  Solution_Direct.resize(nPoint,nVar);
  Sensitivity.resize(nPoint,nDim) = su2double(0.0);

  for (unsigned long iPoint = 0; iPoint < nPoint; ++iPoint)
    for (unsigned long iVar = 0; iVar < nVar; ++iVar)
      Solution(iPoint,iVar) = sol[iVar];
}

void CDiscAdjVariable::Set_External_To_DualTimeDer() {
  assert(External.size() == DualTime_Derivative.size());
  parallelCopy(External.size(), DualTime_Derivative.data(), External.data());
}

void CDiscAdjVariable::AllocateAdjointSolutionExtra(unsigned long nVarExtra) {
  if (nVarExtra == 0) return;
  SolutionExtra.resize(nVarExtra) = su2double(1e-16);
  /*--- These are only for multizone, but since nVarExtra is small we allocate by default. ---*/
  SolutionExtra_BGS_k.resize(nVarExtra) = su2double(1e-16);
  ExternalExtra.resize(nVarExtra) = su2double(0.0);
}
