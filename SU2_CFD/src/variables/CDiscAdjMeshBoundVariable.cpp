/*!
 * \file CDiscAdjMeshVariable.cpp
 * \brief Main subroutines for the discrete adjoint mesh variable structure.
 * \author Ruben Sanchez
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


#include "../../include/variables/CDiscAdjMeshBoundVariable.hpp"


CDiscAdjMeshBoundVariable::CDiscAdjMeshBoundVariable(unsigned long npoint, unsigned long ndim, CConfig *config) {

  nPoint = npoint;
  nVar = ndim;
  nDim = ndim;

  /*--- Allocate the solution array. ---*/
  Solution.resize(nPoint,nVar) = su2double(0.0);

  VertexMap.Reset(nPoint);
}

void CDiscAdjMeshBoundVariable::AllocateBoundaryVariables(CConfig *config) {

  if (VertexMap.GetIsValid()) return; // nothing to do

  /*--- Count number of vertices and build map ---*/

  unsigned long nBoundPt = VertexMap.Build();

  /*--- Allocate ---*/

  bool fsi = false;

  /*--- Initialize Boundary Displacement container to 0.0 ---*/

  Bound_Disp_Sens.resize(nBoundPt,nDim) = su2double(0.0);
  Bound_Disp_Direct.resize(nBoundPt,nDim) = su2double(0.0);

  /*--- Container for the BGS solution at the previous iteration ---*/

  if (fsi) Solution_BGS_k.resize(nBoundPt,nDim) = su2double(0.0);

}

void CDiscAdjMeshBoundVariable::Set_BGSSolution_k() {
  Solution_BGS_k = Bound_Disp_Sens;
}

void CDiscAdjMeshBoundVariable::Restore_BGSSolution_k() {
  Bound_Disp_Sens = Solution_BGS_k;
}
