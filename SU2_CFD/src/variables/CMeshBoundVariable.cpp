/*!
 * \file CMeshBoundVariable.cpp
 * \brief Definition of the boundary variables for mesh motion using a pseudo-elastic approach.
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

#include "../../include/variables/CMeshBoundVariable.hpp"

CMeshBoundVariable::CMeshBoundVariable(unsigned long npoint, unsigned long ndim, CConfig *config) :
  CMeshVariable(npoint, ndim, config) {

  VertexMap.Reset(nPoint);
}

void CMeshBoundVariable::AllocateBoundaryVariables(CConfig *config) {

  if (VertexMap.GetIsValid()) return; // nothing to do

  /*--- Count number of vertices and build map ---*/

  unsigned long nBoundPt = VertexMap.Build();

  /*--- Allocate ---*/

  Boundary_Displacement.resize(nBoundPt,nDim) = su2double(0.0);
  if (config->GetTime_Domain()) Boundary_Velocity.resize(nBoundPt,nDim) = su2double(0.0);
}

void CMeshBoundVariable::Register_BoundDisp() {
  for (unsigned long iVertex = 0; iVertex < Boundary_Displacement.rows(); iVertex++)
    for (unsigned long iVar = 0; iVar < nVar; iVar++)
      AD::RegisterInput(Boundary_Displacement(iVertex,iVar));
}
