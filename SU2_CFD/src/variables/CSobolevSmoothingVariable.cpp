/*!
 * \file CSobolevSmoothingVariable.cpp
 * \brief Definition of the variables for gradient smoothing problems.
 * \author T. Dick
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

#include "../../include/variables/CSobolevSmoothingVariable.hpp"

CSobolevSmoothingVariable::CSobolevSmoothingVariable(unsigned long npoint, unsigned long ndim, const CConfig *config) :
  CVariable(npoint, ndim, config) {

  nDim = ndim;
  nPoint = npoint;

  Sensitivity.resize(nPoint,nDim) = su2double(0.0);

  BoundaryVertexMap.Reset(nPoint);
}

void CSobolevSmoothingVariable::MarkAsBoundaryPoint(unsigned long iPoint) {
  BoundaryVertexMap.SetIsVertex(iPoint, true);
}

bool CSobolevSmoothingVariable::GetIsBoundaryPoint(unsigned long iPoint) const {
  return BoundaryVertexMap.GetIsVertex(iPoint);
}

void CSobolevSmoothingVariable::AllocateBoundaryVariables() {
  if (BoundaryVertexMap.GetIsValid()) return;  // nothing to do

  /*--- Build the map ---*/
  BoundaryVertexMap.Build();
}
