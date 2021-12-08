/*!
 * \file CSobolevSmoothingVariable.cpp
 * \brief Definition of the variables for gradient smoothing problems.
 * \author T. Dick
 * \version 7.2.1 "Blackbird"
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

#include "../../include/variables/CSobolevSmoothingVariable.hpp"

CSobolevSmoothingVariable::CSobolevSmoothingVariable(unsigned long npoint, unsigned long ndim, CConfig *config) :
  CVariable(npoint, ndim, config) {

  nDim = ndim;
  nPoint = npoint;

  Sensitivity.resize(nPoint,nDim) = su2double(0.0);

  boundary_vertex = new bool[nPoint];
  for (unsigned long iPoint=0; iPoint<nPoint; iPoint++) {
    boundary_vertex[iPoint] = false;
  }
  nBoundPoints = 0;
}

CSobolevSmoothingVariable::~CSobolevSmoothingVariable() {
  if (boundary_vertex != NULL) delete [] boundary_vertex;
}

void CSobolevSmoothingVariable::MarkAsBoundaryPoint(unsigned long iPoint) {
  boundary_vertex[iPoint] = true;
  nBoundPoints++;
}

bool CSobolevSmoothingVariable::IsBoundaryPoint(unsigned long iPoint) {
  return boundary_vertex[iPoint];
}

unsigned int CSobolevSmoothingVariable::GetNBoundPoints() {
  return nBoundPoints;
}
