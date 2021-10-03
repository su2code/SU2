/*!
 * \file turb_sources.cpp
 * \brief Implementation of numerics classes for integration of
 *        turbulence source-terms.
 * \author T. Kattmann
 * \version 7.2.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2021, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../../include/numerics/species/species_sources.hpp"

#include "../../../include/numerics/CNumerics.hpp"

CSourceSpecies::CSourceSpecies(unsigned short val_nDim, unsigned short val_nVar, const CConfig* config)
    : CNumerics(val_nDim, val_nVar, config) {
  residual = new su2double[nVar]();
  jacobian = new su2double*[nVar];

  for (unsigned short iVar = 0; iVar < nVar; ++iVar) jacobian[iVar] = new su2double[nVar]();
}

CSourceSpecies::~CSourceSpecies() {
  delete[] residual;
  if (jacobian) {
    for (unsigned short iVar = 0; iVar < nVar; ++iVar) delete[] jacobian[iVar];
    delete[] jacobian;
  }
}

CNumerics::ResidualType<> CSourceSpecies::ComputeResidual(const CConfig* config) {
  residual = new su2double[nVar]();
  jacobian = new su2double*[nVar];

  for (unsigned short iVar = 0; iVar < nVar; ++iVar) jacobian[iVar] = new su2double[nVar]();

  return ResidualType<>(residual, jacobian, nullptr);
}
