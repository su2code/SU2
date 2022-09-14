/*!
 * \file species_sources.cpp
 * \brief Implementation of numerics classes for integration of
 *        species transport source-terms.
 * \author T. Kattmann
 * \version 7.3.1 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2022, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../../include/variables/CEulerVariable.hpp"
#include "../../../include/variables/CIncEulerVariable.hpp"
#include "../../../include/variables/CNEMOEulerVariable.hpp"

CSourceBase_Species::CSourceBase_Species(unsigned short val_nDim, unsigned short val_nVar, const CConfig* config)
    : CNumerics(val_nDim, val_nVar, config) {
  residual = new su2double[nVar]();
  jacobian = new su2double*[nVar];
  for (unsigned short iVar = 0; iVar < nVar; iVar++) {
    jacobian[iVar] = new su2double[nVar]();
  }
}

CSourceBase_Species::~CSourceBase_Species() {
  delete[] residual;
  if (jacobian) {
    for (unsigned short iVar = 0; iVar < nVar; iVar++) {
      delete[] jacobian[iVar];
    }
    delete[] jacobian;
  }
}

template <class T>
CSourceAxisymmetric_Species<T>::CSourceAxisymmetric_Species(unsigned short val_nDim, unsigned short val_nVar,
                                                         const CConfig* config)
    : CSourceBase_Species(val_nDim, val_nVar, config),
    idx(val_nDim, config->GetnSpecies()) {
  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
}

template <class T>
CNumerics::ResidualType<> CSourceAxisymmetric_Species<T>::ComputeResidual(const CConfig* config) {
  for (unsigned short iVar = 0; iVar < nVar; iVar++) residual[iVar] = 0.0;

  if (implicit) {
    for (unsigned short iVar = 0; iVar < nVar; iVar++) {
      for (unsigned short jVar = 0; jVar < nVar; jVar++) jacobian[iVar][jVar] = 0.0;
    }
  }

  return ResidualType<>(residual, jacobian, nullptr);
}

/*--- Explicit instantiations until we don't move this to the hpp. ---*/
template class CSourceAxisymmetric_Species<CEulerVariable::CIndices<unsigned short> >;
template class CSourceAxisymmetric_Species<CIncEulerVariable::CIndices<unsigned short> >;
template class CSourceAxisymmetric_Species<CNEMOEulerVariable::CIndices<unsigned short> >;
