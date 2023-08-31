/*!
 * \file CRadVariable.cpp
 * \brief Definition of the radiation variables
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

#include "../../include/variables/CRadVariable.hpp"

CRadVariable::CRadVariable(unsigned long npoint, unsigned long ndim, unsigned long nvar, CConfig *config)
: CVariable(npoint, ndim, nvar, config) {

  /*--- The first term is the source term */
  /*--- The second term is the Jacobian   */
  Radiative_SourceTerm.resize(nPoint,2) = su2double(0.0);

  /*--- We need to initialize some containers */
  /*--- Gradient related fields ---*/
  Gradient.resize(nPoint,nVar,nDim,0.0);

  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {
    Rmatrix.resize(nPoint,nDim,nDim,0.0);
  }

  Max_Lambda_Visc.resize(nPoint);
  Delta_Time.resize(nPoint);

  /* Volumetric heat source boolean initialization. */
  Vol_HeatSource.resize(nPoint) = false;

}
