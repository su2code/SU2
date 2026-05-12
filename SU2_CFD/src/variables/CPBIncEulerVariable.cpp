/*!
 * \file CPBIncEulerVariable.cpp
 * \brief Definition of the variable classes for pressure based incompressible flow.
 * \author F. Palacios, T. Economon
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

#include "../../include/variables/CPBIncEulerVariable.hpp"
#include "../../include/fluid/CFluidModel.hpp"

CPBIncEulerVariable::CPBIncEulerVariable(su2double pressure, const su2double *velocity, su2double enthalpy,
                                     unsigned long npoint, unsigned long ndim, unsigned long nvar, const CConfig *config)
  : CFlowVariable(npoint, ndim, nvar, ndim + 10,
                  ndim + (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED ? 2 : 4), config),
    indices(ndim, 0) {

}
