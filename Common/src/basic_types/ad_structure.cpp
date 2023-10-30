/*!
 * \file ad_structure.cpp
 * \brief Main subroutines for the algorithmic differentiation (AD) structure.
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

#include "../../include/basic_types/datatype_structure.hpp"

namespace AD {
#ifdef CODI_REVERSE_TYPE
/*--- Initialization of the global variables ---*/

TapePosition StartPosition, EndPosition;
std::vector<TapePosition> TapePositions;

bool PreaccActive = false;
#ifdef HAVE_OPDI
SU2_OMP(threadprivate(PreaccActive))
#endif

bool PreaccEnabled = true;

codi::PreaccumulationHelper<su2double> PreaccHelper;
#ifdef HAVE_OPDI
SU2_OMP(threadprivate(PreaccHelper))
#endif

ExtFuncHelper FuncHelper;

#endif

void Initialize() {
#ifdef CODI_REVERSE_TYPE
  FuncHelper.disableInputPrimalStore();
  FuncHelper.disableOutputPrimalStore();
#endif
}

void Finalize() {}
}  // namespace AD
