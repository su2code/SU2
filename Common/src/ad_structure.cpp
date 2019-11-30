/*!
 * \file ad_structure.cpp
 * \brief Main subroutines for the algorithmic differentiation (AD) structure.
 * \author T. Albring
 * \version 7.0.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation 
 * (http://su2foundation.org)
 *
 * Copyright 2012-2019, SU2 Contributors (cf. AUTHORS.md)
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

#include "../include/datatype_structure.hpp"

namespace AD {
#ifdef CODI_REVERSE_TYPE
  /*--- Initialization of the global variables ---*/

  int adjointVectorPosition = 0;

  std::vector<su2double::GradientData> inputValues;
  std::vector<su2double::GradientData> localInputValues;
  std::vector<su2double*> localOutputValues;

  su2double::TapeType& globalTape = su2double::getGlobalTape();
  su2double::TapeType::Position StartPosition, EndPosition;
  std::vector<su2double::TapeType::Position> TapePositions;

  bool Status = false;
  bool PreaccActive = false;
  bool PreaccEnabled = true;

  codi::PreaccumulationHelper<su2double> PreaccHelper;

  ExtFuncHelper* FuncHelper;

#endif
}
