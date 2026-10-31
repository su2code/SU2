/*!
 * \file transition_data.hpp
 * \brief Diagnostic values shared by the LM transition numerics and solver.
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

#pragma once

#include "../../Common/include/basic_types/datatype_structure.hpp"

/*! \brief Diagnostic values of the two-equation and simplified LM transition models. */
struct TransitionLMData {
  su2double criticalReynolds = 1.0;
  su2double momentumThicknessReynolds = 0.0;
  su2double turbulenceIntensity = 0.0;
  su2double pressureGradient = 0.0;
  su2double streamwiseVelocityGradient = 0.0;
  su2double vorticityReynolds = 0.0;
  su2double production = 0.0;
  su2double destruction = 0.0;
  su2double onset1 = 0.0;
  su2double onset2 = 0.0;
  su2double onset3 = 0.0;
  su2double onset = 0.0;
};
