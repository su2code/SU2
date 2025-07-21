/*!
 * \file SU2_GEO.hpp
 * \brief Headers of the main subroutines of the code SU2_GEO.
 *        The subroutines and functions are in the <i>SU2_GEO.cpp</i> file.
 * \author F. Palacios, T. Economon
 * \version 8.2.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2025, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../Common/include/parallelization/mpi_structure.hpp"

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>

#include "../../Common/include/geometry/CPhysicalGeometry.hpp"
#include "../../Common/include/CConfig.hpp"
#include "../../Common/include/grid_movement/CSurfaceMovement.hpp"
#include "../../Common/include/grid_movement/CFreeFormDefBox.hpp"

using namespace std;
