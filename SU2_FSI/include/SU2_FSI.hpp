/*!
 * \file SU2_FSI.hpp
 * \brief Headers of the main subroutines of the code SU2_FSI.
 *        The subroutines and functions are in the <i>SU2_FSI.cpp</i> file.
 * \author R. Sanchez, F. Palacios, T. Economon
 * \version 3.2.9 "eagle"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
 *                      Dr. Thomas D. Economon (economon@stanford.edu).
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
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

#ifdef HAVE_MPI
  #include "mpi.h"
#endif
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>

#include "../../SU2_CFD/include/solver_structure.hpp"
#include "../../SU2_CFD/include/output_structure.hpp"
#include "../../SU2_CFD/include/integration_structure.hpp"
#include "../../SU2_CFD/include/numerics_structure.hpp"
#include "../../SU2_CFD/include/definition_structure.hpp"
#include "../../SU2_CFD/include/iteration_structure.hpp"

#include "../../Common/include/geometry_structure.hpp"
#include "../../Common/include/config_structure.hpp"
#include "../../Common/include/grid_movement_structure.hpp"

using namespace std;
