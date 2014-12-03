/*!
 * \file SU2_MSH.hpp
 * \brief Headers of the main subroutines of the code SU2_MSH.
 *        The subroutines and functions are in the <i>SU2_MSH.cpp</i> file.
 * \author F. Palacios
 * \version 3.2.5 "eagle"
 *
 * Copyright (C) 2012-2014 SU2 <https://github.com/su2code>.
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

#include "../../Common/include/geometry_structure.hpp"
#include "../../Common/include/config_structure.hpp"
#include "../../Common/include/grid_adaptation_structure.hpp"

using namespace std;
