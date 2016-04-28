/*!
 * \file inlet_funct.hpp
 * \brief Headers of the for the inlet functions.
 * \author D. Manosalvas
 * \version 4.1.2 "Cardinal"
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
 * Copyright (C) 2012-2016 SU2, the open-source CFD code.
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

#include "../../Common/include/mpi_structure.hpp"

#include <cmath>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <complex>

#include "fluid_model.hpp"
#include "numerics_structure.hpp"
#include "variable_structure.hpp"
#include "../../Common/include/gauss_structure.hpp"
#include "../../Common/include/element_structure.hpp"
#include "../../Common/include/geometry_structure.hpp"
#include "../../Common/include/config_structure.hpp"
#include "../../Common/include/matrix_structure.hpp"
#include "../../Common/include/vector_structure.hpp"
#include "../../Common/include/linear_solvers_structure.hpp"
#include "../../Common/include/grid_movement_structure.hpp"

using namespace std;

/*!
 * \brief Transform the real location value to a number between -1 and 1.
 * \return A value between -1 and 1.
 */
su2double ScaleCoordinate(su2double y_max, su2double y_min, su2double y);

/*!
 * \brief Evaluates a polynomial thst represent the inlet velocity profile.
 * \return Normalized velocity profile from 0 to 1.
 */
su2double poly2D(su2double C1, su2double C2, su2double C3, su2double C4, su2double C5, su2double y);


/*!
 * \brief Evaluates a piecewise velocity profile and includes the amplitude.
 * \return Velocity profile with amplitude A.
 */
su2double polydisc(su2double A , su2double y_max, su2double y_min, su2double y);