/*!
 * \file SU2_CFD.hpp
 * \brief Headers of the main subroutines of the code SU2_CFD.
 *        The subroutines and functions are in the <i>SU2_CFD.cpp</i> file.
 * \author F. Palacios, T. Economon
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
#pragma once

#include "../../Common/include/parallelization/mpi_structure.hpp"
#include "../../Common/include/parallelization/omp_structure.hpp"
#include "CLI11.hpp"

#include "drivers/CDriver.hpp"
#include "drivers/CSinglezoneDriver.hpp"
#include "drivers/CMultizoneDriver.hpp"
#include "drivers/CDiscAdjSinglezoneDriver.hpp"
#include "drivers/CDiscAdjMultizoneDriver.hpp"
#include "drivers/CDummyDriver.hpp"
#include "output/COutput.hpp"
#include "../../Common/include/fem/fem_geometry_structure.hpp"
#include "../../Common/include/geometry/CGeometry.hpp"
#include "../../Common/include/CConfig.hpp"
#include "../include/definition_structure.hpp"
#include "../include/interfaces/CInterface.hpp"

using namespace std;
