/*!
 * \file SU2_SOL.hpp
 * \brief Headers of the main subroutines of the code SU2_SOL.
 *        The subroutines and functions are in the <i>SU2_SOL.cpp</i> file.
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

#include "../../SU2_CFD/include/solvers/CBaselineSolver.hpp"
#include "../../SU2_CFD/include/solvers/CBaselineSolver_FEM.hpp"
#include "../../SU2_CFD/include/output/CBaselineOutput.hpp"
#include "../../Common/include/geometry/CPhysicalGeometry.hpp"
#include "../../Common/include/CConfig.hpp"

void WriteFiles(CConfig* config, CGeometry* geometry, CSolver** solver_container, COutput* output,
                unsigned long TimeIter);

using namespace std;
