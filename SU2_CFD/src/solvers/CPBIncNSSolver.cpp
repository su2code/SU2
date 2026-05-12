/*!
 * \file CPBIncNSSolver.cpp
 * \brief Main subroutines for solving pressure based Navier-Stokes incompressible flow.
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

#include "../../include/solvers/CPBIncNSSolver.hpp"
#include "../../include/variables/CPBIncNSVariable.hpp"
#include "../../../Common/include/toolboxes/printing_toolbox.hpp"
#include "../../include/solvers/CFVMFlowSolverBase.inl"

/*--- Explicit instantiation of the parent class of CPBIncEulerSolver,
 *    to spread the compilation over two cpp files. ---*/
template class CFVMFlowSolverBase<CPBIncEulerVariable, ENUM_REGIME::INCOMPRESSIBLE>;


CPBIncNSSolver::CPBIncNSSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh) :
  CPBIncEulerSolver(geometry, config, iMesh, true) {
  SU2_ZONE_SCOPED

}
