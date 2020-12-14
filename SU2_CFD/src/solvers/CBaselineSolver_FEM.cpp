/*!
 * \file CBaselineSolver_FEM.cpp
 * \brief Main subroutines for CBaselineSolver_FEM class.
 * \author F. Palacios, T. Economon
 * \version 7.0.8 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2020, SU2 Contributors (cf. AUTHORS.md)
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


#include "../../include/solvers/CBaselineSolver_FEM.hpp"


CBaselineSolver_FEM::CBaselineSolver_FEM(void) : CSolver() { }

CBaselineSolver_FEM::CBaselineSolver_FEM(CGeometry *geometry, CConfig *config) {

  SU2_MPI::Error("Not implemented yet", CURRENT_FUNCTION);
}

void CBaselineSolver_FEM::SetOutputVariables(CGeometry *geometry, CConfig *config) {

  SU2_MPI::Error("Not implemented yet", CURRENT_FUNCTION);
}

void CBaselineSolver_FEM::LoadRestart(CGeometry **geometry, CSolver ***solver, CConfig *config, int val_iter, bool val_update_geo) {

  SU2_MPI::Error("Not implemented yet", CURRENT_FUNCTION);
}

CBaselineSolver_FEM::~CBaselineSolver_FEM(void) { }
