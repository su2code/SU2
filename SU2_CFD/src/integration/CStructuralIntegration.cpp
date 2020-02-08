/*!
 * \file CStructuralIntegration.cpp
 * \brief Space and time integration for structural problems.
 * \author F. Palacios, T. Economon
 * \version 7.0.1 "Blackbird"
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

#include "../../include/integration/CStructuralIntegration.hpp"


CStructuralIntegration::CStructuralIntegration(CConfig *config) : CIntegration(config) { }

void CStructuralIntegration::Structural_Iteration(CGeometry ****geometry, CSolver *****solver_container,
                                                  CNumerics ******numerics_container, CConfig **config, unsigned short RunTime_EqSystem, unsigned short iZone, unsigned short iInst) {

  unsigned short SolContainer_Position = config[iZone]->GetContainerPosition(RunTime_EqSystem);

  /*--- Preprocessing ---*/

  solver_container[iZone][iInst][MESH_0][SolContainer_Position]->Preprocessing(geometry[iZone][iInst][MESH_0], solver_container[iZone][iInst][MESH_0],
      config[iZone], numerics_container[iZone][iInst][MESH_0][SolContainer_Position], MESH_0, NO_RK_ITER, RunTime_EqSystem, false);


  /*--- Space integration ---*/

  Space_Integration_FEM(geometry[iZone][iInst][MESH_0], solver_container[iZone][iInst][MESH_0], numerics_container[iZone][iInst][MESH_0][SolContainer_Position],
                    config[iZone], RunTime_EqSystem);

  /*--- Time integration ---*/

  Time_Integration_FEM(geometry[iZone][iInst][MESH_0], solver_container[iZone][iInst][MESH_0], numerics_container[iZone][iInst][MESH_0][SolContainer_Position],
                config[iZone], RunTime_EqSystem);

  /*--- Postprocessing ---*/

  solver_container[iZone][iInst][MESH_0][SolContainer_Position]->Postprocessing(geometry[iZone][iInst][MESH_0], solver_container[iZone][iInst][MESH_0],
      config[iZone], numerics_container[iZone][iInst][MESH_0][SolContainer_Position],  MESH_0);

}

