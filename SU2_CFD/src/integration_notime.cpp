/*!
 * \file integration_notime.cpp
 * \brief No time stepping methods for integration a PDE without time.
 * \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
 * \version 2.0.7
 *
 * Stanford University Unstructured (SU2).
 * Copyright (C) 2012-2013 Aerospace Design Laboratory (ADL).
 *
 * SU2 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * SU2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SU2. If not, see <http://www.gnu.org/licenses/>.
 */

#include "../include/integration_structure.hpp"

CPotentialIntegration::CPotentialIntegration(CConfig *config) : CIntegration(config) { }

CPotentialIntegration::~CPotentialIntegration(void) { }

void CPotentialIntegration::SetPotential_Solver(CGeometry ***geometry, CSolver ****solver_container, CNumerics *****numerics_container,
                                                CConfig **config, unsigned short RunTime_EqSystem, unsigned short iMesh, unsigned short iZone) {
	
	unsigned short SolContainer_Position = config[iZone]->GetContainerPosition(RunTime_EqSystem);

	/*--- Do some preprocessing stuff ---*/
	solver_container[iZone][iMesh][SolContainer_Position]->Preprocessing(geometry[iZone][iMesh], solver_container[iZone][iMesh], config[iZone], iMesh, 0, RunTime_EqSystem);

	/*--- Space integration ---*/
	Space_Integration(geometry[iZone][iMesh], solver_container[iZone][iMesh], numerics_container[iZone][iMesh][SolContainer_Position], config[iZone], iMesh, 0, RunTime_EqSystem);

	/*--- Solve the linear system ---*/
	Solving_Linear_System(geometry[iZone][iMesh], solver_container[iZone][iMesh][SolContainer_Position], solver_container[iZone][iMesh], config[iZone], iMesh);

	/*--- Compute the gradient ---*/
	if (config[iZone]->GetKind_Gradient_Method() == GREEN_GAUSS) solver_container[iZone][iMesh][SolContainer_Position]->SetSolution_Gradient_GG(geometry[iZone][iMesh], config[iZone]);
	if (config[iZone]->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) solver_container[iZone][iMesh][SolContainer_Position]->SetSolution_Gradient_LS(geometry[iZone][iMesh], config[iZone]);

}
