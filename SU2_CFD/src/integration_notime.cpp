/*!
 * \file integration_notime.cpp
 * \brief No time stepping methods for integration a PDE without time.
 * \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
 * \version 2.0.5
 *
 * Stanford University Unstructured (SU2) Code
 * Copyright (C) 2012 Aerospace Design Laboratory
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "../include/integration_structure.hpp"

CPotentialIntegration::CPotentialIntegration(CConfig *config) : CIntegration(config) { }

CPotentialIntegration::~CPotentialIntegration(void) { }

void CPotentialIntegration::SetPotential_Solver(CGeometry ***geometry, CSolution ****solution_container, CNumerics *****solver_container,
                                                CConfig **config, unsigned short RunTime_EqSystem, unsigned short iMesh, unsigned short iZone) {
	
	unsigned short SolContainer_Position = config[iZone]->GetContainerPosition(RunTime_EqSystem);

	/*--- Send-Receive boundary conditions ---*/
	solution_container[iZone][iMesh][SolContainer_Position]->MPI_Send_Receive(geometry, solution_container, config, iMesh, iZone);

	/*--- Do some preprocessing stuff ---*/
	solution_container[iZone][iMesh][SolContainer_Position]->Preprocessing(geometry[iZone][iMesh], solution_container[iZone][iMesh], solver_container[iZone][iMesh][SolContainer_Position],config[iZone], iMesh, 0, RunTime_EqSystem);

	/*--- Space integration ---*/
	Space_Integration(geometry[iZone][iMesh], solution_container[iZone][iMesh], solver_container[iZone][iMesh][SolContainer_Position], config[iZone], iMesh, 0, RunTime_EqSystem);

	/*--- Solve the linear system ---*/
	Solving_Linear_System(geometry[iZone][iMesh], solution_container[iZone][iMesh][SolContainer_Position], solution_container[iZone][iMesh], config[iZone], iMesh);

	/*--- Compute the gradient ---*/
	if (config[iZone]->GetKind_Gradient_Method() == GREEN_GAUSS) solution_container[iZone][iMesh][SolContainer_Position]->SetSolution_Gradient_GG(geometry[iZone][iMesh], config[iZone]);
	if (config[iZone]->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) solution_container[iZone][iMesh][SolContainer_Position]->SetSolution_Gradient_LS(geometry[iZone][iMesh], config[iZone]);

}
