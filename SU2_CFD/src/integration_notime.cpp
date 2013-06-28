/*!
 * \file integration_notime.cpp
 * \brief No time stepping methods for integration a PDE without time.
 * \author Current Development: Stanford University.
 *         Original Structure: CADES 1.0 (2009).
 * \version 1.0.
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

void CPotentialIntegration::SetPotential_Solver(CGeometry **geometry, CSolution ***solution_container, 
											 CNumerics ****solver_container, CConfig *config, unsigned short RunTime_EqSystem, unsigned short iMesh) {
	
	unsigned short SolContainer_Position = config->GetContainerPosition(RunTime_EqSystem);
	
	Space_Integration(geometry[iMesh], solution_container[iMesh], solver_container[iMesh][SolContainer_Position], config, iMesh, 0, RunTime_EqSystem);
	Solving_Linear_System(geometry[iMesh], solution_container[iMesh][SolContainer_Position], solution_container[iMesh], config, iMesh);
	
	solution_container[iMesh][SolContainer_Position]->SetSolution_Gradient_GG(geometry[iMesh]);
	
}

CEikonalIntegration::CEikonalIntegration(CConfig *config) : CIntegration(config) { }

CEikonalIntegration::~CEikonalIntegration(void) { }

void CEikonalIntegration::SetEikonal_Solver(CGeometry **geometry, CSolution ***solution_container, CConfig *config) {
	
#ifdef EIKONAL_FEM
	/*--- FEM-solver of Eikonal equation, only 2D ---*/
	solution_container[MESH_0][EIKONAL_SOL]->FEMEikonalSolver(geometry[MESH_0], config);
#else
	/*--- Brute-force computation of distances (non-exact) ---*/
	solution_container[MESH_0][EIKONAL_SOL]->SetDistance(geometry[MESH_0], config);
#endif
}
