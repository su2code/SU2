/*!
 * \file solution_direct_eikonal.cpp
 * \brief Main subrotuines for solving direct problems (Euler, Navier-Stokes, etc.).
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

#include "../include/solution_structure.hpp"

CEikonalSolution::CEikonalSolution(void) : CSolution() { }

CEikonalSolution::~CEikonalSolution(void) { }

CEikonalSolution::CEikonalSolution(CGeometry *geometry, CConfig *config) : CSolution() {
	unsigned long iPoint;
	unsigned short nMarker;
	
	nDim = geometry->GetnDim();
	nVar = 1;
	
	nMarker = config->GetnMarker_All(); 
	node = new CVariable*[geometry->GetnPoint()];
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
		node[iPoint] = new CVariable(nDim, nVar, config);
}

void CEikonalSolution::SetDistance(CGeometry *geometry, CConfig *config) {
	unsigned long iPoint;
	double WallDistance;
	
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
		WallDistance = geometry->node[iPoint]->GetWallDistance();
		node[iPoint]->SetSolution(0, WallDistance);
	}	
}
