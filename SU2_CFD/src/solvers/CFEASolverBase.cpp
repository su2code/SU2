/*!
 * \file CFEASolverBase.cpp
 * \brief Common class template for FEA solvers
 * \author T. Dick
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

#include "../../include/solvers/CFEASolverBase.hpp"
#include <algorithm>

CFEASolverBase::CFEASolverBase(LINEAR_SOLVER_MODE mesh_deform_mode) : CSolver(mesh_deform_mode) {

  nElement = 0;
  nDim = 0;
  nMarker = 0;
  nPoint = 0;
  nPointDomain = 0;

  element_container = new CElement** [MAX_TERMS]();
  for (unsigned short iTerm = 0; iTerm < MAX_TERMS; iTerm++)
    element_container[iTerm] = new CElement* [MAX_FE_KINDS*omp_get_max_threads()]();

}

CFEASolverBase::CFEASolverBase(CGeometry *geometry, CConfig *config, LINEAR_SOLVER_MODE mesh_deform_mode) : CSolver(mesh_deform_mode) {

  nElement      = geometry->GetnElem();
  nDim          = geometry->GetnDim();
  nMarker       = geometry->GetnMarker();
  nPoint        = geometry->GetnPoint();
  nPointDomain  = geometry->GetnPointDomain();

  /*--- Create the baseline element_container ---*/
  element_container = new CElement** [MAX_TERMS]();
  for (unsigned short iTerm = 0; iTerm < MAX_TERMS; iTerm++)
    element_container[iTerm] = new CElement* [MAX_FE_KINDS*omp_get_max_threads()]();

}

CFEASolverBase::~CFEASolverBase() {

  if (element_container != nullptr) {
    for (unsigned int iVar = 0; iVar < MAX_TERMS; iVar++) {
      for (unsigned int jVar = 0; jVar < MAX_FE_KINDS*omp_get_max_threads(); jVar++) {
        delete element_container[iVar][jVar];
      }
      delete [] element_container[iVar];
    }
    delete [] element_container;
  }

}

void CFEASolverBase::CommunicateExtraEliminationVertices(const CGeometry* geometry, vector<unsigned long>& myPoints) {

  /*--- communicate the boundary points ---*/

  const unordered_set<unsigned long> markerPoints(myPoints.begin(), myPoints.end());

  vector<unsigned long> numPoints(size);
  unsigned long num = myPoints.size();
  SU2_MPI::Allgather(&num, 1, MPI_UNSIGNED_LONG, numPoints.data(), 1, MPI_UNSIGNED_LONG, SU2_MPI::GetComm());

  /*--- Global to local map for the halo points of the rank (not covered by the CGeometry map). ---*/
  unordered_map<unsigned long, unsigned long> Global2Local;
  for (auto iPoint = nPointDomain; iPoint < nPoint; ++iPoint) {
    Global2Local[geometry->nodes->GetGlobalIndex(iPoint)] = iPoint;
  }

  /*--- Populate elimination list. ---*/
  ExtraVerticesToEliminate.clear();

  for (int i = 0; i < size; ++i) {
    /*--- Send our point list. ---*/
    if (rank == i) {
      SU2_MPI::Bcast(myPoints.data(), numPoints[i], MPI_UNSIGNED_LONG, rank, SU2_MPI::GetComm());
      continue;
    }

    /*--- Receive point list. ---*/
    vector<unsigned long> theirPoints(numPoints[i]);
    SU2_MPI::Bcast(theirPoints.data(), numPoints[i], MPI_UNSIGNED_LONG, i, SU2_MPI::GetComm());

    for (auto iPointGlobal : theirPoints) {
      /*--- Check if the rank has the point. ---*/
      auto it = Global2Local.find(iPointGlobal);
      if (it == Global2Local.end()) continue;

      /*--- If the point is not covered by this rank's markers, mark it for elimination. ---*/
      if (markerPoints.count(iPointGlobal) == 0)
        ExtraVerticesToEliminate.push_back(it->second);
    }
  }

}
