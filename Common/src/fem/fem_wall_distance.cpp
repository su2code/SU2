/*!
 * \file fem_wall_distance.cpp
 * \brief Main subroutines for computing the wall distance for the FEM solver.
 * \author E. van der Weide
 * \version 7.0.6 "Blackbird"
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

#include "../../include/fem/fem_geometry_structure.hpp"
#include "../../include/adt_structure.hpp"

std::unique_ptr<CADTElemClass> CMeshFEM_DG::ComputeViscousWallADT(const CConfig *config) const {

  return NULL;

}

/*!
 * \brief Set wall distances a specific value
 */
void CMeshFEM_DG::SetWallDistance(su2double val){

}

void CMeshFEM_DG::SetWallDistance(const CConfig *config, CADTElemClass *WallADT){

}
