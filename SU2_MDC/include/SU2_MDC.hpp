/*!
 * \file SU2_MDC.hpp
 * \brief Headers of the main subroutines of the code SU2_MDC.
 *        The subroutines and functions are in the <i>SU2_MDC.cpp</i> file.
 * \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
 * \version 2.0.6
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

#pragma once

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>

#include "../../Common/include/geometry_structure.hpp"
#include "../../Common/include/config_structure.hpp"
#include "../../Common/include/grid_movement_structure.hpp"

using namespace std;

/*!
 * \brief Gets the number of zones in the mesh file for SU2_MDC.
 * \param[in] val_mesh_filename - Name of the file with the grid information.
 * \param[in] val_format - Format of the file with the grid information.
 * \param[in] config - Definition of the particular problem.
 * \return Total number of domains in the grid file.
 */
unsigned short GetnZone(string val_mesh_filename, unsigned short val_format, CConfig *config);

/*!
 * \brief Write a multi-zone mesh file after deforming the mesh using SU2_MDC.
 * \param[in] geometry - Multi-zone physical geometry container.
 * \param[in] config - Definition of the particular problem.
 * \param[in] nZone - Number of zones in the mesh.
 */
void SetMultiZone_MeshFile(CPhysicalGeometry **geometry, CConfig **config, unsigned short nZone);

/*!
 * \brief Write a single-zone mesh file after deforming the mesh using SU2_MDC.
 * \param[in] geometry - Multi-zone physical geometry container.
 * \param[in] config - Definition of the particular problem.
 * \param[in] surface_movement - durface information including FFD box.
 */
void SetSingleZone_MeshFile(CPhysicalGeometry *geometry, CConfig *config, CSurfaceMovement *surface_movement);