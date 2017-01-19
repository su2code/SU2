/*!
 * \file definition_structure.hpp
 * \brief Headers of the main subroutines used by SU2_CFD.
 *        The subroutines and functions are in the <i>definition_structure.cpp</i> file.
 * \author F. Palacios, T. Economon
 * \version 5.0.0 "Raven"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
 *                      Dr. Thomas D. Economon (economon@stanford.edu).
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *                 Prof. Edwin van der Weide's group at the University of Twente.
 *                 Prof. Vincent Terrapon's group at the University of Liege.
 *
 * Copyright (C) 2012-2017 SU2, the open-source CFD code.
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

#pragma once

#include "../../Common/include/mpi_structure.hpp"

#include <ctime>

#include "../../Common/include/gauss_structure.hpp"
#include "../../Common/include/element_structure.hpp"
#include "driver_structure.hpp"
#include "iteration_structure.hpp"
#include "solver_structure.hpp"
#include "integration_structure.hpp"
#include "output_structure.hpp"
#include "numerics_structure.hpp"
#include "transfer_structure.hpp"
#include "../../Common/include/geometry_structure.hpp"
#include "../../Common/include/config_structure.hpp"
#include "../../Common/include/interpolation_structure.hpp"



using namespace std;

/*! 
 * \brief Gets the number of zones in the mesh file.
 * \param[in] val_mesh_filename - Name of the file with the grid information.
 * \param[in] val_format - Format of the file with the grid information.
 * \param[in] config - Definition of the particular problem.
 * \return Total number of zones in the grid file.
 */
unsigned short GetnZone(string val_mesh_filename, unsigned short val_format, CConfig *config);

/*! 
 * \brief Gets the number of dimensions in the mesh file
 * \param[in] val_mesh_filename - Name of the file with the grid information.
 * \param[in] val_format - Format of the file with the grid information.
 * \return Total number of domains in the grid file.
 */
unsigned short GetnDim(string val_mesh_filename, unsigned short val_format);

/*!
 * \brief Performs an analysis of the mesh partitions for distributed memory calculations.
 * \param[in] geometry - Geometrical definition of the problem.
 * \param[in] config - Definition of the particular problem.
 */
void Partition_Analysis(CGeometry *geometry, CConfig *config);
