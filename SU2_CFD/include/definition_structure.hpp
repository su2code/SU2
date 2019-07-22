/*!
 * \file definition_structure.hpp
 * \brief Headers of the main subroutines used by SU2_CFD.
 *        The subroutines and functions are in the <i>definition_structure.cpp</i> file.
 * \author F. Palacios, T. Economon
 * \version 6.2.0 "Falcon"
 *
 * The current SU2 release has been coordinated by the
 * SU2 International Developers Society <www.su2devsociety.org>
 * with selected contributions from the open-source community.
 *
 * The main research teams contributing to the current release are:
 *  - Prof. Juan J. Alonso's group at Stanford University.
 *  - Prof. Piero Colonna's group at Delft University of Technology.
 *  - Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *  - Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *  - Prof. Rafael Palacios' group at Imperial College London.
 *  - Prof. Vincent Terrapon's group at the University of Liege.
 *  - Prof. Edwin van der Weide's group at the University of Twente.
 *  - Lab. of New Concepts in Aeronautics at Tech. Institute of Aeronautics.
 *
 * Copyright 2012-2019, Francisco D. Palacios, Thomas D. Economon,
 *                      Tim Albring, and the SU2 contributors.
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
#include "../../Common/include/fem_geometry_structure.hpp"
#include "../../Common/include/geometry_structure.hpp"
#include "../../Common/include/config_structure.hpp"
#include "../../Common/include/interpolation_structure.hpp"



using namespace std;

/*!
 * \brief Performs an analysis of the mesh partitions for distributed memory calculations.
 * \param[in] geometry - Geometrical definition of the problem.
 * \param[in] config - Definition of the particular problem.
 */
void Partition_Analysis(CGeometry *geometry, CConfig *config);

/*!
 * \brief Performs an analysis of the mesh partitions for distributed memory calculations for the FEM solver.
 * \param[in] geometry - Geometrical definition of the problem.
 * \param[in] config - Definition of the particular problem.
 */
void Partition_Analysis_FEM(CGeometry *geometry, CConfig *config);
