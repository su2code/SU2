/*!
 * \file definition_structure.hpp
 * \brief Headers of the main subroutines used by SU2_CFD.
 *        The subroutines and functions are in the <i>definition_structure.cpp</i> file.
 * \author Current Development: Stanford University.
 *         Original Structure: CADES 1.0 (2009).
 * \version 1.1.
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

#include <ctime>

#ifndef NO_MPI
#include <mpi.h>
#endif

#include "solution_structure.hpp"
#include "integration_structure.hpp"
#include "output_structure.hpp"
#include "numerics_structure.hpp"
#include "../../Common/include/geometry_structure.hpp"
#include "../../Common/include/config_structure.hpp"

using namespace std;

/*! 
 * \brief Gets the number of domains in the mesh file
 * \param[in] config - Definition of the particular problem.
 * \param[in] val_mesh_filename - Name of the file with the grid information.
 * \param[in] val_format - Format of the file with the grid information.
 * \return Total number of domains in the grid file.
 */
unsigned short GetnDomain(CConfig *config, string val_mesh_filename, unsigned short val_format);

/*! 
 * \brief Definition and allocation of all solution classes.
 * \param[in] solution_container - Container vector with all the solutions.
 * \param[in] geometry - Geometrical definition of the problem.
 * \param[in] config - Definition of the particular problem.
 */
void Solution_Definition(CSolution ***solution_container, CGeometry **geometry, CConfig *config);

/*! 
 * \brief Definition and allocation of all integration classes.
 * \param[in] integration_container - Container vector with all the integration methods.
 * \param[in] geometry - Geometrical definition of the problem.
 * \param[in] config - Definition of the particular problem.
 */
void Integration_Definition(CIntegration **integration_container, CGeometry **geometry, CConfig *config);

/*! 
 * \brief Definition and allocation of all solver classes.
 * \param[in] solver_container - Description of the numerical method (the way in which the equations are solved).
 * \param[in] solution_container - Container vector with all the solutions.
 * \param[in] geometry - Geometrical definition of the problem.
 * \param[in] config - Definition of the particular problem.
 */
void Solver_Definition(CNumerics ****solver_container, CSolution ***solution_container, CGeometry **geometry, CConfig *config);

/*! 
 * \brief Do the geometrical preprocessing.
 * \param[in] geometry - Geometrical definition of the problem.
 * \param[in] config - Definition of the particular problem.
 * \param[in] val_nDomain - Total number of domains in the grid file.
 */
void Geometrical_Definition(CGeometry ***geometry, CConfig *config, unsigned short val_nDomain);

/*!
 * \brief Deallocation of the pointers solver_container, solution_container, integration_container, geometry, output, config.
 * \param[in] solver_container - Description of the numerical method (the way in which the equations are solved).
 * \param[in] solution_container - Container vector with all the solutions.
 * \param[in] integration_container - Container vector with all the integration methods.
 * \param[in] output - Pointer to the COutput class.
 * \param[in] geometry - Geometrical definition of the problem.
 * \param[in] config - Definition of the particular problem.
 */
void Solver_Deallocation(CNumerics ****solver_container, CSolution ***solution_container, CIntegration **integration_container,
					  COutput *output, CGeometry **geometry, CConfig *config);

/*! 
 * \brief Deallocation of the pointers geometry.
 * \param[in] geometry - Geometrical definition of the problem.
 * \param[in] config - Definition of the particular problem.
 */
void Geometrical_Deallocation(CGeometry **geometry, CConfig *config);
