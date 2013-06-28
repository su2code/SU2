/*!
 * \file definition_structure.hpp
 * \brief Headers of the main subroutines used by SU2_CFD.
 *        The subroutines and functions are in the <i>definition_structure.cpp</i> file.
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
 * \brief Definition of all the different solver possibilities and allocate the solutions.
 * \param[in] solver_container - Description of the numerical method (the way in which the equations are solved).
 * \param[in] solution_container - Container vector with all the solutions.
 * \param[in] integration_container - Container vector with all the integration methods.
 * \param[in] geometry - Geometrical definition of the problem.
 * \param[in] config - Definition of the particular problem.
 */
void Solver_Definition(CNumerics ****solver_container, CSolution ***solution_container, CIntegration **integration_container, 
					   CGeometry **geometry, CConfig *config);

/*! 
 * \brief Do the geometrical preprocessing.
 * \param[in] geometry - Geometrical definition of the problem.
 * \param[in] config - Definition of the particular problem.
 */
void Geometrical_Definition(CGeometry **geometry, CConfig *config);

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
