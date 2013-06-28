/*!
 * \file solution_baseline.hpp
 * \brief Headers of the main subroutines for solving partial differential equations. 
 *        The subroutines and functions are in the <i>solution_structure.cpp</i>, 
 *        <i>solution_direct.cpp</i>, <i>solution_adjoint.cpp</i>, and 
 *        <i>solution_linearized.cpp</i> files.
 * \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
 * \version 2.0.2
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

#include <cmath>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>
#ifndef NO_MPI
#include <mpi.h>
#endif

#include "variable_structure.hpp"
#include "../../Common/include/geometry_structure.hpp"
#include "../../Common/include/config_structure.hpp"

using namespace std;

/*! 
 * \class CSolution
 * \brief Main class for defining the PDE solution, it requires 
 * a child class for each particular solver (Euler, Navier-Stokes, Plasma, etc.) 
 * \author F. Palacios.
 * \version 2.0.1
 */
class CSolution {
protected:
	unsigned short nVar,					/*!< \brief Number of variables of the problem. */
	nDim;													/*!< \brief Number of dimensions of the problem. */
	unsigned long nPoint;					/*!< \brief Number of points of the computational grid. */

public:
    
	CVariable** node;	/*!< \brief Vector which the define the variables for each problem. */

	/*! 
	 * \brief Constructor of the class. 
	 */
	CSolution(void);

	/*! 
	 * \brief Destructor of the class. 
	 */
	virtual ~CSolution(void);

    /*!
	 * \brief Get the number of variables of the problem.
	 */
	unsigned short GetnVar(void);
};

/*! 
 * \class CBaselineSolution
 * \brief Main class for defining the Euler's flow solver.
 * \ingroup Euler_Equations
 * \author F. Palacios.
 * \version 2.0.1
 */
class CBaselineSolution : public CSolution {
public:

	/*! 
	 * \brief Constructor of the class. 
	 */
	CBaselineSolution(void);

	/*! 
	 * \overload
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CBaselineSolution(CGeometry *geometry, CConfig *config, unsigned short iMesh);

	/*! 
	 * \brief Destructor of the class. 
	 */
	virtual ~CBaselineSolution(void);

};
