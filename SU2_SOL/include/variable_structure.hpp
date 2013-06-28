/*!
 * \file variable_structure.hpp
 * \brief Headers of the main subroutines for storing all the variables for 
 *        each kind of governing equation (direct, adjoint and linearized).
 *        The subroutines and functions are in the <i>variable_structure.cpp</i> file.
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
#include <iostream>

#include "../../Common/include/config_structure.hpp"

using namespace std;

/*! 
 * \class CVariable
 * \brief Main class for defining the variables.
 * \author F. Palacios.
 * \version 2.0.1
 */
class CVariable {
protected:

	double *Solution;		/*!< \brief Solution of the problem. */
	static unsigned short nDim;		/*!< \brief Number of dimension of the problem. */
	unsigned short nVar;		/*!< \brief Number of variables of the problem, 
													 note that this variable cannnot be static, it is possible to 
													 have different number of nVar in the same problem. */

public:

	/*!
	 * \brief Constructor of the class. 
	 */
	CVariable(void);

	/*!
	 * \overload 
	 * \param[in] val_ndim - Number of dimensions of the problem.		 
	 * \param[in] val_nvar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.	 
	 */
	CVariable(unsigned short val_ndim, unsigned short val_nvar, CConfig *config);

	/*!
	 * \brief Destructor of the class. 
	 */
	virtual ~CVariable(void);

	/*!
	 * \brief Set the value of the solution.
	 * \param[in] val_solution - Solution of the problem.
	 */
	void SetSolution(double *val_solution);

	/*!
	 * \overload
	 * \param[in] val_var - Index of the variable.
	 * \param[in] val_solution - Value of the solution for the index <i>val_var</i>.
	 */
	void SetSolution(unsigned short val_var, double val_solution);

	/*!
	 * \brief Get the solution.
	 * \param[in] val_var - Index of the variable.
	 * \return Value of the solution for the index <i>val_var</i>.
	 */
	double GetSolution(unsigned short val_var);

  /*!
	 * \brief Get the solution of the problem.
	 * \return Pointer to the solution vector.
	 */
	double *GetSolution(void);
  
};

/*! 
 * \class CBaselineVariable
 * \brief Main class for defining the variables of the Euler's solver.
 * \ingroup Euler_Equations
 * \author F. Palacios.
 * \version 2.0.1
 */
class CBaselineVariable : public CVariable {
public:

	/*!
	 * \brief Constructor of the class. 
	 */
	CBaselineVariable(void);

	/*!
	 * \overload
	 * \param[in] val_solution - Pointer to the flow value (initialization value).
	 * \param[in] val_ndim - Number of dimensions of the problem.
	 * \param[in] val_nvar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.	 
	 */		
	CBaselineVariable(double *val_solution, unsigned short val_ndim, unsigned short val_nvar, CConfig *config);

	/*!
	 * \brief Destructor of the class. 
	 */		
	virtual ~CBaselineVariable(void);

};
