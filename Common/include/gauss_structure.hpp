/*!
 * \file gauss_structure.hpp
 * \brief Headers of the Finite Element structure (gaussian points)
 *        The subroutines and functions are in the <i>gauss_structure.cpp</i> file.
 * \author R. Sanchez
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

#include "mpi_structure.hpp"

#include <cmath>
#include <iostream>
#include <cstdlib>

#include "config_structure.hpp"
#include "geometry_structure.hpp"

using namespace std;


/*!
 * \class CGaussVariable
 * \brief Main class for defining the gaussian points.
 * \author R. Sanchez
 * \version 5.0.0 "Raven"
 */
class CGaussVariable {
protected:

	su2double **GradNi_Xj,		// Gradient of the shape functions N[i] respect to the reference configuration
	**GradNi_xj;			// Gradient of the shape functions N[i] respect to the current configuration
	su2double *Ni;				// Shape functions N[i] at the gaussian point
	su2double J_X,				// Jacobian of the element evaluated at the current Gauss Point respect to the reference configuration
	J_x;					// Jacobian of the element evaluated at the current Gauss Point respect to the current configuration
	unsigned short iGaussPoint;	// Identifier of the Gauss point considered

public:

	/*!
	 * \brief Constructor of the class.
	 */
	CGaussVariable(void);

  /*!
	 * \overload
	 * \param[in] val_nvar - Number of variables of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	//CGaussVariable(unsigned short val_nvar, CConfig *config);

   /*!
	  * \overload
	  * \param[in] val_iGauss - ID of the Gaussian Point
	  * \param[in] val_nDim - Number of dimensions of the problem.
	  * \param[in] config - Definition of the particular problem.
	*/
	CGaussVariable(unsigned short val_iGauss, unsigned short val_nDim, unsigned short val_nNodes);

	/*!
	 * \brief Destructor of the class.
	 */
	virtual ~CGaussVariable(void);

	void SetGradNi_Xj(su2double val_GradNi_Xj, unsigned short val_iDim, unsigned short val_Ni);

	void SetGradNi_xj(su2double val_GradNi_xj, unsigned short val_iDim, unsigned short val_Ni);

	void SetNi(su2double val_ShapeNi, unsigned short val_Ni);

	void SetJ_X(su2double valJ_X);

	void SetJ_x(su2double valJ_x);

	su2double **GetGradNi_Xj(void);

	su2double GetGradNi_Xj(unsigned short val_Ni, unsigned short val_iDim);

	su2double GetGradNi_xj(unsigned short val_Ni, unsigned short val_iDim);

	su2double GetNi(unsigned short val_Ni);

	su2double GetJ_X(void);

	su2double GetJ_x(void);

	unsigned short Get_iGauss(void);

};


#include "gauss_structure.inl"
