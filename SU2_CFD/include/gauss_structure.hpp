/*!
 * \file gauss_structure.hpp
 * \brief Headers of the finite element structure (gaussian points)
 *        The subroutines and functions are in the <i>gauss_structure.cpp</i> file.
 * \author R. Sanchez
 * \version 4.0.0 "Cardinal"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
 *                      Dr. Thomas D. Economon (economon@stanford.edu).
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *
 * Copyright (C) 2012-2015 SU2, the open-source CFD code.
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

#ifdef HAVE_MPI
  #include "mpi.h"
#endif
#include <cmath>
#include <iostream>
#include <cstdlib>

#include "../../Common/include/config_structure.hpp"
#include "../../Common/include/geometry_structure.hpp"

using namespace std;


/*!
 * \class CGaussVariable
 * \brief Main class for defining the gaussian points.
 * \author R. Sanchez
 * \version 4.0.0 "Cardinal"
 */
class CGaussVariable {
protected:

	double **GradNi_Xj,		// Gradient of the shape functions N[i] respect to the reference configuration
	**GradNi_xj;				// Gradient of the shape functions N[i] respect to the current configuration
	double J_X,				// Jacobian of the element evaluated at the current Gauss Point respect to the reference configuration
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

	void SetGradNi_Xj(double val_GradNi_Xj, unsigned short val_iDim, unsigned short val_Ni);

	void SetGradNi_xj(double val_GradNi_xj, unsigned short val_iDim, unsigned short val_Ni);

	void SetJ_X(double valJ_X);

	void SetJ_x(double valJ_x);

	double **GetGradNi_Xj(void);

	double GetGradNi_Xj(unsigned short val_Ni, unsigned short val_iDim);

	double GetJ_X(void);

	double GetJ_x(void);

	unsigned short Get_iGauss(void);

};


#include "gauss_structure.inl"
