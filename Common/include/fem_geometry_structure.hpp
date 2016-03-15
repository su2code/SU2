/*!
 * \file fem_geometry_structure.hpp
 * \brief Headers of the main subroutines for creating the geometrical structure for the FEM solver.
 *        The subroutines and functions are in the <i>fem_geometry_structure.cpp</i> file.
 * \author E. van der Weide
 * \version 4.1.0 "Cardinal"
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

#include "./geometry_structure.hpp"

using namespace std;

/*!
 * \class CMeshFEM
 * \brief Base class for the FEM solver.
 * \author E. van der Weide
 * \version 4.1.0 "Cardinal"
 */
class CMeshFEM: public CGeometry {
protected:

public:
  
	/*!
	 * \brief Constructor of the class.
	 */
	CMeshFEM(void);

  /*!
	 * \overload
	 * \brief Redistributes the grid over the ranks and creates the halo layer.
  * \param[in] geometry - The linear distributed grid that must be redistributed.
	 * \param[in] config   - Definition of the particular problem.
	 */
  CMeshFEM(CGeometry *geometry, CConfig *config);
  
	/*!
	 * \brief Destructor of the class.
	 */
	virtual ~CMeshFEM(void);
};

/*!
 * \class CMeshFEM_DG
 * \brief Class which contains all the variables for the DG FEM solver.
 * \author E. van der Weide
 * \version 4.1.0 "Cardinal"
 */
class CMeshFEM_DG: public CMeshFEM {

public:

 /*!
  * \brief Constructor of the class.
  */
 CMeshFEM_DG(void);

  /*!
  * \overload
  * \brief Redistributes the grid over the ranks and creates the halo layer.
  * \param[in] geometry - The linear distributed grid that must be redistributed.
  * \param[in] config   - Definition of the particular problem.
  */
  CMeshFEM_DG(CGeometry *geometry, CConfig *config);
 
 /*!
  * \brief Destructor of the class.
  */
 ~CMeshFEM_DG(void);
};

#include "fem_geometry_structure.inl"
