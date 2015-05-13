/*!
 * \file definition_structure.hpp
 * \brief Headers of the main subroutines used by SU2_CFD.
 *        The subroutines and functions are in the <i>definition_structure.cpp</i> file.
 * \author F. Palacios, T. Economon
 * \version 3.2.9 "eagle"
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
#include <ctime>

#include "solver_structure.hpp"
#include "integration_structure.hpp"
#include "output_structure.hpp"
#include "numerics_structure.hpp"
#include "../../Common/include/geometry_structure.hpp"
#include "../../Common/include/config_structure.hpp"

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
 * \brief Definition and allocation of all solution classes.
 * \param[in] solver_container - Container vector with all the solutions.
 * \param[in] geometry - Geometrical definition of the problem.
 * \param[in] config - Definition of the particular problem.
 * \param[in] iZone - Index of the zone.
 */
void Solver_Preprocessing(CSolver ***solver_container, CGeometry **geometry, CConfig *config, unsigned short iZone);

/*! 
 * \brief De-allocation of all solution classes.
 * \param[in] solver_container - Container vector with all the solutions.
 * \param[in] geometry - Geometrical definition of the problem.
 * \param[in] config - Definition of the particular problem.
 * \param[in] iZone - Index of the zone.
 */
void Solver_Postprocessing(CSolver ***solver_container, CGeometry **geometry, CConfig *config, unsigned short iZone);


/*!
 * \brief Definition and allocation of all integration classes.
 * \param[in] integration_container - Container vector with all the integration methods.
 * \param[in] geometry - Geometrical definition of the problem.
 * \param[in] config - Definition of the particular problem.
 * \param[in] iZone - Index of the zone.
 */
void Integration_Preprocessing(CIntegration **integration_container, CGeometry **geometry, CConfig *config, unsigned short iZone);

/*! 
 * \brief Definition and allocation of all solver classes.
 * \param[in] numerics_container - Description of the numerical method (the way in which the equations are solved).
 * \param[in] solver_container - Container vector with all the solutions.
 * \param[in] geometry - Geometrical definition of the problem.
 * \param[in] config - Definition of the particular problem.
 * \param[in] iZone - Index of the zone.
 */
void Numerics_Preprocessing(CNumerics ****numerics_container, CSolver ***solver_container, CGeometry **geometry, CConfig *config, unsigned short iZone);

/*! 
 * \brief De-allocation of all solver classes.
 * \param[in] numerics_container - Description of the numerical method (the way in which the equations are solved).
 * \param[in] solver_container - Container vector with all the solutions.
 * \param[in] geometry - Geometrical definition of the problem.
 * \param[in] config - Definition of the particular problem.
 * \param[in] iZone - Index of the zone.
 */
void Numerics_Postprocessing(CNumerics ****numerics_container, CSolver ***solver_container, CGeometry **geometry, CConfig *config, unsigned short iZone);


/*!
 * \brief Do the geometrical preprocessing.
 * \param[in] geometry - Geometrical definition of the problem.
 * \param[in] config - Definition of the particular problem.
 * \param[in] val_nZone - Total number of zones.
 */
void Geometrical_Preprocessing(CGeometry ***geometry, CConfig **config, unsigned short val_nZone);
