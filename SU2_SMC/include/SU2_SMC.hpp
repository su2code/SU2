/*!
 * \file SU2_SMC.hpp
 * \brief Headers of the main subroutines of the code SU2_SMC.
 *        The subroutines and functions are in the <i>SU2_SMC.cpp</i> file.
 * \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
 * \version 2.0.6
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

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>

#include "../../Common/include/geometry_structure.hpp"
#include "../../Common/include/config_structure.hpp"

using namespace std;

/*!
 * \brief Gets the number of zones in the mesh file for SU2_SMC.
 * \author T. Economon
 * \param[in] val_mesh_filename - Name of the file with the grid information.
 * \param[in] val_format - Format of the file with the grid information.
 * \param[in] config - Definition of the particular problem.
 * \return Total number of domains in the grid file.
 */
unsigned short GetnZone(string val_mesh_filename, unsigned short val_format, CConfig *config);

/*!
 * \brief Gets the number of dimensions in the mesh file for SU2_SMC.
 * \author T. Economon
 * \param[in] val_mesh_filename - Name of the file with the grid information.
 * \param[in] val_format - Format of the file with the grid information.
 * \return Physical dimension of the mesh (2-D or 3-D).
 */
unsigned short GetnDim(string val_mesh_filename, unsigned short val_format);

/*!
 * \brief Perform the geometric preprocessing for SU2_SMC.
 * \author T. Economon
 * \param[in] geometry_container - Geometric definition of the problem.
 * \param[in] config_container - Definition of the particular problem.
 */
void Geometric_Definition(CGeometry *geometry_container, CConfig *config_container);

/*!
 * \brief Match the zones at each sliding interface on a point-by-point basis.
 * \author T. Economon
 * \param[in] geometry_container - Geometric definition of the problem.
 * \param[in] config_container - Definition of the particular problem.
 * \param[in] nZone - Total number of zones in the mesh.
 */
void MatchSliding_Interfaces(CGeometry **geometry_container, CConfig **config_container, unsigned short nZone);

/*!
 * \brief .
 * \author T. Economon
 * \param[in] geometry_container - Geometric definition of the problem.
 * \param[in] config_container - Definition of the particular problem.
 * \param[in] nZone - Total number of zones in the mesh.
 */
//void ConstructSliding_Interfaces();

/*!
 * \brief .
 * \author T. Economon
 * \param[in] geometry_container - Geometric definition of the problem.
 * \param[in] config_container - Definition of the particular problem.
 * \param[in] nZone - Total number of zones in the mesh.
 */
//void CreateNew_Zones();

/*!
 * \brief Write a Tecplot file of the current multi-zone geometry.
 * \author T. Economon
 * \param[in] config_filename - Name of the file where the Tecplot
 *            information is going to be stored.
 */
void SetTecplot(CPoint ***node, unsigned long *nPoint, CPrimalGrid ***elem, 
                unsigned long *nElem, CConfig **config_container, unsigned short nZone, unsigned short nDim);

/*!
 * \brief Write an SU2 mesh file containing the new sliding interfaces built with SU2_SMC.
 * \author T. Economon
 * \param[in] config - Definition of the particular problem.
 * \param[in] val_mesh_out_filename - Name of the output file.
 */
void SetMeshFile(CPoint ***node, unsigned long *nPoint, CPrimalGrid ***elem,
                 unsigned long *nElem, CPrimalGrid**** bound, unsigned long **nElem_Bound,
                 CPrimalGrid**** newBound, unsigned long **nNewElem_Bound,
                 CConfig **config_container, unsigned short nZone, unsigned short nDim);
