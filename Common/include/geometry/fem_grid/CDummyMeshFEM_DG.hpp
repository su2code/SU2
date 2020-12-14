/*!
 * \file CDummyMeshFEM_DG.hpp
 * \brief Class definition for a dummy mesh object for the DG-FEM solver.
 *        The implementations are in the <i>CDummyMeshFEM_DG.cpp</i> file.
 * \author E. van der Weide
 * \version 7.0.8 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2020, SU2 Contributors (cf. AUTHORS.md)
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

#include "CMeshFEM_DG.hpp"

using namespace std;

/*!
 * \class CDummyMeshFEM_DG
 * \brief Class for defining a DG geometry that does not contain any points/elements.
 *        Can be used for initializing other classes that depend on the geometry without
 *        going through the time-consuming mesh initialization and paritioning.
 * \author T. Albring
 */
class CDummyMeshFEM_DG : public CMeshFEM_DG {

public:
  /*!
   * \brief Constructor of the class
   * \param[in] config - Definition of the particular problem.
   */
  CDummyMeshFEM_DG(CConfig *config);

};
