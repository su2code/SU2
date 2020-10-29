/*!
 * \file CVolumeElementFEM_Base.hpp
 * \brief Base class for a volume element for the FEM solver.
 *        The implementations are in the <i>CVolumeElementFEM_Base.cpp</i> file.
 * \author E. van der Weide
 * \version 7.0.7 "Blackbird"
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

#include "../../fem/CFEMStandardElementBase.hpp"

using namespace std;

/*!
 * \class CVolumeElementFEM_Base
 * \brief Base class to store a volume element for the FEM solver.
 * \author E. van der Weide
 * \version 7.0.7 "Blackbird"
 */
class CVolumeElementFEM_Base {
public:
  bool elemIsOwned;             /*!< \brief Whether or not this is an owned element. */
  bool JacIsConsideredConstant; /*!< \brief Whether or not the Jacobian of the transformation
                                            to the standard element is considered constant. */

  int rankOriginal;            /*!< \brief The rank where the original volume is stored. For
                                           the owned volumes, this is simply the current rank. */

  unsigned short VTK_Type;     /*!< \brief Element type using the VTK convention. */
  unsigned short nPolyGrid;    /*!< \brief Polynomial degree for the geometry of the element. */
  unsigned short nPolySol;     /*!< \brief Polynomial degree for the solution of the element. */
  unsigned short nDOFsGrid;    /*!< \brief Number of DOFs for the geometry of the element. */
  unsigned short nDOFsSol;     /*!< \brief Number of DOFs for the solution of the element. */
  unsigned short nFaces;       /*!< \brief Number of faces of the element. */
  unsigned short timeLevel;    /*!< \brief Time level of the element when time accurate local
                                           time stepping is employed. */

  unsigned int factTimeLevel;  /*!< \brief Number of local time steps for this element
                                           compared to the largest time step when time
                                           accurate local time stepping is employed. */

  unsigned long elemIDGlobal;        /*!< \brief Global element ID of this element. */

  vector<bool> JacFacesIsConsideredConstant; /*!< \brief Vector with the booleans whether the Jacobian of the
                                                         transformation to the standard element is constant
                                                         for the faces. */
  vector<bool> ElementOwnsFaces;             /*!< \brief Vector with the booleans whether the element is the
                                                         owner of the faces. */

  vector<unsigned long> nodeIDsGrid; /*!< \brief Vector with the node IDs of the grid for this element. */

  su2double lenScale;                /*!< \brief Length scale of the element. */


  CFEMStandardElementBase *standardElemGrid; /*!< \brief Pointer to the standard element for the grid. */
};
