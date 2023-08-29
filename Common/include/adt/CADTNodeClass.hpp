/*!
 * \file CADTNodeClass.hpp
 * \brief Class for storing the information needed in a node of an ADT.
 * \author E. van der Weide
 * \version 8.0.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2023, SU2 Contributors (cf. AUTHORS.md)
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
#include "../basic_types/datatype_structure.hpp"

/*!
 * \class CADTNodeClass
 * \ingroup ADT
 * \brief  Class for storing the information needed in a node of an ADT.
 * \author E. van der Weide
 */
struct CADTNodeClass {
  bool childrenAreTerminal[2]; /*!< \brief Whether or not the child leaves are terminal. */
  unsigned long children[2];   /*!< \brief Child leaves. If childrenAreTerminal is true the children
                                           contain the point ID's or bounding box ID's. Note that it
                                           is allowed that one child is termimal and the other is not. */
  unsigned long centralNodeID; /*!< \brief ID of a node, which is near the center of the leaf. */

  su2double* xMin; /*!< \brief The minimum coordinates of this leaf. It points to a position in the large
                               vector, which contains the coordinates of all leaves. */
  su2double* xMax; /*!< \brief The maximum coordinates of this leaf. It points to a position in the large
                               vector, which contains the coordinates of all leaves. */
};
