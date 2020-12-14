/*!
 * \file CFEMStandardInternalFaceGrid.hpp
 * \brief Class for the FEM standard element for internal matching faces
 *        for the grid.
 *        The functions are in the <i>CFEMStandardInternalFaceGrid.cpp</i> file.
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

#include "CFEMStandardElementBase.hpp"

/*!
 * \class CFEMStandardInternalFaceGrid
 * \brief Class which defines the variables and methods for the standard
 *        standard element for an internal matching face for the grid.
 * \author E. van der Weide
 * \version 7.0.8 "Blackbird"
 */
class CFEMStandardInternalFaceGrid final {
private:
  CFEMStandardElementBase *elem0;  /*!< \brief Standard element on side 0 of the internal matching face. */
  CFEMStandardElementBase *elem1;  /*!< \brief Standard element on side 1 of the internal matching face. */

public:
  /*!
   * \brief Default constructor of the class, deleted to make sure the
   *        overloaded constructor is always used.
   */
  CFEMStandardInternalFaceGrid() = delete;

  /*!
   * \overload
   * \param[in] val_elem0 - Standard surface elements for side 0.
   * \param[in] val_elem1 - Standard surface elements for side 1
   */
  CFEMStandardInternalFaceGrid(CFEMStandardElementBase *val_elem0,
                               CFEMStandardElementBase *val_elem1);

  /*!
   * \brief Destructor. Nothing to be done.
   */
  ~CFEMStandardInternalFaceGrid() = default;

private:

};
