/*!
 * \file CFEMStandardFaceTriGrid.hpp
 * \brief Class for the FEM surface triangle standard element for the grid.
 *        The functions are in the <i>CFEMStandardFaceTriGrid.cpp</i> file.
 * \author E. van der Weide
 * \version 7.0.6 "Blackbird"
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

#include "CFEMStandardTri.hpp"

/*!
 * \class CFEMStandardFaceTriGrid
 * \brief Class which defines the variables and methods for the surface
 *        triangle standard element for the grid.
 * \author E. van der Weide
 * \version 7.0.6 "Blackbird"
 */
class CFEMStandardFaceTriGrid final: public CFEMStandardTri {

public:
  /*!
   * \brief Default constructor of the class, deleted to make sure the
   *        overloaded constructor is always used.
   */
  CFEMStandardFaceTriGrid() = delete;

  /*!
   * \overload
   * \param[in] val_nPoly      - Polynomial degree of the grid for this element.
   * \param[in] val_orderExact - Polynomial order that must be integrated exactly
   *                             by the integration rule.
   */
  CFEMStandardFaceTriGrid(const unsigned short val_nPoly,
                          const unsigned short val_orderExact);

  /*!
   * \brief Destructor. Nothing to be done.
   */
  ~CFEMStandardFaceTriGrid() = default;

private:

};
