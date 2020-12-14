/*!
 * \file CFEMStandardLine.hpp
 * \brief Base class for the FEM line standard element.
 *        The functions are in the <i>CFEMStandardLine.cpp</i> file.
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
 * \class CFEMStandardLine
 * \brief Base class which defines the variables and methods for the
 *        line standard element.
 * \author E. van der Weide
 * \version 7.0.8 "Blackbird"
 */
class CFEMStandardLine: public CFEMStandardElementBase {

public:
  /*!
   * \brief Default constructor of the class, deleted to make sure the
   *        overloaded constructor is always used.
   */
  CFEMStandardLine() = delete;

  /*!
   * \overload
   * \param[in] val_nPoly      - Polynomial degree of the grid for this element.
   * \param[in] val_orderExact - Polynomial order that must be integrated exactly
   *                             by the integration rule.
   */
  CFEMStandardLine(const unsigned short val_nPoly,
                   const unsigned short val_orderExact);

  /*!
   * \brief Destructor. Nothing to be done.
   */
  virtual ~CFEMStandardLine() = default;

protected:

  vector<passivedouble> rLineDOFsEqui; /*!< \brief 1D parametric coordinates of the grid
                                                   DOFs when equidistant spacing is used. */
  vector<passivedouble> rLineDOFsLGL;  /*!< \brief 1D parametric coordinates of the grid
                                                   DOFs when the LGL distribution is used. */
  vector<passivedouble> rLineInt;      /*!< \brief 1D parametric coordinates of the
                                                   integration points. */
  vector<passivedouble> wLineInt;      /*!< \brief Weights of the 1D integration points. */
};
