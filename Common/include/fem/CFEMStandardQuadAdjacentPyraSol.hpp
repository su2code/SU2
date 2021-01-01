/*!
 * \file CFEMStandardQuadAdjacentPyraSol.hpp
 * \brief Class for the FEM standard surface quadrilateral element
 *        adjacent to a pyramid for the solution.
 *        The functions are in the <i>CFEMStandardQuadAdjacentPyraSol.cpp</i> file.
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

#include "CFEMStandardPyraBase.hpp"
#include "CGemmBase.hpp"

/*!
 * \class CFEMStandardQuadAdjacentPyraSol.
 * \brief Class which defines the variables and methods for the quadrilateral
 *        standard surface element adjacent to a pyramid for the solution.
 * \author E. van der Weide
 * \version 7.0.8 "Blackbird"
 */
class CFEMStandardQuadAdjacentPyraSol final: public CFEMStandardPyraBase {

public:
  /*!
   * \brief Default constructor of the class, deleted to make sure the
   *        overloaded constructor is always used.
   */
  CFEMStandardQuadAdjacentPyraSol() = delete;

  /*!
   * \overload
   * \param[in] val_nPoly       - Polynomial degree of the grid for this element.
   * \param[in] val_orderExact  - Polynomial degree that must be integrated exactly.
   * \param[in] val_faceID_Elem - This is the face ID of the adjacent volume element
   *                              to which this surface element corresponds.
   * \param[in] val_orientation - Orientation of this surface element relative to
   *                              the adjacent volume element.
   * \param[in] val_gemm_1      - Pointer to the gemm type used to carry out the
   *                              first gemm functionality for this standard face.
   * \param[in] val_gemm_2      - Pointer to the gemm type used to carry out the
   *                              second gemm functionality for this standard face.
   */
  CFEMStandardQuadAdjacentPyraSol(const unsigned short val_nPoly,
                                  const unsigned short val_orderExact,
                                  const unsigned short val_faceID_Elem,
                                  const unsigned short val_orientation,
                                  CGemmBase           *val_gemm_1,
                                  CGemmBase           *val_gemm_2);

  /*!
   * \brief Destructor. Nothing to be done.
   */
  ~CFEMStandardQuadAdjacentPyraSol() = default;

private:

};
