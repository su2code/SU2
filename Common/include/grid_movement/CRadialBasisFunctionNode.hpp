/*!
 * \file CRadialBasisFunctionNode.hpp
 * \brief Declaration of the RBF node class that stores nodal information for the RBF interpolation.
 * \author F. van Steen
 * \version 8.0.1 "Harrier"
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

#include "../../../Common/include/geometry/CGeometry.hpp"

/*!
 * \class CRadialBasisFunctionNode
 * \brief Class for defining a Radial Basis Function Node
 * \author F. van Steen
 */

class CRadialBasisFunctionNode{
 protected:
  unsigned long idx;          /*!< \brief Global index. */
  unsigned short marker_idx;  /*!< \brief Marker index. */
  unsigned long vertex_idx;   /*!< \brief Vertex index. */
    

 public:

  /*!
  * \brief Constructor of the class.
  */
  CRadialBasisFunctionNode(unsigned long idx_val, unsigned short marker_val, unsigned long vertex_val);

  /*!
  * \brief Returns global index.
  */
  inline unsigned long GetIndex(void){return idx;}

  /*!
  * \brief Returns vertex index.
  */
  inline unsigned long GetVertex(void){return vertex_idx;}

  /*!
  * \brief Returns marker index.
  */
  inline unsigned short GetMarker(void){return marker_idx;}

};