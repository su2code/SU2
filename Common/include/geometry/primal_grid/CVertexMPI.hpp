/*!
 * \file CVertexMPI.hpp
 * \brief Headers of the main subroutines for storing the primal grid structure.
 *        The subroutines and functions are in the <i>primal_grid_structure.cpp</i> file.
 * \author F. Palacios
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

#include "CPrimalGrid.hpp"

/*! \class CVertexMPIConnectivity
 * Defines the connectivity of a vertex element for parallelization.
 * See CPrimalGridWithConnectivity.
 */
struct CVertexMPIConnectivity {
  // All the index directions should actually have 0 elements, but we cannot declare arrays of length 0.
  enum { nNodes = 1 };
  enum { nFaces = 0 };
  enum { maxNodesFace = 0 };
  enum { VTK_Type = VERTEX };
  static constexpr unsigned short nNodesFace[1] = {0};
  static constexpr unsigned short Faces[1][1] = {{0}};
  static constexpr unsigned short nNeighbor_Nodes[1] = {0};
  static constexpr unsigned short Neighbor_Nodes[1][1] = {{0}};
};

/*!
 * \class CVertexMPI
 * \brief Class for vertex element definition. This kind
 *        of element is used in the parallelization stuff.
 * \author F. Palacios
 */
class CVertexMPI final : public CPrimalGridWithConnectivity<CVertexMPIConnectivity> {
 private:
  /*! \brief Definition of the rotation, translation of the solution at the vertex. */
  unsigned short Rotation_Type;

 public:
  /*!
   * \brief Constructor using the nodes and index.
   * \param[in] val_point - Index of the 1st triangle point read from the grid file.
   */
  CVertexMPI(unsigned long val_point);

  /*!
   * \brief Get the type of rotation/traslation that must be applied.
   * \return Type of the element using VTK nomenclature.
   */
  inline unsigned short GetRotation_Type(void) const override { return Rotation_Type; }

  /*!
   * \brief Set the type of rotation/traslation that must be applied.
   * \param[in] val_rotation_type - Kind of rotation/traslation that must be applied.
   */
  inline void SetRotation_Type(unsigned short val_rotation_type) override { Rotation_Type = val_rotation_type; }

  inline void Change_Orientation() override {}
};
